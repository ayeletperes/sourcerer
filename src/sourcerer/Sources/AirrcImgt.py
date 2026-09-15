"""
AIRR-C germline sets blended with IMGT

A germline reference that takes each locus from the source that curates it best:
the immunoglobulin V, D and J from OGRDB's AIRR-C sets, and everything OGRDB does
not cover -- all of the T-cell receptor, and the immunoglobulin constants that
are not in a published set -- from IMGT. It is the reference nf-core/airrflow
builds for its ``airrc-imgt`` database type.

Rather than reimplement either source, this composes them: it asks the OGRDB
source for the immunoglobulin sets and the IMGT source for the gap, tags each
download with which one it came from, and lets each build its own files back into
one reference tree, so an OGRDB allele lands as ``airrc_...`` and an IMGT allele
as ``imgt_...`` exactly as they would from the sources alone. Amino acid V is not
included, matching what airrflow's airrc-imgt build uses.
"""

# Info
__author__ = 'Ayelet Peres'

# Imports
import logging
from datetime import UTC

# Sourcerer imports
from sourcerer.Reference import KIND_AA, ReferenceReport
from sourcerer.Sources.Base import Query
from sourcerer.Sources.Germline import ReferenceSource
from sourcerer.Sources.Imgt import ImgtSource
from sourcerer.Sources.Ogrdb import OgrdbSource

log = logging.getLogger(__name__)

#: The IMGT immunoglobulin constants to take per species: the ones OGRDB has no
#: set for. Human IGHC comes from OGRDB, so only the light constants are taken;
#: mouse has no constant set at all, so all three heavy and light come from IMGT.
IMGT_IG_CONSTANTS = {'human': ('IGKC', 'IGLC'),
                     'mouse': ('IGHC', 'IGKC', 'IGLC')}

#: Which download a unit came from, recorded so buildReference can hand each unit
#: back to the source that knows how to read it.
VIA = 'via'


class AirrcImgtSource(ReferenceSource):
    """
    OGRDB immunoglobulin sets blended with IMGT for TR and the IG constants.
    """

    name = 'airrc-imgt'
    description = ('AIRR-C immunoglobulin sets blended with IMGT for TR and the '
                   'remaining constants')
    homepage = 'https://ogrdb.airr-community.org/'
    collections = ('human', 'mouse')
    collection_help = {'human': 'Homo sapiens blended reference',
                       'mouse': 'Mus musculus blended reference'}
    license = ('OGRDB data under CC BY 4.0 and IMGT data under the IMGT terms of '
               'use; cite both OGRDB and IMGT')
    citation = OgrdbSource.citation + ImgtSource.citation

    def __init__(self, client, schema=None):
        """
        Arguments:
          client (HttpClient): the shared HTTP client, passed to both sources.
          schema (SourceSchema): the loaded snapshot, or None to load on demand.
        """
        super().__init__(client, schema)
        self._ogrdb = OgrdbSource(client)
        self._imgt = ImgtSource(client)

    def harvestSchema(self):
        """
        Build a snapshot for the blended source.

        The blend takes a whole species at a time and has no filters of its own,
        so the snapshot is just its collections; the drift checks that matter live
        on the imgt and ogrdb snapshots this composes.

        Returns:
          SourceSchema: the snapshot.
        """
        from datetime import datetime

        from sourcerer.Schema import Collection, SourceSchema
        from sourcerer.Version import __version__

        return SourceSchema(
            source=self.name,
            harvested=datetime.now(UTC).strftime('%Y-%m-%dT%H:%M:%SZ'),
            harvested_by='sourcerer %s' % __version__,
            source_urls={'ogrdb': self._ogrdb.name, 'imgt': self._imgt.name},
            collections={sp: Collection(name=sp) for sp in self.collections})

    def searchUnits(self, query):
        """
        Resolve a query to the OGRDB and IMGT files the blend needs.

        Arguments:
          query (Query): the validated request; collection is the species.

        Returns:
          list: DataUnit objects, each tagged with which source produced it.
        """
        species = query.collection

        units = []
        for unit in self._ogrdb.searchUnits(
                Query(collection=species, filters={'locus': '*'})):
            unit.metadata[VIA] = self._ogrdb.name
            units.append(unit)

        for unit in self._imgtGapUnits(species):
            unit.metadata[VIA] = self._imgt.name
            units.append(unit)

        if query.limit is not None:
            units = units[:query.limit]

        return units

    def _imgtGapUnits(self, species):
        """
        Pick the IMGT files that fill what OGRDB does not cover.

        That is every T-cell receptor chain, and the immunoglobulin constants
        without an OGRDB set; the immunoglobulin V, D and J come from OGRDB, and
        amino acid V is not part of this blend.

        Arguments:
          species (str): the species.

        Returns:
          list: DataUnit objects from the IMGT source.
        """
        constants = IMGT_IG_CONSTANTS.get(species, ())
        gap = []
        for unit in self._imgt.searchUnits(
                Query(collection=species,
                      filters={'locus': '*', 'segment': '*'})):
            meta = unit.metadata
            if meta['kind'] == KIND_AA:
                continue
            if meta['locus'].startswith('TR'):
                gap.append(unit)
            elif meta['kind'] == 'constant' and meta['chain'] in constants:
                gap.append(unit)

        return gap

    def buildReference(self, entries, reference_dir):
        """
        Let each source build its own files back into one reference tree.

        Arguments:
          entries (list): (DataUnit, Path) pairs from the fetch step.
          reference_dir (Path): the reference_base root.

        Returns:
          ReferenceReport: the files written by both sources.
        """
        ogrdb_entries = [pair for pair in entries
                         if pair[0].metadata.get(VIA) == self._ogrdb.name]
        imgt_entries = [pair for pair in entries
                        if pair[0].metadata.get(VIA) == self._imgt.name]

        report = ReferenceReport()
        report.written.extend(
            self._ogrdb.buildReference(ogrdb_entries, reference_dir).written)
        report.written.extend(
            self._imgt.buildReference(imgt_entries, reference_dir).written)

        return report
