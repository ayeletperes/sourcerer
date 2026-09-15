"""
OGRDB (AIRR Community germline sets)

OGRDB publishes curated germline sets through a small REST API. Resolving a set
to a download takes three calls -- species to a numeric id, id to the sets it
holds, set to its latest release -- after which the FASTA is fetched twice, once
ungapped and once IMGT-gapped, because a set carries V, D, J and C together and
each segment is taken from the form that suits it.

The segment split is the load-bearing piece and is kept verbatim from airrdb:
V is taken gapped, so Change-O keeps the IMGT numbering it needs; D and J are
taken ungapped; and constant regions are taken gapped. The one ambiguity is the
delta locus, where the diversity segment IGHD and the constant IGHD share a
name; they are told apart by length, since the constant is far longer. Getting
this wrong silently files an allele under the wrong segment, so it is covered by
fixtures.

OGRDB serves only immunoglobulin sets, and only for the species it has curated,
so a TR request or an uncovered locus resolves to nothing rather than an error.
"""

# Info
__author__ = 'Ayelet Peres'

# Imports
import logging
import re
from datetime import UTC
from urllib.parse import quote

# Sourcerer imports
from sourcerer.Exceptions import OgrdbParseError
from sourcerer.Reference import (
    KIND_CONSTANT,
    KIND_VDJ,
    ReferenceReport,
    parseFasta,
)
from sourcerer.Sources.Base import DataUnit
from sourcerer.Sources.Germline import ReferenceSource

log = logging.getLogger(__name__)

#: Endpoint.
API = 'https://ogrdb.airr-community.org/api_v2'

#: Species as OGRDB labels them.
SPECIES_LABEL = {'human': 'Homo sapiens', 'mouse': 'Mus musculus'}

#: The germline sets to fetch per species and locus, and the chains each covers.
#: A locus can need more than one set: mouse splits V and J across strain-specific
#: and all-strain sets. Kept as data so a new set is one line, not new code.
SETS = {
    ('human', 'IGH'): (('IGH_VDJ', ('IGHV', 'IGHD', 'IGHJ')), ('IGHC', ('IGHC',))),
    ('human', 'IGK'): (('IGKappa_VJ', ('IGKV', 'IGKJ')),),
    ('human', 'IGL'): (('IGLambda_VJ', ('IGLV', 'IGLJ')),),
    ('mouse', 'IGH'): (('C57BL/6 IGH', ('IGHV', 'IGHD', 'IGHJ')),),
    ('mouse', 'IGK'): (('C57BL/6J IGKV', ('IGKV',)),
                       ('IGKJ (all strains)', ('IGKJ',))),
    ('mouse', 'IGL'): (('C57BL/6J IGLV', ('IGLV',)),
                       ('IGLJ (all strains)', ('IGLJ',))),
}

#: Loci OGRDB covers, offered as a filter flag.
LOCI = ('IGH', 'IGK', 'IGL')

#: The two forms each set is fetched in.
FORMATS = ('ungapped', 'gapped')

#: A constant region under 100 nucleotides is really the delta diversity segment
#: wearing the same name; see the module docstring.
CONSTANT_MIN_LENGTH = 100


def normalizeVersion(value):
    """
    Render a release version without a trailing '.0'.

    OGRDB reports the version as a number, so an integer release arrives as
    '3.0' where the download URL wants '3'.

    Arguments:
      value: the reported version.

    Returns:
      str: the version as it appears in a download URL.
    """
    text = str(value)

    return text[:-2] if text.endswith('.0') else text


def safeSetName(set_name):
    """
    Make a filesystem-safe token from a set name.

    Set names carry spaces and slashes (``C57BL/6J IGKV``) that must not become
    directory separators in the raw mirror, but the name is never parsed back:
    the real set name travels in the unit metadata.

    Arguments:
      set_name (str): the OGRDB set name.

    Returns:
      str: an identifier-safe token.
    """
    return re.sub(r'[^A-Za-z0-9]+', '_', set_name).strip('_')


def bucketChain(name, sequence):
    """
    Decide which reference chain a germline allele belongs to.

    V, D and J are filed under their four-character chain (``IGHV``); a constant
    allele is filed under its locus constant (``IGHM`` -> ``IGHC``) so every
    isotype of a locus lands in one file, as airrflow expects. Returns None for a
    name too short to classify.

    Arguments:
      name (str): the allele name.
      sequence (str): its sequence, used only to tell the delta segments apart.

    Returns:
      tuple: (chain, kind) or None.
    """
    if len(name) < 4:
        return None

    segment = name[3]
    if segment == 'V':
        return name[:4], KIND_VDJ
    if segment == 'J':
        return name[:4], KIND_VDJ
    if segment == 'D':
        # IGHD is both the diversity segment (short) and the delta constant
        # (long); length is the only thing that separates them.
        if len(sequence.replace('.', '')) < CONSTANT_MIN_LENGTH:
            return name[:4], KIND_VDJ
        return name[:3] + 'C', KIND_CONSTANT

    return name[:3] + 'C', KIND_CONSTANT


def _splitSegments(forms):
    """
    Sort a set's alleles into reference chains, form by form.

    V is taken from the gapped alleles, so its IMGT numbering survives, along
    with the constant regions; D and J are taken from the ungapped alleles. An
    allele that appears in both forms is therefore filed once, from the form its
    segment is read from.

    Arguments:
      forms (dict): 'ungapped' and 'gapped' each mapping allele name to sequence.

    Returns:
      dict: (chain, kind) to a list of (name, sequence) tuples.
    """
    chains = {}
    for name, sequence in forms.get('gapped', {}).items():
        target = bucketChain(name, sequence)
        if target is None:
            continue
        chain, kind = target
        if kind == KIND_CONSTANT or chain[3] == 'V':
            chains.setdefault(target, []).append((name, sequence))

    for name, sequence in forms.get('ungapped', {}).items():
        target = bucketChain(name, sequence)
        if target is None:
            continue
        chain, kind = target
        if kind == KIND_VDJ and chain[3] in ('D', 'J'):
            chains.setdefault(target, []).append((name, sequence))

    return chains


class OgrdbSource(ReferenceSource):
    """
    The OGRDB (AIRR Community) germline reference source.
    """

    name = 'ogrdb'
    aliases = ('airrc',)
    prefix = 'airrc'
    description = 'OGRDB: AIRR Community curated immunoglobulin germline sets'
    homepage = 'https://ogrdb.airr-community.org/'
    collections = ('human', 'mouse')
    collection_help = {'human': 'Homo sapiens curated IG sets',
                       'mouse': 'Mus musculus curated IG sets'}

    license = 'CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)'
    citation = (
        'Lees WD, Busse CE, Corcoran M, et al. OGRDB: a reference database of '
        'inferred immune receptor genes. Nucleic Acids Res. '
        '2020;48(D1):D964-D970. doi:10.1093/nar/gkz822',
    )

    # -- API client -------------------------------------------------------

    def speciesId(self, species_label):
        """
        Resolve a species label to its OGRDB id.

        Arguments:
          species_label (str): the label, e.g. 'Homo sapiens'.

        Returns:
          str: the species id.

        Raises:
          OgrdbParseError: if the species is not listed.
        """
        payload = self.client.get('%s/germline/species' % API).json()
        for item in payload.get('species', []):
            if item.get('label') == species_label:
                return item['id']

        raise OgrdbParseError('OGRDB does not list species %r; the species '
                              'endpoint changed or the species was withdrawn'
                              % species_label)

    def resolveSetId(self, species_id, locus, set_name):
        """
        Resolve a set name to its germline set id.

        Arguments:
          species_id (str): the OGRDB species id.
          locus (str): the locus, e.g. 'IGH'.
          set_name (str): the set name.

        Returns:
          str: the germline set id.

        Raises:
          OgrdbParseError: if the set is not found for the species and locus.
        """
        payload = self.client.get('%s/germline/sets/%s' % (API, species_id)).json()
        for item in payload.get('germline_species', []):
            if (item.get('germline_set_name') == set_name
                    and item.get('locus') == locus):
                return item['germline_set_id']

        raise OgrdbParseError("OGRDB has no set %r for locus %s; the set was "
                              'renamed or withdrawn' % (set_name, locus))

    def latestRelease(self, set_id):
        """
        Read the latest release version and date of a set.

        Arguments:
          set_id (str): the germline set id.

        Returns:
          tuple: (version, release_date) with the date truncated to YYYY-MM-DD.

        Raises:
          OgrdbParseError: if the release payload is not shaped as expected.
        """
        safe = quote(set_id, safe='.')
        payload = self.client.get('%s/germline/set/%s/latest' % (API, safe)).json()
        try:
            record = payload['GermlineSet'][0]

            return normalizeVersion(record['release_version']), \
                record['release_date'][:10]
        except (KeyError, IndexError, TypeError) as error:
            raise OgrdbParseError('OGRDB latest-release payload for %s is not '
                                  'shaped as expected (%s)' % (set_id, error))

    def fastaUrl(self, set_id, version, fmt, human):
        """
        Build a set's FASTA download URL.

        Arguments:
          set_id (str): the germline set id.
          version (str): the release version.
          fmt (str): 'ungapped' or 'gapped'.
          human (bool): whether the species is human, which takes the ``_ex``
            endpoint variant.

        Returns:
          str: the absolute download URL.
        """
        safe = quote(set_id, safe='.')
        suffix = '_ex' if human else ''

        return '%s/germline/set/%s/%s/%s%s' % (API, safe, version, fmt, suffix)

    # -- schema -----------------------------------------------------------

    def harvestSchema(self):
        """
        Contact OGRDB and build a fresh snapshot of the loci it curates.

        The species and sets endpoints are queried, so a refresh both verifies
        the API is answering and records which of the loci this source consumes
        are actually available -- the drift signal that matters for OGRDB.

        Returns:
          SourceSchema: the harvested snapshot.
        """
        from datetime import datetime

        from sourcerer.Schema import Collection, Field, SourceSchema
        from sourcerer.Version import __version__

        collections = {}
        for sp in self.collections:
            species_id = self.speciesId(SPECIES_LABEL[sp])
            payload = self.client.get('%s/germline/sets/%s'
                                      % (API, species_id)).json()
            available = {x.get('locus') for x in payload.get('germline_species', [])}
            loci = tuple(x for x in LOCI if x in available)
            collections[sp] = Collection(name=sp,
                                         fields=(Field(name='locus', values=loci),))

        return SourceSchema(
            source=self.name,
            harvested=datetime.now(UTC).strftime('%Y-%m-%dT%H:%M:%SZ'),
            harvested_by='sourcerer %s' % __version__,
            source_urls={'api': API},
            parse_contracts={'segment_split': 'V,C gapped; D,J ungapped; '
                                              'delta D vs C by length'},
            collections=collections)

    # -- search and build -------------------------------------------------

    def searchUnits(self, query):
        """
        Resolve a query to the germline files to fetch.

        Each set is fetched twice, ungapped and gapped, so both are present when
        the segments are split in buildReference. The set id and latest version
        are resolved here so the download URLs are concrete.

        Arguments:
          query (Query): the validated request; collection is the species and the
            locus filter narrows which sets are fetched.

        Returns:
          list: DataUnit objects, two per set.
        """
        species = query.collection
        label = SPECIES_LABEL[species]
        human = species == 'human'
        locus_filter = query.filters.get('locus', '*')

        species_id = self.speciesId(label)
        units = []
        for (set_species, locus), sets in SETS.items():
            if set_species != species:
                continue
            if locus_filter not in ('*', locus):
                continue

            for set_name, chains in sets:
                set_id = self.resolveSetId(species_id, locus, set_name)
                version, _date = self.latestRelease(set_id)
                for fmt in FORMATS:
                    units.append(DataUnit(
                        unit_id='%s.%s.fasta' % (safeSetName(set_name), fmt),
                        collection=species,
                        url=self.fastaUrl(set_id, version, fmt, human),
                        metadata={'species': species, 'locus': locus,
                                  'set_name': set_name, 'format': fmt,
                                  'set_id': set_id, 'version': version,
                                  'chains': list(chains)}))

        if query.limit is not None:
            units = units[:query.limit]

        return units

    def buildReference(self, entries, reference_dir):
        """
        Split the downloaded sets into per-chain reference FASTAs.

        Arguments:
          entries (list): (DataUnit, Path) pairs from the fetch step.
          reference_dir (Path): the reference_base root.

        Returns:
          ReferenceReport: the files written.
        """
        report = ReferenceReport()
        for species in sorted({unit.metadata['species'] for unit, _ in entries}):
            forms = {fmt: {} for fmt in FORMATS}
            for unit, path in entries:
                if unit.metadata['species'] != species:
                    continue
                fmt = unit.metadata['format']
                for header, sequence in parseFasta(path.read_text()):
                    forms[fmt][header.split()[0]] = sequence

            chains = _splitSegments(forms)
            for (chain, kind), records in sorted(chains.items()):
                written = self.writeChain(reference_dir, species, kind, chain,
                                          records)
                report.written.append((chain, written))
                log.info('%s: %d sequences', written.name, len(records))

        return report
