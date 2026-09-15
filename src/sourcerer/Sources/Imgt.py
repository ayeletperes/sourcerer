"""
IMGT/GENE-DB

IMGT exposes no data API. Its GENE-DB GENElect page takes a numbered query and a
chain and returns an HTML page with the FASTA embedded in its second ``<pre>``
block; the first holds the query echo. A query can fail and still return HTTP
200, so a valid answer is not the status alone but the presence of that second
block with sequence in it, which is what isValidResponse checks and the weekly
API canary relies on.

One germline file per (species, chain) is fetched, matching how airrflow's
bin/fetch_references.sh drives GENElect: query 7.14 for V, D and J nucleotide,
14.1 for constant (7.5 for the mouse kappa and lambda constants, which 14.1 does
not serve), and 7.3 for the translated V. The files are written into the
reference_base with their IMGT headers and gaps intact; see sourcerer.Reference
for why cleaning is deferred to the database build.
"""

# Info
__author__ = 'Ayelet Peres'

# Imports
import logging
from datetime import UTC
from urllib.parse import quote

from bs4 import BeautifulSoup

# Sourcerer imports
from sourcerer.Exceptions import ImgtParseError
from sourcerer.Reference import (
    KIND_AA,
    KIND_CONSTANT,
    KIND_VDJ,
    ReferenceReport,
    parseFasta,
)
from sourcerer.Sources.Base import DataUnit
from sourcerer.Sources.Germline import ReferenceSource

log = logging.getLogger(__name__)

#: Endpoints.
GENELECT = 'https://www.imgt.org/genedb/GENElect'
RELEASE_URL = 'https://www.imgt.org/download/GENE-DB/RELEASE'

#: GENElect query numbers, per chain kind.
Q_VDJ = '7.14'          # V, D and J nucleotide
Q_CONSTANT = '14.1'     # constant nucleotide
Q_CONSTANT_MOUSE = '7.5'  # mouse IGKC and IGLC, which 14.1 does not serve
Q_AA = '7.3'            # translated V

#: Species as GENElect wants them in the query string, and as they appear in the
#: FASTA headers. The query form is pre-encoded so it is not double-escaped.
SPECIES_QUERY = {'human': 'Homo%20sapiens', 'mouse': 'Mus%20musculus'}
SPECIES_LABEL = {'human': 'Homo sapiens', 'mouse': 'Mus musculus'}

#: Chains fetched as V/D/J nucleotide.
VDJ_CHAINS = ('IGHV', 'IGHD', 'IGHJ', 'IGKV', 'IGKJ', 'IGLV', 'IGLJ',
              'TRAV', 'TRAJ', 'TRBV', 'TRBD', 'TRBJ',
              'TRDV', 'TRDD', 'TRDJ', 'TRGV', 'TRGJ')

#: Chains fetched as constant nucleotide.
CONSTANT_CHAINS = ('IGHC', 'IGKC', 'IGLC', 'TRAC', 'TRBC', 'TRGC', 'TRDC')

#: Chains fetched as translated V.
AA_CHAINS = ('IGHV', 'IGKV', 'IGLV', 'TRAV', 'TRBV', 'TRDV', 'TRGV')

#: Loci and segments a search can be narrowed to, offered as filter flags.
LOCI = ('IGH', 'IGK', 'IGL', 'TRA', 'TRB', 'TRG', 'TRD')
SEGMENTS = ('V', 'D', 'J', 'C')


def buildQueryUrl(species, query, chain, label=None):
    """
    Build a GENElect query URL.

    Arguments:
      species (str): the species key, e.g. 'human'.
      query (str): the GENElect query number, e.g. '7.14'.
      chain (str): the chain, e.g. 'IGHV'.
      label (str): an optional IMGTlabel qualifier.

    Returns:
      str: the absolute query URL.
    """
    url = '%s?query=%s+%s&species=%s' % (GENELECT, quote(query), chain,
                                         SPECIES_QUERY[species])
    if label:
        url += '&IMGTlabel=%s' % label

    return url


def isValidResponse(html):
    """
    Report whether a GENElect reply actually carries a germline FASTA.

    GENElect answers a failed query with HTTP 200 and an error page, so a live
    check cannot trust the status code. A real answer has a second ``<pre>``
    block, and that block has sequence in it. Both conditions are required.

    Arguments:
      html (str): the GENElect reply body.

    Returns:
      bool: True if the reply contains a non-empty germline FASTA.
    """
    blocks = BeautifulSoup(html, 'html.parser').find_all('pre')
    if len(blocks) < 2:
        return False

    return '>' in blocks[1].get_text()


def extractFasta(html, species):
    """
    Pull the germline FASTA out of a GENElect reply.

    The FASTA is the second ``<pre>`` block. The species name in the headers has
    its spaces replaced with underscores, as airrflow does, so a header stays one
    whitespace-delimited field.

    Arguments:
      html (str): the GENElect reply body.
      species (str): the species key, for the header rewrite.

    Returns:
      str: the FASTA text.

    Raises:
      ImgtParseError: if the reply has no second ``<pre>`` block, which means the
        query failed or the page layout changed.
    """
    blocks = BeautifulSoup(html, 'html.parser').find_all('pre')
    if len(blocks) < 2:
        raise ImgtParseError(
            'GENElect reply has fewer than two <pre> blocks; the query failed or '
            'the page layout changed')

    text = blocks[1].get_text()
    label = SPECIES_LABEL.get(species)
    if label:
        text = text.replace(label, label.replace(' ', '_'))

    return text


def _chainPlan(species):
    """
    Enumerate every (chain, kind, query, locus, segment) this source fetches.

    Arguments:
      species (str): the species key, used to route the mouse constant queries.

    Returns:
      list: dicts describing one germline file each.
    """
    plan = []
    for chain in VDJ_CHAINS:
        plan.append({'chain': chain, 'kind': KIND_VDJ, 'query': Q_VDJ,
                     'locus': chain[:3], 'segment': chain[3]})
    for chain in CONSTANT_CHAINS:
        query = Q_CONSTANT
        if species == 'mouse' and chain in ('IGKC', 'IGLC'):
            query = Q_CONSTANT_MOUSE
        plan.append({'chain': chain, 'kind': KIND_CONSTANT, 'query': query,
                     'locus': chain[:3], 'segment': 'C'})
    for chain in AA_CHAINS:
        plan.append({'chain': chain, 'kind': KIND_AA, 'query': Q_AA,
                     'locus': chain[:3], 'segment': 'V'})

    return plan


class ImgtSource(ReferenceSource):
    """
    The IMGT/GENE-DB germline reference source.
    """

    name = 'imgt'
    prefix = 'imgt'
    description = 'IMGT/GENE-DB: germline V, D, J and C reference sequences'
    homepage = 'https://www.imgt.org/genedb/'
    collections = ('human', 'mouse')
    collection_help = {'human': 'Homo sapiens germline reference',
                       'mouse': 'Mus musculus germline reference'}

    #: IMGT's reuse terms are not an open-data licence; germline data may be used
    #: for research on condition IMGT is cited. Recorded so a reader of a download
    #: directory sees the obligation without having to consult IMGT separately.
    license = ('IMGT terms of use (https://www.imgt.org/about/termsofuse.php); '
               'cite IMGT, the international ImMunoGeneTics information system')
    citation = (
        'Lefranc MP, Giudicelli V, Duroux P, et al. IMGT, the international '
        'ImMunoGeneTics information system 25 years on. Nucleic Acids Res. '
        '2015;43(Database issue):D413-D422. doi:10.1093/nar/gku1056',
    )

    def harvestSchema(self):
        """
        Contact IMGT for its current release and build a fresh snapshot.

        GENElect has no field-listing endpoint, so the searchable vocabulary is
        the fixed set of loci and segments this source knows how to query. The
        network call to the release file is what turns a refresh into a genuine
        liveness check rather than a rewrite of a constant.

        Returns:
          SourceSchema: the harvested snapshot.
        """
        from datetime import datetime

        from sourcerer.Schema import Collection, Field, SourceSchema
        from sourcerer.Version import __version__

        # A liveness check, not stored: the release tag changes with every IMGT
        # build, and keeping it in the snapshot would make a monthly refresh
        # rewrite a tracked file with no change in the vocabulary.
        log.info('IMGT release: %s', self.fetchRelease() or 'unknown')

        fields = (Field(name='locus', values=LOCI),
                  Field(name='segment', values=SEGMENTS))
        collections = {sp: Collection(name=sp, fields=fields)
                       for sp in self.collections}

        return SourceSchema(
            source=self.name,
            harvested=datetime.now(UTC).strftime('%Y-%m-%dT%H:%M:%SZ'),
            harvested_by='sourcerer %s' % __version__,
            source_urls={'genelect': GENELECT, 'release': RELEASE_URL},
            parse_contracts={'fasta_block': 'second <pre> element',
                             'header': 'pipe-delimited, allele in field 2'},
            collections=collections)

    def fetchRelease(self):
        """
        Return IMGT's current GENE-DB release tag.

        Returns:
          str: the release tag, e.g. '202619-7', or '' if it cannot be read.
        """
        try:
            return self.client.get(RELEASE_URL).text.strip()
        except Exception as error:
            log.warning('could not read the IMGT release tag: %s', error)
            return ''

    def searchUnits(self, query):
        """
        Resolve a query to the germline files to fetch.

        Arguments:
          query (Query): the validated request; collection is the species, and
            the locus and segment filters narrow which chains are fetched.

        Returns:
          list: DataUnit objects, one per germline file.
        """
        species = query.collection
        locus = query.filters.get('locus', '*')
        segment = query.filters.get('segment', '*')

        units = []
        for item in _chainPlan(species):
            if locus not in ('*', item['locus']):
                continue
            if segment not in ('*', item['segment']):
                continue

            url = buildQueryUrl(species, item['query'], item['chain'])
            unit_id = '%s/%s.html' % (item['kind'], item['chain'])
            units.append(DataUnit(
                unit_id=unit_id, collection=species, url=url,
                metadata={'species': species, 'chain': item['chain'],
                          'kind': item['kind'], 'locus': item['locus'],
                          'segment': item['segment'], 'query': item['query']}))

        if query.limit is not None:
            units = units[:query.limit]

        return units

    def buildReference(self, entries, reference_dir):
        """
        Extract each downloaded GENElect page into the reference tree.

        Arguments:
          entries (list): (DataUnit, Path) pairs from the fetch step.
          reference_dir (Path): the reference_base root.

        Returns:
          ReferenceReport: the files written.
        """
        report = ReferenceReport()
        for unit, path in entries:
            html = path.read_text()
            fasta = extractFasta(html, unit.metadata['species'])
            records = parseFasta(fasta)
            written = self.writeChain(reference_dir, unit.metadata['species'],
                                      unit.metadata['kind'],
                                      unit.metadata['chain'], records)
            report.written.append((unit.metadata['chain'], written))
            log.info('%s: %d sequences', written.name, len(records))

        return report
