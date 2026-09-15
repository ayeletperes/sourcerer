"""
Unit tests for the IMGT source
"""

# Info
__author__ = 'Ayelet Peres'

# Imports
import os
import tempfile
import unittest
from pathlib import Path

# Sourcerer imports
from sourcerer.Exceptions import ImgtParseError
from sourcerer.Sources.Base import DataUnit, Query
from sourcerer.Sources.Imgt import (
    ImgtSource,
    buildQueryUrl,
    extractFasta,
    isValidResponse,
)

test_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(test_path, 'data')


def readFixture(name):
    """Read a captured fixture from tests/data."""
    with open(os.path.join(data_path, name)) as handle:
        return handle.read()


class TestQueryUrl(unittest.TestCase):
    """
    Tests for GENElect URL construction
    """

    def test_encodes_query_and_species(self):
        """The query number is escaped and the species is sent pre-encoded."""
        url = buildQueryUrl('human', '7.14', 'IGHV')
        self.assertEqual(
            url,
            'https://www.imgt.org/genedb/GENElect?query=7.14+IGHV'
            '&species=Homo%20sapiens')

    def test_appends_label(self):
        """An IMGTlabel qualifier is appended when given."""
        url = buildQueryUrl('mouse', '8.1', 'IGHV', label='L-PART1+L-PART2')
        self.assertTrue(url.endswith('&IMGTlabel=L-PART1+L-PART2'))


class TestResponseParsing(unittest.TestCase):
    """
    Tests for reading a GENElect reply
    """

    def test_valid_response_has_second_pre_with_fasta(self):
        """A real reply has a second <pre> block carrying a FASTA."""
        self.assertTrue(isValidResponse(readFixture('imgt_ighd.html')))

    def test_error_page_is_not_valid_despite_http_200(self):
        """An error page with a single <pre> is rejected."""
        self.assertFalse(isValidResponse(readFixture('imgt_error.html')))

    def test_extract_returns_fasta_with_underscored_species(self):
        """extractFasta pulls the FASTA and underscores the species name."""
        fasta = extractFasta(readFixture('imgt_ighd.html'), 'human')
        self.assertIn('>X97051|IGHD1-1*01|Homo_sapiens|F', fasta)
        self.assertNotIn('Homo sapiens', fasta)
        self.assertEqual(fasta.count('>'), 3)

    def test_extract_raises_on_error_page(self):
        """A page with no second <pre> is a parse error, not an empty result."""
        with self.assertRaises(ImgtParseError):
            extractFasta(readFixture('imgt_error.html'), 'human')


class TestSearchUnits(unittest.TestCase):
    """
    Tests for enumerating the germline files to fetch
    """

    def setUp(self):
        self.source = ImgtSource(client=None)

    def _query(self, **filters):
        resolved = {'locus': '*', 'segment': '*'}
        resolved.update(filters)
        return Query(collection='human', filters=resolved)

    def test_unfiltered_covers_vdj_constant_and_aa(self):
        """With no filter every VDJ, constant and AA chain is scheduled."""
        units = self.source.searchUnits(self._query())
        # 17 VDJ + 7 constant + 7 AA V.
        self.assertEqual(len(units), 31)
        kinds = {u.metadata['kind'] for u in units}
        self.assertEqual(kinds, {'vdj', 'constant', 'vdj_aa'})

    def test_locus_and_segment_filter(self):
        """--locus IGH --segment V leaves the nucleotide and amino acid V."""
        units = self.source.searchUnits(self._query(locus='IGH', segment='V'))
        self.assertEqual({u.metadata['chain'] for u in units}, {'IGHV'})
        self.assertEqual({u.metadata['kind'] for u in units}, {'vdj', 'vdj_aa'})

    def test_mouse_light_constant_uses_special_query(self):
        """Mouse IGKC and IGLC take query 7.5, which 14.1 does not serve."""
        source = ImgtSource(client=None)
        units = source.searchUnits(
            Query(collection='mouse',
                  filters={'locus': 'IGK', 'segment': 'C'}))
        self.assertEqual(len(units), 1)
        self.assertIn('query=7.5+IGKC', units[0].url)


class TestBuildReference(unittest.TestCase):
    """
    Tests for extracting downloaded pages into the reference tree
    """

    def test_writes_chain_fasta_from_page(self):
        """A downloaded GENElect page becomes one per-chain reference FASTA."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            raw = tmp / 'IGHD.html'
            raw.write_text(readFixture('imgt_ighd.html'))
            unit = DataUnit(
                unit_id='vdj/IGHD.html', collection='human', url='x',
                metadata={'species': 'human', 'chain': 'IGHD',
                          'kind': 'vdj', 'locus': 'IGH', 'segment': 'D'})

            source = ImgtSource(client=None)
            report = source.buildReference([(unit, raw)], tmp / 'reference_base')

            written = (tmp / 'reference_base' / 'human' / 'vdj'
                       / 'imgt_human_IGHD.fasta')
            self.assertTrue(written.exists())
            self.assertEqual(written.read_text().count('>'), 3)
            self.assertEqual(len(report.written), 1)


class TestArchivePin(unittest.TestCase):
    """
    Tests for re-downloading a pinned release from the genedb-releases archive
    """

    class Listing:
        def __init__(self, names):
            self.payload = [{'name': n, 'type': 'dir'} for n in names]

        def json(self):
            return self.payload

    class ListingClient:
        def __init__(self, names):
            self.names = names

        def get(self, url):
            return TestArchivePin.Listing(self.names)

    def test_pinned_search_yields_archive_bulk_units(self):
        """A pinned release resolves to the nucleotide and amino acid bulk files."""
        source = ImgtSource(client=self.ListingClient(
            ['2026-08-03_GENEDB_202631-7']))
        source.pinRelease('202631-7')
        units = source.searchUnits(Query(collection='human',
                                         filters={'locus': '*', 'segment': '*'}))
        self.assertEqual({u.metadata['group'] for u in units}, {'nt', 'aa'})
        self.assertTrue(all(u.metadata['archive'] == 'genedb' for u in units))
        self.assertTrue(all(u.metadata['release'] == '202631-7' for u in units))

    def test_archive_units_can_be_narrowed_to_chains_and_groups(self):
        """A caller wanting part of the reference fetches one bulk file, filtered."""
        source = ImgtSource(client=self.ListingClient(
            ['2026-08-03_GENEDB_202631-7']))
        source.pinRelease('202631-7')
        units = source.archiveUnits('human', chains=('TRAV', 'IGKC'),
                                    groups=('nt',))
        self.assertEqual(len(units), 1)
        self.assertEqual(units[0].metadata['group'], 'nt')
        self.assertEqual(units[0].metadata['chains'], ['IGKC', 'TRAV'])

    def test_a_substituted_release_is_recorded_as_a_substitute(self):
        """The archive lacking the wanted release is recorded, not papered over."""
        import yaml

        source = ImgtSource(client=self.ListingClient(
            ['2026-08-03_GENEDB_202631-7']))
        source.pinRelease('202629-7')          # not archived
        units = source.searchUnits(Query(collection='human',
                                         filters={'locus': '*', 'segment': '*'}))
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            source.writeReferenceMetadata(tmp, units)
            record = yaml.safe_load((tmp / 'IMGT.yaml').read_text())
        entry = record['species']['human']
        self.assertEqual(entry['release'], '202631-7')
        self.assertEqual(entry['requested'], '202629-7')
        self.assertFalse(entry['exact'])

    def test_build_from_archive_splits_the_bulk_into_chains(self):
        """A bulk file is filtered and split into per-chain reference FASTAs."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            bulk = tmp / 'nt.fasta'
            bulk.write_text(
                '>X1|IGHV1-2*02|Homo sapiens|F\nAC.GT\n'
                '>X2|IGHJ4*02|Homo sapiens|F\nACTG\n'
                '>X3|MICA*01|Homo sapiens|F\nGGGG\n')  # dropped: not IG/TR
            unit = DataUnit(
                unit_id='genedb/human_nt.fasta', collection='human', url='x',
                metadata={'species': 'human', 'archive': 'genedb', 'group': 'nt'})

            source = ImgtSource(client=None)
            source.buildReference([(unit, bulk)], tmp / 'reference_base')

            base = tmp / 'reference_base' / 'human'
            self.assertTrue((base / 'vdj' / 'imgt_human_IGHV.fasta').exists())
            self.assertTrue((base / 'vdj' / 'imgt_human_IGHJ.fasta').exists())
            self.assertFalse(list((base / 'constant').glob('*MIC*'))
                             if (base / 'constant').exists() else [])


class TestReleaseCapture(unittest.TestCase):
    """
    Tests for when the recorded GENE-DB release is read
    """

    class RollingClient:
        """A client whose release tag changes, as IMGT's does when it rebuilds."""

        def __init__(self, release):
            self.release = release

        def get(self, url):
            holder = type('R', (), {'text': self.release})

            return holder()

    def test_release_is_the_one_current_when_the_download_started(self):
        """
        A rebuild part way through must not be recorded as the release fetched.

        GENElect serves only the current build, so the tag read after a download
        may name a release whose sequences were never the ones written. Reading
        it when the units are resolved ties the record to what was fetched.
        """
        import yaml

        client = self.RollingClient('202631-7')
        source = ImgtSource(client=client)
        units = source.searchUnits(Query(collection='human',
                                         filters={'locus': 'IGH', 'segment': 'V'}))
        client.release = '202632-7'      # IMGT rebuilds mid-download

        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            source.writeReferenceMetadata(tmp, units)
            record = yaml.safe_load((tmp / 'IMGT.yaml').read_text())

        entry = record['species']['human']
        self.assertEqual(entry['release'], '202631-7')
        self.assertNotIn('requested', entry)    # nothing was pinned

    def test_no_units_writes_no_sidecar(self):
        """A source that fetched nothing has no release to record."""
        source = ImgtSource(client=self.RollingClient('202631-7'))
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(source.writeReferenceMetadata(Path(tmp), []), [])
            self.assertFalse((Path(tmp) / 'IMGT.yaml').exists())


if __name__ == '__main__':
    unittest.main()
