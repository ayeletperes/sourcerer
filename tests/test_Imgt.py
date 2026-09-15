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


if __name__ == '__main__':
    unittest.main()
