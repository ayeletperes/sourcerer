"""
Unit tests for the OGRDB source
"""

# Info
__author__ = 'Ayelet Peres'

# Imports
import os
import tempfile
import unittest
from pathlib import Path

# Sourcerer imports
from sourcerer.Reference import KIND_CONSTANT, KIND_VDJ
from sourcerer.Sources.Base import DataUnit, Query
from sourcerer.Sources.Ogrdb import (
    OgrdbSource,
    bucketChain,
    normalizeVersion,
    safeSetName,
)

test_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(test_path, 'data')


def readFixture(name):
    """Read a captured fixture from tests/data."""
    with open(os.path.join(data_path, name)) as handle:
        return handle.read()


class Canned:
    """A response exposing the .json() the OGRDB client reads."""

    def __init__(self, payload):
        self._payload = payload

    def json(self):
        return self._payload


class StubClient:
    """An OGRDB API client backed by canned payloads, no network."""

    SPECIES = {'species': [{'label': 'Homo sapiens', 'id': '9606'}]}
    SETS = {'germline_species': [
        {'germline_set_name': 'IGH_VDJ', 'locus': 'IGH',
         'germline_set_id': '9606.IGH_VDJ'},
        {'germline_set_name': 'IGHC', 'locus': 'IGH',
         'germline_set_id': '9606.IGHC'},
        {'germline_set_name': 'IGKappa_VJ', 'locus': 'IGK',
         'germline_set_id': '9606.IGK'},
        {'germline_set_name': 'IGLambda_VJ', 'locus': 'IGL',
         'germline_set_id': '9606.IGL'}]}
    LATEST = {'GermlineSet': [
        {'release_version': 2.0, 'release_date': '2024-06-01T00:00:00'}]}

    def get(self, url):
        if '/latest' in url:
            return Canned(self.LATEST)
        if '/germline/sets/' in url:
            return Canned(self.SETS)
        if '/germline/species' in url:
            return Canned(self.SPECIES)
        raise AssertionError('unexpected url %s' % url)


class TestHelpers(unittest.TestCase):
    """
    Tests for the small pure helpers
    """

    def test_normalize_version_strips_trailing_zero(self):
        """An integer release reported as 3.0 becomes 3 for the URL."""
        self.assertEqual(normalizeVersion(3.0), '3')
        self.assertEqual(normalizeVersion('2.1'), '2.1')

    def test_safe_set_name_replaces_separators(self):
        """A set name with spaces and slashes becomes one safe token."""
        self.assertEqual(safeSetName('C57BL/6J IGKV'), 'C57BL_6J_IGKV')


class TestBucketChain(unittest.TestCase):
    """
    Tests for classifying an allele into a reference chain
    """

    def test_v_and_j_use_four_character_chain(self):
        """V and J file under their own four-character chain."""
        self.assertEqual(bucketChain('IGKV1-12*01', 'A' * 300),
                         ('IGKV', KIND_VDJ))
        self.assertEqual(bucketChain('IGKJ1*01', 'A' * 38), ('IGKJ', KIND_VDJ))

    def test_short_ighd_is_the_diversity_segment(self):
        """A short IGHD is the D segment and files under vdj."""
        self.assertEqual(bucketChain('IGHD1-1*01', 'A' * 17), ('IGHD', KIND_VDJ))

    def test_long_ighd_is_the_delta_constant(self):
        """A long IGHD is the delta constant and files under the locus constant."""
        self.assertEqual(bucketChain('IGHD*01', 'A' * 400), ('IGHC', KIND_CONSTANT))

    def test_isotype_constant_files_under_locus_constant(self):
        """A heavy isotype such as IGHM files under IGHC."""
        self.assertEqual(bucketChain('IGHM*01', 'A' * 400), ('IGHC', KIND_CONSTANT))


class TestSearchUnits(unittest.TestCase):
    """
    Tests for resolving a query to downloads
    """

    def test_emits_two_forms_per_set_with_human_ex_endpoint(self):
        """Each set is fetched ungapped and gapped, human via the _ex endpoint."""
        source = OgrdbSource(client=StubClient())
        units = source.searchUnits(
            Query(collection='human', filters={'locus': 'IGK'}))

        self.assertEqual(len(units), 2)
        formats = {u.metadata['format'] for u in units}
        self.assertEqual(formats, {'ungapped', 'gapped'})
        for unit in units:
            self.assertTrue(unit.url.endswith('_ex'))
            self.assertIn('/9606.IGK/2/', unit.url)
            self.assertEqual(unit.metadata['version'], '2')

    def test_locus_filter_narrows_to_one_locus(self):
        """A locus filter fetches only that locus's sets."""
        source = OgrdbSource(client=StubClient())
        units = source.searchUnits(
            Query(collection='human', filters={'locus': 'IGK'}))
        self.assertEqual({u.metadata['locus'] for u in units}, {'IGK'})

    def test_wildcard_covers_every_configured_locus(self):
        """With no locus filter, every immunoglobulin locus is fetched."""
        source = OgrdbSource(client=StubClient())
        units = source.searchUnits(
            Query(collection='human', filters={'locus': '*'}))
        self.assertEqual({u.metadata['locus'] for u in units},
                         {'IGH', 'IGK', 'IGL'})


class TestBuildReference(unittest.TestCase):
    """
    Tests for splitting downloaded sets into per-chain FASTAs
    """

    def _entries(self, tmp):
        entries = []
        for fmt, fixture in (('ungapped', 'ogrdb_igk_ungapped.fasta'),
                             ('gapped', 'ogrdb_igk_gapped.fasta')):
            path = tmp / ('%s.fasta' % fmt)
            path.write_text(readFixture(fixture))
            unit = DataUnit(
                unit_id='IGKappa_VJ.%s.fasta' % fmt, collection='human', url='x',
                metadata={'species': 'human', 'locus': 'IGK',
                          'set_name': 'IGKappa_VJ', 'format': fmt,
                          'chains': ['IGKV', 'IGKJ']})
            entries.append((unit, path))
        return entries

    def test_v_from_gapped_and_j_from_ungapped(self):
        """V keeps its gaps from the gapped form; J comes from the ungapped."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            source = OgrdbSource(client=None)
            source.buildReference(self._entries(tmp), tmp / 'reference_base')

            vdj = tmp / 'reference_base' / 'human' / 'vdj'
            v = (vdj / 'airrc_human_IGKV.fasta').read_text()
            j = (vdj / 'airrc_human_IGKJ.fasta').read_text()

            self.assertIn('IGKV1-12*01', v)
            self.assertIn('.', v)                 # gapped V keeps its IMGT gaps
            self.assertIn('IGKJ1*01', j)
            self.assertNotIn('.', j)              # J is taken ungapped


class TestAlias(unittest.TestCase):
    """
    Tests for the airrc alias
    """

    def test_airrc_resolves_to_ogrdb(self):
        """'airrc' is an alias that resolves to the ogrdb source."""
        from sourcerer.Sources import canonicalName, getSource

        self.assertEqual(canonicalName('airrc'), 'ogrdb')
        self.assertIsInstance(getSource('airrc', client=None), OgrdbSource)

    def test_ogrdb_declares_the_alias(self):
        """The source lists airrc among its aliases."""
        self.assertIn('airrc', OgrdbSource.aliases)


if __name__ == '__main__':
    unittest.main()
