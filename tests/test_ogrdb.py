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


class TestPinningAndMetadata(unittest.TestCase):
    """
    Tests for --from version pinning and the AIRRC.yaml sidecar
    """

    def test_search_uses_latest_by_default(self):
        """With no pin, a set is fetched at its latest version and date."""
        source = OgrdbSource(client=StubClient())
        unit = source.searchUnits(Query(collection='human',
                                        filters={'locus': 'IGK'}))[0]
        self.assertEqual(unit.metadata['version'], '2')
        self.assertEqual(unit.metadata['release_date'], '2024-06-01')

    def test_pin_overrides_the_version_and_date(self):
        """A pinned set is fetched at its recorded version, not the latest."""
        source = OgrdbSource(client=StubClient())
        source.pinSets([{'set': 'IGKappa_VJ', 'version': '1',
                         'release_date': '2020-01-01'}])
        unit = source.searchUnits(Query(collection='human',
                                        filters={'locus': 'IGK'}))[0]
        self.assertEqual(unit.metadata['version'], '1')
        self.assertEqual(unit.metadata['release_date'], '2020-01-01')
        self.assertIn('/9606.IGK/1/', unit.url)

    def test_metadata_has_one_entry_per_set(self):
        """AIRRC.yaml carries one record per set, not one per download form."""
        import yaml

        source = OgrdbSource(client=StubClient())
        units = source.searchUnits(Query(collection='human',
                                         filters={'locus': 'IGH'}))
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            source.writeReferenceMetadata(tmp, units)
            record = yaml.safe_load((tmp / 'AIRRC.yaml').read_text())
            names = {s['set'] for s in record['sets']}
            self.assertEqual(names, {'IGH_VDJ', 'IGHC'})
            self.assertEqual(record['sets'][0]['release_date'], '2024-06-01')

    def test_a_set_missing_from_the_pins_warns_before_taking_latest(self):
        """A part-pinned reference is what --from exists to avoid, so it is said."""
        source = OgrdbSource(client=StubClient())
        source.pinSets([{'set': 'IGH_VDJ', 'version': '1',
                         'release_date': '2020-01-01'}])
        with self.assertLogs('sourcerer.Sources.Ogrdb', level='WARNING') as logged:
            units = source.searchUnits(Query(collection='human',
                                             filters={'locus': 'IGH'}))
        self.assertTrue(any('IGHC' in line for line in logged.output))
        pinned = {u.metadata['set_name']: u.metadata['version'] for u in units}
        self.assertEqual(pinned['IGH_VDJ'], '1')      # pinned
        self.assertNotEqual(pinned['IGHC'], '1')      # fell back to latest

    def test_enable_doi_is_off_until_asked_for(self):
        """DOI resolution scrapes a web UI, so it is opt-in."""
        source = OgrdbSource(client=None)
        self.assertFalse(source._resolve_doi)
        source.enableDoi()
        self.assertTrue(source._resolve_doi)

    def test_pin_release_is_rejected_for_ogrdb(self):
        """OGRDB cannot pin an IMGT release; asking is a clear error."""
        from sourcerer.Exceptions import SourcererError

        with self.assertRaises(SourcererError):
            OgrdbSource(client=None).pinRelease('202631-7')


class TestDoi(unittest.TestCase):
    """
    Tests for the opt-in Zenodo DOI resolver
    """

    class Text:
        def __init__(self, text):
            self.text = text

    class DoiClient:
        PAGE = ('<input name="csrf_token" type="hidden" value="TOK">')
        TABLE = ('<table><tr>'
                 '<td>IGH_VDJ</td><td>10</td><td>2026-05-27</td>'
                 '<td><a href="https://doi.org/10.5281/zenodo.20409587">x</a></td>'
                 '</tr></table>')

        def __init__(self):
            self.posted = None

        def get(self, url):
            return TestDoi.Text(self.PAGE)

        def post(self, url, data=None):
            self.posted = data
            return TestDoi.Text(self.TABLE)

    def test_resolves_doi_and_zenodo_id_from_the_matching_row(self):
        """The row matching set, version and date yields the DOI and record id."""
        source = OgrdbSource(client=self.DoiClient())
        found = source.resolveDoi('Homo sapiens', 'IGH_VDJ', '10', '2026-05-27')
        self.assertEqual(found['doi'], '10.5281/zenodo.20409587')
        self.assertEqual(found['zenodo_record_id'], '20409587')
        self.assertIn('20409587', found['zenodo_url'])

    def test_no_matching_row_returns_empty_not_error(self):
        """A version the page does not list resolves to nothing, not a crash."""
        source = OgrdbSource(client=self.DoiClient())
        self.assertEqual(source.resolveDoi('Homo sapiens', 'IGH_VDJ', '99',
                                           '2026-05-27'), {})


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
