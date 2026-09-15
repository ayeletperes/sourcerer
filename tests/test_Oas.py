"""
Unit tests for the OAS source module
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import collections
import gzip
import os
import unittest

# Sourcerer imports
from sourcerer.Exceptions import OasParseError
from sourcerer.Sources import Oas

test_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(test_path, 'data')


def readFixture(name):
    """Read a fixture, transparently decompressing a gzipped one."""
    path = os.path.join(data_path, name)
    if name.endswith('.gz'):
        with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as handle:
            return handle.read()

    with open(path, encoding='utf-8', errors='replace') as handle:
        return handle.read()


class TestFormSchema(unittest.TestCase):
    """
    Tests for harvesting the searchable vocabulary from a search form
    """

    def setUp(self):
        self.paired = Oas.parseFormSchema(readFixture('oas_paired_form.html'),
                                          'paired')
        self.unpaired = Oas.parseFormSchema(readFixture('oas_unpaired_form.html'),
                                            'unpaired')

    def names(self, fields):
        return [x['name'] for x in fields]

    def test_paired_fields(self):
        self.assertEqual(self.names(self.paired),
                         ['Species', 'Age', 'BSource', 'BType', 'Vaccine',
                          'Disease', 'Subject', 'Longitudinal'])

    def test_unpaired_fields(self):
        self.assertEqual(self.names(self.unpaired),
                         ['Species', 'BSource', 'BType', 'Longitudinal', 'Age',
                          'Disease', 'Subject', 'Vaccine', 'Chain', 'Isotype',
                          'Primer'])

    def test_collections_have_different_fields(self):
        """
        The two collections are not interchangeable.

        The predecessor tool hardcoded paired as having Isotype and lacking
        BSource and BType, which is the exact inverse of what the form offers.
        Paired has no Isotype at all and unpaired has a Primer field that paired
        does not.
        """
        paired, unpaired = set(self.names(self.paired)), set(self.names(self.unpaired))
        self.assertNotIn('Isotype', paired)
        self.assertIn('Isotype', unpaired)
        self.assertIn('Primer', unpaired)
        self.assertNotIn('Primer', paired)
        self.assertTrue({'BSource', 'BType', 'Subject'} <= paired)

    def test_wildcard_is_separated_from_the_vocabulary(self):
        species = next(x for x in self.paired if x['name'] == 'Species')
        self.assertEqual(species['wildcard'], '*')
        self.assertNotIn('*', species['values'])
        self.assertIn('human', species['values'])

    def test_presence_only_fields_are_flagged(self):
        """
        Paired Age, Subject and Longitudinal filter on presence, not on value.

        Recording defined and undefined as a vocabulary would make the builder
        offer them as real choices and would make any comparison against the
        unpaired forms look like a mass deletion of values.
        """
        for name in ('Age', 'Subject', 'Longitudinal'):
            found = next(x for x in self.paired if x['name'] == name)
            self.assertTrue(found['pseudo_values'], name)
            self.assertEqual(found['values'], [])

        species = next(x for x in self.paired if x['name'] == 'Species')
        self.assertFalse(species['pseudo_values'])

    def test_escaped_commas_are_restored(self):
        """Option labels escape embedded commas, which must be undone."""
        btype = next(x for x in self.paired if x['name'] == 'BType')
        commas = [x for x in btype['values'] if ',' in x]
        self.assertTrue(commas)
        self.assertFalse(any('\\,' in x for x in btype['values']))

    def test_missing_form_raises(self):
        with self.assertRaises(OasParseError):
            Oas.parseFormSchema('<html><body>nothing here</body></html>', 'paired')


class TestSearchReply(unittest.TestCase):
    """
    Tests for reading a paired search reply
    """

    @classmethod
    def setUpClass(cls):
        cls.html = readFixture('oas_paired_search_all.html.gz')

    def test_totals(self):
        totals = Oas.parseSearchTotals(self.html)
        self.assertEqual(totals['studies'], 21)
        self.assertGreater(totals['sequences'], 3000000)

    def test_download_urls(self):
        urls = Oas.parseDownloadUrls(self.html)
        self.assertEqual(len(urls), 610)
        self.assertTrue(all(x.startswith(Oas.DOWNLOAD_BASE) for x in urls))

    def test_both_directory_layouts_are_present(self):
        """
        Paired data is not one uniform path shape.

        The majority of units sit under csv_paired/ rather than csv/, so anything
        that assumed a single layout would mishandle most of the collection.
        """
        urls = Oas.parseDownloadUrls(self.html)
        segments = collections.Counter(x.split('/')[-2] for x in urls)
        self.assertEqual(set(segments), {'csv', 'csv_paired'})
        self.assertGreater(segments['csv_paired'], segments['csv'])

    def test_missing_script_raises(self):
        with self.assertRaises(OasParseError):
            Oas.parseDownloadUrls('<html>no script</html>')

    def test_table_rows_use_canonical_field_names(self):
        """
        Result columns are renamed to the spelling the rest of the tool uses.

        OAS calls the same concept Organism here and Species on the form and in
        the data unit metadata.
        """
        rows = Oas.parseSearchTable(self.html)
        self.assertGreater(len(rows), 100)

        first = rows[0]
        self.assertIn('unit_id', first)
        self.assertIn('Species', first)
        self.assertNotIn('Organism', first)
        self.assertIn('Subject', first)
        self.assertNotIn('Individual', first)

    def test_table_rows_match_download_urls(self):
        """Every table row identifies a unit that also appears in the script."""
        rows = Oas.parseSearchTable(self.html)
        urls = Oas.parseDownloadUrls(self.html)
        from_urls = {Oas.unitIdFromUrl(x)[1] for x in urls}

        self.assertTrue({x['unit_id'] for x in rows} <= from_urls)


class TestIdentifiers(unittest.TestCase):
    """
    Tests for the opaque identifier and URL rules
    """

    def test_unit_id_round_trip(self):
        url = (Oas.DOWNLOAD_BASE +
               'paired/Alsoiussi_2020/csv/SRR11528761_paired.csv.gz')
        collection, unit_id = Oas.unitIdFromUrl(url)

        self.assertEqual(collection, 'paired')
        self.assertEqual(unit_id, 'Alsoiussi_2020/csv/SRR11528761_paired.csv.gz')
        self.assertEqual(Oas.urlFromUnitId(collection, unit_id), url)

    def test_unit_id_handles_the_other_layout(self):
        """
        A unit whose filename carries no run accession round trips unchanged.

        This is why identifiers are opaque: there is nothing to parse out of
        '1_S1__1_Paired_All.csv.gz'.
        """
        url = (Oas.DOWNLOAD_BASE +
               'paired/Phad_2022/csv_paired/1_S1__1_Paired_All.csv.gz')
        collection, unit_id = Oas.unitIdFromUrl(url)

        self.assertEqual(unit_id, 'Phad_2022/csv_paired/1_S1__1_Paired_All.csv.gz')
        self.assertEqual(Oas.urlFromUnitId(collection, unit_id), url)

    def test_unknown_collection_raises(self):
        with self.assertRaises(OasParseError):
            Oas.unitIdFromUrl('https://example.org/webapps/ngsdb/other/x.csv.gz')

    def test_catalog_key_maps_to_a_download_url(self):
        key = ('/vols/naga-datasets/oas/unpaired/Banerjee_2017/csv/'
               'SRR5060321_Heavy_Bulk.csv.gz')
        self.assertEqual(
            Oas.urlFromCatalogKey(key),
            Oas.DOWNLOAD_BASE + 'unpaired/Banerjee_2017/csv/'
            'SRR5060321_Heavy_Bulk.csv.gz')

    def test_unexpected_catalog_key_raises(self):
        with self.assertRaises(OasParseError):
            Oas.urlFromCatalogKey('/some/other/mount/unpaired/x.csv.gz')

    def test_unit_stem_is_unique_per_unit(self):
        """Two units sharing a filename must not share an identifier prefix."""
        one = Oas.unitStem('StudyA/csv/SRR1_paired.csv.gz')
        two = Oas.unitStem('StudyB/csv/SRR1_paired.csv.gz')

        self.assertNotEqual(one, two)


class TestHelpers(unittest.TestCase):
    """
    Tests for the small field level mappings
    """

    def test_locus_comes_from_v_call(self):
        """
        The file's own locus column cannot distinguish kappa from lambda.

        OAS writes single letters, and 'L' covers both IGL and, in its own
        encoding, anything light. v_call is unambiguous.
        """
        self.assertEqual(Oas.deriveLocus('IGKV1-39*01', 'K'), 'IGK')
        self.assertEqual(Oas.deriveLocus('IGLV1-47*01', 'L'), 'IGL')
        self.assertEqual(Oas.deriveLocus('IGHV3-73*02', 'H'), 'IGH')
        self.assertEqual(Oas.deriveLocus('TRBV20-1*01', ''), 'TRB')

    def test_locus_falls_back_to_the_letter(self):
        self.assertEqual(Oas.deriveLocus('', 'H'), 'IGH')
        self.assertEqual(Oas.deriveLocus('', ''), '')

    def test_sentinel_isotypes_do_not_become_calls(self):
        """
        Bulk and All mean 'not isotype resolved', not a constant region call.

        Copying them into c_call invents a measurement the experiment never made.
        """
        self.assertEqual(Oas.isotypeToCall('Bulk'), '')
        self.assertEqual(Oas.isotypeToCall('All'), '')
        self.assertEqual(Oas.isotypeToCall(''), '')
        self.assertEqual(Oas.isotypeToCall('IGHG'), 'IGHG')

    def test_boolean_spellings(self):
        self.assertEqual(Oas.toAirrBool('T'), 'T')
        self.assertEqual(Oas.toAirrBool('true'), 'T')
        self.assertEqual(Oas.toAirrBool('F'), 'F')
        self.assertEqual(Oas.toAirrBool('0'), 'F')
        self.assertEqual(Oas.toAirrBool(''), '')

    def test_null_tokens(self):
        self.assertTrue(Oas.isNull('no'))
        self.assertTrue(Oas.isNull('None'))
        self.assertFalse(Oas.isNull('PBMC'))

    def test_cell_barcode_drops_the_contig(self):
        self.assertEqual(Oas.cellBarcode('AAACCTGAGTCAATAG-1_contig_2'),
                         'AAACCTGAGTCAATAG-1')
        self.assertEqual(Oas.cellBarcode('AAACCTGAGTCAATAG-1_contig_11'),
                         'AAACCTGAGTCAATAG-1')
        self.assertEqual(Oas.cellBarcode(''), '')

    def test_cell_barcode_leaves_an_unrecognized_shape_alone(self):
        """
        An identifier that is not barcode_contig is returned whole.

        Truncating it on a guess would invent a cell grouping that the file does
        not support.
        """
        self.assertEqual(Oas.cellBarcode('read_00417'), 'read_00417')


class StubDetailClient:
    """
    A client that answers detail page requests from a canned body.

    Arguments:
      body (str): the HTML to return.
      fail (bool): raise instead of answering, to exercise the failure path.
    """

    def __init__(self, body='', fail=False):
        self.body = body
        self.fail = fail
        self.urls = []

    def get(self, url, **kwargs):
        """Record the request and return a response-like object."""
        self.urls.append(url)
        if self.fail:
            raise OSError('detail page unavailable')

        return collections.namedtuple('Response', 'text')(self.body)


class TestCatalogEnrichment(unittest.TestCase):
    """
    Tests for filling in the fields only a unit's detail page carries
    """

    def setUp(self):
        self.detail = readFixture('oas_dataunit_paired_detail.html')

    def makeRows(self):
        """Two units, one already read and one never attempted."""
        return [{'unit_id': 'Study_A/csv/one_paired.csv.gz', 'collection': 'paired',
                 'BSource': 'PBMC', 'BType': 'Memory-B-Cells',
                 'detail_status': 'ok'},
                {'unit_id': 'Study_B/csv_paired/two.csv.gz', 'collection': 'paired',
                 'BSource': '', 'BType': '', 'detail_status': ''}]

    def test_auto_skips_units_already_read(self):
        """
        The default pass costs one request per unit that still needs one.

        A monthly refresh that re-read all 610 detail pages would be both slow
        and impolite to a host that gives the data away.
        """
        rows = self.makeRows()
        client = StubDetailClient(self.detail)
        source = Oas.OasSource(client)

        source.enrichCatalog(rows)

        self.assertEqual(len(client.urls), 1)
        self.assertIn('two.csv.gz', client.urls[0])

    def test_force_rereads_every_unit(self):
        """
        --refresh-details all is the only recovery from a detail layout change.

        Once a unit is marked ok it is never selected again, so without a way to
        override that, a page whose layout changed would keep its stale values
        forever.
        """
        rows = self.makeRows()
        client = StubDetailClient(self.detail)
        source = Oas.OasSource(client)

        source.enrichCatalog(rows, force=True)

        self.assertEqual(len(client.urls), 2)

    def test_a_failed_page_keeps_earlier_values(self):
        """
        Enrichment degrades to stale-but-correct, never to silently emptied.

        A transient failure must not blank BSource and BType, and must leave the
        unit eligible for another attempt rather than writing it off.
        """
        rows = self.makeRows()
        source = Oas.OasSource(StubDetailClient(fail=True))

        enriched = source.enrichCatalog(rows, force=True)

        self.assertEqual(enriched, 0)
        self.assertEqual(rows[0]['BSource'], 'PBMC')
        self.assertEqual(rows[0]['detail_status'], 'failed')
        self.assertTrue(all(Oas.needsDetail(x) for x in rows))

    def test_an_unparseable_page_is_an_error_not_an_empty_result(self):
        """
        A scraper that returns {} on failure produces confidently wrong output.

        Without this the caller would mark the unit enriched, leaving BSource and
        BType blank forever with nothing recording that parsing broke.
        """
        with self.assertRaises(OasParseError):
            Oas.parseDetailPage('<html><body><p>Service unavailable</p></body></html>')

        rows = self.makeRows()
        source = Oas.OasSource(StubDetailClient('<html><body></body></html>'))

        source.enrichCatalog(rows, force=True)

        self.assertTrue(all(x['detail_status'] == 'failed' for x in rows))
        self.assertEqual(rows[0]['BSource'], 'PBMC')

    def test_limit_caps_the_number_of_fetches(self):
        rows = self.makeRows()
        client = StubDetailClient(self.detail)
        source = Oas.OasSource(client)

        source.enrichCatalog(rows, limit=1, force=True)

        self.assertEqual(len(client.urls), 1)


if __name__ == '__main__':
    unittest.main()
