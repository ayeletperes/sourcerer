"""
Unit tests for the OAS source module
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import collections
import csv
import gzip
import os
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest import mock

# Sourcerer imports
from sourcerer.Cli import getArgParser
from sourcerer.Exceptions import OasParseError
from sourcerer.Http import HttpClient
from sourcerer.Sources import Oas
from sourcerer.Sources.Base import DataUnit
from tests.FakeHttp import FakeResponse, FakeSession

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

    def test_a_value_with_an_escaped_comma_is_unescaped(self):
        """
        Detail pages escape a comma the same way the search form does.

        Left escaped, a value like BType's 'Plasmablasts\\, Memory B cells and
        activated T cells' can never match --btype's validated (unescaped)
        filter: Catalog.filterCatalog does an exact-value comparison, so the
        search would silently return zero hits -- the failure mode the
        snapshot's validated filters exist to prevent. Built inline, as a
        synthetic page, rather than a second committed detail-page fixture,
        since the point being pinned is entirely in this one cell.
        """
        page = ('<html><body><table>'
               '<tr><td>BType</td>'
               '<td>Plasmablasts\\, Memory B cells and activated T cells</td></tr>'
               '</table></body></html>')

        found = Oas.parseDetailPage(page)

        self.assertEqual(found['BType'],
                         'Plasmablasts, Memory B cells and activated T cells')

        rows = self.makeRows()
        source = Oas.OasSource(StubDetailClient(page))
        source.enrichCatalog(rows, force=True)

        self.assertEqual(rows[0]['BType'],
                         'Plasmablasts, Memory B cells and activated T cells')


class TestCatalogFingerprint(unittest.TestCase):
    """
    Tests for condensing the unpaired catalog into a fingerprint
    """

    def test_fingerprint_names_its_own_collection(self):
        """
        The fingerprint records which collection it covers rather than
        leaving a reader (Drift.compareFingerprints, Drift.findAnomalies) to
        assume one. Today that is always 'unpaired' -- OAS publishes no
        machine readable paired index -- but the label still has to come
        from the catalog rows actually fingerprinted, not be hardcoded
        downstream.
        """
        payload = {
            '/vols/naga-datasets/oas/unpaired/Study_A/csv/one.csv.gz':
                {'Species': 'human'},
            '/vols/naga-datasets/oas/unpaired/Study_B/csv/two.csv.gz':
                {'Species': 'mouse_BALB/c'},
        }

        fingerprint = Oas.buildFingerprint(b'raw', {}, payload)

        self.assertEqual(fingerprint['collection'], 'unpaired')

    def test_a_mixed_collection_payload_leaves_the_label_unset(self):
        """
        A single fixed collection is what today's fingerprint is inherently
        about; a payload that somehow named more than one is a fact worth
        surfacing as an unlabeled (source wide) finding rather than guessing
        which collection the drift is really about.
        """
        payload = {
            '/vols/naga-datasets/oas/unpaired/Study_A/csv/one.csv.gz':
                {'Species': 'human'},
            '/vols/naga-datasets/oas/paired/Study_B/csv/two.csv.gz':
                {'Species': 'human'},
        }

        fingerprint = Oas.buildFingerprint(b'raw', {}, payload)

        self.assertIsNone(fingerprint['collection'])


class TestKnownFields(unittest.TestCase):
    """
    The contract between the snapshot and the code

    The monthly refresh updates the snapshot mechanically; this test is what
    turns an unmapped upstream field into a red CI run on that refresh PR.
    When it fails, either map the new field to its AIRR or samplesheet column
    or record an explicit decision to ignore it -- silence is the one option
    the predecessor tool took, and it shipped a stale field list for years.
    """

    def test_every_snapshot_field_is_understood(self):
        from sourcerer.Schema import loadSchema

        schema = loadSchema('oas')
        for name in schema.collection_names:
            for item in schema.getCollection(name).fields:
                self.assertIn(
                    item.name, Oas.KNOWN_FIELDS,
                    "Unmapped OAS field '%s' in collection '%s': add it to "
                    'Oas.KNOWN_FIELDS with its AIRR/samplesheet mapping or an '
                    'explicit note that it is a search filter only'
                    % (item.name, name))


def makeUnit(unit_id, collection='paired', **metadata):
    """Build a data unit carrying the given OAS metadata."""
    return DataUnit(unit_id=unit_id, collection=collection,
                    url='https://example.invalid/%s' % unit_id,
                    metadata=metadata)


class TestSamplesheetRow(unittest.TestCase):
    """
    Tests for mapping OAS metadata to airrflow samplesheet columns
    """

    def test_species_is_lowercased_and_defaults_to_human(self):
        """Species missing entirely still gets a usable value."""
        self.assertEqual(
            Oas.samplesheetRow(makeUnit('A_2020/x.csv.gz'))['species'], 'human')
        self.assertEqual(
            Oas.samplesheetRow(
                makeUnit('A_2020/x.csv.gz', Species='mouse_C57BL/6'))['species'],
            'mouse_c57bl/6')

    def test_subject_no_is_preserved_rather_than_falling_back_to_study(self):
        """
        OAS's own "no" for Subject is kept as is, not replaced by the study name.

        Substituting the study would falsely tell airrflow that every
        otherwise-unidentified unit in that study is the same subject, pooling
        unrelated individuals into one clonal group.
        """
        row = Oas.samplesheetRow(makeUnit('Corinaldesi_2024/csv_paired/a.csv.gz',
                                          study='Corinaldesi_2024', Subject='no'))

        self.assertEqual(row['subject_id'], 'no')

    def test_subject_absent_still_falls_back_to_study(self):
        """
        With no Subject value at all, the study name is still the best guess.

        Unlike an explicit null token such as "no", an absent value carries no
        information of its own to preserve.
        """
        row = Oas.samplesheetRow(makeUnit('Corinaldesi_2024/csv_paired/a.csv.gz',
                                          study='Corinaldesi_2024'))

        self.assertEqual(row['subject_id'], 'Corinaldesi_2024')

    def test_longitudinal_is_carried_through(self):
        """A real Longitudinal value from OAS reaches its own column."""
        row = Oas.samplesheetRow(makeUnit('A_2020/csv/a.csv.gz', Longitudinal='yes'))

        self.assertEqual(row['longitudinal'], 'yes')

    def test_longitudinal_absent_or_no_becomes_na(self):
        """
        Like Age, Longitudinal is a presence flag: "no" and absent both mean
        the design carries no longitudinal information, so both collapse to
        the same 'NA' placeholder airrflow expects.
        """
        self.assertEqual(
            Oas.samplesheetRow(
                makeUnit('A_2020/csv/a.csv.gz', Longitudinal='no'))['longitudinal'],
            'NA')
        self.assertEqual(
            Oas.samplesheetRow(makeUnit('B_2020/csv/b.csv.gz'))['longitudinal'],
            'NA')

    def test_single_cell_is_keyed_on_the_paired_collection(self):
        """
        Only paired data is single cell, driven by the unit's own collection
        rather than hardcoded -- the R implementation assumed TRUE because it
        only ever handled paired.
        """
        self.assertEqual(
            Oas.samplesheetRow(
                makeUnit('A/a.csv.gz', collection='paired'))['single_cell'],
            'TRUE')
        self.assertEqual(
            Oas.samplesheetRow(
                makeUnit('A/a.csv.gz', collection='unpaired'))['single_cell'],
            'FALSE')

    def test_tissue_defaults_to_unknown(self):
        """airrflow requires a tissue value, so a missing BSource gets one."""
        self.assertEqual(
            Oas.samplesheetRow(makeUnit('A/a.csv.gz'))['tissue'], 'unknown')
        self.assertEqual(
            Oas.samplesheetRow(makeUnit('A/a.csv.gz', BSource='PBMC'))['tissue'],
            'PBMC')


class TestCountUnresolvedSubjects(unittest.TestCase):
    """
    Tests for the unresolved-subject count `handleDownload` warns from
    """

    def test_counts_only_null_sentinel_subjects(self):
        """A mix of real, missing, and sentinel Subject values counts right."""
        entries = [
            (makeUnit('A_2020/x.csv.gz', Subject='Donor-1'), Path('a')),
            (makeUnit('B_2020/y.csv.gz', Subject='no'), Path('b')),
            (makeUnit('C_2020/z.csv.gz', Subject='None'), Path('c')),
            (makeUnit('D_2020/w.csv.gz'), Path('d')),
        ]

        self.assertEqual(Oas.countUnresolvedSubjects(entries), 3)

    def test_zero_when_every_unit_has_a_subject(self):
        """Nothing to warn about when OAS recorded a subject for every unit."""
        entries = [(makeUnit('A_2020/x.csv.gz', Subject='Donor-1'), Path('a')),
                  (makeUnit('B_2020/y.csv.gz', Subject='Donor-2'), Path('b'))]

        self.assertEqual(Oas.countUnresolvedSubjects(entries), 0)


#: A samplesheet header narrow enough for verify's own tests: it only reads
#: sample_id, sample_name and subject_id, so the rest is set dressing.
VERIFY_COLUMNS = ('sample_id', 'filename', 'subject_id', 'species', 'sample_name')

#: One esearch/esummary/efetch round trip resolving SRR1 to BL-110_VDJ, the
#: same fixture shape test_Ncbi.py exercises in isolation; this only checks
#: that handleOasVerify wires it into the evidence TSV and --apply correctly.
NCBI_ROUTES = {
    'esearch': FakeResponse(200, b'<eSearchResult><IdList><Id>1</Id>'
                                 b'</IdList></eSearchResult>'),
    'esummary': FakeResponse(200,
        b'<eSummaryResult><DocSum><Id>1</Id>'
        b'<Item Name="ExpXml" Type="String">'
        b'&lt;Summary&gt;&lt;Title&gt;GSM1: BL-110_VDJ&lt;/Title&gt;&lt;/Summary&gt;'
        b'&lt;Biosample&gt;SAMN1&lt;/Biosample&gt;</Item>'
        b'<Item Name="Runs" Type="String">'
        b'&lt;Run acc="SRR1" total_spots="1"/&gt;</Item>'
        b'</DocSum></eSummaryResult>'),
    'efetch': FakeResponse(200,
        b'<BioSampleSet><BioSample accession="SAMN1">'
        b'<Ids><Id db="BioSample">SAMN1</Id></Ids>'
        b'<Description><Title>BL-110_VDJ</Title></Description>'
        b'</BioSample></BioSampleSet>'),
}


class TestHandleOasVerify(unittest.TestCase):
    """
    Tests for the verify command's evidence report
    """

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        self.samplesheet = self.tmp / 'samplesheet_airrflow_airr.tsv'

    def writeSamplesheet(self, rows):
        """Write a samplesheet with just the columns verify needs."""
        with open(self.samplesheet, 'w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(VERIFY_COLUMNS),
                                    delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)

    def runVerify(self, extra_argv=()):
        """Parse a real commandline and run the verify handler against a fake NCBI."""
        argv = ['oas', 'verify', str(self.samplesheet)] + list(extra_argv)
        args = getArgParser().parse_args(argv)

        fake_client = HttpClient(delay=0, backoff=0,
                                 session=FakeSession(lambda method, url, headers, i: next(
                                     response for substring, response in NCBI_ROUTES.items()
                                     if substring in url)))
        with mock.patch('sourcerer.Sources.Oas.HttpClient', return_value=fake_client):
            return Oas.handleOasVerify(args)

    def readReport(self, path=None):
        """Read back the evidence report as a list of dicts."""
        path = path or self.samplesheet.with_name(
            self.samplesheet.stem + '.ncbi_evidence' + self.samplesheet.suffix)
        with open(path, newline='') as handle:
            return list(csv.DictReader(handle, delimiter='\t'))

    def test_a_real_subject_id_is_looked_up_too_and_compared(self):
        """
        A row that already has a subject_id is still cross-referenced.

        OAS recording a subject is not proof it is correct -- a typo, a
        short code reused across studies, or a pooled run naming several
        donors under one value are all real failure modes -- so verify
        looks the accession up regardless, and ncbi_sample_name always carries
        NCBI's own text rather than a copy of subject_id. Here NCBI's
        BL-110_VDJ has nothing in common with 'Donor-2', so subject_check
        reports 'differs'.
        """
        self.writeSamplesheet([{'sample_id': 'ssr_1', 'filename': 'x.tsv',
                               'subject_id': 'Donor-2', 'species': 'human',
                               'sample_name': 'Study/csv_paired/SRR1_1_Paired_All.csv.gz'}])

        self.assertEqual(self.runVerify(), 0)

        rows = self.readReport()
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['status'], 'ok')
        self.assertEqual(rows[0]['ncbi_sample_name'], 'BL-110_VDJ')
        self.assertEqual(rows[0]['ncbi_subject_suggested'], 'BL-110')
        self.assertEqual(rows[0]['subject_check'], 'differs')

    def test_a_real_subject_id_that_matches_ncbi_agrees(self):
        """subject_check reports 'agrees' when subject_id is NCBI's own text."""
        self.writeSamplesheet([{'sample_id': 'ssr_1', 'filename': 'x.tsv',
                               'subject_id': 'BL-110', 'species': 'human',
                               'sample_name': 'Study/csv_paired/SRR1_1_Paired_All.csv.gz'}])

        self.assertEqual(self.runVerify(), 0)

        rows = self.readReport()
        self.assertEqual(rows[0]['subject_check'], 'agrees')

    def test_a_pooled_subject_id_is_never_looked_up_as_a_single_subject(self):
        """
        OAS's own Subject field can itself name a pool of donors.

        'donor 21; 22; 23 and 24' is exactly the shape OAS's paired catalog
        uses for a 10x hashed/pooled run; subject_check must recognize this
        from subject_id alone, the same way it recognizes NCBI's own pooled
        text, rather than reporting a false 'differs'.
        """
        self.writeSamplesheet([{'sample_id': 'ssr_1', 'filename': 'x.tsv',
                               'subject_id': 'donor 21; 22; 23 and 24',
                               'species': 'human',
                               'sample_name': 'Study/csv_paired/SRR1_1_Paired_All.csv.gz'}])

        self.assertEqual(self.runVerify(), 0)

        rows = self.readReport()
        self.assertEqual(rows[0]['subject_check'], 'pooled')

    def test_no_accession_leaves_subject_check_unresolved_only_when_null(self):
        """
        A row whose sample_name carries no accession, but does have a real
        subject_id, is neither 'unresolved' (that's for a null subject_id)
        nor comparable -- it is 'unverified'.
        """
        self.writeSamplesheet([{'sample_id': 'ssr_1', 'filename': 'x.tsv',
                               'subject_id': 'Donor-2', 'species': 'human',
                               'sample_name': 'not-an-accession'}])

        empty_client = HttpClient(delay=0, backoff=0, session=FakeSession(
            lambda *a: (_ for _ in ()).throw(AssertionError('no network call was expected'))))
        args = getArgParser().parse_args(['oas', 'verify', str(self.samplesheet)])
        with mock.patch('sourcerer.Sources.Oas.HttpClient', return_value=empty_client):
            self.assertEqual(Oas.handleOasVerify(args), 0)

        rows = self.readReport()
        self.assertEqual(rows[0]['status'], 'no_accession')
        self.assertEqual(rows[0]['subject_check'], 'unverified')
        self.assertEqual(rows[0]['ncbi_sample_name'], '')

    def test_resolved_row_gets_both_the_raw_name_and_a_suggested_subject(self):
        """
        A resolved row's report names its BioSample, a check link, NCBI's raw
        sample name and the subject suggested from it -- no flag needed to
        choose between them.
        """
        self.writeSamplesheet([{'sample_id': 'ssr_1', 'filename': 'x.tsv',
                               'subject_id': 'no', 'species': 'human',
                               'sample_name': 'Study/csv_paired/SRR1_1_Paired_All.csv.gz'}])

        self.assertEqual(self.runVerify(), 0)

        rows = self.readReport()
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['status'], 'ok')
        self.assertEqual(rows[0]['biosample_accession'], 'SAMN1')
        self.assertEqual(rows[0]['ncbi_sample_name'], 'BL-110_VDJ')
        self.assertEqual(rows[0]['ncbi_subject_suggested'], 'BL-110')
        self.assertEqual(rows[0]['subject_check'], 'unresolved')
        self.assertIn('SAMN1', rows[0]['biosample_url'])

    def test_report_carries_every_input_column(self):
        """
        The report is a superset of the input, not a separate NCBI-only
        file: airrflow-required columns absent from NCBI_EVIDENCE_COLUMNS
        (filename, species, ...) must survive untouched, in their original
        position, so the report can be used as airrflow --input directly.
        """
        self.writeSamplesheet([{'sample_id': 'ssr_1', 'filename': 'fasta/x.fasta',
                               'subject_id': 'no', 'species': 'human',
                               'sample_name': 'Study/csv_paired/SRR1_1_Paired_All.csv.gz'}])

        self.assertEqual(self.runVerify(), 0)

        with open(self.samplesheet.with_name(
                self.samplesheet.stem + '.ncbi_evidence' + self.samplesheet.suffix),
                newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            fields = reader.fieldnames
            row = next(reader)

        self.assertEqual(fields, list(VERIFY_COLUMNS) + list(Oas.NCBI_EVIDENCE_COLUMNS))
        self.assertEqual(row['filename'], 'fasta/x.fasta')
        self.assertEqual(row['species'], 'human')

    def test_default_report_path_moves_the_extension_rather_than_appending_it(self):
        """
        The default path is <stem>.ncbi_evidence<ext>, not
        <stem><ext>.ncbi_evidence.tsv -- one extension, not two.
        """
        self.writeSamplesheet([{'sample_id': 'ssr_1', 'filename': 'x.tsv',
                               'subject_id': 'Donor-2', 'species': 'human',
                               'sample_name': 'Study/csv_paired/SRR1_1_Paired_All.csv.gz'}])

        self.assertEqual(self.runVerify(), 0)

        self.assertTrue((self.tmp / 'samplesheet_airrflow_airr.ncbi_evidence.tsv').exists())
        self.assertFalse(Path(str(self.samplesheet) + '.ncbi_evidence.tsv').exists())

    def test_out_overrides_the_default_report_path(self):
        """--out sends the report somewhere other than the default sidecar path."""
        self.writeSamplesheet([{'sample_id': 'ssr_1', 'filename': 'x.tsv',
                               'subject_id': 'no', 'species': 'human',
                               'sample_name': 'Study/csv_paired/SRR1_1_Paired_All.csv.gz'}])
        out = self.tmp / 'report.tsv'

        self.assertEqual(self.runVerify(['--out', str(out)]), 0)

        self.assertTrue(out.exists())
        self.assertFalse(self.samplesheet.with_name(
            self.samplesheet.stem + '.ncbi_evidence' + self.samplesheet.suffix).exists())
        self.assertEqual(self.readReport(out)[0]['ncbi_sample_name'], 'BL-110_VDJ')

    def test_missing_required_column_is_reported_by_name(self):
        """A samplesheet missing a column verify needs names it in the error."""
        with open(self.samplesheet, 'w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=['sample_id', 'filename'],
                                    delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerow({'sample_id': 'ssr_1', 'filename': 'x.tsv'})
        args = getArgParser().parse_args(['oas', 'verify', str(self.samplesheet)])

        with self.assertRaises(Exception) as raised:
            Oas.handleOasVerify(args)
        self.assertIn('subject_id', str(raised.exception))
        self.assertIn('sample_name', str(raised.exception))

    def test_same_subject_id_in_two_studies_is_warned_about(self):
        """
        A short subject_id reused across studies is a real collision.

        airrflow keys a subject on subject_id alone, so two studies sharing
        'Donor-2' would otherwise merge silently; this is worth a log
        warning independent of subject_check, which only compares each row
        against NCBI and cannot see across rows.
        """
        columns = list(VERIFY_COLUMNS) + ['study']
        with open(self.samplesheet, 'w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t',
                                    lineterminator='\n')
            writer.writeheader()
            writer.writerow({'sample_id': 'ssr_1', 'filename': 'a.tsv',
                             'subject_id': 'Donor-2', 'species': 'human',
                             'sample_name': 'StudyA/x.csv.gz', 'study': 'StudyA'})
            writer.writerow({'sample_id': 'ssr_2', 'filename': 'b.tsv',
                             'subject_id': 'Donor-2', 'species': 'human',
                             'sample_name': 'StudyB/y.csv.gz', 'study': 'StudyB'})

        with self.assertLogs('sourcerer', level='WARNING') as logs:
            self.assertEqual(self.runVerify(), 0)

        self.assertTrue(any('Donor-2' in message for message in logs.output))
        self.assertTrue(any('StudyA' in message and 'StudyB' in message
                            for message in logs.output))


if __name__ == '__main__':
    unittest.main()
