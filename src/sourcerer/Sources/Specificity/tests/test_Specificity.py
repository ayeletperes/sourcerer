"""
Unit tests for specificity sources (IEDB) and the Range-header paginator
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import json
import tempfile
import unittest
import zipfile
from io import BytesIO
from pathlib import Path

import pandas

# Sourcerer imports
from sourcerer.Sources.Specificity import REGISTRY, getSpecificitySource
from sourcerer.Sources.Specificity.Iedb import (_CHAIN_DIRECT, _CHAIN_EXTRA,
                                                _CHAIN_PREFERRED,
                                                _SHARED_COLUMNS, API_BASE,
                                                BULK_URL, IedbSource)
from sourcerer.Sources.Specificity.Paginate import pageByRange
from sourcerer.Exceptions import HttpError, IedbParseError
from sourcerer.Http import HttpClient
from sourcerer.Sources.Base import DataUnit, Query
from tests.FakeHttp import FakeResponse, FakeSession, rangeHandler


def blankReceptorRow():
    """
    A receptor CSV row, as a {(category, field): value} dict, every value ''.

    Built from the module's own column mapping tables rather than a hand
    written list, so a fixture row always has exactly the columns the mapping
    code reads, however that mapping grows.
    """
    row = {key: '' for key in _SHARED_COLUMNS}
    for chain_label in ('Chain 1', 'Chain 2'):
        for column in ('Type', 'Nucleotide Sequence', 'Protein Sequence'):
            row[(chain_label, column)] = ''
        for curated, calculated in _CHAIN_PREFERRED.values():
            row[(chain_label, curated)] = ''
            row[(chain_label, calculated)] = ''
        for column in _CHAIN_DIRECT.values():
            row[(chain_label, column)] = ''
        for column in _CHAIN_EXTRA.values():
            row[(chain_label, column)] = ''

    return row


def makeReceptorFrame(rows):
    """
    Build a receptor chunk: a DataFrame with a (category, field) MultiIndex.

    Arguments:
      rows (list): one dict per row, as returned by blankReceptorRow(),
        optionally overridden with the fields a test cares about.

    Returns:
      pandas.DataFrame: columns in the same shape readUnit produces.
    """
    columns = pandas.MultiIndex.from_tuples(list(rows[0].keys()))

    return pandas.DataFrame([[row[c] for c in columns] for row in rows],
                            columns=columns)


def makeClient(handler):
    """Build an HttpClient with no politeness delay over a scripted session."""
    return HttpClient(delay=0, backoff=0, session=FakeSession(handler))


def itemRangeHandler(pages):
    """
    Build a handler serving JSON pages via item Range headers, PostgREST style.

    Arguments:
      pages (list): one list of records per page, in the order they are served.
        A request past the last page gets a 416, the way IEDB's API ends
        pagination.

    Returns:
      callable: a handler for FakeSession.
    """
    total = sum(len(page) for page in pages)

    def handler(method, url, headers, index):
        if index >= len(pages):
            return FakeResponse(416, b'')

        page = pages[index]
        start = sum(len(p) for p in pages[:index])
        end = start + len(page) - 1
        body = json.dumps(page).encode('utf-8')

        return FakeResponse(206, body, {
            'Content-Range': 'items %d-%d/%d' % (start, end, total)})

    return handler


class TestPageByRange(unittest.TestCase):
    """
    Tests for the PostgREST Range-header paginator
    """

    def test_pages_until_content_range_total_is_reached(self):
        client = makeClient(itemRangeHandler([[{'id': 1}, {'id': 2}], [{'id': 3}]]))
        records = [r for batch in pageByRange(client, 'https://x/y', page_size=2)
                  for r in batch]

        self.assertEqual(records, [{'id': 1}, {'id': 2}, {'id': 3}])

    def test_a_short_final_page_also_ends_pagination(self):
        """A page shorter than page_size ends iteration even with no total."""
        def handler(method, url, headers, index):
            if index == 0:
                return FakeResponse(206, json.dumps([{'id': 1}]).encode())
            raise AssertionError('paginator should have stopped after one page')

        client = makeClient(handler)
        records = [r for batch in pageByRange(client, 'https://x/y', page_size=5)
                  for r in batch]

        self.assertEqual(records, [{'id': 1}])

    def test_416_past_the_last_row_ends_pagination(self):
        client = makeClient(itemRangeHandler([[{'id': 1}]]))
        records = [r for batch in pageByRange(client, 'https://x/y', page_size=1)
                  for r in batch]

        self.assertEqual(records, [{'id': 1}])

    def test_unexpected_status_raises(self):
        client = makeClient(lambda method, url, headers, index: FakeResponse(404))

        with self.assertRaises(HttpError):
            list(pageByRange(client, 'https://x/y'))


def makeZip(members):
    """
    Build ZIP bytes containing the given filename to text-content mapping.
    """
    buffer = BytesIO()
    with zipfile.ZipFile(buffer, 'w') as archive:
        for name, content in members.items():
            archive.writestr(name, content)

    return buffer.getvalue()


class TestIedbSource(unittest.TestCase):
    """
    Tests for the IEDB specificity source
    """

    def test_registered(self):
        self.assertIs(REGISTRY['iedb'], IedbSource)
        client = makeClient(lambda *a: FakeResponse(200))
        self.assertIsInstance(getSpecificitySource('iedb', client), IedbSource)

    def test_unknown_source_raises(self):
        client = makeClient(lambda *a: FakeResponse(200))
        with self.assertRaises(KeyError):
            getSpecificitySource('not-a-real-db', client)

    def test_search_units_returns_one_unit_per_table(self):
        source = IedbSource(makeClient(lambda *a: FakeResponse(200)))

        for table in source.collections:
            units = source.searchUnits(Query(collection=table, filters={}))
            self.assertEqual(len(units), 1)
            self.assertEqual(units[0].collection, table)

        bcr = source.searchUnits(Query(collection='bcr', filters={}))[0]
        self.assertEqual(bcr.url, BULK_URL)
        bcell = source.searchUnits(
            Query(collection='bcell', filters={'qualitative_measure': 'Positive'}))[0]
        self.assertEqual(bcell.url, API_BASE + '/bcell_search')
        self.assertEqual(bcell.metadata, {'qualitative_measure': 'Positive'})

    def test_fetch_unit_extracts_the_right_member_from_the_shared_zip(self):
        """bcr and tcr are both members of one shared ZIP, extracted by name."""
        zip_bytes = makeZip({
            'doc/bcr_full_v3.csv': 'a,b\n1,2\n',
            'doc/tcr_full_v3.csv': 'a,b\n3,4\n',
        })
        client = makeClient(rangeHandler(zip_bytes))
        source = IedbSource(client)
        unit = DataUnit(unit_id='bcr_full_v3.csv', collection='bcr', url=BULK_URL)

        with tempfile.TemporaryDirectory() as outdir:
            result = source.fetchUnit(unit, Path(outdir))
            self.assertEqual(result.path.read_text(), 'a,b\n1,2\n')
            # The ZIP itself is cached alongside the extracted table, so a
            # second table sharing it does not re-download it.
            self.assertTrue((Path(outdir) / 'receptor_full_v3.zip').exists())

    def test_fetch_unit_raises_when_the_zip_member_is_missing(self):
        zip_bytes = makeZip({'doc/bcr_full_v3.csv': 'a,b\n1,2\n'})
        client = makeClient(rangeHandler(zip_bytes))
        source = IedbSource(client)
        unit = DataUnit(unit_id='tcr_full_v3.csv', collection='tcr', url=BULK_URL)

        with tempfile.TemporaryDirectory() as outdir:
            with self.assertRaises(IedbParseError):
                source.fetchUnit(unit, Path(outdir))

    def test_fetch_unit_pages_the_api_and_writes_json(self):
        client = makeClient(itemRangeHandler([[{'receptor_group_id': 1,
                                               'bcell_id': 2}]]))
        source = IedbSource(client)
        unit = DataUnit(unit_id='bcr_to_bcell.json', collection='bcr_to_bcell',
                       url=API_BASE + '/bcr_to_bcell')

        with tempfile.TemporaryDirectory() as outdir:
            result = source.fetchUnit(unit, Path(outdir))
            written = json.loads(result.path.read_text())
            self.assertEqual(written, [{'receptor_group_id': 1, 'bcell_id': 2}])

    def test_normalize_chunk_adds_provenance_columns(self):
        source = IedbSource(makeClient(lambda *a: FakeResponse(200)))
        unit = DataUnit(unit_id='bcr_to_bcell.json', collection='bcr_to_bcell',
                       url=API_BASE + '/bcr_to_bcell', metadata={})
        chunk = pandas.DataFrame([{'receptor_group_id': 1, 'bcell_id': 2}])
        report = {'rows_in': 0, 'rows_out': 0}

        frame = source.normalizeChunk({}, chunk, unit, 0, report)

        self.assertEqual(frame['sourcerer_source'].tolist(), ['iedb'])
        self.assertEqual(frame['sourcerer_collection'].tolist(), ['bcr_to_bcell'])
        self.assertEqual(frame['sourcerer_unit_id'].tolist(), ['bcr_to_bcell.json'])
        self.assertEqual(len(frame['sourcerer_row_hash'].iloc[0]), 12)
        self.assertEqual(report['rows_in'], 1)
        self.assertEqual(report['rows_out'], 1)

    def test_normalize_chunk_applies_the_qualitative_measure_filter(self):
        """
        bcell_search has no server side filter for qualitative_measure, so a
        non wildcard value is applied client side during normalization.
        """
        source = IedbSource(makeClient(lambda *a: FakeResponse(200)))
        unit = DataUnit(unit_id='bcell_search.json', collection='bcell',
                       url=API_BASE + '/bcell_search',
                       metadata={'qualitative_measure': 'Positive'})
        chunk = pandas.DataFrame([
            {'qualitative_measure': 'Positive'},
            {'qualitative_measure': 'Negative'},
        ])
        report = {'rows_in': 0, 'rows_out': 0}

        frame = source.normalizeChunk({}, chunk, unit, 0, report)

        self.assertEqual(frame['qualitative_measure'].tolist(), ['Positive'])
        self.assertEqual(report['rows_in'], 2)
        self.assertEqual(report['rows_out'], 1)

    def test_normalize_chunk_keeps_everything_for_the_wildcard(self):
        source = IedbSource(makeClient(lambda *a: FakeResponse(200)))
        unit = DataUnit(unit_id='bcell_search.json', collection='bcell',
                       url=API_BASE + '/bcell_search',
                       metadata={'qualitative_measure': '*'})
        chunk = pandas.DataFrame([
            {'qualitative_measure': 'Positive'},
            {'qualitative_measure': 'Negative'},
        ])

        frame = source.normalizeChunk({}, chunk, unit, 0,
                                      {'rows_in': 0, 'rows_out': 0})

        self.assertEqual(len(frame), 2)

    def test_read_unit_bulk_csv_honors_the_two_header_rows(self):
        """
        The receptor CSVs carry a category row over a field name row; reading
        them with a single header row would take category labels ('Chain 1',
        'Chain 2', ...) as column names and the real field names as a data row.
        """
        source = IedbSource(makeClient(lambda *a: FakeResponse(200)))
        unit = DataUnit(unit_id='bcr_full_v3.csv', collection='bcr', url=BULK_URL)

        row = blankReceptorRow()
        row[('Receptor', 'IEDB Receptor ID')] = '1'
        row[('Chain 1', 'Type')] = 'heavy'
        row[('Chain 1', 'Nucleotide Sequence')] = 'ACGT'
        frame_in = makeReceptorFrame([row])

        with tempfile.TemporaryDirectory() as outdir:
            path = Path(outdir) / 'bcr_full_v3.csv'
            frame_in.to_csv(path, index=False)
            metadata, chunks = source.readUnit(path, unit)
            frame = pandas.concat(list(chunks), ignore_index=True)

        self.assertEqual(metadata, {})
        self.assertEqual(frame[('Chain 1', 'Type')].tolist(), ['heavy'])
        self.assertEqual(frame[('Chain 1', 'Nucleotide Sequence')].tolist(),
                         ['ACGT'])

    def test_read_unit_api_json(self):
        source = IedbSource(makeClient(lambda *a: FakeResponse(200)))
        unit = DataUnit(unit_id='bcr_to_bcell.json', collection='bcr_to_bcell',
                       url=API_BASE + '/bcr_to_bcell')

        with tempfile.TemporaryDirectory() as outdir:
            path = Path(outdir) / 'bcr_to_bcell.json'
            path.write_text(json.dumps([{'receptor_group_id': 1, 'bcell_id': 2}]))
            metadata, chunks = source.readUnit(path, unit)
            frame = pandas.concat(list(chunks), ignore_index=True)

        self.assertEqual(metadata, {})
        self.assertEqual(frame.to_dict('records'),
                         [{'receptor_group_id': 1, 'bcell_id': 2}])


class TestIedbReceptorAirrMapping(unittest.TestCase):
    """
    Tests for mapping a receptor CSV row to per-chain AIRR records
    """

    def test_bcr_and_tcr_are_airr_collections(self):
        self.assertEqual(IedbSource.airr_collections, frozenset(['bcr', 'tcr']))
        self.assertEqual(IedbSource.airr_collections & set(IedbSource.collections),
                         frozenset(['bcr', 'tcr']))

    def _normalize(self, rows, collection='bcr'):
        source = IedbSource(makeClient(lambda *a: FakeResponse(200)))
        unit = DataUnit(unit_id='bcr_full_v3.csv', collection=collection, url=BULK_URL)
        report = {'rows_in': 0, 'rows_out': 0}

        frame = source.normalizeChunk({}, makeReceptorFrame(rows), unit, 0, report)

        return frame, report

    def test_a_two_chain_receptor_becomes_two_rows_sharing_a_cell_id(self):
        row = blankReceptorRow()
        row[('Receptor', 'IEDB Receptor ID')] = '42'
        row[('Chain 1', 'Type')] = 'heavy'
        row[('Chain 1', 'Nucleotide Sequence')] = 'ACGT'
        row[('Chain 2', 'Type')] = 'kappa_light'
        row[('Chain 2', 'Nucleotide Sequence')] = 'TTGG'

        frame, report = self._normalize([row])

        self.assertEqual(frame['sequence_id'].tolist(), ['42_1', '42_2'])
        self.assertEqual(frame['cell_id'].tolist(), ['42', '42'])
        self.assertEqual(frame['locus'].tolist(), ['IGH', 'IGK'])
        self.assertEqual(frame['sequence'].tolist(), ['ACGT', 'TTGG'])
        self.assertEqual(report['rows_in'], 1)
        self.assertEqual(report['rows_out'], 2)

    def test_a_single_chain_receptor_produces_only_the_present_chain(self):
        """A missing second chain must not become a fabricated empty row."""
        row = blankReceptorRow()
        row[('Receptor', 'IEDB Receptor ID')] = '7'
        row[('Chain 1', 'Type')] = 'alpha'
        row[('Chain 1', 'Protein Sequence')] = 'EVQ'

        frame, report = self._normalize([row], collection='tcr')

        self.assertEqual(len(frame), 1)
        self.assertEqual(frame['sequence_id'].tolist(), ['7_1'])
        self.assertEqual(frame['locus'].tolist(), ['TRA'])
        self.assertEqual(report['rows_out'], 1)

    def test_curated_value_is_preferred_over_calculated(self):
        row = blankReceptorRow()
        row[('Receptor', 'IEDB Receptor ID')] = '1'
        row[('Chain 1', 'Type')] = 'heavy'
        row[('Chain 1', 'Curated V Gene')] = 'IGHV1-2*02'
        row[('Chain 1', 'Calculated V Gene')] = 'IGHV1-2*01'

        frame, _ = self._normalize([row])

        self.assertEqual(frame['v_call'].tolist(), ['IGHV1-2*02'])

    def test_calculated_value_is_used_when_curated_is_blank(self):
        row = blankReceptorRow()
        row[('Receptor', 'IEDB Receptor ID')] = '1'
        row[('Chain 1', 'Type')] = 'heavy'
        row[('Chain 1', 'Calculated V Gene')] = 'IGHV1-2*01'

        frame, _ = self._normalize([row])

        self.assertEqual(frame['v_call'].tolist(), ['IGHV1-2*01'])

    def test_an_unresolvable_chain_type_leaves_locus_blank_but_keeps_the_row(self):
        row = blankReceptorRow()
        row[('Receptor', 'IEDB Receptor ID')] = '1'
        row[('Chain 1', 'Type')] = 'IgNAR'
        row[('Chain 1', 'Protein Sequence')] = 'EVQ'

        frame, _ = self._normalize([row])

        self.assertEqual(len(frame), 1)
        self.assertEqual(frame['locus'].tolist(), [''])
        self.assertEqual(frame['iedb_chain_type'].tolist(), ['IgNAR'])

    def test_junction_aa_length_is_derived_and_left_blank_when_absent(self):
        row1 = blankReceptorRow()
        row1[('Receptor', 'IEDB Receptor ID')] = '1'
        row1[('Chain 1', 'Type')] = 'heavy'
        row1[('Chain 1', 'Junction Calculated')] = 'CARDYW'

        row2 = blankReceptorRow()
        row2[('Receptor', 'IEDB Receptor ID')] = '2'
        row2[('Chain 1', 'Type')] = 'heavy'

        frame, _ = self._normalize([row1, row2])

        self.assertEqual(frame['junction_aa_length'].tolist(), ['6', ''])

    def test_shared_receptor_and_epitope_context_is_carried_on_every_chain(self):
        row = blankReceptorRow()
        row[('Receptor', 'IEDB Receptor ID')] = '9'
        row[('Receptor', 'Group IRI')] = 'https://www.iedb.org/receptor/9'
        row[('Epitope', 'Name')] = 'ACE2 loop'
        row[('Chain 1', 'Type')] = 'heavy'
        row[('Chain 2', 'Type')] = 'kappa_light'

        frame, _ = self._normalize([row])

        self.assertEqual(frame['receptor_group_id'].tolist(),
                         ['https://www.iedb.org/receptor/9'] * 2)
        self.assertEqual(frame['epitope_name'].tolist(), ['ACE2 loop'] * 2)

    def test_writes_as_a_header_and_row_valid_rearrangement(self):
        """
        The mapped frame is what `sourcerer specificity iedb download bcr`
        actually writes: Convert.writeAirr fills in the required fields the
        source data does not carry (rev_comp, cigars, ...) as blank columns,
        the same way it does for every other source.
        """
        from sourcerer import Convert

        row = blankReceptorRow()
        row[('Receptor', 'IEDB Receptor ID')] = '1'
        row[('Chain 1', 'Type')] = 'heavy'
        row[('Chain 1', 'Nucleotide Sequence')] = 'ACGT'
        row[('Chain 1', 'Junction Calculated')] = 'CARDYW'

        frame, _ = self._normalize([row])

        with tempfile.TemporaryDirectory() as outdir:
            validation = Convert.writeAirr(iter([frame]), Path(outdir) / 'bcr.tsv')

        self.assertTrue(validation['header_valid'])
        self.assertEqual(validation['rows_invalid'], 0)
