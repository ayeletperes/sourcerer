"""
Unit tests for reading OAS data units and converting them to AIRR
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import gzip
import io
import os
import tempfile
import unittest
import zlib
from pathlib import Path

import pandas

# Sourcerer imports
from sourcerer import Convert
from sourcerer.Exceptions import OasParseError
from sourcerer.Gzip import decompressPrefix
from sourcerer.Sources import Oas

try:
    import changeo.Receptor  # noqa: F401
    import presto.Annotation  # noqa: F401
    HAS_IMMCANTATION_TOOLS = True
except ImportError:
    # Neither is a dependency of sourcerer; they are the downstream consumers of
    # the FASTA headers, so the round trip is verified only where present.
    HAS_IMMCANTATION_TOOLS = False

test_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(test_path, 'data')

#: unit id, fixture, collection
PAIRED_CSV = ('Alsoiussi_2020/csv/SRR11528761_paired.csv.gz',
              'SRR11528761_paired.head.csv.gz', 'paired')
PAIRED_CSV_PAIRED = ('Phad_2022/csv_paired/1_S1__1_Paired_All.csv.gz',
                     '1_S1__1_Paired_All.head.csv.gz', 'paired')
UNPAIRED = ('Banerjee_2017/csv/SRR5060321_Heavy_Bulk.csv.gz',
            'SRR5060321_Heavy_Bulk.head.csv.gz', 'unpaired')


def convert(case, chunksize=50000, prefix_ids=False):
    """Convert a fixture end to end and return (frame, report, metadata)."""
    unit_id, fixture, collection = case
    metadata, chunks = Oas.readDataUnit(os.path.join(data_path, fixture),
                                        chunksize=chunksize)
    report = Oas.newReport()
    frames, offset = [], 0
    for chunk in chunks:
        frames.append(Oas.normalizeChunk(metadata, chunk, unit_id, collection,
                                         offset, report,
                                         prefix_ids=prefix_ids))
        offset += len(chunk)

    return pandas.concat(frames, ignore_index=True), report, metadata


class TestReadDataUnit(unittest.TestCase):
    """
    Tests for the data unit reader
    """

    def test_reads_metadata_and_records(self):
        metadata, chunks = Oas.readDataUnit(
            os.path.join(data_path, PAIRED_CSV[1]))
        frame = pandas.concat(list(chunks), ignore_index=True)

        self.assertEqual(metadata['Run'], 'SRR11528761')
        self.assertEqual(metadata['Species'], 'mouse_C57BL/6')
        self.assertEqual(len(frame.columns), 180)

    def test_reads_a_multi_member_gzip_stream(self):
        """
        Data units are two gzip members: metadata, then CSV.

        A decoder that stops at the first member reports success while returning
        only the metadata, so this asserts the records are actually reached.
        """
        raw = Path(data_path, PAIRED_CSV[1]).read_bytes()

        # One member only, as a naive decoder would see it.
        naive = zlib.decompressobj(31).decompress(raw)
        self.assertLess(len(naive), 1000)

        # Multi member aware.
        self.assertGreater(len(decompressPrefix(raw)), 10000)

    def test_metadata_line_may_span_physical_lines(self):
        """
        The metadata is a quoted CSV field and may contain newlines.

        Reading one physical line would truncate it, so the reader consumes a
        whole CSV record instead.
        """
        payload = io.BytesIO()
        with gzip.GzipFile(fileobj=payload, mode='wb', mtime=0) as handle:
            handle.write(b'"{""Run"": ""X"",\n ""Species"": ""human""}"\n'
                         b'sequence_id,sequence,v_call\n'
                         b's1,ACGT,IGHV1-2*02\n')

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'multiline.csv.gz'
            path.write_bytes(payload.getvalue())
            metadata, chunks = Oas.readDataUnit(path)
            frame = pandas.concat(list(chunks), ignore_index=True)

        self.assertEqual(metadata['Run'], 'X')
        self.assertEqual(metadata['Species'], 'human')
        self.assertEqual(len(frame), 1)

    def test_non_json_first_record_raises(self):
        """A changed layout is reported, never swallowed into an empty result."""
        payload = io.BytesIO()
        with gzip.GzipFile(fileobj=payload, mode='wb', mtime=0) as handle:
            handle.write(b'not json at all\nsequence_id\ns1\n')

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'bad.csv.gz'
            path.write_bytes(payload.getvalue())
            with self.assertRaises(OasParseError):
                Oas.readDataUnit(path)


class TestPairedPivot(unittest.TestCase):
    """
    Tests for the wide to long conversion of paired units
    """

    def test_two_rows_per_cell(self):
        frame, report, _ = convert(PAIRED_CSV)

        self.assertEqual(report['rows_out'], 2 * report['rows_in'])
        self.assertEqual(len(frame), 38)
        self.assertEqual(frame['cell_id'].nunique(), 19)
        self.assertEqual(frame.groupby('cell_id').size().unique().tolist(), [2])

    def test_the_two_layouts_have_different_columns(self):
        """
        Paired units do not share one schema.

        The csv_paired layout carries Redundancy, c_region and per chain Isotype
        that the csv layout does not, so the stems have to be discovered per file.
        """
        narrow, _, _ = convert(PAIRED_CSV)
        wide, _, _ = convert(PAIRED_CSV_PAIRED)

        self.assertNotIn('c_region', narrow.columns)
        self.assertIn('c_region', wide.columns)

    def test_disagreeing_barcodes_still_make_one_cell(self):
        """
        A cell is resolved once per row, not once per chain.

        Deriving it from each chain's own identifier would split the cell in two
        the moment the two columns disagreed, and nothing downstream could tell
        that from a cell that genuinely had a single chain. The disagreement is
        counted instead, so it is visible rather than silent.
        """
        chunk = pandas.DataFrame(
            {'sequence_id_heavy': ['AAAC-1_contig_1', 'TTTT-1_contig_1'],
             'sequence_id_light': ['AAAC-1_contig_2', 'GGGG-1_contig_2'],
             'sequence_heavy': ['A', 'A'], 'sequence_light': ['C', 'C'],
             'v_call_heavy': ['IGHV1-2*01', 'IGHV1-2*01'],
             'v_call_light': ['IGKV1-5*01', 'IGKV1-5*01']})
        report = Oas.newReport()

        frame = Oas.normalizeChunk({}, chunk, 'S/csv/u.csv.gz', 'paired',
                                   0, report)

        self.assertEqual(report['cell_barcode_mismatch'], 1)
        self.assertEqual(frame['cell_id'].nunique(), 2)
        self.assertEqual(frame.groupby('cell_id').size().unique().tolist(), [2])
        # The heavy chain's barcode wins, and the light row joins it.
        self.assertEqual(sorted(set(frame['cell_id'])), ['AAAC-1', 'TTTT-1'])

    def test_rows_without_identifiers_fall_back_to_the_row_index(self):
        """
        A layout with no sequence_id still has to produce usable identifiers.
        """
        chunk = pandas.DataFrame({'sequence_heavy': ['A'], 'sequence_light': ['C'],
                                  'v_call_heavy': ['IGHV1-2*01'],
                                  'v_call_light': ['IGKV1-5*01']})

        frame = Oas.normalizeChunk({}, chunk, 'S/csv/u.csv.gz', 'paired', 7,
                                   Oas.newReport())

        self.assertEqual(frame['cell_id'].tolist(),
                         ['S_csv_u_cell_000000007'] * 2)
        self.assertEqual(frame['sequence_id'].tolist(),
                         ['S_csv_u_cell_000000007_heavy',
                          'S_csv_u_cell_000000007_light'])

    def test_asymmetric_chains_raise(self):
        frame = pandas.DataFrame({'sequence_heavy': ['A'], 'v_call_light': ['B']})
        with self.assertRaises(OasParseError):
            Oas.splitChains(frame)

    def test_unsuffixed_columns_raise(self):
        frame = pandas.DataFrame({'sequence_heavy': ['A'], 'sequence_light': ['B'],
                                  'stray': ['C']})
        with self.assertRaises(OasParseError):
            Oas.splitChains(frame)


class TestChunkInvariance(unittest.TestCase):
    """
    Tests that chunking cannot change the output
    """

    def test_identical_across_chunk_sizes(self):
        """
        Converting in chunks equals converting whole, byte for byte.

        Identifiers are derived from each chunk's global offset rather than from
        a row's position within its chunk, and this is what proves it.
        """
        for case in (PAIRED_CSV, PAIRED_CSV_PAIRED, UNPAIRED):
            whole, _, _ = convert(case, chunksize=100000)
            in_threes, _, _ = convert(case, chunksize=3)
            in_sevens, _, _ = convert(case, chunksize=7)

            pandas.testing.assert_frame_equal(whole, in_threes)
            pandas.testing.assert_frame_equal(whole, in_sevens)

    def test_identifiers_are_stable_and_unique(self):
        frame, _, _ = convert(PAIRED_CSV, chunksize=5)
        again, _, _ = convert(PAIRED_CSV, chunksize=50000)

        self.assertEqual(frame['sequence_id'].tolist(),
                         again['sequence_id'].tolist())
        self.assertEqual(frame['sequence_id'].nunique(), len(frame))


class TestFieldMapping(unittest.TestCase):
    """
    Tests for the individual field decisions
    """

    def test_locus_distinguishes_kappa_from_lambda(self):
        frame, _, _ = convert(PAIRED_CSV_PAIRED)
        self.assertTrue({'IGH', 'IGK', 'IGL'} <= set(frame['locus']))

    def test_c_call_only_where_a_real_isotype_exists(self):
        """
        Per chain isotypes become c_call; sentinels stay empty.

        Light chains in this fixture are all recorded as Bulk, so they must have
        no constant region call at all.
        """
        frame, _, _ = convert(PAIRED_CSV_PAIRED)
        heavy = frame[frame['locus'] == 'IGH']
        light = frame[frame['locus'] != 'IGH']

        self.assertTrue(set(heavy['c_call']) & {'IGHA', 'IGHG', 'IGHM', 'IGHD'})
        self.assertEqual(set(light['c_call']), {''})

    def test_unit_level_isotype_sentinel_is_not_copied(self):
        """The csv layout has no per chain Isotype and a unit Isotype of 'All'."""
        frame, _, metadata = convert(PAIRED_CSV)

        self.assertEqual(metadata['Isotype'], 'All')
        self.assertEqual(set(frame['c_call']), {''})

    def test_duplicate_count_from_redundancy(self):
        """
        Redundancy becomes duplicate_count.

        Values are text throughout: records stream straight to TSV, so holding a
        numeric dtype would only reintroduce float formatting on output.
        """
        frame, _, _ = convert(UNPAIRED)
        self.assertEqual(frame['duplicate_count'].iloc[0], '7')
        self.assertEqual(frame['duplicate_count'].iloc[1], '2223')

    def test_duplicate_count_defaults_when_absent(self):
        frame, report, _ = convert(PAIRED_CSV)

        self.assertEqual(set(frame['duplicate_count']), {'1'})
        self.assertEqual(report['missing_duplicate_count'], len(frame))

    def test_paired_keeps_the_source_identifier(self):
        """
        The 10x barcode and contig are the identifier, not a row counter.

        They are what joins a row back to the file it came from, and a counter
        carries none of that. Unique within a unit, which is what one output file
        per unit needs.
        """
        frame, _, _ = convert(PAIRED_CSV)

        self.assertTrue(frame['sequence_id'].str.contains('_contig_').all())
        self.assertTrue(frame['sequence_id'].is_unique)

    def test_cell_id_is_the_barcode_without_the_contig(self):
        """
        Both chains of a cell carry the same barcode, so both land on one cell.
        """
        frame, _, _ = convert(PAIRED_CSV)

        for _, row in frame.iterrows():
            self.assertEqual(row['cell_id'],
                             row['sequence_id'].rsplit('_contig_', 1)[0])
        self.assertEqual(frame.groupby('cell_id').size().unique().tolist(), [2])

    def test_original_is_recorded_only_when_the_identifier_was_rewritten(self):
        """
        A provenance column repeating sequence_id in every row is noise.

        A value here means the identifier was changed, which is the only case
        where a user needs the original to join back on.
        """
        plain, _, _ = convert(PAIRED_CSV)
        self.assertEqual(set(plain['sourcerer_original_sequence_id']), {''})

        prefixed, _, _ = convert(PAIRED_CSV, prefix_ids=True)
        self.assertTrue(
            prefixed['sourcerer_original_sequence_id'].str.contains(
                '_contig_').all())

    def test_prefixing_namespaces_identifiers_for_merged_output(self):
        """
        10x barcodes come from a fixed whitelist and recur in every unit.

        Combining units without a prefix would merge unrelated cells silently,
        so anything writing several units into one file must prefix.
        """
        plain, _, _ = convert(PAIRED_CSV)
        prefixed, _, _ = convert(PAIRED_CSV, prefix_ids=True)

        stem = Oas.unitStem(PAIRED_CSV[0])
        self.assertTrue(prefixed['sequence_id'].str.startswith(stem).all())
        self.assertTrue(prefixed['cell_id'].str.startswith(stem).all())
        self.assertEqual(prefixed['cell_id'].nunique(),
                         plain['cell_id'].nunique())

    def test_unpaired_has_no_source_identifier(self):
        frame, _, _ = convert(UNPAIRED)
        self.assertEqual(set(frame['sourcerer_original_sequence_id']), {''})
        self.assertTrue(frame['sequence_id'].is_unique)

    def test_consumed_columns_are_dropped(self):
        frame, _, _ = convert(PAIRED_CSV_PAIRED)

        self.assertNotIn('Redundancy', frame.columns)
        self.assertNotIn('Isotype', frame.columns)
        self.assertFalse([x for x in frame.columns if x.startswith('_')])


class TestAirrOutput(unittest.TestCase):
    """
    Tests for the AIRR writer
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def writeCase(self, case, **kwargs):
        frame, _, _ = convert(case)
        out = Path(self.tmp.name) / 'out.tsv'
        validation = Convert.writeAirr([frame], out, **kwargs)

        return out, validation

    def test_output_validates(self):
        """
        Every fixture converts to a file the airr package accepts.

        Both paired layouts and unpaired are checked, because they have different
        column sets and only one of them was covered by the original R pipeline.
        """
        import airr

        for case in (PAIRED_CSV, PAIRED_CSV_PAIRED, UNPAIRED):
            out, validation = self.writeCase(case)

            self.assertTrue(validation['header_valid'], case[0])
            self.assertEqual(validation['rows_invalid'], 0, case[0])
            self.assertTrue(airr.validate_rearrangement(str(out)), case[0])

    def test_integers_are_not_written_in_float_form(self):
        """
        OAS writes whole numbers as '1.0', which AIRR rejects.

        Without coercion every integer field of every row fails validation, which
        would bury any genuine problem in noise.
        """
        out, _ = self.writeCase(UNPAIRED)
        frame = pandas.read_csv(out, sep='\t', dtype=str, na_filter=False)

        self.assertFalse(frame['v_alignment_start'].str.contains(r'\.').any())

    def test_strict_drops_non_airr_columns(self):
        out, _ = self.writeCase(PAIRED_CSV, strict=True)
        frame = pandas.read_csv(out, sep='\t', dtype=str, na_filter=False)

        self.assertNotIn('ANARCI_status', frame.columns)
        self.assertNotIn('sourcerer_unit_id', frame.columns)

    def test_extras_are_kept_by_default(self):
        out, _ = self.writeCase(PAIRED_CSV)
        frame = pandas.read_csv(out, sep='\t', dtype=str, na_filter=False)

        self.assertIn('sourcerer_unit_id', frame.columns)
        self.assertIn('sourcerer_row_hash', frame.columns)


class TestFastaOutput(unittest.TestCase):
    """
    Tests for the FASTA writer
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def test_missing_annotation_is_omitted_not_written_as_a_placeholder(self):
        """
        An annotation with no value is left out of the header entirely.

        MakeDb.py parses these annotations back into the rearrangement table, so
        a placeholder would land in the output as a literal 'NA' where AIRR wants
        an empty value. pRESTO reads the header as a dictionary, so omitting the
        key is unambiguous rather than a misaligned field.
        """
        frame, _, _ = convert(PAIRED_CSV)
        frame.loc[frame.index[0], 'cell_id'] = ''
        out = Path(self.tmp.name) / 'out.fasta'
        written = Convert.writeFasta([frame], out)

        headers = [x for x in out.read_text().splitlines() if x.startswith('>')]
        self.assertEqual(written, len(headers))
        self.assertNotIn('=NA', headers[0])
        self.assertNotIn('|', headers[0])

    def test_only_cell_id_is_carried(self):
        """
        Fields IgBLAST derives for itself are not carried.

        MakeDb.py applies header annotations on top of the record it parsed from
        the aligner, so carrying c_call or locus would overwrite an alignment
        based call with the source's own. cell_id is the one thing IgBLAST
        cannot recover from the sequence.
        """
        self.assertEqual(Convert.FASTA_ANNOTATIONS, ('cell_id',))

        frame, _, _ = convert(PAIRED_CSV)
        out = Path(self.tmp.name) / 'out.fasta'
        Convert.writeFasta([frame], out)

        headers = [x for x in out.read_text().splitlines() if x.startswith('>')]
        for header in headers:
            self.assertNotIn('c_call=', header)
            self.assertNotIn('locus=', header)

    def test_header_is_presto_annotation_format(self):
        """
        The header is an identifier followed by pipe separated key=value pairs.

        Keys are the AIRR column names, so a header and the rearrangement TSV
        can never disagree about what a field is called.
        """
        frame, _, _ = convert(PAIRED_CSV)
        out = Path(self.tmp.name) / 'out.fasta'
        Convert.writeFasta([frame], out)

        header = out.read_text().splitlines()[0]
        identifier, *annotations = header.lstrip('>').split('|')
        record = frame.iloc[0]

        self.assertEqual(identifier, record['sequence_id'])
        self.assertEqual([x.split('=', 1)[0] for x in annotations],
                         list(Convert.FASTA_ANNOTATIONS))
        self.assertEqual(annotations[0], 'cell_id=%s' % record['cell_id'])

    @unittest.skipUnless(HAS_IMMCANTATION_TOOLS,
                         'presto and changeo are not installed')
    def test_header_round_trips_to_the_airr_cell_id_column(self):
        """
        The header survives the trip through MakeDb.py back into AIRR.

        This is the whole reason the annotation exists: airrflow's assembled
        mode converts FASTA with IgBLAST, and cell pairing is preserved only if
        Change-O maps the header key onto the cell_id column. Asserting the
        header text alone would not catch a key that parses but maps elsewhere.
        """
        from changeo.Receptor import AIRRSchema, ChangeoSchema
        from presto.Annotation import parseAnnotation

        frame, _, _ = convert(PAIRED_CSV_PAIRED)
        out = Path(self.tmp.name) / 'out.fasta'
        Convert.writeFasta([frame], out)

        header = next(x for x in out.read_text().splitlines()
                      if x.startswith('>'))
        parsed = parseAnnotation(header.lstrip('>'))
        parsed.pop('ID')
        columns = {AIRRSchema.fromReceptor(ChangeoSchema.toReceptor(k)): v
                   for k, v in parsed.items()}

        self.assertEqual(list(columns), ['cell_id'])
        self.assertIn(columns['cell_id'], set(frame['cell_id']))


if __name__ == '__main__':
    unittest.main()
