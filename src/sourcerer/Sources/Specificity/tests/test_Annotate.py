"""
Unit tests for antigen specificity annotation (Annotate.py, AnnotateCli.py)
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import io
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path
from unittest import mock

import pandas

# Sourcerer imports
from sourcerer.Cli import getArgParser
from sourcerer.Exceptions import AnnotationError
from sourcerer.Sources.Specificity.Annotate import (runExactAnnotation,
                                                    runFuzzyAnnotation,
                                                    runLevenshteinAnnotation,
                                                    runPairedExactAnnotation)
from sourcerer.Sources.Specificity.AnnotateCli import handleAnnotate


def frame(rows):
    return pandas.DataFrame(rows)


def fakeCompairrRun(pairs_rows):
    """
    Build a subprocess.run replacement that writes `pairs_rows` to the `-p`
    path instead of actually invoking compairr, so runExactAnnotation can be
    tested without the real binary installed.
    """
    import subprocess

    def _run(cmd, check=True, capture_output=True, text=True):
        pairs_path = Path(cmd[cmd.index('-p') + 1])
        pandas.DataFrame(pairs_rows).to_csv(pairs_path, sep='\t', index=False)
        return subprocess.CompletedProcess(cmd, 0)

    return _run


def noCompairrCallAllowed(cmd, **kwargs):
    raise AssertionError('compairr should not have been invoked')


def fakeCompairrExact(cmd, check=True, capture_output=True, text=True):
    """
    A subprocess.run replacement that actually performs exact (junction_aa,
    plus v_call/j_call when present) matching between the two AIRR TSVs
    _writeAirrTsv wrote, instead of returning one canned pairs table.

    Needed wherever a test drives more than one compairr invocation with
    genuinely different content per call (e.g. runPairedExactAnnotation,
    which calls runExactAnnotation once for heavy rows and once for light) --
    fakeCompairrRun's single fixed pairs table can't distinguish between them.
    """
    import subprocess

    query_path, reference_path = Path(cmd[2]), Path(cmd[3])
    pairs_path = Path(cmd[cmd.index('-p') + 1])

    query = pandas.read_csv(query_path, sep='\t', dtype=str)
    reference = pandas.read_csv(reference_path, sep='\t', dtype=str)
    use_genes = 'v_call' in query.columns

    pairs = []
    for _, q in query.iterrows():
        for _, r in reference.iterrows():
            if q['junction_aa'] != r['junction_aa']:
                continue
            if use_genes and (q['v_call'] != r['v_call'] or q['j_call'] != r['j_call']):
                continue
            pairs.append({'sequence_id_1': q['sequence_id'],
                          'sequence_id_2': r['sequence_id'], 'distance': 0})

    pandas.DataFrame(pairs, columns=['sequence_id_1', 'sequence_id_2', 'distance']
                    ).to_csv(pairs_path, sep='\t', index=False)

    return subprocess.CompletedProcess(cmd, 0)


class TestRunExactAnnotation(unittest.TestCase):
    """
    Tests for the CompAIRR-backed exact annotation
    """

    def test_raises_when_compairr_is_not_found(self):
        queries = frame([{'cdr3': 'CARDYW'}])
        reference = frame([{'id': 'r1', 'cdr3': 'CARDYW'}])

        with mock.patch('shutil.which', return_value=None):
            with self.assertRaises(AnnotationError):
                runExactAnnotation(queries, reference, 'cdr3', id_col='id')

    def test_clean_rows_go_through_compairr(self):
        queries = frame([{'cdr3': 'CARDYW'}, {'cdr3': 'CQQYNS'}])
        reference = frame([{'id': 'r1', 'cdr3': 'CARDYW'}])
        pairs = [{'sequence_id_1': '0', 'sequence_id_2': 'r1', 'distance': 0}]

        with mock.patch('shutil.which', return_value='/usr/bin/compairr'), \
             mock.patch('subprocess.run', side_effect=fakeCompairrRun(pairs)):
            result = runExactAnnotation(queries, reference, 'cdr3', id_col='id')

        self.assertEqual(result['hit_ids'].tolist(), ['r1', ''])
        self.assertEqual(result['n_hits_total'].tolist(), [1, 0])

    def test_x_containing_query_rows_skip_compairr_entirely(self):
        """
        A query set that is all 'X'/blank never has a clean-vs-clean bulk to
        run, so compairr must not be invoked at all -- only the pure-Python
        fallback handles it.
        """
        queries = frame([{'cdr3': 'CARXYW'}])
        reference = frame([{'id': 'r1', 'cdr3': 'CARDYW'}])

        with mock.patch('shutil.which', return_value='/usr/bin/compairr'), \
             mock.patch('subprocess.run', side_effect=noCompairrCallAllowed):
            result = runExactAnnotation(queries, reference, 'cdr3', id_col='id')

        self.assertEqual(result['hit_ids'].tolist(), ['r1'])

    def test_x_containing_reference_rows_are_still_reachable(self):
        """
        A clean query can still hit an 'X'-containing reference row, via the
        clean-queries-vs-X-reference fallback pass.
        """
        queries = frame([{'cdr3': 'CARDYW'}])
        reference = frame([{'id': 'r1', 'cdr3': 'CARXYW'}])

        with mock.patch('shutil.which', return_value='/usr/bin/compairr'), \
             mock.patch('subprocess.run', side_effect=noCompairrCallAllowed):
            result = runExactAnnotation(queries, reference, 'cdr3', id_col='id')

        self.assertEqual(result['hit_ids'].tolist(), ['r1'])

    def test_gene_gating_is_forwarded_to_the_compairr_command(self):
        queries = frame([{'cdr3': 'CARDYW', 'v_call': 'IGHV1-2'}])
        reference = frame([{'id': 'r1', 'cdr3': 'CARDYW', 'v_call': 'IGHV1-2'}])
        pairs = [{'sequence_id_1': '0', 'sequence_id_2': 'r1', 'distance': 0}]

        seen_cmds = []

        def _run(cmd, **kwargs):
            seen_cmds.append(cmd)
            return fakeCompairrRun(pairs)(cmd, **kwargs)

        with mock.patch('shutil.which', return_value='/usr/bin/compairr'), \
             mock.patch('subprocess.run', side_effect=_run):
            runExactAnnotation(queries, reference, 'cdr3', id_col='id',
                               vgene_col='v_call')

        # use_genes=True -> no '-g' flag appended.
        self.assertNotIn('-g', seen_cmds[0])

    def test_no_gene_columns_passes_the_g_flag(self):
        queries = frame([{'cdr3': 'CARDYW'}])
        reference = frame([{'id': 'r1', 'cdr3': 'CARDYW'}])
        pairs = [{'sequence_id_1': '0', 'sequence_id_2': 'r1', 'distance': 0}]

        seen_cmds = []

        def _run(cmd, **kwargs):
            seen_cmds.append(cmd)
            return fakeCompairrRun(pairs)(cmd, **kwargs)

        with mock.patch('shutil.which', return_value='/usr/bin/compairr'), \
             mock.patch('subprocess.run', side_effect=_run):
            runExactAnnotation(queries, reference, 'cdr3', id_col='id')

        self.assertIn('-g', seen_cmds[0])


class TestRunFuzzyAnnotation(unittest.TestCase):
    """
    Tests for substitution-tolerant (Hamming) annotation
    """

    def test_within_budget_hits(self):
        queries = frame([{'cdr3': 'CARDYW'}])
        reference = frame([{'id': 'r1', 'cdr3': 'CARQYW'}])  # 1 substitution

        result = runFuzzyAnnotation(queries, reference, ['cdr3'], max_mismatches=1,
                                    id_col='id')

        self.assertEqual(result['hit_ids'].tolist(), ['r1'])
        self.assertEqual(result['min_mismatches'].tolist(), [1.0])

    def test_over_budget_does_not_hit(self):
        queries = frame([{'cdr3': 'CARDYW'}])
        reference = frame([{'id': 'r1', 'cdr3': 'CQRQYW'}])  # 2 substitutions

        result = runFuzzyAnnotation(queries, reference, ['cdr3'], max_mismatches=1,
                                    id_col='id')

        self.assertEqual(result['n_hits_total'].tolist(), [0])
        self.assertTrue(pandas.isna(result['min_mismatches'].iloc[0]))

    def test_exact_cols_gate_candidates_before_fuzzy_comparison(self):
        queries = frame([{'cdr3': 'CARDYW', 'v_call': 'IGHV1-2'}])
        reference = frame([
            {'id': 'r1', 'cdr3': 'CARQYW', 'v_call': 'IGHV3-9'},  # wrong V-gene
        ])

        result = runFuzzyAnnotation(queries, reference, ['cdr3'], max_mismatches=1,
                                    id_col='id', exact_cols=['v_call'])

        self.assertEqual(result['n_hits_total'].tolist(), [0])

    def test_threshold_is_fraction_scales_with_query_length(self):
        # length 10, 15% -> floor(1.5) == 1 substitution allowed
        queries = frame([{'cdr3': 'AAAAAAAAAA'}])
        reference = frame([{'id': 'close', 'cdr3': 'AAAAAAAAAB'},
                           {'id': 'far', 'cdr3': 'AABBAAAAAB'}])

        result = runFuzzyAnnotation(queries, reference, ['cdr3'], max_mismatches=0.15,
                                    id_col='id', threshold_is_fraction=True)

        self.assertEqual(result['hit_ids'].tolist(), ['close'])

    def test_mismatch_positions_single_column_are_bare(self):
        queries = frame([{'cdr3': 'AAAA'}])
        reference = frame([{'id': 'r1', 'cdr3': 'ABAB'}])

        result = runFuzzyAnnotation(queries, reference, ['cdr3'], max_mismatches=2,
                                    id_col='id')

        self.assertEqual(result['mismatch_positions'].tolist(), ['1,3'])

    def test_mismatch_positions_multi_column_are_column_prefixed(self):
        queries = frame([{'a': 'AA', 'b': 'AA'}])
        reference = frame([{'id': 'r1', 'a': 'AB', 'b': 'BA'}])

        result = runFuzzyAnnotation(queries, reference, ['a', 'b'], max_mismatches=2,
                                    id_col='id')

        self.assertEqual(result['mismatch_positions'].tolist(), ['0:1,1:0'])

    def test_best_hit_per_id_is_kept_when_several_reference_rows_share_it(self):
        queries = frame([{'cdr3': 'AAAA'}])
        reference = frame([
            {'id': 'r1', 'cdr3': 'AAAB'},  # distance 1
            {'id': 'r1', 'cdr3': 'ABBB'},  # distance 3, same id, worse
        ])

        result = runFuzzyAnnotation(queries, reference, ['cdr3'], max_mismatches=3,
                                    id_col='id')

        self.assertEqual(result['n_hits_total'].tolist(), [1])
        self.assertEqual(result['min_mismatches'].tolist(), [1.0])


class TestRunPairedExactAnnotation(unittest.TestCase):
    """
    Tests for the heavy+light reconciliation (exact V-gene, J-gene and CDR3
    agreement, both chains against the same reference antibody)
    """

    def _antibodies(self):
        """
        5 query antibodies, one per classification, against 2 reference
        antibodies -- see the inline comments for which reference row each
        is built to agree or disagree with.
        """
        queries = frame([
            # q1: heavy and light both agree with refA -> both
            {'cell_id': 'q1', 'locus': 'IGH', 'cdr3_aa': 'CARDYW',
             'v_call': 'IGHV1', 'j_call': 'IGHJ1'},
            {'cell_id': 'q1', 'locus': 'IGL', 'cdr3_aa': 'CQQYNS',
             'v_call': 'IGKV1', 'j_call': 'IGKJ1'},
            # q2: heavy agrees with refA, light agrees with nothing -> heavy_only
            {'cell_id': 'q2', 'locus': 'IGH', 'cdr3_aa': 'CARDYW',
             'v_call': 'IGHV1', 'j_call': 'IGHJ1'},
            {'cell_id': 'q2', 'locus': 'IGL', 'cdr3_aa': 'CDIFFER',
             'v_call': 'IGKV9', 'j_call': 'IGKJ9'},
            # q3: heavy agrees with nothing, light agrees with refA -> light_only
            {'cell_id': 'q3', 'locus': 'IGH', 'cdr3_aa': 'CDIFF2',
             'v_call': 'IGHV9', 'j_call': 'IGHJ9'},
            {'cell_id': 'q3', 'locus': 'IGL', 'cdr3_aa': 'CQQYNS',
             'v_call': 'IGKV1', 'j_call': 'IGKJ1'},
            # q4: neither chain agrees with anything -> neither
            {'cell_id': 'q4', 'locus': 'IGH', 'cdr3_aa': 'CDIFF3',
             'v_call': 'IGHV8', 'j_call': 'IGHJ8'},
            {'cell_id': 'q4', 'locus': 'IGL', 'cdr3_aa': 'CDIFF4',
             'v_call': 'IGKV8', 'j_call': 'IGKJ8'},
            # q5: heavy agrees with refA, light agrees with refB -> conflicting
            {'cell_id': 'q5', 'locus': 'IGH', 'cdr3_aa': 'CARDYW',
             'v_call': 'IGHV1', 'j_call': 'IGHJ1'},
            {'cell_id': 'q5', 'locus': 'IGL', 'cdr3_aa': 'CZZZZZZ',
             'v_call': 'IGKV7', 'j_call': 'IGKJ7'},
        ])
        reference = frame([
            {'cell_id': 'refA', 'locus': 'IGH', 'cdr3_aa': 'CARDYW',
             'v_call': 'IGHV1', 'j_call': 'IGHJ1'},
            {'cell_id': 'refA', 'locus': 'IGL', 'cdr3_aa': 'CQQYNS',
             'v_call': 'IGKV1', 'j_call': 'IGKJ1'},
            {'cell_id': 'refB', 'locus': 'IGH', 'cdr3_aa': 'CWWWWWW',
             'v_call': 'IGHV6', 'j_call': 'IGHJ6'},
            {'cell_id': 'refB', 'locus': 'IGL', 'cdr3_aa': 'CZZZZZZ',
             'v_call': 'IGKV7', 'j_call': 'IGKJ7'},
        ])
        return queries, reference

    def test_classifications(self):
        queries, reference = self._antibodies()

        with mock.patch('shutil.which', return_value='/usr/bin/compairr'), \
             mock.patch('subprocess.run', side_effect=fakeCompairrExact):
            result = runPairedExactAnnotation(queries, reference)

        by_cell_id = result.set_index('cell_id')['classification'].to_dict()
        self.assertEqual(by_cell_id, {
            'q1': 'both', 'q2': 'heavy_only', 'q3': 'light_only',
            'q4': 'neither', 'q5': 'conflicting',
        })

    def test_both_hit_ids_is_the_intersection_not_the_union(self):
        queries, reference = self._antibodies()

        with mock.patch('shutil.which', return_value='/usr/bin/compairr'), \
             mock.patch('subprocess.run', side_effect=fakeCompairrExact):
            result = runPairedExactAnnotation(queries, reference)

        row = result.set_index('cell_id').loc['q1']
        self.assertEqual(row['heavy_hit_ids'], 'refA')
        self.assertEqual(row['light_hit_ids'], 'refA')
        self.assertEqual(row['both_hit_ids'], 'refA')

        row5 = result.set_index('cell_id').loc['q5']
        self.assertEqual(row5['heavy_hit_ids'], 'refA')
        self.assertEqual(row5['light_hit_ids'], 'refB')
        self.assertEqual(row5['both_hit_ids'], '')


class TestRunLevenshteinAnnotation(unittest.TestCase):
    """
    Tests for indel-tolerant (Levenshtein) annotation
    """

    def test_a_substitution_within_budget_hits(self):
        queries = frame([{'cdr3': 'CARDYW'}])
        reference = frame([{'id': 'r1', 'cdr3': 'CARQYW'}])  # 1 substitution

        result = runLevenshteinAnnotation(queries, reference, 'cdr3', max_edits=1,
                                          id_col='id')

        self.assertEqual(result['hit_ids'].tolist(), ['r1'])
        self.assertEqual(result['min_mismatches'].tolist(), [1.0])

    def test_an_insertion_within_budget_hits(self):
        """Unlike Hamming, Levenshtein tolerates a length difference."""
        queries = frame([{'cdr3': 'CARDYW'}])
        reference = frame([{'id': 'r1', 'cdr3': 'CARDXYW'}])  # 1 insertion

        result = runLevenshteinAnnotation(queries, reference, 'cdr3', max_edits=1,
                                          id_col='id')

        self.assertEqual(result['hit_ids'].tolist(), ['r1'])
        self.assertEqual(result['min_mismatches'].tolist(), [1.0])

    def test_over_budget_does_not_hit(self):
        queries = frame([{'cdr3': 'CARDYW'}])
        reference = frame([{'id': 'r1', 'cdr3': 'CQRQYW'}])  # 2 substitutions

        result = runLevenshteinAnnotation(queries, reference, 'cdr3', max_edits=1,
                                          id_col='id')

        self.assertEqual(result['n_hits_total'].tolist(), [0])
        self.assertTrue(pandas.isna(result['min_mismatches'].iloc[0]))

    def test_exact_cols_gate_candidates_first(self):
        queries = frame([{'cdr3': 'CARDYW', 'v_call': 'IGHV1-2'}])
        reference = frame([
            {'id': 'r1', 'cdr3': 'CARQYW', 'v_call': 'IGHV3-9'},  # wrong V-gene
        ])

        result = runLevenshteinAnnotation(queries, reference, 'cdr3', max_edits=1,
                                          id_col='id', exact_cols=['v_call'])

        self.assertEqual(result['n_hits_total'].tolist(), [0])

    def test_best_hit_per_id_is_kept(self):
        queries = frame([{'cdr3': 'AAAA'}])
        reference = frame([
            {'id': 'r1', 'cdr3': 'AAAB'},   # distance 1
            {'id': 'r1', 'cdr3': 'AABBB'},  # distance 3, same id, worse
        ])

        result = runLevenshteinAnnotation(queries, reference, 'cdr3', max_edits=3,
                                          id_col='id')

        self.assertEqual(result['n_hits_total'].tolist(), [1])
        self.assertEqual(result['min_mismatches'].tolist(), [1.0])


class TestAnnotateCliParser(unittest.TestCase):
    """
    Tests for `annotate`'s parser wiring into the top level CLI
    """

    def test_annotate_is_a_top_level_subcommand(self):
        parser = getArgParser()
        with mock.patch('sys.stdout', new_callable=io.StringIO) as out:
            with self.assertRaises(SystemExit):
                parser.parse_args(['annotate', '--help'])

        text = out.getvalue()
        self.assertIn('--query', text)
        self.assertIn('--reference', text)
        self.assertIn('--compare-cols', text)
        self.assertIn('--exact', text)
        self.assertIn('--fuzzy', text)
        self.assertIn('--identity', text)
        self.assertIn('--levenshtein', text)
        self.assertIn('--paired', text)
        self.assertIn('--vgene-col', text)
        self.assertIn('--compairr-bin', text)

    def test_exact_fuzzy_paired_are_mutually_exclusive(self):
        parser = getArgParser()
        with mock.patch('sys.stderr', new_callable=io.StringIO):
            with self.assertRaises(SystemExit):
                parser.parse_args(['annotate', '--query', 'q.tsv',
                                   '--reference', 'r.tsv', '--fuzzy', '2',
                                   '--paired', '-o', 'out.tsv'])

    def test_required_flags_are_enforced(self):
        parser = getArgParser()
        with mock.patch('sys.stderr', new_callable=io.StringIO):
            with self.assertRaises(SystemExit):
                parser.parse_args(['annotate'])


class TestHandleAnnotate(unittest.TestCase):
    """
    Tests for the end-to-end file-in, file-out handler
    """

    def _writeTsv(self, path, rows):
        pandas.DataFrame(rows).to_csv(path, sep='\t', index=False)

    def _baseArgs(self, **overrides):
        defaults = dict(compare_cols=['cdr3_aa'], fuzzy=None, identity=None,
                        levenshtein=None, paired=False, id_col='sequence_id',
                        vgene_col=None, jgene_col=None, compairr_bin=None,
                        threads=1, locus=None)
        defaults.update(overrides)
        return Namespace(**defaults)

    def test_exact_annotation_writes_a_result_tsv(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            query_path = tmp / 'query.tsv'
            reference_path = tmp / 'reference.tsv'
            out_path = tmp / 'out.tsv'

            self._writeTsv(query_path, [{'sequence_id': 'q1', 'cdr3_aa': 'CARDYW'}])
            self._writeTsv(reference_path,
                           [{'sequence_id': 'ref1', 'cdr3_aa': 'CARDYW'}])

            pairs = [{'sequence_id_1': '0', 'sequence_id_2': 'ref1', 'distance': 0}]
            with mock.patch('shutil.which', return_value='/usr/bin/compairr'), \
                 mock.patch('subprocess.run', side_effect=fakeCompairrRun(pairs)):
                status = handleAnnotate(self._baseArgs(
                    query=query_path, reference=reference_path, out=out_path))

            self.assertEqual(status, 0)
            written = pandas.read_csv(out_path, sep='\t')
            self.assertEqual(written['hit_ids'].tolist(), ['ref1'])

    def test_exact_annotation_rejects_more_than_one_compare_col(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            query_path = tmp / 'query.tsv'
            reference_path = tmp / 'reference.tsv'

            self._writeTsv(query_path, [{'sequence_id': 'q1', 'cdr3_aa': 'CARDYW',
                                        'cdr3b_aa': 'X'}])
            self._writeTsv(reference_path, [{'sequence_id': 'ref1', 'cdr3_aa': 'CARDYW',
                                            'cdr3b_aa': 'X'}])

            with self.assertRaises(AnnotationError):
                handleAnnotate(self._baseArgs(
                    query=query_path, reference=reference_path,
                    compare_cols=['cdr3_aa', 'cdr3b_aa'],
                    out=tmp / 'out.tsv'))

    def test_fuzzy_annotation_still_works_without_compairr(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            query_path = tmp / 'query.tsv'
            reference_path = tmp / 'reference.tsv'
            out_path = tmp / 'out.tsv'

            self._writeTsv(query_path, [{'sequence_id': 'q1', 'locus': 'IGH',
                                        'cdr3_aa': 'CARDYW', 'v_call': 'IGHV1',
                                        'j_call': 'IGHJ1'}])
            self._writeTsv(reference_path, [{'sequence_id': 'ref1', 'locus': 'IGH',
                                            'cdr3_aa': 'CARQYW', 'v_call': 'IGHV1',
                                            'j_call': 'IGHJ1'}])

            status = handleAnnotate(self._baseArgs(
                query=query_path, reference=reference_path, out=out_path,
                fuzzy=1))

            self.assertEqual(status, 0)
            written = pandas.read_csv(out_path, sep='\t')
            self.assertEqual(written['hit_ids'].tolist(), ['ref1'])

    def test_fuzzy_drops_light_chain_rows_regardless_of_locus_flag(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            query_path = tmp / 'query.tsv'
            reference_path = tmp / 'reference.tsv'
            out_path = tmp / 'out.tsv'

            self._writeTsv(query_path, [
                {'sequence_id': 'q1', 'locus': 'IGH', 'cdr3_aa': 'CARDYW',
                 'v_call': 'IGHV1', 'j_call': 'IGHJ1'},
                {'sequence_id': 'q2', 'locus': 'IGL', 'cdr3_aa': 'CQQYNS',
                 'v_call': 'IGKV1', 'j_call': 'IGKJ1'},
            ])
            self._writeTsv(reference_path, [
                {'sequence_id': 'ref1', 'locus': 'IGH', 'cdr3_aa': 'CARQYW',
                 'v_call': 'IGHV1', 'j_call': 'IGHJ1'},
                {'sequence_id': 'ref2', 'locus': 'IGL', 'cdr3_aa': 'CQQYNS',
                 'v_call': 'IGKV1', 'j_call': 'IGKJ1'},
            ])

            # --locus left unset: fuzzy still only ever sees the IGH rows.
            handleAnnotate(self._baseArgs(
                query=query_path, reference=reference_path, out=out_path,
                fuzzy=1))

            written = pandas.read_csv(out_path, sep='\t')
            self.assertEqual(len(written), 1)
            self.assertEqual(written['hit_ids'].tolist(), ['ref1'])

    def test_fuzzy_with_a_conflicting_locus_raises(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            query_path = tmp / 'query.tsv'
            reference_path = tmp / 'reference.tsv'

            self._writeTsv(query_path, [{'sequence_id': 'q1', 'locus': 'IGL',
                                        'cdr3_aa': 'CQQYNS'}])
            self._writeTsv(reference_path, [{'sequence_id': 'ref1', 'locus': 'IGL',
                                            'cdr3_aa': 'CQQYNS'}])

            with self.assertRaises(AnnotationError):
                handleAnnotate(self._baseArgs(
                    query=query_path, reference=reference_path,
                    out=tmp / 'out.tsv', fuzzy=1, locus='IGL'))

    def test_locus_filters_both_sides_before_annotating(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            query_path = tmp / 'query.tsv'
            reference_path = tmp / 'reference.tsv'
            out_path = tmp / 'out.tsv'

            self._writeTsv(query_path, [
                {'sequence_id': 'q1', 'locus': 'IGH', 'cdr3_aa': 'CARDYW'},
                {'sequence_id': 'q2', 'locus': 'IGL', 'cdr3_aa': 'CQQYNS'},
            ])
            self._writeTsv(reference_path, [
                {'sequence_id': 'ref1', 'locus': 'IGH', 'cdr3_aa': 'CARDYW'},
                {'sequence_id': 'ref2', 'locus': 'IGL', 'cdr3_aa': 'CQQYNS'},
            ])

            pairs = [{'sequence_id_1': '0', 'sequence_id_2': 'ref1', 'distance': 0}]
            with mock.patch('shutil.which', return_value='/usr/bin/compairr'), \
                 mock.patch('subprocess.run', side_effect=fakeCompairrRun(pairs)):
                handleAnnotate(self._baseArgs(
                    query=query_path, reference=reference_path, out=out_path,
                    locus='IGH'))

            written = pandas.read_csv(out_path, sep='\t')
            # Only the IGH row of the query survives the filter.
            self.assertEqual(len(written), 1)
            self.assertEqual(written['hit_ids'].tolist(), ['ref1'])

    def test_fuzzy_always_gates_on_v_and_j_gene(self):
        """V/J-gene gating under --fuzzy is automatic, not opt-in via a flag
        -- a CDR3 within budget still misses if the V-gene disagrees."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            query_path = tmp / 'query.tsv'
            reference_path = tmp / 'reference.tsv'
            out_path = tmp / 'out.tsv'

            self._writeTsv(query_path, [{'sequence_id': 'q1', 'locus': 'IGH',
                                        'cdr3_aa': 'CARDYW', 'v_call': 'IGHV1-2',
                                        'j_call': 'IGHJ1'}])
            self._writeTsv(reference_path, [{'sequence_id': 'ref1', 'locus': 'IGH',
                                            'cdr3_aa': 'CARQYW', 'v_call': 'IGHV3-9',
                                            'j_call': 'IGHJ1'}])

            handleAnnotate(self._baseArgs(
                query=query_path, reference=reference_path, out=out_path,
                fuzzy=1))

            written = pandas.read_csv(out_path, sep='\t')
            self.assertEqual(written['n_hits_total'].tolist(), [0])

    def test_identity_translates_percent_to_a_mismatch_fraction(self):
        """length 10, 85% identity -> up to floor(1.5) == 1 substitution."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            query_path = tmp / 'query.tsv'
            reference_path = tmp / 'reference.tsv'
            out_path = tmp / 'out.tsv'

            self._writeTsv(query_path, [{'sequence_id': 'q1', 'locus': 'IGH',
                                        'cdr3_aa': 'AAAAAAAAAA', 'v_call': 'IGHV1',
                                        'j_call': 'IGHJ1'}])
            self._writeTsv(reference_path, [
                {'sequence_id': 'close', 'locus': 'IGH', 'cdr3_aa': 'AAAAAAAAAB',
                 'v_call': 'IGHV1', 'j_call': 'IGHJ1'},
                {'sequence_id': 'far', 'locus': 'IGH', 'cdr3_aa': 'AABBAAAAAB',
                 'v_call': 'IGHV1', 'j_call': 'IGHJ1'},
            ])

            handleAnnotate(self._baseArgs(
                query=query_path, reference=reference_path, out=out_path,
                identity=85))

            written = pandas.read_csv(out_path, sep='\t')
            self.assertEqual(written['hit_ids'].tolist(), ['close'])

    def test_identity_out_of_range_raises(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            query_path = tmp / 'query.tsv'
            reference_path = tmp / 'reference.tsv'

            self._writeTsv(query_path, [{'sequence_id': 'q1', 'locus': 'IGH',
                                        'cdr3_aa': 'CARDYW', 'v_call': 'IGHV1',
                                        'j_call': 'IGHJ1'}])
            self._writeTsv(reference_path, [{'sequence_id': 'ref1', 'locus': 'IGH',
                                            'cdr3_aa': 'CARDYW', 'v_call': 'IGHV1',
                                            'j_call': 'IGHJ1'}])

            with self.assertRaises(AnnotationError):
                handleAnnotate(self._baseArgs(
                    query=query_path, reference=reference_path,
                    out=tmp / 'out.tsv', identity=100))

    def test_levenshtein_tolerates_an_indel(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            query_path = tmp / 'query.tsv'
            reference_path = tmp / 'reference.tsv'
            out_path = tmp / 'out.tsv'

            self._writeTsv(query_path, [{'sequence_id': 'q1', 'locus': 'IGH',
                                        'cdr3_aa': 'CARDYW', 'v_call': 'IGHV1',
                                        'j_call': 'IGHJ1'}])
            self._writeTsv(reference_path, [{'sequence_id': 'ref1', 'locus': 'IGH',
                                            'cdr3_aa': 'CARDXYW', 'v_call': 'IGHV1',
                                            'j_call': 'IGHJ1'}])  # 1 insertion

            status = handleAnnotate(self._baseArgs(
                query=query_path, reference=reference_path, out=out_path,
                levenshtein=1))

            self.assertEqual(status, 0)
            written = pandas.read_csv(out_path, sep='\t')
            self.assertEqual(written['hit_ids'].tolist(), ['ref1'])

    def test_levenshtein_rejects_more_than_one_compare_col(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            query_path = tmp / 'query.tsv'
            reference_path = tmp / 'reference.tsv'

            self._writeTsv(query_path, [{'sequence_id': 'q1', 'locus': 'IGH',
                                        'cdr3_aa': 'CARDYW', 'cdr1_aa': 'X',
                                        'v_call': 'IGHV1', 'j_call': 'IGHJ1'}])
            self._writeTsv(reference_path, [{'sequence_id': 'ref1', 'locus': 'IGH',
                                            'cdr3_aa': 'CARDYW', 'cdr1_aa': 'X',
                                            'v_call': 'IGHV1', 'j_call': 'IGHJ1'}])

            with self.assertRaises(AnnotationError):
                handleAnnotate(self._baseArgs(
                    query=query_path, reference=reference_path,
                    compare_cols=['cdr3_aa', 'cdr1_aa'],
                    out=tmp / 'out.tsv', levenshtein=1))

    def test_paired_mode_writes_a_classification_tsv(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            query_path = tmp / 'query.tsv'
            reference_path = tmp / 'reference.tsv'
            out_path = tmp / 'out.tsv'

            self._writeTsv(query_path, [
                {'cell_id': 'q1', 'locus': 'IGH', 'cdr3_aa': 'CARDYW',
                 'v_call': 'IGHV1', 'j_call': 'IGHJ1'},
                {'cell_id': 'q1', 'locus': 'IGL', 'cdr3_aa': 'CQQYNS',
                 'v_call': 'IGKV1', 'j_call': 'IGKJ1'},
            ])
            self._writeTsv(reference_path, [
                {'cell_id': 'ref1', 'locus': 'IGH', 'cdr3_aa': 'CARDYW',
                 'v_call': 'IGHV1', 'j_call': 'IGHJ1'},
                {'cell_id': 'ref1', 'locus': 'IGL', 'cdr3_aa': 'CQQYNS',
                 'v_call': 'IGKV1', 'j_call': 'IGKJ1'},
            ])

            with mock.patch('shutil.which', return_value='/usr/bin/compairr'), \
                 mock.patch('subprocess.run', side_effect=fakeCompairrExact):
                status = handleAnnotate(self._baseArgs(
                    query=query_path, reference=reference_path, out=out_path,
                    paired=True, compare_cols=None))

            self.assertEqual(status, 0)
            written = pandas.read_csv(out_path, sep='\t')
            self.assertEqual(written['classification'].tolist(), ['both'])
            self.assertEqual(written['both_hit_ids'].tolist(), ['ref1'])
