"""
Unit tests for drift detection
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import json
import os
import unittest

# Sourcerer imports
from sourcerer import Drift
from sourcerer.Drift import Finding, Snapshot
from sourcerer.Exceptions import ProbeIncompleteError
from sourcerer.Schema import Collection, Field, SourceSchema
from sourcerer.Sources import Oas

test_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(test_path, 'data')


def makeSchema(paired_fields=None, unpaired_fields=None):
    """Build a small two collection schema from field tuples."""
    def build(fields):
        return tuple(Field(name=n, values=tuple(v)) for n, v in fields)

    return SourceSchema(
        source='oas',
        collections={
            'paired': Collection(name='paired', fields=build(
                paired_fields or [('Species', ['human'])])),
            'unpaired': Collection(name='unpaired', fields=build(
                unpaired_fields or [('Species', ['human'])]))})


class TestSeverity(unittest.TestCase):
    """
    Tests for the ordered severity threshold
    """

    def test_overall_level_is_the_maximum(self):
        """A run's level is the worst finding, not the last or the most common."""
        findings = [Finding('additive', 'unit', 'x'),
                    Finding('structural', 'field', 'y'),
                    Finding('additive', 'unit', 'z'),
                    Finding('anomaly', 'partial-key', 'w')]

        self.assertEqual(Drift.overallLevel(findings), 'structural')
        self.assertEqual(Drift.overallLevel([]), 'none')

    def test_exit_code_reflects_the_threshold(self):
        findings = [Finding('removed', 'value', 'x')]

        self.assertEqual(Drift.exitCode(findings, 'never'), 0)
        self.assertEqual(Drift.exitCode(findings, 'structural'), 0)
        self.assertEqual(Drift.exitCode(findings, 'removed'), 25)
        self.assertEqual(Drift.exitCode(findings, 'additive'), 25)
        self.assertEqual(Drift.exitCode([], 'additive'), 0)

    def test_every_level_has_a_distinct_ascending_code(self):
        codes = [Drift.LEVEL_CODES[x] for x in Drift.LEVELS]

        self.assertEqual(codes, sorted(set(codes)))


class TestSchemaComparison(unittest.TestCase):
    """
    Tests for form schema drift classification
    """

    def test_a_new_field_is_structural(self):
        old = makeSchema(unpaired_fields=[('Species', ['human'])])
        new = makeSchema(unpaired_fields=[('Species', ['human']),
                                          ('Primer', ['p1'])])
        findings = Drift.compareSchemas(old, new)

        self.assertEqual([x.level for x in findings], ['structural'])
        self.assertIn("'Primer' appeared", findings[0].message)

    def test_a_lost_field_is_structural(self):
        old = makeSchema(paired_fields=[('Species', ['human']),
                                        ('Isotype', ['IGHG'])])
        new = makeSchema(paired_fields=[('Species', ['human'])])
        findings = Drift.compareSchemas(old, new)

        self.assertEqual([x.level for x in findings], ['structural'])
        self.assertIn("'Isotype' disappeared", findings[0].message)

    def test_vocabulary_gain_and_loss_split_levels(self):
        old = makeSchema(paired_fields=[('Disease', ['None', 'HIV'])])
        new = makeSchema(paired_fields=[('Disease', ['None', 'Dengue'])])
        findings = Drift.compareSchemas(old, new)
        levels = sorted(x.level for x in findings)

        self.assertEqual(levels, ['additive', 'removed'])

    def test_a_parse_contract_change_is_structural(self):
        old = makeSchema()
        new = SourceSchema(source='oas', collections=dict(old.collections),
                           parse_contracts={'count_regex': 'changed'})
        findings = Drift.compareSchemas(old, new)

        self.assertTrue(findings)
        self.assertTrue(all(x.level == 'structural' for x in findings))


class TestCatalogComparison(unittest.TestCase):
    """
    Tests for catalog drift classification
    """

    def test_new_units_are_additive_and_grouped_by_study(self):
        old = [{'unit_id': 'A_2020/csv/x.csv.gz', 'study': 'A_2020'}]
        new = old + [{'unit_id': 'B_2026/csv/y.csv.gz', 'study': 'B_2026'},
                     {'unit_id': 'B_2026/csv/z.csv.gz', 'study': 'B_2026'}]
        findings = Drift.compareCatalogs('paired', old, new)

        self.assertEqual(len(findings), 1)
        self.assertEqual(findings[0].level, 'additive')
        self.assertIn('2 new unit(s) in B_2026', findings[0].message)

    def test_lost_units_are_removed(self):
        old = [{'unit_id': 'A_2020/csv/x.csv.gz', 'study': 'A_2020'}]
        findings = Drift.compareCatalogs('paired', old, [])

        self.assertEqual([x.level for x in findings], ['removed'])


class TestContractComparison(unittest.TestCase):
    """
    Tests for data contract drift classification
    """

    def probe(self, **overrides):
        entry = {'unit_id': 'A_2020/csv/x.csv.gz', 'dir_segment': 'csv',
                 'gzip_members': 2, 'n_columns': 180}
        entry.update(overrides)
        return entry

    def wrap(self, segments, probes):
        return {'collections': {'paired': {
            'path_layouts': {'observed_dir_segments': segments,
                             'observed_filename_patterns': {}},
            'probe_units': probes}}}

    def test_a_novel_path_layout_is_additive(self):
        """Paths are opaque, so a layout the code has never seen costs nothing."""
        old = self.wrap({'csv': 10}, [self.probe()])
        new = self.wrap({'csv': 10, 'csv_paired': 5}, [self.probe()])
        findings = Drift.compareContracts(old, new)

        self.assertEqual([x.level for x in findings], ['additive'])
        self.assertIn("csv_paired", findings[0].message)

    def test_a_fact_deviation_is_structural(self):
        old = self.wrap({'csv': 10}, [self.probe(n_columns=178)])
        new = self.wrap({'csv': 10}, [self.probe(n_columns=180)])
        findings = Drift.compareContracts(old, new)

        self.assertEqual([x.level for x in findings], ['structural'])
        self.assertIn('n_columns', findings[0].message)

    def test_a_replaced_pin_is_removed_not_structural(self):
        """
        Two different units may legitimately differ in row level facts, so a
        replaced pin must not read as a format change.
        """
        old = self.wrap({'csv': 10}, [self.probe(n_columns=178)])
        new = self.wrap({'csv': 10},
                        [self.probe(unit_id='B_2026/csv/y.csv.gz')])
        findings = Drift.compareContracts(old, new)

        self.assertEqual([x.level for x in findings], ['removed'])


class TestFingerprintComparison(unittest.TestCase):
    """
    Tests for catalog fingerprint drift classification
    """

    def test_growth_is_additive_and_shrinkage_is_removed(self):
        old = {'n_units': 10, 'key_counts': {}, 'value_types': {}}
        grown = {'n_units': 12, 'key_counts': {}, 'value_types': {}}
        shrunk = {'n_units': 8, 'key_counts': {}, 'value_types': {}}

        self.assertEqual([x.level for x in Drift.compareFingerprints(old, grown)],
                         ['additive'])
        self.assertEqual([x.level for x in Drift.compareFingerprints(old, shrunk)],
                         ['removed'])

    def test_a_key_on_every_unit_is_structural_on_a_few_additive(self):
        old = {'n_units': 10, 'key_counts': {'Run': 10}, 'value_types': {}}
        every = {'n_units': 10, 'key_counts': {'Run': 10, 'New': 10},
                 'value_types': {}}
        a_few = {'n_units': 10, 'key_counts': {'Run': 10, 'New': 2},
                 'value_types': {}}

        self.assertEqual([x.level for x in Drift.compareFingerprints(old, every)],
                         ['structural'])
        self.assertEqual([x.level for x in Drift.compareFingerprints(old, a_few)],
                         ['additive'])

    def test_a_value_type_change_is_structural(self):
        old = {'n_units': 10, 'key_counts': {'Total': 10},
               'value_types': {'Total': ['int']}}
        new = {'n_units': 10, 'key_counts': {'Total': 10},
               'value_types': {'Total': ['str']}}
        findings = Drift.compareFingerprints(old, new)

        self.assertEqual([x.level for x in findings], ['structural'])


class TestNewSnapshotChecks(unittest.TestCase):
    """
    Tests for the checks that inspect the fresh snapshot alone
    """

    def test_a_partial_key_is_an_anomaly(self):
        """The Organism-on-one-unit case: neither absent nor universal."""
        snapshot = Snapshot(fingerprint={
            'n_units': 15631, 'key_counts': {'Species': 15631, 'Organism': 1}})
        findings = Drift.findAnomalies(snapshot)

        self.assertEqual([x.level for x in findings], ['anomaly'])
        self.assertIn("'Organism' is present on 1 of 15631", findings[0].message)

    def test_an_unresolvable_unit_id_is_structural(self):
        """
        Paths are opaque, but they must still be preservable as local relative
        paths; one that escapes upward or collides breaks the mirror.
        """
        snapshot = Snapshot(catalogs={'paired': [
            {'unit_id': '../escape.csv.gz'},
            {'unit_id': 'A_2020/csv/x.csv.gz'},
            {'unit_id': 'A_2020/csv//x.csv.gz'}]})
        findings = Drift.checkPathSafety(snapshot)
        levels = [x.level for x in findings]

        self.assertEqual(levels, ['structural', 'structural'])
        self.assertIn('escape', findings[0].message)
        self.assertIn('collide', findings[1].message)

    def test_a_dead_url_is_structural(self):
        class DeadClient:
            def probeAlive(self, url):
                return False

        snapshot = Snapshot(catalogs={'paired': [
            {'unit_id': 'A_2020/csv/x.csv.gz',
             'url': 'https://example.org/x.csv.gz'}]})
        findings = Drift.probeUrls(snapshot, DeadClient())

        self.assertEqual([x.level for x in findings], ['structural'])
        self.assertIn('does not resolve', findings[0].message)


class TestHarvestFailureIsNotDrift(unittest.TestCase):
    """
    Tests for the probe-incomplete condition
    """

    def test_an_exhausted_probe_propagates_as_a_failure(self):
        """
        A probe that cannot complete must raise rather than classify: a slow or
        truncated response is a harvest failure, and reporting it as drift
        would open a pull request claiming the format changed when it did not.
        """
        class CappedClient:
            def readRanges(self, url, is_complete, **kwargs):
                raise ProbeIncompleteError('cap reached')

        source = Oas.OasSource(CappedClient())
        rows = [{'unit_id': 'A_2020/csv/x.csv.gz', 'dir_segment': 'csv',
                 'url': 'https://example.org/x.csv.gz',
                 'n_unique_sequences': '5'}]

        with self.assertRaises(ProbeIncompleteError):
            source.harvestContracts({'paired': rows}, existing=None)


class TestStaleSnapshotFixture(unittest.TestCase):
    """
    Tests against the deliberately stale schema_prev snapshot

    This is the regression the project exists to catch: the predecessor tool
    hardcoded the paired field list, OAS moved on, and nothing noticed. The
    stale fixture records that old state; the drift check must classify the
    difference against the current packaged snapshot as structural on both the
    form axis and the file-format axis.
    """

    @classmethod
    def setUpClass(cls):
        cls.old = Drift.loadSnapshotDir(
            'oas', path=os.path.join(data_path, 'schema_prev'))
        cls.new = Drift.loadSnapshotDir('oas')
        cls.findings = Drift.checkDrift(cls.old, cls.new, client=None)

    def test_overall_level_is_structural(self):
        self.assertEqual(Drift.overallLevel(self.findings), 'structural')

    def structural(self, category, collection):
        return [x for x in self.findings
                if x.level == 'structural' and x.category == category
                and x.collection == collection]

    def test_form_axis_new_and_lost_fields_are_structural(self):
        paired = {x.message for x in self.structural('field', 'paired')}
        unpaired = {x.message for x in self.structural('field', 'unpaired')}

        self.assertTrue(any("'Primer' appeared" in x for x in unpaired))
        self.assertTrue(any("'BSource' appeared" in x for x in paired))
        self.assertTrue(any("'Isotype' disappeared" in x for x in paired))

    def test_file_format_axis_is_structural(self):
        contract = self.structural('data-contract', 'paired')

        self.assertTrue(any('n_columns' in x.message for x in contract))

    def test_the_novel_directory_layout_is_additive(self):
        layout = [x for x in self.findings
                  if x.category == 'layout' and 'csv_paired' in x.message]

        self.assertTrue(layout)
        self.assertTrue(all(x.level == 'additive' for x in layout))


class TestReporting(unittest.TestCase):
    """
    Tests for report assembly and rendering
    """

    def setUp(self):
        self.findings = [Finding('structural', 'field', "field 'X' appeared",
                                 'paired'),
                         Finding('additive', 'unit', '3 new unit(s) in A_2026')]

    def test_report_is_json_serializable_and_complete(self):
        report = Drift.buildReport('oas', 'git:HEAD', self.findings)
        parsed = json.loads(json.dumps(report))

        self.assertEqual(parsed['overall_level'], 'structural')
        self.assertEqual(len(parsed['findings']), 2)
        self.assertEqual(parsed['findings'][0]['collection'], 'paired')

    def test_markdown_orders_worst_first(self):
        report = Drift.buildReport('oas', 'git:HEAD', self.findings)
        text = Drift.renderMarkdown(report)

        self.assertIn('**Overall level: structural**', text)
        self.assertLess(text.index('## structural'), text.index('## additive'))

    def test_markdown_says_so_when_quiet(self):
        report = Drift.buildReport('oas', 'git:HEAD', [])
        text = Drift.renderMarkdown(report)

        self.assertIn('No drift detected', text)


if __name__ == '__main__':
    unittest.main()
