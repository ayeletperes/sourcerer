"""
Unit tests for airrflow samplesheet generation
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import csv
import shutil
import tempfile
import unittest
from pathlib import Path

# Sourcerer imports
from sourcerer.Airrflow import (
    SAMPLESHEET_COLUMNS,
    buildSamplesheet,
    loadSamplesheet,
    targetLocus,
)
from sourcerer.Exceptions import SourcererError
from sourcerer.Sources.Base import DataUnit


def makeUnit(unit_id, **metadata):
    """Build a paired data unit carrying the given source metadata."""
    return DataUnit(unit_id=unit_id, collection='paired',
                    url='https://example.invalid/%s' % unit_id,
                    metadata=metadata)


def readRows(path):
    """Read a samplesheet back as a list of dicts."""
    with open(path, newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


class TestTargetLocus(unittest.TestCase):
    """
    Tests for the receptor class collapse
    """

    def test_collapses_ig_loci(self):
        """Immunoglobulin loci collapse to IG for pcr_target_locus."""
        self.assertEqual(targetLocus(['IGH', 'IGK']), 'IG')

    def test_rejects_mixed_receptor_classes(self):
        """
        A unit mixing IG and TR has no single pcr_target_locus.

        Picking either one would silently mislabel half the data, so this is an
        error rather than a choice.
        """
        with self.assertRaises(ValueError):
            targetLocus(['IGH', 'TRB'])


class TestSamplesheetMerge(unittest.TestCase):
    """
    Tests for accumulating a samplesheet across several downloads
    """

    def setUp(self):
        self.outdir = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.outdir, ignore_errors=True)
        self.sheet = self.outdir / 'samplesheet_airrflow_fasta.tsv'

    def write(self, units, loci=None):
        """Run the builder over the given units, writing to the shared sheet."""
        entries = [(unit, self.outdir / 'fasta' / ('%s.fasta' % unit.unit_id))
                   for unit in units]
        loci = loci or {unit.unit_id: {'IGH', 'IGK'} for unit in units}

        return buildSamplesheet(entries, self.sheet, 'paired',
                                root=self.outdir, loci=loci)

    def test_second_download_appends_rather_than_replacing(self):
        """
        A second download into the same outdir keeps the first one's rows.

        Assembling a dataset over several filtered downloads is normal usage.
        Rewriting the sheet from only the current run would leave the earlier
        run's converted files on disk with nothing describing them.
        """
        self.write([makeUnit('A_2020/csv/a.csv.gz', Species='mouse_C57BL/6')])
        self.write([makeUnit('B_2024/csv_paired/b.csv.gz', Species='human')])

        rows = readRows(self.sheet)
        self.assertEqual([x['sample_name'] for x in rows],
                         ['A_2020/csv/a.csv.gz', 'B_2024/csv_paired/b.csv.gz'])
        self.assertEqual([x['sample_id'] for x in rows], ['ssr_1', 'ssr_2'])
        self.assertEqual([x['species'] for x in rows], ['mouse_c57bl/6', 'human'])

    def test_existing_sample_ids_never_renumber(self):
        """
        A unit keeps its sample_id when new units are added around it.

        Users reference sample_id in downstream airrflow configuration, so
        renumbering an existing sample when the sheet grows would silently
        repoint that configuration at different data.
        """
        first = makeUnit('Z_2020/csv/z.csv.gz')
        self.write([first])
        self.write([makeUnit('A_2019/csv/a.csv.gz'), first])

        rows = {x['sample_name']: x['sample_id'] for x in readRows(self.sheet)}
        self.assertEqual(rows['Z_2020/csv/z.csv.gz'], 'ssr_1')
        self.assertEqual(rows['A_2019/csv/a.csv.gz'], 'ssr_2')

    def test_rerunning_a_unit_does_not_duplicate_it(self):
        """Re-downloading a unit already described updates its row in place."""
        unit = makeUnit('A_2020/csv/a.csv.gz', Species='human')
        self.write([unit])
        self.write([unit])

        self.assertEqual(len(readRows(self.sheet)), 1)

    def test_rerunning_fills_blanks_but_keeps_edits(self):
        """
        A repeat run fills empty fields without overwriting existing values.

        sex is not derivable from OAS at all, so a hand-edited value is the only
        way it is ever populated and must survive. Conversely tissue starts empty
        and should pick up a value once detail page enrichment supplies one.
        """
        unit_id = 'A_2020/csv/a.csv.gz'
        self.write([makeUnit(unit_id)])

        rows = readRows(self.sheet)
        rows[0]['sex'] = 'female'
        with open(self.sheet, 'w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(SAMPLESHEET_COLUMNS),
                                    delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)

        self.write([makeUnit(unit_id, BSource='PBMC')])

        merged = readRows(self.sheet)[0]
        self.assertEqual(merged['sex'], 'female')
        self.assertEqual(merged['tissue'], 'PBMC')

    def test_refuses_to_overwrite_a_foreign_file(self):
        """
        A file that is not a sourcerer samplesheet is never rewritten.

        Merging depends on the column layout, and a hand-built samplesheet that
        happens to occupy the expected path represents work that cannot be
        regenerated.
        """
        self.sheet.write_text('sample\tfile\nfoo\tbar.fasta\n')

        with self.assertRaises(SourcererError):
            self.write([makeUnit('A_2020/csv/a.csv.gz')])

    def test_missing_file_loads_as_empty(self):
        """A samplesheet that does not exist yet reads as no rows."""
        self.assertEqual(loadSamplesheet(self.outdir / 'absent.tsv'), [])

    def test_subject_fallback_names_the_assigned_sample_id(self):
        """
        A unit with no subject falls back to a name derived from its sample_id.

        The fallback is applied after the merge assigns identifiers, so it has to
        agree with the id the row actually ended up with.
        """
        self.write([makeUnit('A_2020/csv/a.csv.gz')])
        self.write([makeUnit('B_2020/csv/b.csv.gz')])

        rows = readRows(self.sheet)
        self.assertEqual(rows[1]['sample_id'], 'ssr_2')
        self.assertEqual(rows[1]['subject_id'], 'ssr_2_subj')


if __name__ == '__main__':
    unittest.main()
