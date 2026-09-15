"""
Unit tests for the download provenance record
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import tempfile
import unittest
from pathlib import Path

import yaml

# Sourcerer imports
from sourcerer import Provenance
from sourcerer.Exceptions import SourcererError
from sourcerer.Sources.Base import DataUnit, DownloadResult


def makeUnit(unit_id='study/a.csv.gz'):
    """Build a DataUnit and a matching DownloadResult."""
    unit = DataUnit(unit_id=unit_id, collection='paired',
                    url='https://example.org/%s' % unit_id, n_sequences=7)
    result = DownloadResult(unit=unit, path=Path('/out/raw/paired') / unit_id,
                            sha256='a' * 64, size_bytes=11)

    return unit, result


class TestUnitRecord(unittest.TestCase):
    """
    Tests for describing one downloaded unit
    """

    def test_paths_are_recorded_relative_to_the_output_root(self):
        """
        Moving or renaming the download directory must not invalidate the record.
        """
        unit, result = makeUnit()
        record = Provenance.buildUnitRecord(
            unit, result, Path('/out'), {'fasta': Path('/out/fasta/a.fasta')})

        self.assertEqual(record['raw'], 'raw/paired/study/a.csv.gz')
        self.assertEqual(record['outputs'], {'fasta': 'fasta/a.fasta'})

    def test_a_path_outside_the_root_is_kept_whole(self):
        """A path that cannot be relativized is recorded as it is, not mangled."""
        unit, result = makeUnit()
        record = Provenance.buildUnitRecord(unit, result, Path('/elsewhere'))

        self.assertEqual(record['raw'], '/out/raw/paired/study/a.csv.gz')

    def test_the_digest_is_recorded(self):
        """
        The digest is the point of the file.

        It is computed during every download and was previously only logged, so
        nothing on disk recorded what the raw mirror was verified against.
        """
        unit, result = makeUnit()
        record = Provenance.buildUnitRecord(unit, result, Path('/out'))

        self.assertEqual(record['sha256'], 'a' * 64)
        self.assertEqual(record['size_bytes'], 11)
        self.assertEqual(record['n_sequences'], 7)


class TestMerge(unittest.TestCase):
    """
    Tests for accumulating across several downloads
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.out = Path(self.tmp.name)

    def write(self, units, formats, filters=None):
        return Provenance.writeDownloadMetadata(
            self.out, 'oas', 'paired', filters or {}, None, formats, units)

    def read(self):
        return yaml.safe_load((self.out / Provenance.DOWNLOAD_METADATA).read_text())

    def test_a_second_download_does_not_orphan_the_first(self):
        """
        Assembling a dataset over several downloads is the normal case.

        Rewriting the file from one run's units alone would drop every unit
        fetched earlier, leaving files on disk that nothing describes.
        """
        first, first_result = makeUnit('study/a.csv.gz')
        second, second_result = makeUnit('study/b.csv.gz')

        self.write([Provenance.buildUnitRecord(first, first_result, self.out)],
                   ['raw'])
        self.write([Provenance.buildUnitRecord(second, second_result, self.out)],
                   ['raw'])

        record = self.read()
        self.assertEqual([x['unit_id'] for x in record['units']],
                         ['study/a.csv.gz', 'study/b.csv.gz'])
        self.assertEqual(len(record['runs']), 2)

    def test_converting_again_accumulates_outputs(self):
        """
        A unit converted to a second format keeps the first one's entry.

        The earlier output is still on disk, so unrecording it would describe a
        directory that does not match reality.
        """
        unit, result = makeUnit()

        self.write([Provenance.buildUnitRecord(unit, result, self.out,
                                               {'airr': self.out / 'airr/a.tsv'})],
                   ['airr'])
        self.write([Provenance.buildUnitRecord(unit, result, self.out,
                                               {'fasta': self.out / 'fasta/a.fasta'})],
                   ['fasta'])

        record = self.read()
        self.assertEqual(len(record['units']), 1)
        self.assertEqual(record['units'][0]['outputs'],
                         {'airr': 'airr/a.tsv', 'fasta': 'fasta/a.fasta'})

    def test_every_run_is_recorded(self):
        """Each run appends its own entry, so the history is not overwritten."""
        unit, result = makeUnit()
        self.write([Provenance.buildUnitRecord(unit, result, self.out)],
                   ['raw'], filters={'Species': 'human'})
        self.write([Provenance.buildUnitRecord(unit, result, self.out)],
                   ['fasta'], filters={'Species': 'mouse'})

        runs = self.read()['runs']
        self.assertEqual([x['filters'] for x in runs],
                         [{'Species': 'human'}, {'Species': 'mouse'}])
        self.assertEqual([x['formats'] for x in runs], [['raw'], ['fasta']])

    def test_header_keys_come_first(self):
        """
        The identifying keys stay at the top however the file was last written.

        A reader opening the file should see what it is before a unit list that
        may run to hundreds of entries.
        """
        unit, result = makeUnit()
        self.write([Provenance.buildUnitRecord(unit, result, self.out)], ['raw'])
        self.write([Provenance.buildUnitRecord(unit, result, self.out)], ['raw'])

        text = (self.out / Provenance.DOWNLOAD_METADATA).read_text()
        self.assertTrue(text.startswith('sourcerer_metadata_version:'))

    def test_a_foreign_file_is_not_overwritten(self):
        """
        A file sourcerer did not write is never rewritten.

        The name is a plausible one for a user to have chosen themselves, and
        clobbering hand-written provenance would destroy work that cannot be
        regenerated.
        """
        path = self.out / Provenance.DOWNLOAD_METADATA
        path.write_text('notes: downloaded these by hand\n')

        unit, result = makeUnit()
        with self.assertRaises(SourcererError):
            self.write([Provenance.buildUnitRecord(unit, result, self.out)],
                       ['raw'])

        self.assertIn('by hand', path.read_text())


class TestNaming(unittest.TestCase):
    """
    Tests for what the file is called
    """

    def test_the_name_claims_no_pipeline(self):
        """
        Raw source files are not airrflow input in either mode.

        Naming this file after a samplesheet or after airrflow would invite
        feeding a directory of .csv.gz to a pipeline that cannot read it.
        """
        self.assertNotIn('airrflow', Provenance.DOWNLOAD_METADATA)
        self.assertNotIn('samplesheet', Provenance.DOWNLOAD_METADATA)
        self.assertTrue(Provenance.DOWNLOAD_METADATA.endswith('.yml'))


if __name__ == '__main__':
    unittest.main()
