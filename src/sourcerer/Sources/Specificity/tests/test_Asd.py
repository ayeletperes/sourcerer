"""
Unit tests for the ASD specificity source
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

import pandas

# Sourcerer imports
from sourcerer.Schema import loadSchema
from sourcerer.Sources.Specificity import REGISTRY, getSpecificitySource
from sourcerer.Sources.Specificity.Asd import (COLLECTION, DRIVE_FOLDER_URL,
                                               AsdSource)
from sourcerer.Exceptions import AsdFetchError, AsdParseError
from sourcerer.Http import HttpClient
from sourcerer.Sources.Base import DataUnit, Query
from tests.FakeHttp import FakeResponse, FakeSession


def makeClient(handler):
    """Build an HttpClient with no politeness delay over a scripted session."""
    return HttpClient(delay=0, backoff=0, session=FakeSession(handler))


class FakeGdown(types.ModuleType):
    """
    A stand-in for the gdown module, driven by a callback.
    """

    def __init__(self, make_files):
        """
        Arguments:
          make_files (callable): called as make_files(Path) once
            download_folder is invoked, to populate the destination directory.
        """
        super().__init__('gdown')
        self._make_files = make_files
        self.calls = []

    def download_folder(self, url, output=None, quiet=True, use_cookies=False):
        """Record the call and populate output via the configured callback."""
        self.calls.append({'url': url, 'output': output})
        self._make_files(Path(output))


def writeParquet(path, records):
    """Write a small parquet file at path from a list of dict records."""
    pandas.DataFrame.from_records(records).to_parquet(path)


class TestAsdSource(unittest.TestCase):
    """
    Tests for the ASD specificity source
    """

    def test_registered(self):
        self.assertIs(REGISTRY['asd'], AsdSource)
        client = makeClient(lambda *a: FakeResponse(200))
        self.assertIsInstance(getSpecificitySource('asd', client), AsdSource)

    def test_packaged_schema_offers_a_dataset_field(self):
        schema = loadSchema('asd')
        collection = schema.getCollection(COLLECTION)
        self.assertIn('dataset', collection.field_names)
        self.assertIn('hiv', collection.getField('dataset').values)

    def test_search_units_returns_one_unit_carrying_the_filters(self):
        source = AsdSource(makeClient(lambda *a: FakeResponse(200)))

        units = source.searchUnits(
            Query(collection=COLLECTION, filters={'dataset': 'hiv'}))

        self.assertEqual(len(units), 1)
        self.assertEqual(units[0].unit_id, COLLECTION)
        self.assertEqual(units[0].url, DRIVE_FOLDER_URL)
        self.assertEqual(units[0].metadata, {'dataset': 'hiv'})

    def test_fetch_unit_downloads_the_folder_with_gdown(self):
        def make_files(path):
            path.mkdir(parents=True, exist_ok=True)
            (path / 'part-0.snappy.parquet').write_bytes(b'PAR1fake-parquet-bytes')

        fake_gdown = FakeGdown(make_files)
        source = AsdSource(makeClient(lambda *a: FakeResponse(200)))
        unit = DataUnit(unit_id=COLLECTION, collection=COLLECTION,
                        url=DRIVE_FOLDER_URL)

        with tempfile.TemporaryDirectory() as outdir:
            with mock.patch.dict(sys.modules, {'gdown': fake_gdown}):
                result = source.fetchUnit(unit, Path(outdir))

            self.assertEqual(len(fake_gdown.calls), 1)
            self.assertEqual(fake_gdown.calls[0]['url'], DRIVE_FOLDER_URL)
            self.assertTrue(result.path.is_dir())
            self.assertTrue(any(result.path.glob('*.parquet')))
            self.assertEqual(len(result.sha256), 64)
            self.assertGreater(result.size_bytes, 0)

    def test_fetch_unit_skips_a_directory_that_already_has_parquet_files(self):
        source = AsdSource(makeClient(lambda *a: FakeResponse(200)))
        unit = DataUnit(unit_id=COLLECTION, collection=COLLECTION,
                        url=DRIVE_FOLDER_URL)

        def unexpected_call(*args, **kwargs):
            raise AssertionError('gdown should not be called when resuming')

        with tempfile.TemporaryDirectory() as outdir:
            dest = Path(outdir) / COLLECTION / COLLECTION
            dest.mkdir(parents=True)
            (dest / 'part-0.snappy.parquet').write_bytes(b'PAR1already-here')

            with mock.patch.dict(sys.modules,
                                 {'gdown': types.SimpleNamespace(
                                     download_folder=unexpected_call)}):
                result = source.fetchUnit(unit, Path(outdir), resume=True)

            self.assertEqual(result.path, dest)

    def test_fetch_unit_raises_when_gdown_is_not_installed(self):
        source = AsdSource(makeClient(lambda *a: FakeResponse(200)))
        unit = DataUnit(unit_id=COLLECTION, collection=COLLECTION,
                        url=DRIVE_FOLDER_URL)

        with tempfile.TemporaryDirectory() as outdir:
            # Setting a module to None in sys.modules makes the next `import
            # gdown` raise ImportError, without requiring gdown to actually be
            # absent from the environment running the tests.
            with mock.patch.dict(sys.modules, {'gdown': None}):
                with self.assertRaises(AsdFetchError):
                    source.fetchUnit(unit, Path(outdir))

    def test_fetch_unit_raises_when_the_folder_has_no_parquet_files(self):
        fake_gdown = FakeGdown(lambda path: path.mkdir(parents=True, exist_ok=True))
        source = AsdSource(makeClient(lambda *a: FakeResponse(200)))
        unit = DataUnit(unit_id=COLLECTION, collection=COLLECTION,
                        url=DRIVE_FOLDER_URL)

        with tempfile.TemporaryDirectory() as outdir:
            with mock.patch.dict(sys.modules, {'gdown': fake_gdown}):
                with self.assertRaises(AsdFetchError):
                    source.fetchUnit(unit, Path(outdir))

    def test_read_unit_reads_every_partition_as_a_chunk(self):
        source = AsdSource(makeClient(lambda *a: FakeResponse(200)))
        unit = DataUnit(unit_id=COLLECTION, collection=COLLECTION,
                        url=DRIVE_FOLDER_URL)

        with tempfile.TemporaryDirectory() as outdir:
            path = Path(outdir)
            writeParquet(path / 'part-0.parquet',
                        [{'dataset': 'hiv', 'heavy_sequence': 'QVQ'}])
            writeParquet(path / 'part-1.parquet',
                        [{'dataset': 'covid-19', 'heavy_sequence': 'EVQ'}])

            metadata, chunks = source.readUnit(path, unit)
            frame = pandas.concat(list(chunks), ignore_index=True)

        self.assertEqual(metadata, {})
        self.assertEqual(sorted(frame['dataset'].tolist()), ['covid-19', 'hiv'])

    def test_read_unit_raises_when_no_parquet_files_are_present(self):
        source = AsdSource(makeClient(lambda *a: FakeResponse(200)))
        unit = DataUnit(unit_id=COLLECTION, collection=COLLECTION,
                        url=DRIVE_FOLDER_URL)

        with tempfile.TemporaryDirectory() as outdir:
            with self.assertRaises(AsdParseError):
                source.readUnit(Path(outdir), unit)

    def test_normalize_chunk_serializes_nested_metadata_to_json(self):
        source = AsdSource(makeClient(lambda *a: FakeResponse(200)))
        unit = DataUnit(unit_id=COLLECTION, collection=COLLECTION,
                        url=DRIVE_FOLDER_URL, metadata={})
        chunk = pandas.DataFrame([
            {'dataset': 'hiv', 'metadata': {'target_pdb': '5ovw'}},
        ])
        report = {'rows_in': 0, 'rows_out': 0}

        frame = source.normalizeChunk({}, chunk, unit, 0, report)

        self.assertEqual(frame['metadata'].iloc[0], '{"target_pdb": "5ovw"}')
        self.assertEqual(frame['sourcerer_source'].tolist(), ['asd'])
        self.assertEqual(frame['sourcerer_collection'].tolist(), [COLLECTION])
        self.assertEqual(len(frame['sourcerer_row_hash'].iloc[0]), 12)
        self.assertEqual(report['rows_in'], 1)
        self.assertEqual(report['rows_out'], 1)

    def test_normalize_chunk_applies_the_dataset_filter(self):
        source = AsdSource(makeClient(lambda *a: FakeResponse(200)))
        unit = DataUnit(unit_id=COLLECTION, collection=COLLECTION,
                        url=DRIVE_FOLDER_URL, metadata={'dataset': 'hiv'})
        chunk = pandas.DataFrame([
            {'dataset': 'hiv', 'metadata': None},
            {'dataset': 'covid-19', 'metadata': None},
        ])
        report = {'rows_in': 0, 'rows_out': 0}

        frame = source.normalizeChunk({}, chunk, unit, 0, report)

        self.assertEqual(frame['dataset'].tolist(), ['hiv'])
        self.assertEqual(report['rows_in'], 2)
        self.assertEqual(report['rows_out'], 1)

    def test_normalize_chunk_keeps_everything_for_the_wildcard(self):
        source = AsdSource(makeClient(lambda *a: FakeResponse(200)))
        unit = DataUnit(unit_id=COLLECTION, collection=COLLECTION,
                        url=DRIVE_FOLDER_URL, metadata={'dataset': '*'})
        chunk = pandas.DataFrame([
            {'dataset': 'hiv', 'metadata': None},
            {'dataset': 'covid-19', 'metadata': None},
        ])

        frame = source.normalizeChunk({}, chunk, unit, 0,
                                      {'rows_in': 0, 'rows_out': 0})

        self.assertEqual(len(frame), 2)


if __name__ == '__main__':
    unittest.main()
