"""
Unit tests for the CoV-AbDab specificity source
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import tempfile
import unittest
from pathlib import Path

import pandas

# Sourcerer imports
from sourcerer.Schema import loadSchema
from sourcerer.Sources.Specificity import REGISTRY, getSpecificitySource
from sourcerer.Sources.Specificity.CovAbDab import (COLLECTION, HOME_URL,
                                                     CovAbDabSource,
                                                     findCsvUrl)
from sourcerer.Exceptions import CovAbDabParseError
from sourcerer.Http import HttpClient
from sourcerer.Sources.Base import DataUnit, Query
from tests.FakeHttp import FakeResponse, FakeSession, rangeHandler


def makeClient(handler):
    """Build an HttpClient with no politeness delay over a scripted session."""
    return HttpClient(delay=0, backoff=0, session=FakeSession(handler))


HOME_HTML = """
<html><body>
<a href="/webapps/covabdab/static/downloads/CoV-AbDab_080224.csv">Download</a>
</body></html>
"""


class TestFindCsvUrl(unittest.TestCase):
    """
    Tests for scraping the current bulk CSV link off the homepage
    """

    def test_resolves_a_relative_download_link_to_an_absolute_url(self):
        client = makeClient(lambda *a: FakeResponse(200, HOME_HTML.encode()))

        url = findCsvUrl(client)

        self.assertEqual(
            url, 'https://opig.stats.ox.ac.uk/webapps/covabdab/static/'
                'downloads/CoV-AbDab_080224.csv')

    def test_raises_when_no_csv_link_is_present(self):
        client = makeClient(
            lambda *a: FakeResponse(200, b'<html><body>no links here</body></html>'))

        with self.assertRaises(CovAbDabParseError):
            findCsvUrl(client)


class TestCovAbDabSource(unittest.TestCase):
    """
    Tests for the CoV-AbDab specificity source
    """

    def test_registered(self):
        self.assertIs(REGISTRY['covabdab'], CovAbDabSource)
        client = makeClient(lambda *a: FakeResponse(200))
        self.assertIsInstance(getSpecificitySource('covabdab', client),
                              CovAbDabSource)

    def test_packaged_schema_loads_and_offers_the_antibodies_collection(self):
        schema = loadSchema('covabdab')
        self.assertEqual(schema.collection_names, (COLLECTION,))

    def test_search_units_scrapes_the_current_csv_url(self):
        client = makeClient(lambda *a: FakeResponse(200, HOME_HTML.encode()))
        source = CovAbDabSource(client)

        units = source.searchUnits(Query(collection=COLLECTION, filters={}))

        self.assertEqual(len(units), 1)
        self.assertEqual(units[0].collection, COLLECTION)
        self.assertEqual(units[0].unit_id, 'CoV-AbDab_080224.csv')
        self.assertTrue(units[0].url.startswith(HOME_URL.rstrip('/')) or
                        units[0].url.startswith('https://opig.stats.ox.ac.uk'))

    def test_fetch_unit_downloads_the_resolved_csv(self):
        """fetchUnit falls back to SourceBase's plain URL download."""
        body = b'a,b\n1,2\n'
        client = makeClient(rangeHandler(body))
        source = CovAbDabSource(client)
        unit = DataUnit(unit_id='CoV-AbDab_080224.csv', collection=COLLECTION,
                        url='https://opig.stats.ox.ac.uk/x.csv')

        with tempfile.TemporaryDirectory() as outdir:
            result = source.fetchUnit(unit, Path(outdir))
            self.assertEqual(result.path.read_bytes(), body)

    def test_read_unit_strips_the_bom_and_parses_the_csv(self):
        source = CovAbDabSource(makeClient(lambda *a: FakeResponse(200)))
        unit = DataUnit(unit_id='x.csv', collection=COLLECTION, url='https://x/x.csv')

        with tempfile.TemporaryDirectory() as outdir:
            path = Path(outdir) / 'x.csv'
            path.write_bytes('﻿Name,Ab or Nb\nCurtis_3548,Ab\n'.encode('utf-8-sig'))
            metadata, chunks = source.readUnit(path, unit)
            frame = pandas.concat(list(chunks), ignore_index=True)

        self.assertEqual(metadata, {})
        self.assertEqual(list(frame.columns), ['Name', 'Ab or Nb'])
        self.assertEqual(frame.to_dict('records'),
                         [{'Name': 'Curtis_3548', 'Ab or Nb': 'Ab'}])

    def test_normalize_chunk_adds_provenance_columns(self):
        source = CovAbDabSource(makeClient(lambda *a: FakeResponse(200)))
        unit = DataUnit(unit_id='x.csv', collection=COLLECTION, url='https://x/x.csv')
        chunk = pandas.DataFrame([{'Name': 'Curtis_3548', 'Ab or Nb': 'Ab'}])
        report = {'rows_in': 0, 'rows_out': 0}

        frame = source.normalizeChunk({}, chunk, unit, 0, report)

        self.assertEqual(frame['sourcerer_source'].tolist(), ['covabdab'])
        self.assertEqual(frame['sourcerer_collection'].tolist(), [COLLECTION])
        self.assertEqual(frame['sourcerer_unit_id'].tolist(), ['x.csv'])
        self.assertEqual(len(frame['sourcerer_row_hash'].iloc[0]), 12)
        self.assertEqual(report['rows_in'], 1)
        self.assertEqual(report['rows_out'], 1)

    def test_normalize_chunk_handles_an_empty_frame(self):
        source = CovAbDabSource(makeClient(lambda *a: FakeResponse(200)))
        unit = DataUnit(unit_id='x.csv', collection=COLLECTION, url='https://x/x.csv')
        chunk = pandas.DataFrame(columns=['Name', 'Ab or Nb'])

        frame = source.normalizeChunk({}, chunk, unit, 0,
                                      {'rows_in': 0, 'rows_out': 0})

        self.assertEqual(len(frame), 0)


if __name__ == '__main__':
    unittest.main()
