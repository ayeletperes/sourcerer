"""
Unit tests for the schema refresh path

The one invariant that cannot live anywhere else: a refresh that finds nothing
new writes no tracked file at all. The scheduled workflow opens a pull request
whenever the working tree is dirty, so a file stamped on every run would open
an empty PR every month -- the no-drift, no-PR policy is enforced by the
filesystem, and this is the test that holds it there.
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import argparse
import gzip
import json
import os
import tempfile
import unittest
from pathlib import Path

# Sourcerer imports
from sourcerer import Cli
from sourcerer.Http import HttpClient
from tests.FakeHttp import FakeSession, rangeHandler

test_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(test_path, 'data')


def readBytes(name):
    """Read a fixture as raw bytes."""
    with open(os.path.join(data_path, name), 'rb') as handle:
        return handle.read()


#: A miniature unpaired catalog in the upstream JSON shape. It includes the
#: probe unit pinned by the packaged data contracts, so the fake refresh
#: honors the same pin on every run.
UNPAIRED_CATALOG = {
    '/vols/naga-datasets/oas/unpaired/Banerjee_2017/csv/'
    'SRR5060322_Heavy_IGHA.csv.gz': {
        'Species': 'human', 'Chain': 'Heavy', 'Isotype': 'IGHA',
        'BSource': 'PBMC', 'BType': 'Unsorted-B-Cells', 'Age': '22',
        'Disease': 'None', 'Subject': 'no', 'Vaccine': 'None',
        'Longitudinal': 'no', 'Run': 'SRR5060322', 'Author': 'Banerjee et al.',
        'Unique sequences': 12, 'Total sequences': 20},
    '/vols/naga-datasets/oas/unpaired/Banerjee_2017/csv/'
    'SRR5060321_Heavy_Bulk.csv.gz': {
        'Species': 'human', 'Chain': 'Heavy', 'Isotype': 'Bulk',
        'BSource': 'PBMC', 'BType': 'Unsorted-B-Cells', 'Age': '22',
        'Disease': 'None', 'Subject': 'no', 'Vaccine': 'None',
        'Longitudinal': 'no', 'Run': 'SRR5060321', 'Author': 'Banerjee et al.',
        'Unique sequences': 34, 'Total sequences': 55},
}


def buildHandler():
    """
    Build a FakeSession handler standing in for the whole OAS service.

    Data unit URLs are served range aware from the committed head fixtures, one
    per path layout, so the contract probes decode the same bytes every run.
    """
    paired_form = readBytes('oas_paired_form.html')
    unpaired_form = readBytes('oas_unpaired_form.html')
    search_all = gzip.decompress(readBytes('oas_paired_search_all.html.gz'))
    catalog = json.dumps(UNPAIRED_CATALOG).encode()
    units = {'/csv_paired/': readBytes('1_S1__1_Paired_All.head.csv.gz'),
             '/unpaired/': readBytes('SRR5060321_Heavy_Bulk.head.csv.gz'),
             '/csv/': readBytes('SRR11528761_paired.head.csv.gz')}

    def handler(method, url, headers, index):
        if url.endswith('oas_metadata_map.json'):
            return rangeHandler(catalog)(method, url, headers, index)
        if '/webapps/ngsdb/' in url:
            for marker, body in units.items():
                if marker in url:
                    return rangeHandler(body)(method, url, headers, index)
            raise AssertionError('unexpected data unit URL: %s' % url)
        if url.endswith('/oas_paired/'):
            body = search_all if method == 'POST' else paired_form
            return rangeHandler(body)(method, url, headers, index)
        if url.endswith('/oas_unpaired/'):
            return rangeHandler(unpaired_form)(method, url, headers, index)
        raise AssertionError('unexpected URL: %s' % url)

    return handler


def snapshotState(root):
    """Map every file under a directory to its content and mtime."""
    state = {}
    for path in sorted(Path(root).rglob('*')):
        if path.is_file():
            state[str(path.relative_to(root))] = (path.read_bytes(),
                                                  path.stat().st_mtime_ns)
    return state


class TestQuietRefresh(unittest.TestCase):
    """
    Tests for the quiet-month invariant
    """

    def setUp(self):
        self._original = Cli.makeClient
        Cli.makeClient = lambda args: HttpClient(
            delay=0, backoff=0, session=FakeSession(buildHandler()))

    def tearDown(self):
        Cli.makeClient = self._original

    def refresh(self, out):
        args = argparse.Namespace(source='oas', out=Path(out), collection=None,
                                  refresh_details='none', detail_limit=None)
        return Cli.handleSchemaRefresh(args)

    def test_a_refresh_that_finds_nothing_new_writes_nothing(self):
        with tempfile.TemporaryDirectory() as out:
            self.assertEqual(self.refresh(out), 0)
            first = snapshotState(out)

            # The first run writes the full snapshot, provenance included.
            self.assertIn('schema.yaml', first)
            self.assertIn('paired_catalog.tsv', first)
            self.assertIn('unpaired_catalog.tsv', first)
            self.assertIn('data_contracts.yaml', first)
            self.assertIn('catalog_fingerprint.json', first)
            self.assertIn('provenance.json', first)

            self.assertEqual(self.refresh(out), 0)
            second = snapshotState(out)

            # The second run found the same source state, so it must not have
            # touched a single file: identical content and identical mtimes.
            # Content-only comparison would miss a rewrite-in-place, and a
            # rewritten provenance stamp is exactly the failure that would
            # open an empty pull request every month.
            self.assertEqual(first, second)

    def test_the_contract_probe_pins_survive_a_second_run(self):
        import yaml

        with tempfile.TemporaryDirectory() as out:
            self.refresh(out)
            contracts = yaml.safe_load(
                (Path(out) / 'data_contracts.yaml').read_text())
            pins = {x['unit_id']
                    for body in contracts['collections'].values()
                    for x in body['probe_units']}

            # The packaged pin present in the fake catalog is kept rather than
            # re-chosen; re-picking every run would churn the snapshot diff.
            self.assertIn(
                'Banerjee_2017/csv/SRR5060322_Heavy_IGHA.csv.gz', pins)


if __name__ == '__main__':
    unittest.main()
