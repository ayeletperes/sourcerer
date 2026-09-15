"""
Live API checks for the reference sources

These contact IMGT and OGRDB for real, so they are skipped unless SOURCERER_LIVE
is set in the environment. The weekly check-apis workflow sets it and runs this
module on its own; a red run there is the early warning that an upstream API
changed shape before a user's download hits the same failure.

The endpoint constants come from the source modules, so this tests exactly what
production calls rather than a second copy of the URLs.
"""

# Info
__author__ = 'Ayelet Peres'

# Imports
import os
import unittest

# Sourcerer imports
from sourcerer.Http import HttpClient
from sourcerer.Sources.Imgt import (
    ImgtSource,
    buildQueryUrl,
    extractFasta,
    isValidResponse,
)
from sourcerer.Sources.Ogrdb import OgrdbSource

LIVE = os.environ.get('SOURCERER_LIVE')


@unittest.skipUnless(LIVE, 'set SOURCERER_LIVE=1 to contact IMGT and OGRDB')
class TestImgtLive(unittest.TestCase):
    """
    Live checks against IMGT/GENE-DB
    """

    def setUp(self):
        self.source = ImgtSource(client=HttpClient())

    def test_genelect_returns_a_germline_fasta(self):
        """A GENElect query returns a page with a real FASTA in it."""
        url = buildQueryUrl('human', '7.14', 'IGHD')
        html = self.source.client.get(url).text
        self.assertTrue(isValidResponse(html),
                        'GENElect no longer returns a second <pre> with a FASTA')
        self.assertIn('>', extractFasta(html, 'human'))

    def test_release_tag_is_readable(self):
        """The GENE-DB release tag is still published and non-empty."""
        self.assertTrue(self.source.fetchRelease(),
                        'the IMGT release tag could not be read')


@unittest.skipUnless(LIVE, 'set SOURCERER_LIVE=1 to contact IMGT and OGRDB')
class TestOgrdbLive(unittest.TestCase):
    """
    Live checks against OGRDB
    """

    def setUp(self):
        self.source = OgrdbSource(client=HttpClient())

    def test_harvest_schema_sees_the_consumed_loci(self):
        """species and sets resolve, and human still exposes IGH, IGK and IGL."""
        schema = self.source.harvestSchema()
        human = schema.getCollection('human').getField('locus')
        for locus in ('IGH', 'IGK', 'IGL'):
            self.assertIn(locus, human.values,
                          'OGRDB no longer lists %s for human' % locus)

    def test_set_resolves_to_a_downloadable_fasta(self):
        """A set resolves to a release whose FASTA download is non-empty."""
        from sourcerer.Sources.Base import Query

        units = self.source.searchUnits(
            Query(collection='human', filters={'locus': 'IGK'}))
        self.assertTrue(units, 'no OGRDB units resolved for human IGK')

        body = self.source.client.get(units[0].url).text
        self.assertIn('>', body, 'the OGRDB FASTA download was empty')


if __name__ == '__main__':
    unittest.main()
