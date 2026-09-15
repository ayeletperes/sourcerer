"""
Unit tests for the IMGT historical-release archive
"""

# Info
__author__ = 'Ayelet Peres'

# Imports
import unittest

# Sourcerer imports
from sourcerer import Genedb
from sourcerer.Reference import KIND_AA, KIND_CONSTANT, KIND_VDJ


class FakeResponse:
    """A stand-in exposing the .json() the contents listing needs."""

    def __init__(self, payload):
        self._payload = payload

    def json(self):
        return self._payload


class ListingClient:
    """A client whose only call is the contents API, returning canned entries."""

    def __init__(self, names):
        self.payload = [{'name': name, 'type': 'dir'} for name in names]

    def get(self, url):
        return FakeResponse(self.payload)


class TestReleaseKey(unittest.TestCase):
    """
    Tests for ordering release tags
    """

    def test_orders_by_build_then_revision(self):
        """A YYYYWW-R tag sorts by build then revision, as a date would."""
        self.assertLess(Genedb.releaseKey('202630-7'), Genedb.releaseKey('202631-7'))
        self.assertLess(Genedb.releaseKey('202631-1'), Genedb.releaseKey('202631-2'))

    def test_unparseable_tag_is_lowest(self):
        """A tag that is not a release does not crash the ordering."""
        self.assertEqual(Genedb.releaseKey('nonsense'), (0, 0))


class TestResolveRelease(unittest.TestCase):
    """
    Tests for finding a release in the archive
    """

    NAMES = ('2026-07-27_GENEDB_202630-7', '2026-08-03_GENEDB_202631-7',
             '2026-08-10_GENEDB_202632-7')

    def test_exact_match_is_preferred(self):
        """A release the archive holds is returned exactly, flagged exact."""
        client = ListingClient(self.NAMES)
        dirname, tag, exact = Genedb.resolveRelease(client, '202631-7')
        self.assertEqual(tag, '202631-7')
        self.assertEqual(dirname, '2026-08-03_GENEDB_202631-7')
        self.assertTrue(exact)

    def test_missing_release_falls_back_to_closest(self):
        """A release not archived resolves to the nearest, flagged inexact."""
        client = ListingClient(self.NAMES)
        _dir, tag, exact = Genedb.resolveRelease(client, '202629-7')
        self.assertEqual(tag, '202630-7')
        self.assertFalse(exact)


class TestBulkUrl(unittest.TestCase):
    """
    Tests for building the raw archive URL
    """

    def test_encodes_the_plus_in_the_filename(self):
        """The + in the functionality suffix is percent-encoded for the raw host."""
        url = Genedb.bulkUrl('2026-08-03_GENEDB_202631-7', 'nt')
        self.assertIn('%2B', url)
        self.assertNotIn('+', url)


class TestSelectChains(unittest.TestCase):
    """
    Tests for filtering a bulk FASTA into one species' chains
    """

    #: A minimal bulk file: a human V, a human J, a long delta constant that must
    #: file under IGHC, a short delta D that must stay a D segment, a mouse V from
    #: one strain and the same allele from a second strain, and another species.
    BULK = (
        '>X02897|IGHV1-2*02|Homo sapiens|F\nAC.GT\n'
        '>J00256|IGHJ4*02|Homo sapiens|F\nACTACTGG\n'
        '>D78345|IGHD*01|Homo sapiens|F\n%s\n'
        '>Z12345|IGHD3-10*01|Homo sapiens|F\nGTATTAC\n'
        '>M11\t|IGHV1-1*01|Mus musculus_C57BL/6|F\nAAAA\n'
        '>M22|IGHV1-1*01|Mus musculus_BALB/c|F\nCCCC\n'
        '>K99|IGHV9*01|Bos taurus|F\nGGGG\n'
    ) % ('A' * 120)

    def test_human_v_j_and_delta_constant_are_bucketed(self):
        """V and J file under their chain; a long IGHD files under IGHC."""
        chains = Genedb.selectChains(self.BULK, 'human', is_aa=False)
        self.assertIn(('IGHV', KIND_VDJ), chains)
        self.assertIn(('IGHJ', KIND_VDJ), chains)
        self.assertIn(('IGHC', KIND_CONSTANT), chains)
        self.assertIn(('IGHD', KIND_VDJ), chains)  # the short one stays a D

    def test_species_filter_excludes_other_species(self):
        """Only the requested species is kept; another is dropped."""
        chains = Genedb.selectChains(self.BULK, 'human', is_aa=False)
        headers = [h for records in chains.values() for h, _ in records]
        self.assertFalse(any('Bos taurus' in h for h in headers))

    def test_mouse_matches_every_strain_but_keeps_one_per_allele(self):
        """A gene shared across strains is kept once, from its first occurrence."""
        chains = Genedb.selectChains(self.BULK, 'mouse', is_aa=False)
        records = chains[('IGHV', KIND_VDJ)]
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0][1], 'AAAA')  # C57BL/6 came first

    def test_species_binomial_is_underscored_in_the_header(self):
        """The header species is underscored to match the live source."""
        chains = Genedb.selectChains(self.BULK, 'human', is_aa=False)
        header = chains[('IGHV', KIND_VDJ)][0][0]
        self.assertIn('Homo_sapiens', header)

    def test_amino_acid_keeps_only_v(self):
        """From the amino acid set only V is taken, as translated V."""
        chains = Genedb.selectChains(self.BULK, 'human', is_aa=True)
        self.assertEqual(set(chains), {('IGHV', KIND_AA)})


if __name__ == '__main__':
    unittest.main()
