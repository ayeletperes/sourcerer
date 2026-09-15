"""
Unit tests for the airrc-imgt blended source
"""

# Info
__author__ = 'Ayelet Peres'

# Imports
import os
import tempfile
import unittest
from pathlib import Path

# Sourcerer imports
from sourcerer.Sources.AirrcImgt import AirrcImgtSource
from sourcerer.Sources.Base import DataUnit, Query
from tests.test_ogrdb import StubClient

test_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(test_path, 'data')


def readFixture(name):
    """Read a captured fixture from tests/data."""
    with open(os.path.join(data_path, name)) as handle:
        return handle.read()


class TestSearchUnits(unittest.TestCase):
    """
    Tests for composing the OGRDB and IMGT halves of the blend
    """

    def setUp(self):
        self.source = AirrcImgtSource(client=StubClient())
        self.units = self.source.searchUnits(Query(collection='human'))
        self.imgt = [u for u in self.units if u.metadata['via'] == 'imgt']
        self.ogrdb = [u for u in self.units if u.metadata['via'] == 'ogrdb']

    def test_both_sources_contribute(self):
        """The blend draws from OGRDB and IMGT, each tagged with its origin."""
        self.assertTrue(self.ogrdb)
        self.assertTrue(self.imgt)
        self.assertEqual({u.metadata['via'] for u in self.units},
                         {'ogrdb', 'imgt'})

    def test_immunoglobulin_vdj_comes_only_from_ogrdb(self):
        """No IMGT unit supplies an immunoglobulin V, D or J: those are OGRDB's."""
        ig_vdj = [u for u in self.imgt
                  if u.metadata['kind'] == 'vdj'
                  and not u.metadata['locus'].startswith('TR')]
        self.assertEqual(ig_vdj, [])

    def test_imgt_fills_tr_and_the_light_constants(self):
        """IMGT supplies all TR, and the IG constants OGRDB has no set for."""
        loci = {u.metadata['locus'] for u in self.imgt}
        self.assertEqual(loci & {'TRA', 'TRB', 'TRG', 'TRD'},
                         {'TRA', 'TRB', 'TRG', 'TRD'})
        ig_constants = {u.metadata['chain'] for u in self.imgt
                        if u.metadata['kind'] == 'constant'
                        and not u.metadata['locus'].startswith('TR')}
        # Human IGHC is OGRDB's; IMGT provides only the light constants.
        self.assertEqual(ig_constants, {'IGKC', 'IGLC'})

    def test_no_amino_acid_in_the_blend(self):
        """Amino acid V is not part of the airrc-imgt blend."""
        self.assertNotIn('vdj_aa', {u.metadata['kind'] for u in self.imgt})

    def test_mouse_takes_all_ig_constants_from_imgt(self):
        """Mouse has no OGRDB constant set, so all three come from IMGT."""
        # The IMGT gap is pure IMGT selection, so it needs no OGRDB call.
        gap = self.source._imgtGapUnits('mouse')
        ig_constants = {u.metadata['chain'] for u in gap
                        if u.metadata['kind'] == 'constant'
                        and not u.metadata['locus'].startswith('TR')}
        self.assertEqual(ig_constants, {'IGHC', 'IGKC', 'IGLC'})


class TestBuildReference(unittest.TestCase):
    """
    Tests for routing each unit back to the source that produced it
    """

    def test_each_source_writes_its_own_prefix(self):
        """OGRDB units land as airrc_ files and IMGT units as imgt_ files."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            source = AirrcImgtSource(client=StubClient())

            ogrdb_paths = {}
            for fmt in ('ungapped', 'gapped'):
                path = tmp / ('igk_%s.fasta' % fmt)
                path.write_text(readFixture('ogrdb_igk_%s.fasta' % fmt))
                ogrdb_paths[fmt] = path
            ogrdb_entries = [
                (DataUnit(unit_id='IGKappa_VJ.%s.fasta' % fmt, collection='human',
                          url='x', metadata={'species': 'human', 'locus': 'IGK',
                                             'set_name': 'IGKappa_VJ', 'format': fmt,
                                             'via': 'ogrdb'}), path)
                for fmt, path in ogrdb_paths.items()]

            imgt_page = tmp / 'IGHD.html'
            imgt_page.write_text(readFixture('imgt_ighd.html'))
            imgt_entries = [
                (DataUnit(unit_id='vdj/IGHD.html', collection='human', url='x',
                          metadata={'species': 'human', 'chain': 'IGHD',
                                    'kind': 'vdj', 'locus': 'IGH', 'segment': 'D',
                                    'via': 'imgt'}), imgt_page)]

            source.buildReference(ogrdb_entries + imgt_entries,
                                  tmp / 'reference_base')

            vdj = tmp / 'reference_base' / 'human' / 'vdj'
            self.assertTrue((vdj / 'airrc_human_IGKV.fasta').exists())
            self.assertTrue((vdj / 'imgt_human_IGHD.fasta').exists())


class TestPinnedBlend(unittest.TestCase):
    """
    Tests for the blend under a --from pin, where IMGT comes from the archive
    """

    class Listing:
        def __init__(self, names):
            self.payload = [{'name': n, 'type': 'dir'} for n in names]

        def json(self):
            return self.payload

    class ListingClient:
        def __init__(self, names):
            self.names = names

        def get(self, url):
            return TestPinnedBlend.Listing(self.names)

    def pinnedSource(self):
        """A blend whose IMGT half is pinned to an archived release."""
        source = AirrcImgtSource(client=StubClient())
        source._imgt.client = self.ListingClient(['2026-08-03_GENEDB_202631-7'])
        source.pinRelease('202631-7')

        return source

    def test_pinned_gap_resolves_to_the_nucleotide_bulk_file(self):
        """The gap becomes one archive unit, not a per-chain GENElect fetch."""
        gap = self.pinnedSource()._imgtGapUnits('human')
        self.assertEqual(len(gap), 1)
        self.assertEqual(gap[0].metadata['group'], 'nt')
        self.assertEqual(gap[0].metadata['archive'], 'genedb')

    def test_pinned_gap_keeps_only_what_ogrdb_does_not_cover(self):
        """The chain filter carries the gap: TR and the light constants only."""
        gap = self.pinnedSource()._imgtGapUnits('human')
        chains = set(gap[0].metadata['chains'])
        self.assertEqual({c for c in chains if c.startswith('IG')},
                         {'IGKC', 'IGLC'})
        self.assertIn('TRAV', chains)
        # The immunoglobulin V, D and J are OGRDB's; taking them from the bulk
        # file would overwrite the sets the blend exists to use.
        self.assertFalse({'IGHV', 'IGHD', 'IGHJ', 'IGKV', 'IGLV'} & chains)

    def test_pinned_build_writes_only_the_gap_chains(self):
        """Building the bulk file writes the gap and leaves OGRDB's chains alone."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            bulk = tmp / 'nt.fasta'
            bulk.write_text(
                '>X1|IGHV1-2*02|Homo sapiens|F\nAC.GT\n'
                '>X2|TRAV1-1*01|Homo sapiens|F\nACTG\n')
            source = self.pinnedSource()
            unit = source._imgtGapUnits('human')[0]
            unit.metadata['via'] = 'imgt'

            source.buildReference([(unit, bulk)], tmp / 'reference_base')

            vdj = tmp / 'reference_base' / 'human' / 'vdj'
            self.assertTrue((vdj / 'imgt_human_TRAV.fasta').exists())
            self.assertFalse((vdj / 'imgt_human_IGHV.fasta').exists())

    def test_gap_chains_exclude_the_translated_v(self):
        """TRAV is wanted as nucleotide only, never as the amino acid file."""
        gap = AirrcImgtSource(client=StubClient())._imgtGapChains('human')
        self.assertIn(('TRAV', 'vdj'), gap)
        self.assertNotIn(('TRAV', 'vdj_aa'), gap)


if __name__ == '__main__':
    unittest.main()
