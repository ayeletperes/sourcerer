"""
Unit tests for the germline reference output layer
"""

# Info
__author__ = 'Ayelet Peres'

# Imports
import io
import shutil
import tarfile
import tempfile
import unittest
from pathlib import Path

# Sourcerer imports
from sourcerer import Reference
from sourcerer.Exceptions import SourcererError

HAS_MAKEBLASTDB = shutil.which('makeblastdb') is not None


class Canned:
    """A stand-in response exposing just .text for the directory-index mirror."""

    def __init__(self, text=''):
        self.text = text


class EmptyIndexClient:
    """An HTTP client whose directory listings are empty, so nothing mirrors."""

    def get(self, url):
        return Canned('<html><body></body></html>')

    def fetch(self, url, dest, **kwargs):
        raise AssertionError('an empty index should not trigger a fetch')


class TestFastaHelpers(unittest.TestCase):
    """
    Tests for the FASTA parse, name and clean helpers
    """

    def test_parse_round_trips_records(self):
        """parseFasta reads headers and collapses wrapped sequence lines."""
        text = '>a\nACGT\nACGT\n>b\nTTTT\n'
        self.assertEqual(Reference.parseFasta(text),
                         [('a', 'ACGTACGT'), ('b', 'TTTT')])

    def test_parse_ignores_blank_lines(self):
        """Blank lines between records are not mistaken for sequence."""
        self.assertEqual(Reference.parseFasta('\n>a\n\nACGT\n\n'),
                         [('a', 'ACGT')])

    def test_allele_name_from_imgt_pipe_header(self):
        """The allele name is the second pipe field of an IMGT header."""
        header = 'X02897|IGHV1-2*02|Homo sapiens|F|V-REGION'
        self.assertEqual(Reference.alleleName(header), 'IGHV1-2*02')

    def test_allele_name_from_plain_header(self):
        """A header with no pipe yields its first whitespace-delimited token."""
        self.assertEqual(Reference.alleleName('IGKV1-12*01 extra'), 'IGKV1-12*01')

    def test_clean_degaps_uppercases_and_dedups(self):
        """cleanForBlast removes gaps, uppercases, and drops repeated names."""
        records = [('X1|IGHV1-2*02|H', 'ac.gt'),
                   ('X2|IGHV1-2*02|H', 'aaaa'),   # duplicate name, dropped
                   ('IGHV3*01', 'gg..cc')]
        self.assertEqual(
            Reference.cleanForBlast(records),
            [('IGHV1-2*02', 'ACGT'), ('IGHV3*01', 'GGCC')])


class TestReferencePaths(unittest.TestCase):
    """
    Tests for where a chain's FASTA lands in the reference tree
    """

    def test_vdj_path(self):
        """A V/D/J chain lands under vdj/ with the source prefix."""
        path = Reference.referenceFastaPath('/ref', 'imgt', 'human',
                                            Reference.KIND_VDJ, 'IGHV')
        self.assertEqual(path, Path('/ref/human/vdj/imgt_human_IGHV.fasta'))

    def test_constant_path(self):
        """A constant chain lands under constant/."""
        path = Reference.referenceFastaPath('/ref', 'airrc', 'human',
                                            Reference.KIND_CONSTANT, 'IGHC')
        self.assertEqual(path, Path('/ref/human/constant/airrc_human_IGHC.fasta'))

    def test_amino_acid_path_is_disambiguated(self):
        """Amino acid V carries an aa_ tag so it cannot collide with nucleotide V."""
        path = Reference.referenceFastaPath('/ref', 'imgt', 'mouse',
                                            Reference.KIND_AA, 'IGHV')
        self.assertEqual(path,
                         Path('/ref/mouse/vdj_aa/imgt_aa_mouse_IGHV.fasta'))


class TestExtractTar(unittest.TestCase):
    """
    Tests for the path-traversal guard on archive extraction
    """

    def test_rejects_member_escaping_destination(self):
        """A member pointing outside the destination is refused, not written."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            archive = tmp / 'evil.tar'
            with tarfile.open(archive, 'w') as tar:
                info = tarfile.TarInfo('../escaped.txt')
                info.size = 3
                tar.addfile(info, io.BytesIO(b'bad'))

            dest = tmp / 'out'
            dest.mkdir()
            with self.assertRaises(SourcererError):
                Reference.extractTar(archive, dest)
            self.assertFalse((tmp / 'escaped.txt').exists())


@unittest.skipUnless(HAS_MAKEBLASTDB, 'makeblastdb not on PATH')
class TestBuildIgblastBase(unittest.TestCase):
    """
    Tests for aggregating a reference_base into BLAST databases
    """

    def _referenceBase(self, root):
        vdj = root / 'human' / 'vdj'
        vdj.mkdir(parents=True)
        (vdj / 'imgt_human_IGHV.fasta').write_text('>IGHV1-2*02\nACGTACGTACGT\n')
        (vdj / 'imgt_human_IGHJ.fasta').write_text('>IGHJ1*01\nTTTTGGGGCCCC\n')

    def test_builds_present_and_skips_absent(self):
        """Only chains with sequences build a database; the rest are skipped."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            reference = tmp / 'reference_base'
            self._referenceBase(reference)

            report = Reference.buildIgblastBase(reference, tmp / 'igblast_base',
                                                EmptyIndexClient(),
                                                species=['human'])

            self.assertIn('human_ig_v', report.built)
            self.assertIn('human_ig_j', report.built)
            self.assertIn('human_ig_d', report.skipped_empty)
            self.assertIn('human_tr_v', report.skipped_empty)
            self.assertTrue((tmp / 'igblast_base' / 'database'
                             / 'human_ig_v.nsq').exists())

    def test_missing_makeblastdb_is_a_clear_error(self):
        """runMakeblastdb names the missing binary rather than failing obscurely."""
        original = shutil.which
        try:
            shutil.which = lambda name: None
            with tempfile.TemporaryDirectory() as tmp:
                fasta = Path(tmp) / 'x.fasta'
                fasta.write_text('>a\nACGT\n')
                with self.assertRaises(SourcererError) as caught:
                    Reference.runMakeblastdb(fasta, Path(tmp) / 'x', 'nucl')
                self.assertIn('makeblastdb', str(caught.exception))
        finally:
            shutil.which = original


class TestParseReferenceName(unittest.TestCase):
    """
    Tests for reading (species, chain, is_aa) from a filename
    """

    def test_flat_name(self):
        """A bare species_chain name is recognised."""
        self.assertEqual(Reference.parseReferenceName('human_IGHV.fasta'),
                         ('human', 'IGHV', False))

    def test_prefix_is_ignored(self):
        """A source prefix such as imgt_ or airrc_ is allowed and ignored."""
        self.assertEqual(Reference.parseReferenceName('imgt_human_IGHV.fasta'),
                         ('human', 'IGHV', False))
        self.assertEqual(Reference.parseReferenceName('airrc_mouse_IGKC.fasta'),
                         ('mouse', 'IGKC', False))

    def test_aa_marker(self):
        """aa_ marks a translated V, with or without a prefix."""
        self.assertEqual(Reference.parseReferenceName('aa_human_IGHV.fasta'),
                         ('human', 'IGHV', True))
        self.assertEqual(Reference.parseReferenceName('imgt_aa_human_IGHV.fasta'),
                         ('human', 'IGHV', True))

    def test_rejects_unknown_species_or_chain_or_extension(self):
        """A name that is not species_knownchain.fasta is not recognised."""
        self.assertIsNone(Reference.parseReferenceName('rat_IGHV.fasta'))
        self.assertIsNone(Reference.parseReferenceName('human_XYZ.fasta'))
        self.assertIsNone(Reference.parseReferenceName('human_IGHV.txt'))
        self.assertIsNone(Reference.parseReferenceName('notes.fasta'))


class TestPlanReference(unittest.TestCase):
    """
    Tests for planning an IgBLAST build from a reference folder
    """

    def _write(self, path, name, body):
        path.mkdir(parents=True, exist_ok=True)
        (path / name).write_text(body)

    def test_flat_and_nested_layouts_plan_the_same(self):
        """The same files plan identically whether flat or nested."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            flat, nested = tmp / 'flat', tmp / 'nested'
            self._write(flat, 'human_IGHV.fasta', '>IGHV1-2*02\nACGT\n')
            self._write(flat, 'human_IGHJ.fasta', '>IGHJ1*01\nTTGG\n')
            self._write(nested / 'human' / 'vdj', 'imgt_human_IGHV.fasta',
                        '>IGHV1-2*02\nACGT\n')
            self._write(nested / 'human' / 'vdj', 'imgt_human_IGHJ.fasta',
                        '>IGHJ1*01\nTTGG\n')

            built_flat = {b for b, _t, _r in
                          Reference.planReference(flat).databases}
            built_nested = {b for b, _t, _r in
                            Reference.planReference(nested).databases}
            self.assertEqual(built_flat, {'human_ig_v', 'human_ig_j'})
            self.assertEqual(built_flat, built_nested)

    def test_reports_empty_unrecognized_and_duplicates(self):
        """The plan surfaces gaps, unknown names, and dropped duplicate names."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            self._write(tmp, 'human_IGHV.fasta',
                        '>IGHV1-2*02\nACGT\n>IGHV1-2*02\nAAAA\n')  # duplicate name
            self._write(tmp, 'notes.fasta', '>x\nACGT\n')          # unrecognized

            plan = Reference.planReference(tmp)
            self.assertEqual(plan.found_species, ['human'])
            self.assertEqual(plan.duplicates.get('human_ig_v'), 1)
            self.assertIn('human_ig_d', plan.empty)
            self.assertEqual([p.name for p in plan.unrecognized], ['notes.fasta'])
            self.assertTrue(plan.ok)

    def test_species_filter(self):
        """--species narrows the plan to the requested species."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            self._write(tmp, 'human_IGHV.fasta', '>IGHV1-2*02\nACGT\n')
            self._write(tmp, 'mouse_IGHV.fasta', '>IGHV1*01\nACGT\n')

            built = {b for b, _t, _r in
                     Reference.planReference(tmp, species=['human']).databases}
            self.assertEqual(built, {'human_ig_v'})


if __name__ == '__main__':
    unittest.main()
