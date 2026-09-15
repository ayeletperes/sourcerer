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

import yaml

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


class TestMetadata(unittest.TestCase):
    """
    Tests for the reference_base provenance sidecars
    """

    def test_imgt_metadata_records_the_release(self):
        """IMGT.yaml carries the release, the datum a download date does not."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            path = Reference.writeImgtMetadata(tmp, ['human'], '202631-7',
                                               '2026-08-14', 'sourcerer 0.1.0')
            self.assertEqual(path.name, 'IMGT.yaml')
            pins = Reference.loadReferencePins(tmp)
            self.assertEqual(pins['imgt']['species']['human']['release'],
                             '202631-7')
            self.assertIsNone(pins['airrc'])

    def test_airrc_metadata_records_each_set_version(self):
        """AIRRC.yaml carries a version per set, so --from can re-pin them."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            sets = [{'species': 'human', 'locus': 'IGH', 'set': 'IGH_VDJ',
                     'version': '3', 'release_date': '2024-01-01'}]
            Reference.writeAirrcMetadata(tmp, sets, '2026-08-14', 'sourcerer 0.1.0')
            pins = Reference.loadReferencePins(tmp)
            self.assertEqual(pins['airrc']['sets'][0]['version'], '3')

    def test_load_pins_reads_both_from_a_blended_reference(self):
        """A reference_base with both sidecars pins both sources at once."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            Reference.writeImgtMetadata(tmp, ['human'], '202631-7', '2026-08-14',
                                        'x')
            Reference.writeAirrcMetadata(tmp, [{'set': 'IGH_VDJ', 'version': '3'}],
                                         '2026-08-14', 'x')
            pins = Reference.loadReferencePins(tmp)
            self.assertIsNotNone(pins['imgt'])
            self.assertIsNotNone(pins['airrc'])

    def test_load_pins_rejects_a_folder_without_sidecars(self):
        """--from a folder that holds neither sidecar is an error, not a no-op."""
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(SourcererError):
                Reference.loadReferencePins(Path(tmp))

    def test_build_metadata_records_the_build_and_carries_provenance(self):
        """reference build writes a build record and carries the IMGT sidecar in."""
        with tempfile.TemporaryDirectory() as src, \
                tempfile.TemporaryDirectory() as out:
            src, out = Path(src), Path(out)
            Reference.writeImgtMetadata(src, ['human'], '202631-7', '2026-08-14',
                                        'x')
            report = Reference.ReferenceReport(built=['human_ig_v'],
                                               skipped_empty=['human_tr_d'])

            written = Reference.writeBuildMetadata(out, report, src, '2026-08-14',
                                                   'sourcerer 0.1.0')
            self.assertIn(out / 'IMGT.yaml', written)  # carried forward
            record = yaml.safe_load((out / 'sourcerer_build.yaml').read_text())
            self.assertEqual(record['databases'], ['human_ig_v'])
            # Recorded relative to the built directory, so moving the pair keeps
            # the link; resolving it from there must land back on the source.
            self.assertFalse(Path(record['built_from']).is_absolute())
            self.assertEqual((out / record['built_from']).resolve(),
                             src.resolve())


class TestReferenceMap(unittest.TestCase):
    """
    Tests for the manifest that names files the naming rule cannot place
    """

    MANIFEST = ('# file\tspecies\tchain\n'
                '\n'
                'IGH_VDJ_V.fasta\thuman\tIGHV\n'
                'C57BL-6_IGH_V.fasta\tmouse\tIGHV\n'
                'translated.fasta\thuman\tIGHV\taa\n')

    def writeManifest(self, tmp, text=None):
        path = Path(tmp) / 'manifest.tsv'
        path.write_text(self.MANIFEST if text is None else text)
        return path

    def test_reads_species_chain_and_the_aa_marker(self):
        """Three columns place a file; a fourth marks it as translated V."""
        with tempfile.TemporaryDirectory() as tmp:
            mapping = Reference.loadReferenceMap(self.writeManifest(tmp))
        self.assertEqual(mapping['IGH_VDJ_V.fasta'], ('human', 'IGHV', False))
        self.assertEqual(mapping['C57BL-6_IGH_V.fasta'], ('mouse', 'IGHV', False))
        self.assertEqual(mapping['translated.fasta'], ('human', 'IGHV', True))

    def test_unknown_species_or_chain_is_an_error(self):
        """A manifest exists to be explicit, so a bad row fails rather than skips."""
        with tempfile.TemporaryDirectory() as tmp:
            for bad in ('a.fasta\tmartian\tIGHV\n', 'a.fasta\thuman\tIGXV\n',
                        'a.fasta\thuman\n'):
                path = self.writeManifest(tmp, bad)
                with self.assertRaises(SourcererError):
                    Reference.loadReferenceMap(path)

    def test_manifest_places_a_file_the_naming_rule_skips(self):
        """A file with no species or chain in its name is found via the manifest."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / 'IGH_VDJ_V.fasta').write_text('>IGHV1-2*02\nACGT\n')
            files, unrecognized = Reference.discoverReference(tmp)
            self.assertEqual(files, [])
            self.assertEqual(len(unrecognized), 1)

            mapping = Reference.loadReferenceMap(self.writeManifest(tmp))
            files, unrecognized = Reference.discoverReference(tmp, mapping=mapping)
            self.assertEqual(unrecognized, [])
            self.assertEqual(files[0][:3], ('human', 'IGHV', False))

    def test_manifest_overrides_a_name_that_would_parse(self):
        """The manifest wins, so it can correct a misread name, not only place one."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / 'imgt_human_IGHV.fasta').write_text('>IGHV1-2*02\nACGT\n')
            path = self.writeManifest(
                tmp, 'imgt_human_IGHV.fasta\tmouse\tIGKV\n')
            mapping = Reference.loadReferenceMap(path)
            files, _ = Reference.discoverReference(tmp, mapping=mapping)
        self.assertEqual(files[0][:3], ('mouse', 'IGKV', False))

    def test_other_fasta_extensions_are_read(self):
        """A manifest is no use if the file it names is never looked at."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / 'human_IGHV.fa').write_text('>IGHV1-2*02\nACGT\n')
            (tmp / 'human_IGHJ.fna').write_text('>IGHJ4*02\nACGT\n')
            files, unrecognized = Reference.discoverReference(tmp)
        self.assertEqual(unrecognized, [])
        self.assertEqual({f[1] for f in files}, {'IGHV', 'IGHJ'})


class TestShortenName(unittest.TestCase):
    """
    Tests for fitting allele names inside the BLAST identifier limit
    """

    #: A VDJbase-style novel allele name, the only kind that reaches the limit.
    LONG = 'IGHV1-18*01_a12g_t45c_g78a_c101t_a134g_t167c_g200a_c233t_a266g'

    def test_a_normal_name_is_untouched(self):
        """Nothing IMGT or OGRDB publishes is anywhere near the limit."""
        for name in ('IGHV1-2*02', 'IGHV3/OR16-13*01', 'IGKJ0-4JXG*00'):
            self.assertEqual(Reference.shortenName(name), name)

    def test_a_long_name_is_cut_to_exactly_the_limit(self):
        """makeblastdb refuses the whole database over 50, so it must fit."""
        short = Reference.shortenName(self.LONG)
        self.assertGreater(len(self.LONG), Reference.BLAST_NAME_LIMIT)
        self.assertEqual(len(short), Reference.BLAST_NAME_LIMIT)
        self.assertTrue(short.startswith('IGHV1-18*01_'))

    def test_shortening_is_deterministic(self):
        """The aux and ndm files are keyed by name, so it cannot drift per run."""
        self.assertEqual(Reference.shortenName(self.LONG),
                         Reference.shortenName(self.LONG))

    def test_names_sharing_a_head_stay_distinct(self):
        """Two novel alleles off the same gene differ only in the discarded tail."""
        other = self.LONG[:-1] + 'c'
        self.assertNotEqual(Reference.shortenName(self.LONG),
                            Reference.shortenName(other))

    def test_records_report_what_changed(self):
        """The mapping back to the original is what makes a v_call traceable."""
        records, renamed = Reference.shortenForBlast(
            [('IGHV1-2*02', 'ACGT'), (self.LONG, 'TTTT')])
        self.assertEqual(renamed, {self.LONG: Reference.shortenName(self.LONG)})
        self.assertEqual(records[0], ('IGHV1-2*02', 'ACGT'))

    def test_the_plan_carries_and_reports_the_renames(self):
        """A build that had to shorten names says so before it is run."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / 'human_IGHV.fasta').write_text('>%s\nACGT\n' % self.LONG)
            plan = Reference.planReference(tmp)
        self.assertEqual(list(plan.renamed), [self.LONG])
        self.assertIn('shortened', plan.summary())

    def test_the_log_maps_every_shortened_name_back(self):
        """The log is the only record of what a shortened v_call really was."""
        with tempfile.TemporaryDirectory() as tmp:
            path = Reference.writeRenameLog(
                Path(tmp), {self.LONG: Reference.shortenName(self.LONG)})
            rows = [line.split('\t') for line in path.read_text().splitlines()
                    if not line.startswith('#')]
        self.assertEqual(rows, [[self.LONG, Reference.shortenName(self.LONG)]])

    def test_no_log_is_written_when_nothing_was_shortened(self):
        """The usual case writes no file: the caller only logs when it renamed."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / 'human_IGHV.fasta').write_text('>IGHV1-2*02\nACGT\n')
            plan = Reference.planReference(tmp)
            self.assertEqual(plan.renamed, {})
            self.assertFalse((tmp / Reference.RENAME_LOG).exists())

    @unittest.skipUnless(HAS_MAKEBLASTDB, 'makeblastdb is not installed')
    def test_makeblastdb_accepts_the_shortened_name(self):
        """
        The point of the exercise: the long name fails, the shortened one builds.

        makeblastdb -parse_seqids rejects an identifier over 50 characters
        outright, taking the whole database with it.
        """
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            for label, name in (('long', self.LONG),
                                ('short', Reference.shortenName(self.LONG))):
                fasta = tmp / ('%s.fasta' % label)
                fasta.write_text('>%s\nACGTACGTACGTACGTACGTACGT\n' % name)
                if label == 'long':
                    with self.assertRaises(SourcererError):
                        Reference.runMakeblastdb(fasta, tmp / label, 'nucl')
                else:
                    Reference.runMakeblastdb(fasta, tmp / label, 'nucl')
                    self.assertTrue((tmp / 'short.nin').exists())


class TestAuxCoverage(unittest.TestCase):
    """
    Tests for spotting J alleles the mirrored NCBI auxiliary file cannot cover
    """

    def build(self, tmp, j_names):
        tmp = Path(tmp)
        (tmp / 'optional_file').mkdir(parents=True, exist_ok=True)
        (tmp / 'optional_file' / 'mouse_gl.aux').write_text(
            '#name\tj_codon_frame\tchain_type\tj_cdr3_end\textra_bps\n'
            'IGHJ1*01\t0\tJH\t17\t1\n')
        ref = tmp / 'ref'
        ref.mkdir(exist_ok=True)
        (ref / 'mouse_IGHJ.fasta').write_text(
            ''.join('>%s\nACGT\n' % name for name in j_names))

        return Reference.planReference(ref), tmp

    def test_ogrdb_mouse_j_names_are_reported_as_uncovered(self):
        """OGRDB names its mouse J alleles in a way NCBI's file never lists."""
        with tempfile.TemporaryDirectory() as tmp:
            plan, out = self.build(tmp, ['IGHJ0-32C2*00', 'IGHJ0-7IA7*00'])
            missing = Reference.checkAuxCoverage(out, plan)
        self.assertEqual(missing, {'mouse': ['IGHJ0-32C2*00', 'IGHJ0-7IA7*00']})

    def test_a_covered_reference_reports_nothing(self):
        """A J allele NCBI does name needs no warning."""
        with tempfile.TemporaryDirectory() as tmp:
            plan, out = self.build(tmp, ['IGHJ1*01'])
            self.assertEqual(Reference.checkAuxCoverage(out, plan), {})

    def test_no_auxiliary_file_is_not_a_finding(self):
        """Without a mirrored file there is nothing to check against."""
        with tempfile.TemporaryDirectory() as tmp:
            plan, out = self.build(tmp, ['IGHJ0-32C2*00'])
            (Path(out) / 'optional_file' / 'mouse_gl.aux').unlink()
            self.assertEqual(Reference.checkAuxCoverage(out, plan), {})


class TestDescribeReference(unittest.TestCase):
    """
    Tests for reporting what a reference folder is and where it came from
    """

    def build(self, tmp):
        root = Path(tmp) / 'reference_base'
        (root / 'human' / 'vdj').mkdir(parents=True)
        (root / 'human' / 'vdj' / 'imgt_human_IGHV.fasta').write_text(
            '>X1|IGHV1-2*02|Homo_sapiens|F\nACGT\n')
        return root

    def test_reports_the_release_and_the_contents(self):
        """The sidecar and the FASTAs present are both reported."""
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build(tmp)
            Reference.writeImgtMetadata(root, ['human'], '202631-7',
                                        '2026-08-24', 'sourcerer test')
            text = Reference.describeReference(root)
        self.assertIn('202631-7', text)
        self.assertIn('IGHV', text)
        self.assertIn('1 FASTA file(s)', text)

    def test_accepts_the_directory_holding_a_reference_base(self):
        """Pointing at the download root finds the reference_base inside it."""
        with tempfile.TemporaryDirectory() as tmp:
            root = self.build(tmp)
            Reference.writeImgtMetadata(root, ['human'], '202631-7',
                                        '2026-08-24', 'sourcerer test')
            text = Reference.describeReference(Path(tmp))
        self.assertIn('202631-7', text)

    def test_a_folder_without_sidecars_says_so(self):
        """A hand-assembled folder is described, not rejected."""
        with tempfile.TemporaryDirectory() as tmp:
            text = Reference.describeReference(self.build(tmp))
        self.assertIn('no provenance sidecars', text)
        self.assertIn('IGHV', text)


class TestInexactRelease(unittest.TestCase):
    """
    Tests for recording a release the archive did not hold exactly
    """

    def test_plain_download_records_no_request(self):
        """A download of the current release has nothing to have asked for."""
        with tempfile.TemporaryDirectory() as tmp:
            path = Reference.writeImgtMetadata(Path(tmp), ['human'], '202631-7',
                                               '2026-08-24', 'x')
            record = yaml.safe_load(path.read_text())
        entry = record['species']['human']
        self.assertEqual(entry['release'], '202631-7')
        self.assertNotIn('requested', entry)
        self.assertNotIn('exact', entry)

    def test_substituted_release_records_what_was_asked_for(self):
        """A neighbouring release is recorded as a substitute, not as the original."""
        with tempfile.TemporaryDirectory() as tmp:
            path = Reference.writeImgtMetadata(
                Path(tmp), ['human'], '202630-7', '2026-08-24', 'x',
                requested='202629-7', exact=False)
            record = yaml.safe_load(path.read_text())
        entry = record['species']['human']
        self.assertEqual(entry['release'], '202630-7')
        self.assertEqual(entry['requested'], '202629-7')
        self.assertFalse(entry['exact'])


class TestMetadataMerge(unittest.TestCase):
    """
    Tests that a second species downloaded into one reference_base is kept
    """

    def test_a_second_species_does_not_erase_the_first(self):
        """A download writes one species; the folder holds every one written."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            Reference.writeImgtMetadata(tmp, ['human'], '202631-7', '2026-08-24', 'x')
            path = Reference.writeImgtMetadata(tmp, ['mouse'], '202638-7',
                                               '2026-09-24', 'x')
            record = yaml.safe_load(path.read_text())
        self.assertEqual(sorted(record['species']), ['human', 'mouse'])
        # Recorded per species: the two were downloaded from different builds.
        self.assertEqual(record['species']['human']['release'], '202631-7')
        self.assertEqual(record['species']['mouse']['release'], '202638-7')

    def test_re_downloading_a_species_replaces_only_its_entry(self):
        """Fetching human again updates human and leaves mouse untouched."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            Reference.writeImgtMetadata(tmp, ['human'], '202631-7', '2026-08-24', 'x')
            Reference.writeImgtMetadata(tmp, ['mouse'], '202638-7', '2026-09-24', 'x')
            path = Reference.writeImgtMetadata(tmp, ['human'], '202640-1',
                                               '2026-10-01', 'x')
            record = yaml.safe_load(path.read_text())
        self.assertEqual(record['species']['human']['release'], '202640-1')
        self.assertEqual(record['species']['mouse']['release'], '202638-7')

    def test_airrc_sets_from_both_species_are_kept(self):
        """The mouse download joins the human sets rather than replacing them."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            Reference.writeAirrcMetadata(tmp, [{'species': 'human', 'locus': 'IGH',
                                                'set': 'IGH_VDJ', 'version': '10'}],
                                         '2026-08-24', 'x')
            path = Reference.writeAirrcMetadata(
                tmp, [{'species': 'mouse', 'locus': 'IGK',
                       'set': 'IGKJ (all strains)', 'version': '1'}],
                '2026-09-24', 'x')
            record = yaml.safe_load(path.read_text())
        self.assertEqual([(s['species'], s['set']) for s in record['sets']],
                         [('human', 'IGH_VDJ'), ('mouse', 'IGKJ (all strains)')])

    def test_re_downloading_a_set_updates_it_in_place(self):
        """A newer version of a set replaces its entry, it does not duplicate it."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            Reference.writeAirrcMetadata(tmp, [{'species': 'human', 'locus': 'IGH',
                                                'set': 'IGH_VDJ', 'version': '9'}],
                                         '2026-08-24', 'x')
            path = Reference.writeAirrcMetadata(
                tmp, [{'species': 'human', 'locus': 'IGH', 'set': 'IGH_VDJ',
                       'version': '10'}], '2026-09-24', 'x')
            record = yaml.safe_load(path.read_text())
        self.assertEqual(len(record['sets']), 1)
        self.assertEqual(record['sets'][0]['version'], '10')


class TestDiff(unittest.TestCase):
    """
    Tests for comparing two reference folders allele by allele
    """

    def _write(self, root, name, text):
        path = Path(root) / name
        path.write_text(text)

    def test_identical_references_are_reported_same(self):
        """Two folders with the same alleles and sequences compare equal."""
        with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
            self._write(a, 'human_IGHV.fasta', '>IGHV1-2*02\nACGTACGT\n')
            self._write(b, 'human_IGHV.fasta', '>IGHV1-2*02\nACGTACGT\n')
            diff = Reference.diffReference(a, b)
            self.assertTrue(diff.same)
            self.assertEqual(diff.chains[0].identical, 1)

    def test_gaps_and_case_do_not_count_as_a_difference(self):
        """A gapped, lower-case copy of an allele is not a false change."""
        with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
            self._write(a, 'human_IGHV.fasta', '>IGHV1-2*02\nAC..GTACGT\n')
            self._write(b, 'human_IGHV.fasta', '>IGHV1-2*02\nacgtacgt\n')
            self.assertTrue(Reference.diffReference(a, b).same)

    def test_added_removed_and_changed_are_classified(self):
        """Each kind of divergence lands in the right bucket."""
        with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
            self._write(a, 'human_IGHV.fasta',
                        '>IGHV1-2*02\nACGT\n>IGHV3-7*01\nTTTT\n')
            self._write(b, 'human_IGHV.fasta',
                        '>IGHV1-2*02\nACGG\n>IGHV4-4*01\nGGGG\n')
            diff = Reference.diffReference(a, b)
            self.assertFalse(diff.same)
            chain = diff.chains[0]
            self.assertEqual(chain.changed, ['IGHV1-2*02'])
            self.assertEqual(chain.added, ['IGHV4-4*01'])
            self.assertEqual(chain.removed, ['IGHV3-7*01'])

    def test_imgt_pipe_headers_align_with_bare_names(self):
        """An IMGT pipe header and a bare OGRDB name for one allele line up."""
        with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
            self._write(a, 'imgt_human_IGHV.fasta',
                        '>X02897|IGHV1-2*02|Homo_sapiens|F\nACGT\n')
            self._write(b, 'airrc_human_IGHV.fasta', '>IGHV1-2*02\nACGT\n')
            self.assertTrue(Reference.diffReference(a, b).same)


if __name__ == '__main__':
    unittest.main()
