"""
Unit tests for the commandline interface
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import contextlib
import io
import shutil
import tempfile
import unittest
from argparse import ArgumentParser, Namespace, _SubParsersAction
from pathlib import Path
from unittest import mock

import pandas

# Sourcerer imports
from sourcerer import Reference
from sourcerer.Cli import (
    applyPins,
    getArgParser,
    handleDownload,
    handleReferenceDiff,
    handleReferenceDownload,
    handleReferenceShow,
    loadMap,
)
from sourcerer.Exceptions import SourcererError
from sourcerer.Sources.Base import DataUnit, DownloadResult, Query, SourceBase
from sourcerer.Sources.Imgt import ImgtSource
from sourcerer.Sources.Oas import OasSource, newReport
from sourcerer.Sources.Ogrdb import OgrdbSource


class TestArgParser(unittest.TestCase):
    """
    Tests for parser construction
    """

    def test_returns_parser(self):
        """getArgParser returns an ArgumentParser, as autoprogram requires."""
        self.assertIsInstance(getArgParser(), ArgumentParser)

    def test_builds_without_network(self):
        """
        Parser construction performs no network I/O.

        Filter arguments are generated from the packaged snapshot, never from the
        live source. Sphinx builds the docs by calling this function, so a network
        call here would make the documentation build depend on a remote host being
        up. Patching the socket module makes any attempt fail loudly.
        """
        with mock.patch('socket.socket', side_effect=AssertionError('network access')):
            parser = getArgParser()

        self.assertIsInstance(parser, ArgumentParser)

    def test_version_flag(self):
        """--version exits zero rather than falling through to the subcommand check."""
        parser = getArgParser()
        with self.assertRaises(SystemExit) as raised:
            parser.parse_args(['--version'])

        self.assertEqual(raised.exception.code, 0)

    def test_verbose_and_quiet_are_exclusive(self):
        """-v and -q cannot be combined."""
        parser = getArgParser()
        with self.assertRaises(SystemExit):
            parser.parse_args(['-v', '-q'])

    def test_action_help_lists_the_collections(self):
        """
        `sourcerer oas download --help` names the collections it accepts.

        argparse only lists a subparser that was given a help string, so leaving
        it off left the positional section of every action's help empty and the
        user with no way to discover paired and unpaired short of reading the
        source.
        """
        parser = getArgParser()
        with mock.patch('sys.stdout', new_callable=io.StringIO) as out:
            with self.assertRaises(SystemExit):
                parser.parse_args(['oas', 'download', '--help'])

        text = out.getvalue()
        for collection in OasSource.collections:
            self.assertIn(collection, text)
            self.assertIn(OasSource.collection_help[collection], text)

    def test_collection_is_required(self):
        """
        Omitting the collection is a parse error, not a runtime one.

        The named metavar matters here: argparse reports a missing argument by
        its metavar, so the '' used at the other levels would produce an error
        message naming nothing at all.
        """
        parser = getArgParser()
        with mock.patch('sys.stderr', new_callable=io.StringIO) as err:
            with self.assertRaises(SystemExit):
                parser.parse_args(['oas', 'download'])

        self.assertIn('COLLECTION', err.getvalue())

    def test_every_argument_is_documented(self):
        """
        No argument anywhere in the tree is left without help text.

        An undocumented flag is invisible in `--help` and in the generated
        Sphinx page, so it may as well not exist for anyone who did not write
        it. Walking the tree catches the ones added later too.
        """
        undocumented = []

        def walk(parser, path):
            for action in parser._actions:
                if isinstance(action, _SubParsersAction):
                    for name, sub in action.choices.items():
                        walk(sub, '%s %s' % (path, name))
                elif not action.help and action.dest != 'help':
                    flag = '/'.join(action.option_strings) or action.dest
                    undocumented.append('%s %s' % (path, flag))

        walk(getArgParser(), 'sourcerer')

        self.assertEqual(undocumented, [])

    def test_defaults_are_not_stated_twice(self):
        """
        No help string spells out a default the formatter already appends.

        CommonHelpFormatter inherits ArgumentDefaultsHelpFormatter, which adds
        '(default: ...)' on its own, so writing one by hand produced lines
        ending in two of them that disagreed with each other.
        """
        doubled = []

        def walk(parser, path):
            for action in parser._actions:
                if isinstance(action, _SubParsersAction):
                    for name, sub in action.choices.items():
                        walk(sub, '%s %s' % (path, name))
                elif action.help and 'default:' in action.help:
                    flag = '/'.join(action.option_strings) or action.dest
                    doubled.append('%s %s' % (path, flag))

        walk(getArgParser(), 'sourcerer')

        self.assertEqual(doubled, [])


class StubSource(SourceBase):
    """
    A source that serves one unit from memory.

    Only the seams handleDownload actually touches are real: the network and the
    gzip reader are replaced, so the test exercises the command's bookkeeping
    rather than OAS parsing, which test_Oas covers.
    """

    name = 'oas'
    description = 'stub'
    collections = ('paired', 'unpaired')

    unit = DataUnit(unit_id='Study_2020/csv_paired/x_1_Paired_All.csv.gz',
                    collection='paired',
                    url='https://example.invalid/x.csv.gz',
                    metadata={'Species': 'human'}, n_sequences=2)

    def harvestSchema(self):
        raise NotImplementedError

    def searchUnits(self, query):
        return [self.unit]

    def readUnit(self, path, unit):
        raise NotImplementedError

    def normalizeChunk(self, metadata, chunk, unit, offset, report):
        raise NotImplementedError

    def validateQuery(self, collection, filters):
        return Query(collection=collection, filters=filters)

    def fetchUnit(self, unit, outdir, resume=True):
        path = Path(outdir) / unit.unit_id
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b'raw')

        return DownloadResult(unit=unit, path=path, sha256='0' * 64,
                              size_bytes=3)

    def convertUnit(self, path, unit, chunksize=50000):
        frame = pandas.DataFrame(
            {'sequence_id': ['a', 'b'], 'cell_id': ['c1', 'c1'],
             'sequence': ['ACGT', 'TGCA'], 'locus': ['IGH', 'IGK'],
             'c_call': ['IGHM', '']})
        report = newReport()
        report['rows_in'] = 1
        report['rows_out'] = 2
        report['loci'] = {'IGH', 'IGK'}

        return {'Species': 'human'}, iter([frame]), report


class TestHandleDownload(unittest.TestCase):
    """
    Tests for the download command's output bookkeeping
    """

    def setUp(self):
        self.outdir = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.outdir, ignore_errors=True)

    def runDownload(self, *formats):
        """Parse a real commandline and run the download handler against the stub."""
        argv = ['oas', 'download', 'paired', '--outdir', str(self.outdir)]
        for value in formats:
            argv += ['--format', value]

        args = getArgParser().parse_args(argv)
        args.source = 'oas'

        with mock.patch('sourcerer.Cli.getSource', return_value=StubSource(None)):
            return handleDownload(args)

    def test_converted_format_alone_still_records_the_raw_mirror(self):
        """
        Asking only for a converted format does not break the raw bookkeeping.

        The raw mirror is written unconditionally because conversion reads from
        it, so the results dict needs a 'raw' bucket even when the user never
        asked for that format. Keying the dict on the requested formats alone
        raised KeyError on the first unit.
        """
        self.assertEqual(self.runDownload('fasta'), 0)

        self.assertTrue(list(self.outdir.glob('fasta/*.fasta')))
        self.assertTrue((self.outdir / 'samplesheet_airrflow_fasta.tsv').exists())
        self.assertTrue(list(self.outdir.rglob('raw/**/*.csv.gz')))

    def test_raw_only_writes_no_samplesheet(self):
        """Raw OAS files are not an airrflow input, so no samplesheet is written."""
        self.assertEqual(self.runDownload('raw'), 0)

        self.assertEqual(list(self.outdir.glob('samplesheet_*')), [])

    def test_both_formats_write_one_samplesheet_each(self):
        """
        Each converted format gets its own samplesheet.

        airrflow's filename column names exactly one file per sample, so a single
        merged sheet could not describe both outputs.
        """
        self.assertEqual(self.runDownload('airr', 'fasta'), 0)

        self.assertTrue((self.outdir / 'samplesheet_airrflow_airr.tsv').exists())
        self.assertTrue((self.outdir / 'samplesheet_airrflow_fasta.tsv').exists())


def makeReference(root, chain='IGHV', records=(('IGHV1-2*02', 'ACGT'),)):
    """Write a minimal reference_base and return its root."""
    root = Path(root)
    (root / 'human' / 'vdj').mkdir(parents=True, exist_ok=True)
    path = root / 'human' / 'vdj' / ('imgt_human_%s.fasta' % chain)
    path.write_text(''.join('>%s\n%s\n' % record for record in records))

    return root


class TestApplyPins(unittest.TestCase):
    """
    Tests for routing a --from reference's pins to the source that can use them
    """

    def writePins(self, tmp, imgt=True, airrc=True):
        root = Path(tmp)
        if imgt:
            Reference.writeImgtMetadata(root, ['human'], '202631-7',
                                        '2026-08-24', 'x')
            Reference.writeImgtMetadata(root, ['mouse'], '202638-7',
                                        '2026-09-24', 'x')
        if airrc:
            Reference.writeAirrcMetadata(
                root, [{'species': 'human', 'locus': 'IGH', 'set': 'IGH_VDJ',
                        'version': '9', 'release_date': '2024-10-12'}],
                '2026-08-24', 'x')
        return root

    def test_imgt_takes_the_release(self):
        """An imgt source is pinned to the release the reference records."""
        with tempfile.TemporaryDirectory() as tmp:
            source = ImgtSource(client=None)
            applyPins(source, self.writePins(tmp), 'human')
        self.assertEqual(source.release, '202631-7')

    def test_the_release_pinned_is_the_one_for_that_species(self):
        """A reference holding both species must not pin mouse to human's."""
        with tempfile.TemporaryDirectory() as tmp:
            source = ImgtSource(client=None)
            applyPins(source, self.writePins(tmp), 'mouse')
        self.assertEqual(source.release, '202638-7')

    def test_ogrdb_takes_the_set_versions(self):
        """An ogrdb source is pinned to each set version, keyed by set name."""
        with tempfile.TemporaryDirectory() as tmp:
            source = OgrdbSource(client=None)
            applyPins(source, self.writePins(tmp), 'human')
        self.assertEqual(source._pins['IGH_VDJ']['version'], '9')

    def test_ogrdb_ignores_a_release_it_cannot_re_download(self):
        """An IMGT-only reference gives an ogrdb source nothing, and says so."""
        with tempfile.TemporaryDirectory() as tmp:
            source = OgrdbSource(client=None)
            with self.assertRaises(SourcererError):
                applyPins(source, self.writePins(tmp, airrc=False), 'human')

    def test_a_folder_with_no_sidecars_is_an_error(self):
        """--from pointed at an ordinary folder fails rather than fetching latest."""
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(SourcererError):
                applyPins(ImgtSource(client=None), Path(tmp), 'human')


class TestReferenceDiffCommand(unittest.TestCase):
    """
    Tests for the exit code `reference diff` reports
    """

    def run_diff(self, a, b, map_file=None):
        args = mock.Mock(reference_a=a, reference_b=b, species=None,
                         map_file=map_file)
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            code = handleReferenceDiff(args)
        return code, out.getvalue()

    def test_identical_references_exit_zero(self):
        """Nothing changed is a success, and says so."""
        with tempfile.TemporaryDirectory() as tmp:
            a = makeReference(Path(tmp) / 'a')
            b = makeReference(Path(tmp) / 'b')
            code, text = self.run_diff(a, b)
        self.assertEqual(code, 0)
        self.assertIn('identical', text)

    def test_differing_references_exit_non_zero(self):
        """A difference is a non-zero exit, so a re-download can be checked in CI."""
        with tempfile.TemporaryDirectory() as tmp:
            a = makeReference(Path(tmp) / 'a')
            b = makeReference(Path(tmp) / 'b',
                              records=(('IGHV1-2*02', 'TTTT'),))
            code, text = self.run_diff(a, b)
        self.assertEqual(code, 1)
        self.assertIn('changed', text)

    def test_a_missing_folder_is_an_error(self):
        """A path that is not a folder fails rather than comparing nothing."""
        with tempfile.TemporaryDirectory() as tmp:
            a = makeReference(Path(tmp) / 'a')
            with self.assertRaises(SourcererError):
                self.run_diff(a, Path(tmp) / 'nope')


class TestReferenceShowCommand(unittest.TestCase):
    """
    Tests for reporting what a reference folder is
    """

    def test_reports_the_release_and_contents(self):
        """show prints the pinned release and the chains the folder holds."""
        with tempfile.TemporaryDirectory() as tmp:
            root = makeReference(Path(tmp) / 'reference_base')
            Reference.writeImgtMetadata(root, ['human'], '202631-7',
                                        '2026-08-24', 'x')
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                code = handleReferenceShow(mock.Mock(folder=root, map_file=None))
        self.assertEqual(code, 0)
        self.assertIn('202631-7', out.getvalue())
        self.assertIn('IGHV', out.getvalue())


class TestLoadMap(unittest.TestCase):
    """
    Tests for reading the --map manifest off the commandline
    """

    def test_no_manifest_is_none(self):
        """--map is optional; without it nothing is declared."""
        self.assertIsNone(loadMap(mock.Mock(map_file=None)))

    def test_a_missing_manifest_is_an_error(self):
        """A manifest that is not there fails rather than being ignored."""
        with self.assertRaises(SourcererError):
            loadMap(mock.Mock(map_file=Path('/nonexistent/manifest.tsv')))

    def test_a_manifest_is_read(self):
        """The declared files come back keyed by name."""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'manifest.tsv'
            path.write_text('IGH_VDJ_V.fasta\thuman\tIGHV\n')
            mapping = loadMap(mock.Mock(map_file=path))
        self.assertEqual(mapping['IGH_VDJ_V.fasta'], ('human', 'IGHV', False))


class TestMultiSpeciesDownload(unittest.TestCase):
    """
    Tests for `download all`, which puts every species in one reference_base
    """

    def collections(self, source, action):
        """The collection names a source offers for an action."""
        parser = getArgParser()
        for entry in parser._actions:
            if not hasattr(entry, '_name_parser_map'):
                continue
            for act in entry._name_parser_map[source]._actions:
                if not hasattr(act, '_name_parser_map'):
                    continue
                for leaf in act._name_parser_map[action]._actions:
                    if hasattr(leaf, '_name_parser_map'):
                        return set(leaf._name_parser_map)
        return set()

    def test_germline_downloads_offer_all(self):
        """Every germline source can fetch its species into one folder."""
        for source in ('imgt', 'ogrdb', 'airrc-imgt'):
            self.assertIn('all', self.collections(source, 'download'))

    def test_search_does_not(self):
        """Searching two species at once would merge two unrelated hit lists."""
        self.assertNotIn('all', self.collections('imgt', 'search'))

    def test_oas_does_not(self):
        """paired and unpaired are not species; there is nothing to combine."""
        self.assertNotIn('all', self.collections('oas', 'download'))

    def test_every_species_is_fetched_into_one_folder(self):
        """`all` expands to each species, sharing one reference_base."""
        class FakeSource:
            output = 'reference'
            name = 'fake'
            collections = ('human', 'mouse')

            def __init__(self):
                self.searched = []

            def validateQuery(self, collection, filters):
                return Query(collection=collection, filters=filters)

            def searchUnits(self, query):
                self.searched.append(query.collection)
                return []

        source = FakeSource()
        with tempfile.TemporaryDirectory() as tmp:
            args = Namespace(collection='all', source='fake', outdir=Path(tmp),
                             from_ref=None, resolve_doi=False, dry_run=True,
                             limit=None, no_resume=False, igblast=False,
                             igblast_out=None, compare=None)
            self.assertEqual(handleReferenceDownload(args, source), 0)

        self.assertEqual(source.searched, ['human', 'mouse'])


if __name__ == '__main__':
    unittest.main()
