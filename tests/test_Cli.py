"""
Unit tests for the commandline interface
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import csv
import io
import shutil
import tempfile
import unittest
from argparse import ArgumentParser, _SubParsersAction
from pathlib import Path
from unittest import mock

import pandas

# Sourcerer imports
from sourcerer.Cli import (
    NCBI_EVIDENCE_COLUMNS,
    getArgParser,
    handleDownload,
    handleOasVerify,
)
from sourcerer.Http import HttpClient
from sourcerer.Sources.Base import DataUnit, DownloadResult, Query, SourceBase
from sourcerer.Sources.Oas import OasSource, newReport
from tests.FakeHttp import FakeResponse, FakeSession


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


#: A samplesheet header narrow enough for verify's own tests: it only reads
#: sample_id, sample_name and subject_id, so the rest is set dressing.
VERIFY_COLUMNS = ('sample_id', 'filename', 'subject_id', 'species', 'sample_name')

#: One esearch/esummary/efetch round trip resolving SRR1 to BL-110_VDJ, the
#: same fixture shape test_Ncbi.py exercises in isolation; this only checks
#: that handleOasVerify wires it into the evidence TSV and --apply correctly.
NCBI_ROUTES = {
    'esearch': FakeResponse(200, b'<eSearchResult><IdList><Id>1</Id>'
                                 b'</IdList></eSearchResult>'),
    'esummary': FakeResponse(200,
        b'<eSummaryResult><DocSum><Id>1</Id>'
        b'<Item Name="ExpXml" Type="String">'
        b'&lt;Summary&gt;&lt;Title&gt;GSM1: BL-110_VDJ&lt;/Title&gt;&lt;/Summary&gt;'
        b'&lt;Biosample&gt;SAMN1&lt;/Biosample&gt;</Item>'
        b'<Item Name="Runs" Type="String">'
        b'&lt;Run acc="SRR1" total_spots="1"/&gt;</Item>'
        b'</DocSum></eSummaryResult>'),
    'efetch': FakeResponse(200,
        b'<BioSampleSet><BioSample accession="SAMN1">'
        b'<Ids><Id db="BioSample">SAMN1</Id></Ids>'
        b'<Description><Title>BL-110_VDJ</Title></Description>'
        b'</BioSample></BioSampleSet>'),
}


class TestHandleOasVerify(unittest.TestCase):
    """
    Tests for the verify command's evidence report
    """

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        self.samplesheet = self.tmp / 'samplesheet_airrflow_airr.tsv'

    def writeSamplesheet(self, rows):
        """Write a samplesheet with just the columns verify needs."""
        with open(self.samplesheet, 'w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(VERIFY_COLUMNS),
                                    delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)

    def runVerify(self, extra_argv=()):
        """Parse a real commandline and run the verify handler against a fake NCBI."""
        argv = ['oas', 'verify', str(self.samplesheet)] + list(extra_argv)
        args = getArgParser().parse_args(argv)

        fake_client = HttpClient(delay=0, backoff=0,
                                 session=FakeSession(lambda method, url, headers, i: next(
                                     response for substring, response in NCBI_ROUTES.items()
                                     if substring in url)))
        with mock.patch('sourcerer.Cli.HttpClient', return_value=fake_client):
            return handleOasVerify(args)

    def readReport(self, path=None):
        """Read back the evidence report as a list of dicts."""
        path = path or self.samplesheet.with_name(
            self.samplesheet.stem + '.ncbi_evidence' + self.samplesheet.suffix)
        with open(path, newline='') as handle:
            return list(csv.DictReader(handle, delimiter='\t'))

    def test_a_real_subject_id_is_never_looked_up(self):
        """
        A row that already has a subject does not touch the network at all.

        HttpClient is patched to a client whose session raises on any call,
        so a lookup attempt would fail loudly rather than silently succeed --
        the row still appearing correctly in the report is the assertion.
        """
        self.writeSamplesheet([{'sample_id': 'ssr_1', 'filename': 'x.tsv',
                               'subject_id': 'Donor-2', 'species': 'human',
                               'sample_name': 'Study/csv_paired/SRR1_1_Paired_All.csv.gz'}])

        empty_client = HttpClient(delay=0, backoff=0, session=FakeSession(
            lambda *a: (_ for _ in ()).throw(AssertionError('no network call was expected'))))
        args = getArgParser().parse_args(['oas', 'verify', str(self.samplesheet)])
        with mock.patch('sourcerer.Cli.HttpClient', return_value=empty_client):
            self.assertEqual(handleOasVerify(args), 0)

        rows = self.readReport()
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['status'], 'passthrough')
        self.assertEqual(rows[0]['biosample_id'], 'Donor-2')
        self.assertEqual(rows[0]['biosample_id_suggested'], 'Donor-2')

    def test_resolved_row_gets_both_a_raw_and_a_suggested_biosample_id(self):
        """
        A resolved row's report names its BioSample, a check link, and both
        forms of biosample_id -- no flag needed to choose between them.
        """
        self.writeSamplesheet([{'sample_id': 'ssr_1', 'filename': 'x.tsv',
                               'subject_id': 'no', 'species': 'human',
                               'sample_name': 'Study/csv_paired/SRR1_1_Paired_All.csv.gz'}])

        self.assertEqual(self.runVerify(), 0)

        rows = self.readReport()
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['status'], 'ok')
        self.assertEqual(rows[0]['biosample_accession'], 'SAMN1')
        self.assertEqual(rows[0]['biosample_id'], 'BL-110_VDJ')
        self.assertEqual(rows[0]['biosample_id_suggested'], 'BL-110')
        self.assertIn('SAMN1', rows[0]['biosample_url'])

    def test_report_carries_every_input_column(self):
        """
        The report is a superset of the input, not a separate NCBI-only
        file: airrflow-required columns absent from NCBI_EVIDENCE_COLUMNS
        (filename, species, ...) must survive untouched, in their original
        position, so the report can be used as airrflow --input directly.
        """
        self.writeSamplesheet([{'sample_id': 'ssr_1', 'filename': 'fasta/x.fasta',
                               'subject_id': 'no', 'species': 'human',
                               'sample_name': 'Study/csv_paired/SRR1_1_Paired_All.csv.gz'}])

        self.assertEqual(self.runVerify(), 0)

        with open(self.samplesheet.with_name(
                self.samplesheet.stem + '.ncbi_evidence' + self.samplesheet.suffix),
                newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            fields = reader.fieldnames
            row = next(reader)

        self.assertEqual(fields, list(VERIFY_COLUMNS) + list(NCBI_EVIDENCE_COLUMNS))
        self.assertEqual(row['filename'], 'fasta/x.fasta')
        self.assertEqual(row['species'], 'human')

    def test_default_report_path_moves_the_extension_rather_than_appending_it(self):
        """
        The default path is <stem>.ncbi_evidence<ext>, not
        <stem><ext>.ncbi_evidence.tsv -- one extension, not two.
        """
        self.writeSamplesheet([{'sample_id': 'ssr_1', 'filename': 'x.tsv',
                               'subject_id': 'Donor-2', 'species': 'human',
                               'sample_name': 'Study/csv_paired/SRR1_1_Paired_All.csv.gz'}])

        self.assertEqual(self.runVerify(), 0)

        self.assertTrue((self.tmp / 'samplesheet_airrflow_airr.ncbi_evidence.tsv').exists())
        self.assertFalse(Path(str(self.samplesheet) + '.ncbi_evidence.tsv').exists())

    def test_out_overrides_the_default_report_path(self):
        """--out sends the report somewhere other than the default sidecar path."""
        self.writeSamplesheet([{'sample_id': 'ssr_1', 'filename': 'x.tsv',
                               'subject_id': 'no', 'species': 'human',
                               'sample_name': 'Study/csv_paired/SRR1_1_Paired_All.csv.gz'}])
        out = self.tmp / 'report.tsv'

        self.assertEqual(self.runVerify(['--out', str(out)]), 0)

        self.assertTrue(out.exists())
        self.assertFalse(self.samplesheet.with_name(
            self.samplesheet.stem + '.ncbi_evidence' + self.samplesheet.suffix).exists())
        self.assertEqual(self.readReport(out)[0]['biosample_id'], 'BL-110_VDJ')

    def test_missing_required_column_is_reported_by_name(self):
        """A samplesheet missing a column verify needs names it in the error."""
        with open(self.samplesheet, 'w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=['sample_id', 'filename'],
                                    delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerow({'sample_id': 'ssr_1', 'filename': 'x.tsv'})
        args = getArgParser().parse_args(['oas', 'verify', str(self.samplesheet)])

        with self.assertRaises(Exception) as raised:
            handleOasVerify(args)
        self.assertIn('subject_id', str(raised.exception))
        self.assertIn('sample_name', str(raised.exception))


if __name__ == '__main__':
    unittest.main()
