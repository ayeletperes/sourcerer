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
                    metadata={'Species': 'human', 'Subject': 'Donor-1'},
                    n_sequences=2)

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

    def test_no_format_flag_defaults_to_airr_not_raw_alone(self):
        """
        Omitting --format writes something airrflow can use immediately.

        The previous default (raw alone) left a user with nothing to run
        airrflow on until they reran the command with --format; defaulting
        to airr means the first download already produces a samplesheet.
        The raw mirror is still written either way, since conversion reads
        from it.
        """
        self.assertEqual(self.runDownload(), 0)

        self.assertTrue((self.outdir / 'samplesheet_airrflow_airr.tsv').exists())
        self.assertTrue(list(self.outdir.rglob('raw/**/*.csv.gz')))
        self.assertFalse((self.outdir / 'samplesheet_airrflow_fasta.tsv').exists())
        self.assertEqual(list(self.outdir.glob('fasta/*.fasta')), [])

    def test_unresolved_subjects_are_warned_about(self):
        """
        A batch where OAS recorded no subject at all warns, naming the
        samplesheet to run `sourcerer oas verify` against.

        This is the most likely silent-wrong-analysis outcome download can
        produce: every such row gets the same OAS sentinel as subject_id,
        so airrflow would treat them all as one subject unless the user
        runs verify first.
        """
        unit = DataUnit(unit_id=StubSource.unit.unit_id, collection='paired',
                        url=StubSource.unit.url,
                        metadata={'Species': 'human', 'Subject': 'no'},
                        n_sequences=2)
        with mock.patch.object(StubSource, 'unit', unit):
            with self.assertLogs('sourcerer', level='WARNING') as logs:
                self.assertEqual(self.runDownload(), 0)

        self.assertTrue(any('no subject recorded in OAS' in m for m in logs.output))
        self.assertTrue(any('oas verify' in m for m in logs.output))

    def test_resolved_subjects_are_not_warned_about(self):
        """A batch where every unit already has a real subject stays quiet."""
        with self.assertNoLogs('sourcerer', level='WARNING'):
            self.assertEqual(self.runDownload(), 0)

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

    def test_a_real_subject_id_is_looked_up_too_and_compared(self):
        """
        A row that already has a subject_id is still cross-referenced.

        OAS recording a subject is not proof it is correct -- a typo, a
        short code reused across studies, or a pooled run naming several
        donors under one value are all real failure modes -- so verify
        looks the accession up regardless, and biosample_id always carries
        NCBI's own text rather than a copy of subject_id. Here NCBI's
        BL-110_VDJ has nothing in common with 'Donor-2', so subject_check
        reports 'differs'.
        """
        self.writeSamplesheet([{'sample_id': 'ssr_1', 'filename': 'x.tsv',
                               'subject_id': 'Donor-2', 'species': 'human',
                               'sample_name': 'Study/csv_paired/SRR1_1_Paired_All.csv.gz'}])

        self.assertEqual(self.runVerify(), 0)

        rows = self.readReport()
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['status'], 'ok')
        self.assertEqual(rows[0]['biosample_id'], 'BL-110_VDJ')
        self.assertEqual(rows[0]['biosample_id_suggested'], 'BL-110')
        self.assertEqual(rows[0]['subject_check'], 'differs')

    def test_a_real_subject_id_that_matches_ncbi_agrees(self):
        """subject_check reports 'agrees' when subject_id is NCBI's own text."""
        self.writeSamplesheet([{'sample_id': 'ssr_1', 'filename': 'x.tsv',
                               'subject_id': 'BL-110', 'species': 'human',
                               'sample_name': 'Study/csv_paired/SRR1_1_Paired_All.csv.gz'}])

        self.assertEqual(self.runVerify(), 0)

        rows = self.readReport()
        self.assertEqual(rows[0]['subject_check'], 'agrees')

    def test_a_pooled_subject_id_is_never_looked_up_as_a_single_subject(self):
        """
        OAS's own Subject field can itself name a pool of donors.

        'donor 21; 22; 23 and 24' is exactly the shape OAS's paired catalog
        uses for a 10x hashed/pooled run; subject_check must recognize this
        from subject_id alone, the same way it recognizes NCBI's own pooled
        text, rather than reporting a false 'differs'.
        """
        self.writeSamplesheet([{'sample_id': 'ssr_1', 'filename': 'x.tsv',
                               'subject_id': 'donor 21; 22; 23 and 24',
                               'species': 'human',
                               'sample_name': 'Study/csv_paired/SRR1_1_Paired_All.csv.gz'}])

        self.assertEqual(self.runVerify(), 0)

        rows = self.readReport()
        self.assertEqual(rows[0]['subject_check'], 'pooled')

    def test_no_accession_leaves_subject_check_unresolved_only_when_null(self):
        """
        A row whose sample_name carries no accession, but does have a real
        subject_id, is neither 'unresolved' (that's for a null subject_id)
        nor comparable -- it is 'unverified'.
        """
        self.writeSamplesheet([{'sample_id': 'ssr_1', 'filename': 'x.tsv',
                               'subject_id': 'Donor-2', 'species': 'human',
                               'sample_name': 'not-an-accession'}])

        empty_client = HttpClient(delay=0, backoff=0, session=FakeSession(
            lambda *a: (_ for _ in ()).throw(AssertionError('no network call was expected'))))
        args = getArgParser().parse_args(['oas', 'verify', str(self.samplesheet)])
        with mock.patch('sourcerer.Cli.HttpClient', return_value=empty_client):
            self.assertEqual(handleOasVerify(args), 0)

        rows = self.readReport()
        self.assertEqual(rows[0]['status'], 'no_accession')
        self.assertEqual(rows[0]['subject_check'], 'unverified')
        self.assertEqual(rows[0]['biosample_id'], '')

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
        self.assertEqual(rows[0]['subject_check'], 'unresolved')
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

    def test_same_subject_id_in_two_studies_is_warned_about(self):
        """
        A short subject_id reused across studies is a real collision.

        airrflow keys a subject on subject_id alone, so two studies sharing
        'Donor-2' would otherwise merge silently; this is worth a log
        warning independent of subject_check, which only compares each row
        against NCBI and cannot see across rows.
        """
        columns = list(VERIFY_COLUMNS) + ['study']
        with open(self.samplesheet, 'w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t',
                                    lineterminator='\n')
            writer.writeheader()
            writer.writerow({'sample_id': 'ssr_1', 'filename': 'a.tsv',
                             'subject_id': 'Donor-2', 'species': 'human',
                             'sample_name': 'StudyA/x.csv.gz', 'study': 'StudyA'})
            writer.writerow({'sample_id': 'ssr_2', 'filename': 'b.tsv',
                             'subject_id': 'Donor-2', 'species': 'human',
                             'sample_name': 'StudyB/y.csv.gz', 'study': 'StudyB'})

        with self.assertLogs('sourcerer', level='WARNING') as logs:
            self.assertEqual(self.runVerify(), 0)

        self.assertTrue(any('Donor-2' in message for message in logs.output))
        self.assertTrue(any('StudyA' in message and 'StudyB' in message
                            for message in logs.output))


if __name__ == '__main__':
    unittest.main()
