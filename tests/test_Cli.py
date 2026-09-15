"""
Unit tests for the commandline interface
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import contextlib
import csv
import io
import shutil
import tempfile
import unittest
from argparse import ArgumentParser, Namespace, _SubParsersAction
from pathlib import Path
from unittest import mock

import pandas
import yaml

# Sourcerer imports
from sourcerer import Provenance, Reference, Schema
from sourcerer.Cli import (
    NCBI_EVIDENCE_COLUMNS,
    applyPins,
    formatUnitTable,
    getArgParser,
    handleDownload,
    handleOasVerify,
    handleReferenceDiff,
    handleReferenceDownload,
    handleReferenceShow,
    loadMap,
)
from sourcerer.Exceptions import SourcererError
from sourcerer.Http import HttpClient
from sourcerer.Sources.Base import DataUnit, DownloadResult, Query, SourceBase
from sourcerer.Sources.Imgt import ImgtSource
from sourcerer.Sources.Oas import OasSource, newReport
from sourcerer.Sources.Ogrdb import OgrdbSource
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

    def test_presence_flags_spell_out_their_tokens(self):
        """
        A presence-only flag's usage names the tokens it takes.

        `--subject VALUE` said nothing about what VALUE could be; the user only
        learnt 'defined'/'undefined' from the rejection message.
        """
        parser = getArgParser()
        with mock.patch('sys.stdout', new_callable=io.StringIO) as out:
            with self.assertRaises(SystemExit):
                parser.parse_args(['oas', 'search', 'paired', '--help'])

        text = out.getvalue()
        self.assertIn('--subject {*,defined,undefined}', text)
        self.assertIn('--species VALUE', text)


class TestFormatUnitTable(unittest.TestCase):
    """
    Tests for the search/dry-run stdout table
    """

    def test_shows_metadata_beside_each_hit(self):
        """
        The table carries enough metadata to pick units without --out.

        unit_id and a count alone left the user unable to tell one donor's
        run from another's without first saving a catalog.
        """
        units = [DataUnit(unit_id='A_2020/csv/a.csv.gz', collection='paired',
                          url='u', n_sequences=12,
                          metadata={'Species': 'human', 'Disease': 'None',
                                    'Subject': 'Donor-1', 'BSource': 'PBMC'}),
                 DataUnit(unit_id='B_2021/csv/b.csv.gz', collection='paired',
                          url='u', n_sequences=None, metadata={})]

        lines = formatUnitTable(units).split('\n')

        self.assertEqual(lines[0].split(),
                         ['unit_id', 'n_unique_sequences', 'Species', 'Disease',
                          'Subject', 'BSource'])
        self.assertEqual(lines[1].split(),
                         ['A_2020/csv/a.csv.gz', '12', 'human', 'None',
                          'Donor-1', 'PBMC'])
        # A unit with no metadata still lines up under the same header.
        self.assertTrue(lines[2].startswith('B_2021/csv/b.csv.gz'))
        self.assertEqual(len(lines), 3)


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

    def test_both_formats_convert_each_unit_once(self):
        """
        Requesting two formats reads and normalizes each unit once, not twice.

        Conversion is the expensive step (a multi-GB gzip read plus pandas
        normalization), and each format used to pull its own pass through it.
        Both writers are fed from the same chunk stream now, so airr+fasta
        costs one conversion, and both files still come out complete.
        """
        with mock.patch.object(StubSource, 'convertUnit',
                               autospec=True,
                               side_effect=StubSource.convertUnit) as convert:
            self.assertEqual(self.runDownload('airr', 'fasta'), 0)

        self.assertEqual(convert.call_count, 1)
        airr = next(self.outdir.glob('airr/*.tsv')).read_text().splitlines()
        fasta = next(self.outdir.glob('fasta/*.fasta')).read_text().splitlines()
        self.assertEqual(len(airr), 3)
        self.assertEqual(len(fasta), 4)


class TwoUnitStubSource(SourceBase):
    """
    A source serving two units with distinct conversion reports.

    StubSource above always serves exactly one unit with a fixed report, which
    cannot tell "recorded once" from "summed across the whole run" apart. This
    exercises that: each unit's convertUnit returns different counters, so a
    correct run total has to actually add them rather than just reflect
    whichever unit happened to run last.
    """

    name = 'oas'
    description = 'stub'
    collections = ('paired',)

    units = (
        DataUnit(unit_id='Study_2020/csv/one.csv.gz', collection='paired',
                url='https://example.invalid/one.csv.gz',
                metadata={'Species': 'human', 'Subject': 'Donor-1'},
                n_sequences=2),
        DataUnit(unit_id='Study_2020/csv/two.csv.gz', collection='paired',
                url='https://example.invalid/two.csv.gz',
                metadata={'Species': 'human', 'Subject': 'Donor-2'},
                n_sequences=3),
    )

    def harvestSchema(self):
        raise NotImplementedError

    def searchUnits(self, query):
        return list(self.units)

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

        return DownloadResult(unit=unit, path=path, sha256='0' * 64, size_bytes=3)

    def convertUnit(self, path, unit, chunksize=50000):
        frame = pandas.DataFrame(
            {'sequence_id': ['a'], 'cell_id': ['c1'], 'sequence': ['ACGT'],
             'locus': ['IGH'], 'c_call': ['IGHM']})
        report = newReport()
        if unit.unit_id.endswith('one.csv.gz'):
            report.update(rows_in=2, rows_out=2, missing_c_call=1, loci={'IGH'})
        else:
            report.update(rows_in=3, rows_out=3, missing_c_call=2, loci={'IGK'})

        return {'Species': 'human'}, iter([frame]), report


class TestConversionReportProvenance(unittest.TestCase):
    """
    Tests for the run level conversion_report and schema_fingerprint that
    handleDownload writes into download_metadata.yml
    """

    def setUp(self):
        self.outdir = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.outdir, ignore_errors=True)

    def runDownload(self, *formats):
        argv = ['oas', 'download', 'paired', '--outdir', str(self.outdir)]
        for value in formats:
            argv += ['--format', value]

        args = getArgParser().parse_args(argv)
        args.source = 'oas'

        with mock.patch('sourcerer.Cli.getSource', return_value=TwoUnitStubSource(None)):
            return handleDownload(args)

    def lastRun(self):
        text = (self.outdir / Provenance.DOWNLOAD_METADATA).read_text()
        return yaml.safe_load(text)['runs'][-1]

    def test_conversion_counters_sum_across_every_unit_in_the_run(self):
        """
        A run converting two units must report the run as a whole, not
        whichever unit's counters happened to be computed last.
        """
        self.assertEqual(self.runDownload('airr'), 0)

        report = self.lastRun()['conversion_report']
        self.assertEqual(report['rows_in'], 5)
        self.assertEqual(report['rows_out'], 5)
        self.assertEqual(report['missing_c_call'], 3)
        self.assertEqual(sorted(report['loci']), ['IGH', 'IGK'])

    def test_a_raw_only_run_records_no_conversion_report(self):
        """Nothing was converted, so there is nothing to report."""
        self.assertEqual(self.runDownload('raw'), 0)

        self.assertNotIn('conversion_report', self.lastRun())

    def test_the_run_records_a_schema_fingerprint(self):
        """
        Ties the download to the exact snapshot content it was resolved
        against, not only a harvest date and tool version.
        """
        self.assertEqual(self.runDownload('airr'), 0)

        self.assertEqual(self.lastRun()['schema_fingerprint'],
                         Schema.fingerprint('oas'))


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
        looks the accession up regardless, and ncbi_sample_name always carries
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
        self.assertEqual(rows[0]['ncbi_sample_name'], 'BL-110_VDJ')
        self.assertEqual(rows[0]['ncbi_subject_suggested'], 'BL-110')
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
        self.assertEqual(rows[0]['ncbi_sample_name'], '')

    def test_resolved_row_gets_both_the_raw_name_and_a_suggested_subject(self):
        """
        A resolved row's report names its BioSample, a check link, NCBI's raw
        sample name and the subject suggested from it -- no flag needed to
        choose between them.
        """
        self.writeSamplesheet([{'sample_id': 'ssr_1', 'filename': 'x.tsv',
                               'subject_id': 'no', 'species': 'human',
                               'sample_name': 'Study/csv_paired/SRR1_1_Paired_All.csv.gz'}])

        self.assertEqual(self.runVerify(), 0)

        rows = self.readReport()
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['status'], 'ok')
        self.assertEqual(rows[0]['biosample_accession'], 'SAMN1')
        self.assertEqual(rows[0]['ncbi_sample_name'], 'BL-110_VDJ')
        self.assertEqual(rows[0]['ncbi_subject_suggested'], 'BL-110')
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
        self.assertEqual(self.readReport(out)[0]['ncbi_sample_name'], 'BL-110_VDJ')

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
