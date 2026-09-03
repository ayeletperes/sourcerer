"""
sourcerer commandline interface

Download data from online immune repertoire databases and format it for
Immcantation.
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import csv
import logging
import os
import sys
from argparse import ArgumentParser
from pathlib import Path

# Sourcerer imports
from sourcerer import Catalog, Convert, Ncbi, Provenance, Reference
from sourcerer.Airrflow import buildSamplesheet
from sourcerer.Commandline import CommonHelpFormatter, setupLogging
from sourcerer.Exceptions import SourcererError
from sourcerer.Http import HttpClient
from sourcerer.Schema import loadSchema, saveSchema
from sourcerer.Sources import ALIASES, REGISTRY, canonicalName, getSource
from sourcerer.Sources.Oas import isNull
from sourcerer.Version import __date__, __version__

log = logging.getLogger('sourcerer')

#: Output formats the download and convert subcommands can produce.
FORMATS = ('raw', 'airr', 'fasta')

#: Above this many values, a filter flag's help lists only a sample instead of
#: everything, and points at `schema show` for the rest. Enumerated fields are
#: categorical (species, disease, ...), so in practice this rarely bites; it
#: exists so one field having many values can't blow up every --help screen.
VALUE_LIST_CAP = 20


def loadSchemaQuietly(name):
    """
    Load a packaged snapshot, returning None instead of raising.

    Parser construction must work on a checkout that has no snapshot yet,
    otherwise `sourcerer schema refresh` could never be run to create one.

    Arguments:
      name (str): the source name.

    Returns:
      SourceSchema: the snapshot, or None.
    """
    try:
        return loadSchema(name)
    except Exception:
        return None


def addFilterArgs(parser, schema, source, collection):
    """
    Generate one commandline flag per searchable field.

    Flags come from the stored snapshot, never from a hardcoded list, so the day
    the source adds a field its flag appears with no code change. This reads only
    packaged data and performs no network access: it runs at documentation build
    time as well as at runtime.

    Arguments:
      parser (ArgumentParser): the subparser to add to.
      schema (SourceSchema): the snapshot, or None if none is installed.
      source (str): the source name, for the `schema show` pointer on overflow.
      collection (str): which collection's fields to add.
    """
    if schema is None or collection not in schema.collections:
        return

    for item in schema.getCollection(collection).fields:
        if item.pseudo_values:
            summary = 'filter on whether %s is recorded' % item.name
        elif len(item.values) <= VALUE_LIST_CAP:
            summary = '%d values: %s' % (len(item.values), ', '.join(item.values))
        else:
            # Overflow only: today's fields (species, disease, ...) all stay well
            # under the cap. Once `sourcerer build` (the interactive command
            # builder, see plan phase 6) exists, point there instead.
            shown = ', '.join(item.values[:VALUE_LIST_CAP])
            summary = ('%d values, e.g. %s, ... run `sourcerer schema show '
                       '--source %s --collection %s --field %s` for the full list'
                       % (len(item.values), shown, source, collection, item.name))

        parser.add_argument(item.flag, dest='filter_%s' % item.name,
                            metavar='VALUE', default=None, help=summary)


def collectFilters(args):
    """
    Gather the generated filter flags the user actually supplied.

    Arguments:
      args (Namespace): parsed arguments.

    Returns:
      dict: field name to value.
    """
    return {k[len('filter_'):]: v for k, v in vars(args).items()
            if k.startswith('filter_') and v is not None}


def getArgParser():
    """
    Build the top level argument parser.

    Defined as a function returning an ArgumentParser so that
    sphinxcontrib-autoprogram can document the tool, per Immcantation's
    CONTRIBUTING.md.

    Returns:
      argparse.ArgumentParser: the top level parser.
    """
    parser = ArgumentParser(prog='sourcerer', description=__doc__,
                            formatter_class=CommonHelpFormatter)
    # NB: %(prog)s is expanded by argparse, so it must not be part of the string
    # being %-formatted here.
    parser.add_argument('--version', action='version',
                        version='%(prog)s:' + ' %s %s' % (__version__, __date__))

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-v', '--verbose', action='store_true',
                       help='report debug level status messages')
    group.add_argument('-q', '--quiet', action='store_true',
                       help='report errors only')

    commands = parser.add_subparsers(title='subcommands', dest='command',
                                     metavar='')

    sources = commands.add_parser(
        'sources', help='list available sources',
        description='List every data source sourcerer knows how to fetch '
                    'from, along with a one-line description and its homepage.',
        formatter_class=CommonHelpFormatter)
    sources.add_subparsers(dest='action', metavar='').add_parser(
        'list', help='list the sources sourcerer knows about',
        description='List every data source sourcerer knows how to fetch '
                    'from, along with a one-line description and its homepage.',
        formatter_class=CommonHelpFormatter)

    _addSchemaParser(commands)
    _addReferenceParser(commands)
    for name, source in sorted(REGISTRY.items()):
        _addSourceParser(commands, name, source)

    return parser


def _addReferenceParser(commands):
    """Add the reference subcommand: validate and build from a reference folder."""
    reference = commands.add_parser(
        'reference', help='validate a germline reference folder and build '
                          'IgBLAST databases from it',
        description='Check that a folder of germline FASTAs is in a format '
                    'airrflow can use, and build the IgBLAST databases from it. '
                    'Files are recognised by name, in any directory layout: the '
                    'species and chain, with an optional source prefix and an '
                    'optional aa marker for translated V, as in human_IGHV.fasta '
                    'or imgt_human_IGHV.fasta.',
        formatter_class=CommonHelpFormatter)
    actions = reference.add_subparsers(dest='action', metavar='ACTION',
                                       required=True)

    build = actions.add_parser(
        'build', help='validate a reference folder and build IgBLAST databases',
        description='Validate the reference folder and build the IgBLAST '
                    'databases from it. With --check, only validate and report, '
                    'building nothing (and needing no makeblastdb).',
        formatter_class=CommonHelpFormatter)
    build.add_argument('folder', type=Path,
                       help='a reference_base tree or a flat folder of germline '
                            'FASTAs named <species>_<CHAIN>.fasta')
    build.add_argument('--out', type=Path, default=None,
                       help='directory to write igblast_base into; required '
                            'unless --check')
    build.add_argument('--check', action='store_true',
                       help='validate the folder and report what would build, '
                            'without building anything')
    build.add_argument('--species', nargs='+', choices=list(Reference.SPECIES),
                       default=None,
                       help='limit to these species; default is every species '
                            'found in the folder')


def _addSchemaParser(commands):
    """Add the schema subcommand tree."""
    schema = commands.add_parser(
        'schema', help='inspect and refresh snapshots',
        description='Inspect the schema snapshot checked into the package, '
                    'or re-harvest it from the remote source when the '
                    'upstream API changes (new organisms, fields, values, ...).',
        formatter_class=CommonHelpFormatter)
    actions = schema.add_subparsers(dest='action', metavar='')

    show = actions.add_parser(
        'show', help='print a stored snapshot',
        description='Print the schema snapshot stored for a source: its '
                    'collections and how many fields each has, one '
                    'collection\'s fields and how many values each accepts, '
                    'or every value a single field accepts.',
        formatter_class=CommonHelpFormatter)
    show.add_argument('--source', required=True, choices=sorted(REGISTRY) + sorted(ALIASES),
                      help='which source (or alias) to read the snapshot of')
    show.add_argument('--collection', default=None,
                      help='list this collection\'s fields; without it, print '
                           'one summary line per collection')
    show.add_argument('--field', default=None,
                      help='print every value this field accepts, one per line; '
                           'needs --collection')

    refresh = actions.add_parser(
        'refresh', help='re-harvest a snapshot',
        description='Contact a remote source, re-harvest its schema and '
                    'catalogs, and write the result over the packaged '
                    'snapshot (or alongside it, with --out). Previously '
                    'fetched detail-page enrichment is carried forward '
                    'unless --refresh-details asks to redo it.',
        formatter_class=CommonHelpFormatter)
    refresh.add_argument('--source', required=True, choices=sorted(REGISTRY) + sorted(ALIASES),
                         help='which source to contact and re-harvest')
    refresh.add_argument('--out', default=None, type=Path,
                         help='directory to write into; without it the packaged '
                              'snapshot is rewritten in place')
    refresh.add_argument('--collection', action='append', default=None,
                         help='limit to one collection; repeatable')
    refresh.add_argument('--refresh-details', default='auto',
                         choices=['auto', 'all', 'none'],
                         help='fetch per unit detail pages for fields the search '
                              'results omit; auto fetches only units not already '
                              'read, all re-reads every unit')
    refresh.add_argument('--detail-limit', type=int, default=None,
                         help='stop after this many detail pages')


def _addSourceParser(commands, name, source):
    """Add one source's subcommand tree, with a level per collection."""
    schema = loadSchemaQuietly(name)

    parser = commands.add_parser(
        name, aliases=list(source.aliases), help=source.description,
        description='%s\n\nHomepage: %s' % (source.description, source.homepage),
        formatter_class=CommonHelpFormatter)
    actions = parser.add_subparsers(dest='action', metavar='ACTION',
                                    required=True)

    action_help = {
        'search': ('list matching data units',
                   'Search a collection for data units matching the given '
                   'filters and print (or save with --out) a summary of what '
                   'matched. Nothing is downloaded.'),
        'download': ('download matching data units',
                    'Search a collection for data units matching the given '
                    'filters, download each one, and optionally convert it '
                    'to AIRR and/or FASTA, writing a samplesheet for each '
                    'converted format.'),
    }
    for action, (helptext, description) in action_help.items():
        action_parser = actions.add_parser(action, help=helptext,
                                           description=description,
                                           formatter_class=CommonHelpFormatter)
        # The collection is a subcommand rather than a flag because the two
        # collections have genuinely different field sets. argparse cannot vary
        # options by a flag's value, so a flag would force a union and make
        # --help describe fields that do not apply.
        #
        # NB: a named metavar, not the '' used elsewhere, because argparse names
        # the missing argument by its metavar and an empty one produces a
        # required-argument error that says nothing.
        collections = action_parser.add_subparsers(dest='collection',
                                                   metavar='COLLECTION',
                                                   required=True)
        for collection in source.collections:
            # Passing help is what makes argparse list the collection at all.
            collection_help = source.collection_help.get(collection)
            leaf = collections.add_parser(
                collection, help=collection_help,
                description='%s the %s collection (%s), optionally narrowed '
                            'down with the filter flags below.'
                            % (action.capitalize(), collection, collection_help),
                formatter_class=CommonHelpFormatter)
            addFilterArgs(leaf, schema, name, collection)
            leaf.add_argument('--limit', type=int, default=None,
                              help='stop after this many units')
            if action == 'search':
                leaf.add_argument('-o', '--out', type=Path, default=None,
                                  help='write the hits to a TSV file')
            else:
                leaf.add_argument('--outdir', type=Path, required=True,
                                  help='directory to write into')
                leaf.add_argument('--dry-run', action='store_true',
                                  help='report what would be fetched, then stop')
                leaf.add_argument('--no-resume', action='store_true',
                                  help='re-download in full rather than '
                                       'continuing a partly fetched file')
                if source.output == 'reference':
                    # Reference sources build a germline reference_base rather
                    # than converting to AIRR, so they take the igblast options
                    # instead of the format and AIRR-strictness ones.
                    leaf.add_argument('--igblast', action='store_true',
                                      help='also build the IgBLAST databases '
                                           'from the reference; needs makeblastdb '
                                           'on PATH')
                    leaf.add_argument('--igblast-out', type=Path, default=None,
                                      help='where to write igblast_base; '
                                           'defaults to <outdir>/igblast_base')
                else:
                    leaf.add_argument('--format', action='append', dest='formats',
                                      choices=FORMATS,
                                      help='what to write, repeatable to write '
                                           'several; raw mirrors the source files '
                                           'untouched and is always written because '
                                           'the others are converted from it, so '
                                           'omitting this writes raw alone')
                    leaf.add_argument('--strict-airr', action='store_true',
                                      help='drop columns the AIRR schema does '
                                           'not define')

    if name == 'oas':
        _addOasVerifyAction(actions)


def _addOasVerifyAction(actions):
    """
    Add oas verify: cross-reference a samplesheet's unresolved subjects
    against NCBI.

    Not a search/download action: it takes no COLLECTION or filter flags,
    since it reads a samplesheet already on disk rather than querying OAS.
    Every row with subject_id 'no' or 'None' (OAS's own null sentinels; see
    Sources.Oas.isNull) names an SRA run or GEO sample accession in its
    sample_name column, and that accession's NCBI BioSample record usually
    names the sample plainly enough for a human to read off the subject.

    Deliberately two flags, not the larger surface an earlier version of this
    command had (--apply, --use-suggestion, --evidence-out, --limit): the
    report is the one thing this command produces, always with both a raw and
    a suggested biosample_id, so there was nothing left for a flag to switch
    between. It carries every column the input samplesheet had, evidence
    columns appended, so it is a drop-in airrflow input rather than a
    side file -- see buildEvidenceRow and NCBI_EVIDENCE_COLUMNS.
    """
    verify = actions.add_parser(
        'verify', help='cross-reference unresolved subjects against NCBI',
        description='Write an evidence report with one row per samplesheet '
                    'row, every input column carried through unchanged, plus '
                    'subject_id copied through unchanged for a row that '
                    'already had one, or NCBI\'s BioSample record for a row '
                    'OAS left unresolved (subject_id \'no\' or \'None\'), '
                    'looked up from the run/sample accession in its '
                    'sample_name column. biosample_id carries NCBI\'s raw '
                    'sample name -- never guessed at further, since some '
                    'studies\' names need study-specific reading to turn '
                    'into a subject id (see the module docstring in '
                    'Ncbi.py); biosample_id_suggested carries the same value '
                    'with the handful of generic patterns (a trailing locus '
                    'or visit suffix) stripped, the ones safe to normalize '
                    'regardless of study. A pooled/multi-donor run (10x cell '
                    'hashing, e.g.) gets an AMBIGUOUS_POOLED marker naming '
                    'the donor codes in both columns instead of a guessed '
                    'single subject. Because every input column survives, '
                    'the report can be pointed at directly as airrflow '
                    '--input once subject_id (still passed through raw) is '
                    'filled in for any row that needs it.',
        formatter_class=CommonHelpFormatter)
    verify.add_argument('samplesheet', type=Path,
                        help='an airrflow samplesheet to read; only needs '
                             'sample_id, sample_name and subject_id columns, '
                             'so a hand-edited sheet is fine too, but every '
                             'column it has is carried through to the report'
                             )
    verify.add_argument('--out', type=Path, default=None,
                        help='where to write the evidence report; defaults '
                             'to <samplesheet-name>.ncbi_evidence<ext>, e.g. '
                             'samplesheet_airrflow_fasta.ncbi_evidence.tsv '
                             'for samplesheet_airrflow_fasta.tsv')
    verify.add_argument('--ncbi-api-key', default=None,
                        help='an NCBI API key, raising the polite request '
                             'rate from 3/s to 10/s; falls back to the '
                             'NCBI_API_KEY environment variable')


def makeClient(args):
    """Build the shared HTTP client."""
    return HttpClient()


def handleSources(args):
    """List registered sources."""
    for name, source in sorted(REGISTRY.items()):
        print('%-10s %s' % (name, source.description))
        if source.aliases:
            print('%-10s alias: %s' % ('', ', '.join(source.aliases)))
        print('%-10s %s' % ('', source.homepage))
        if source.license:
            print('%-10s license: %s' % ('', source.license))
        for paper in source.citation:
            print('%-10s cite: %s' % ('', paper))

    return 0


def handleSchemaShow(args):
    """Print a stored snapshot."""
    args.source = canonicalName(args.source)
    schema = loadSchema(args.source)

    if args.collection is None:
        print('source:    %s' % schema.source)
        print('harvested: %s by %s' % (schema.harvested, schema.harvested_by))
        for name in schema.collection_names:
            collection = schema.getCollection(name)
            print('  %s: %d fields (%s)'
                  % (name, len(collection.fields),
                     ', '.join(collection.field_names)))
        return 0

    collection = schema.getCollection(args.collection)
    if args.field is None:
        for item in collection.fields:
            kind = 'presence only' if item.pseudo_values else '%d values' % len(item.values)
            print('%-14s %s' % (item.name, kind))
        return 0

    item = collection.getField(args.field)
    if item is None:
        raise SourcererError("no field '%s' in %s %s"
                             % (args.field, args.source, args.collection))
    for value in item.values:
        print(value)

    return 0


def handleSchemaRefresh(args):
    """Re-harvest a snapshot and its catalogs."""
    args.source = canonicalName(args.source)
    client = makeClient(args)
    source = getSource(args.source, client)

    log.info('harvesting %s search schema', args.source)
    schema = source.harvestSchema()

    out = args.out
    if out is None:
        from importlib import resources
        out = Path(str(resources.files('sourcerer').joinpath(
            'data/schemas', args.source)))

    written, changed = saveSchema(schema, out)
    log.info('%s %s', 'wrote' if changed else 'unchanged, left alone:', written)

    wanted = args.collection or list(source.collections)
    for collection in wanted:
        log.info('harvesting %s %s catalog', args.source, collection)
        rows = source.harvestCatalog(collection, schema=schema)

        path = out / ('%s_catalog.tsv' % collection)
        rows = Catalog.mergeEnrichment(Catalog.loadCatalog(path), rows)

        if args.refresh_details != 'none':
            force = args.refresh_details == 'all'
            pending = rows if force else [x for x in rows if Catalog.needsDetail(x)]
            if pending:
                log.info('enriching %d %s units from detail pages',
                         len(pending), collection)
                source.enrichCatalog(rows, limit=args.detail_limit, force=force)

        Catalog.saveCatalog(rows, path)
        log.info('wrote %s (%d units)', path, len(rows))

    return 0


def handleSearch(args):
    """List data units matching a query."""
    client = makeClient(args)
    source = getSource(args.source, client)

    query = source.validateQuery(args.collection, collectFilters(args))
    if args.limit is not None:
        query = type(query)(collection=query.collection, filters=query.filters,
                            limit=args.limit)

    units = source.searchUnits(query)
    total = sum(x.n_sequences or 0 for x in units)
    log.info('%d data units, %s sequences', len(units), format(total, ','))

    rows = [{'unit_id': x.unit_id, 'collection': x.collection,
             'n_unique_sequences': x.n_sequences or '', 'url': x.url,
             **{k: v for k, v in x.metadata.items()
                if k in Catalog.CATALOG_COLUMNS}}
            for x in units]

    if args.out is not None:
        Catalog.saveCatalog(rows, args.out)
        log.info('wrote %s', args.out)
    else:
        for unit in units:
            print('%-64s %10s' % (unit.unit_id, unit.n_sequences or ''))

    return 0


def handleReference(args):
    """Validate a reference folder and, unless --check, build its IgBLAST base."""
    if not args.folder.is_dir():
        raise SourcererError('no such reference folder: %s' % args.folder)

    plan = Reference.planReference(args.folder, species=args.species)
    print(plan.summary())

    if not plan.ok:
        raise SourcererError('no databases can be built from %s; check the file '
                             'names against <species>_<CHAIN>.fasta' % args.folder)

    if args.check:
        return 0

    if args.out is None:
        raise SourcererError('--out is required to build; pass --check to only '
                             'validate the folder')

    Reference.buildFromPlan(plan, args.out, makeClient(args))
    log.info('wrote %s', args.out)

    return 0


def handleReferenceDownload(args, source):
    """Download germline sets and build an airrflow reference_base."""
    query = source.validateQuery(args.collection, collectFilters(args))
    if args.limit is not None:
        query = type(query)(collection=query.collection, filters=query.filters,
                            limit=args.limit)

    units = source.searchUnits(query)
    log.info('%d germline files for %s %s',
             len(units), args.source, args.collection)

    if args.dry_run:
        for unit in units:
            print('%-48s %s' % (unit.unit_id, unit.url))
        log.info('dry run: nothing downloaded')
        return 0

    outdir = Path(args.outdir)
    raw_dir = outdir / 'raw'

    entries, provenance = [], []
    for unit in units:
        result = source.fetchUnit(unit, raw_dir, resume=not args.no_resume)
        entries.append((unit, result.path))
        provenance.append(Provenance.buildUnitRecord(unit, result, outdir, {}))

    reference_dir = outdir / 'reference_base'
    source.buildReference(entries, reference_dir).logSummary()
    log.info('wrote %s', reference_dir)

    formats = ['reference']
    if args.igblast:
        igblast_out = args.igblast_out or (outdir / 'igblast_base')
        Reference.buildIgblastBase(reference_dir, igblast_out, source.client,
                                   species=[args.collection])
        log.info('wrote %s', igblast_out)
        formats.append('igblast')

    record = Provenance.writeDownloadMetadata(
        outdir, args.source, args.collection, collectFilters(args), args.limit,
        formats, provenance, schema=source.schema, license=source.license,
        citation=source.citation)
    log.info('wrote %s', record)

    return 0


def handleDownload(args):
    """Download and optionally convert matching data units."""
    client = makeClient(args)
    source = getSource(args.source, client)

    # Reference sources (germline sets) build a reference_base instead of
    # converting repertoires to AIRR and writing a samplesheet.
    if source.output == 'reference':
        return handleReferenceDownload(args, source)

    query = source.validateQuery(args.collection, collectFilters(args))
    if args.limit is not None:
        query = type(query)(collection=query.collection, filters=query.filters,
                            limit=args.limit)

    units = source.searchUnits(query)
    formats = args.formats or ['raw']
    total = sum(x.n_sequences or 0 for x in units)

    log.info('%d data units, %s sequences, formats: %s',
             len(units), format(total, ','), ', '.join(formats))

    if args.dry_run:
        for unit in units:
            print('%-64s %10s' % (unit.unit_id, unit.n_sequences or ''))
        log.info('dry run: nothing downloaded')
        return 0

    outdir = Path(args.outdir)
    raw_dir = outdir / 'raw'
    # The raw mirror is always written, whether or not it was requested, because
    # converting reads from it. So 'raw' always has a bucket here even when the
    # user asked only for airr or fasta.
    written = {x: [] for x in set(formats) | {'raw'}}
    loci = {}
    provenance = []

    for unit in units:
        result = source.fetchUnit(unit, raw_dir, resume=not args.no_resume)
        written['raw'].append((unit, result.path))
        outputs = {}

        stem = unit.unit_id.replace('/', '_').replace('.csv.gz', '')
        if 'airr' in formats:
            _, chunks, report = source.convertUnit(result.path, unit)
            dest = outdir / 'airr' / ('%s.tsv' % stem)
            validation = Convert.writeAirr(chunks, dest, strict=args.strict_airr)
            Convert.writeValidationReport(validation, dest)
            log.info('%s: %d rows, %d invalid, %d rows in',
                     dest.name, validation['rows_checked'],
                     validation['rows_invalid'], report['rows_in'])
            loci[unit.unit_id] = report['loci']
            written['airr'].append((unit, dest))
            outputs['airr'] = dest

        if 'fasta' in formats:
            _, chunks, report = source.convertUnit(result.path, unit)
            dest = outdir / 'fasta' / ('%s.fasta' % stem)
            Convert.writeFasta(chunks, dest)
            loci.setdefault(unit.unit_id, report['loci'])
            written['fasta'].append((unit, dest))
            outputs['fasta'] = dest

        provenance.append(
            Provenance.buildUnitRecord(unit, result, outdir, outputs))

    # A samplesheet is a derived artifact of a data format, so one is written per
    # converted format rather than one ambiguous sheet naming a single file.
    for fmt in ('airr', 'fasta'):
        if written.get(fmt):
            sheet = outdir / ('samplesheet_airrflow_%s.tsv' % fmt)
            buildSamplesheet(written[fmt], sheet, args.collection, outdir,
                             loci=loci)
            log.info('wrote %s', sheet)

    # Written for every run, including raw-only ones: the raw mirror is the part
    # of the output that cannot be regenerated from anything else here.
    record = Provenance.writeDownloadMetadata(
        outdir, args.source, args.collection, collectFilters(args), args.limit,
        formats, provenance, schema=source.schema, license=source.license,
        citation=source.citation)
    log.info('wrote %s', record)

    return 0


#: The columns `sourcerer oas verify` adds, appended after whatever columns
#: the input samplesheet already had (never inserted among them) -- so the
#: report is the input samplesheet plus evidence, and can be used in its
#: place as airrflow input, rather than a separate file missing the columns
#: airrflow actually reads (`filename`, `species`, `pcr_target_locus`, ...).
NCBI_EVIDENCE_COLUMNS = ('biosample_id', 'biosample_id_suggested', 'status',
                         'accession', 'biosample_accession', 'biosample_url',
                         'pooled_codes')


def readSamplesheetRows(path):
    """
    Read a samplesheet leniently, for verify rather than for the download merge.

    Unlike Airrflow.loadSamplesheet, this accepts any TSV that carries the
    columns verify actually needs, in any order, alongside whatever else a
    hand-edited sheet has picked up. Every column present is kept, not just
    the three verify reads, so the row it came from can be written back out
    whole.

    Arguments:
      path (Path): the samplesheet to read.

    Returns:
      tuple: (fields (list of str), rows (list of dict)), in file order.

    Raises:
      SourcererError: if a required column is missing.
    """
    with open(path, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = reader.fieldnames or []
        missing = {'sample_id', 'sample_name', 'subject_id'} - set(fields)
        if missing:
            raise SourcererError(
                '%s is missing column(s) %s that oas verify needs'
                % (path, ', '.join(sorted(missing))))
        return list(fields), [dict(row) for row in reader]


def buildEvidenceRow(row, found):
    """
    Build one row of the verify evidence report.

    Every column already on `row` passes through verbatim -- this only adds
    NCBI_EVIDENCE_COLUMNS, it never edits or drops what was already on the
    samplesheet. A row that already has a real subject_id is a straight
    passthrough for those too: there is nothing to look up, and 'use that
    one' -- the whole reason biosample_id exists -- means biosample_id and
    biosample_id_suggested both simply are it. A pooled/multi-donor run gets
    the same AMBIGUOUS_POOLED marker in both columns, since there is no
    single subject to suggest either. Anything else takes biosample_id from
    NCBI's raw sample name and biosample_id_suggested from the
    generic-heuristic strip of it (see Ncbi.suggestSubject) -- both always
    written, so joining this report back onto a samplesheet never needs a
    flag to decide which one it gets.

    Arguments:
      row (dict): the samplesheet row.
      found (Ncbi.Evidence): the NCBI lookup for this row's accession, or
        None if subject_id was already real (no lookup was attempted) or no
        accession could be read from sample_name.

    Returns:
      dict: row, plus NCBI_EVIDENCE_COLUMNS.
    """
    subject_id = row.get('subject_id', '')

    if not isNull(subject_id):
        evidence = {'biosample_id': subject_id, 'biosample_id_suggested': subject_id,
                   'status': 'passthrough', 'accession': '', 'biosample_accession': '',
                   'biosample_url': '', 'pooled_codes': ''}
    elif found is None:
        evidence = {'biosample_id': '', 'biosample_id_suggested': '',
                   'status': 'no_accession', 'accession': '', 'biosample_accession': '',
                   'biosample_url': '', 'pooled_codes': ''}
    else:
        pooled_codes = ';'.join(found.pooled_codes)
        if found.status == 'pooled':
            biosample_id = biosample_id_suggested = 'AMBIGUOUS_POOLED:%s' % pooled_codes
        elif found.status == 'ok':
            biosample_id = found.sample_name
            biosample_id_suggested = found.suggested_subject
        else:
            biosample_id = biosample_id_suggested = ''
        evidence = {'biosample_id': biosample_id,
                   'biosample_id_suggested': biosample_id_suggested,
                   'status': found.status, 'accession': found.accession,
                   'biosample_accession': found.biosample_accession,
                   'biosample_url': found.url, 'pooled_codes': pooled_codes}

    return {**row, **evidence}


def handleOasVerify(args):
    """Cross-reference a samplesheet's unresolved subjects against NCBI."""
    fields, rows = readSamplesheetRows(args.samplesheet)
    pending = [row for row in rows if isNull(row.get('subject_id'))]

    accession_by_sample = {}
    for row in pending:
        accession = Ncbi.accessionFromText(row.get('sample_name', ''))
        if accession is not None:
            accession_by_sample[row['sample_id']] = accession

    accessions = set(accession_by_sample.values())
    api_key = args.ncbi_api_key or os.environ.get('NCBI_API_KEY')
    delay = Ncbi.KEYED_DELAY if api_key else Ncbi.DEFAULT_DELAY
    client = HttpClient(delay=delay)
    evidence = (Ncbi.gatherEvidence(client, accessions, api_key=api_key)
               if accessions else {})

    counts = {}
    evidence_rows = []
    for row in rows:
        accession = accession_by_sample.get(row['sample_id'])
        found = evidence.get(accession) if accession is not None else None
        evidence_row = buildEvidenceRow(row, found)
        counts[evidence_row['status']] = counts.get(evidence_row['status'], 0) + 1
        evidence_rows.append(evidence_row)

    # Input columns first, in their own order, then whichever evidence columns
    # were not already among them -- so a samplesheet round-tripped through
    # verify keeps every field it walked in with.
    out_fields = fields + [c for c in NCBI_EVIDENCE_COLUMNS if c not in fields]
    # <name>.ncbi_evidence<ext>, not <name><ext>.ncbi_evidence.tsv: the input's
    # own extension moves after 'ncbi_evidence' rather than getting a second
    # one appended after it.
    out = args.out or args.samplesheet.with_name(
        args.samplesheet.stem + '.ncbi_evidence' + args.samplesheet.suffix)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=out_fields,
                                delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(evidence_rows)

    log.info('%d rows: %s', len(rows),
             ', '.join('%d %s' % (n, status) for status, n in sorted(counts.items())))
    log.info('wrote %s', out)

    return 0


def main():
    """
    Parse the commandline and dispatch to the selected subcommand.

    Returns:
      int: process exit status.
    """
    parser = getArgParser()
    args = parser.parse_args()

    setupLogging(verbose=args.verbose, quiet=args.quiet)

    if args.command is None:
        parser.print_help(sys.stderr)
        return 1

    try:
        if args.command == 'sources':
            return handleSources(args)

        if args.command == 'schema':
            if args.action == 'show':
                return handleSchemaShow(args)
            if args.action == 'refresh':
                return handleSchemaRefresh(args)
            parser.parse_args([args.command, '--help'])

        if args.command == 'reference':
            return handleReference(args)

        source_name = canonicalName(args.command) if args.command else None
        if source_name in REGISTRY:
            # The action and collection levels are required subparsers, so
            # argparse has already rejected a commandline missing either. The
            # command may be an alias (e.g. 'airrc'); resolve it to the canonical
            # source so schema and provenance use one name.
            args.source = source_name
            if args.action == 'search':
                return handleSearch(args)
            if args.action == 'download':
                return handleDownload(args)
            if args.action == 'verify':
                return handleOasVerify(args)

        parser.print_help(sys.stderr)
        return 1
    except SourcererError as error:
        log.error('%s', error)
        return 1


if __name__ == '__main__':
    sys.exit(main())
