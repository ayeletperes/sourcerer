"""
sourcerer commandline interface

Download data from online immune repertoire databases and format it for
Immcantation.
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import logging
import sys
from argparse import ArgumentParser, Namespace
from pathlib import Path

# Sourcerer imports
from sourcerer import Catalog, Convert, Provenance
from sourcerer.Airrflow import buildSamplesheet
from sourcerer.Commandline import CommonHelpFormatter, setupLogging
from sourcerer.Exceptions import SourcererError
from sourcerer.Http import HttpClient
from sourcerer.Schema import loadSchema, saveSchema
from sourcerer.Sources import REGISTRY, getSource
from sourcerer.Sources.Specificity import REGISTRY as SPECIFICITY_REGISTRY
from sourcerer.Sources.Specificity import getSpecificitySource
from sourcerer.Sources.Specificity.AnnotateCli import (addAnnotateParser,
                                                       handleAnnotate)
from sourcerer.Sources.Specificity.BcrCli import addBcrParser, handleBcr
from sourcerer.Sources.Specificity.Provenance import writeSpecificityMetadata
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
    for name, source in sorted(REGISTRY.items()):
        _addSourceParser(commands, name, source)
    _addSpecificityGroupParser(commands)
    addAnnotateParser(commands)

    return parser


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
    show.add_argument('--source', required=True, choices=sorted(REGISTRY),
                      help='which source to read the snapshot of')
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
    refresh.add_argument('--source', required=True, choices=sorted(REGISTRY),
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
        name, help=source.description,
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
                                  help='directory to write into; one '
                                       'subdirectory per format, plus a '
                                       'samplesheet for each converted format')
                leaf.add_argument('--format', action='append', dest='formats',
                                  choices=FORMATS,
                                  help='what to write, repeatable to write '
                                       'several; raw mirrors the source files '
                                       'untouched and is always written because '
                                       'the others are converted from it, so '
                                       'omitting this writes raw alone')
                leaf.add_argument('--dry-run', action='store_true',
                                  help='report what would be fetched, then stop')
                leaf.add_argument('--no-resume', action='store_true',
                                  help='re-download in full rather than '
                                       'continuing a partly fetched file')
                leaf.add_argument('--strict-airr', action='store_true',
                                  help='drop columns the AIRR schema does not define')


def _addSpecificityParser(commands, name, source):
    """
    Add one specificity database's subcommand tree, with a level per table.

    A trimmed copy of `_addSourceParser`: no --format/samplesheet flags, since
    those are repertoire specific, and one extra synthetic 'all' table that
    loops every table the database offers.
    """
    schema = loadSchemaQuietly(name)

    parser = commands.add_parser(
        name, help=source.description,
        description='%s\n\nHomepage: %s' % (source.description, source.homepage),
        formatter_class=CommonHelpFormatter)
    actions = parser.add_subparsers(dest='action', metavar='ACTION',
                                    required=True)

    action_help = {
        'search': ('list what a table contains',
                   'Resolve a table, optionally narrowed by filters, and '
                   'print (or save with --out) a summary of it. Nothing is '
                   'downloaded.'),
        'download': ('download a table',
                    'Fetch a table, optionally narrowed by filters, and write '
                    'it as a raw mirror plus one normalized TSV: a real AIRR '
                    'rearrangement TSV for a table listed in the database\'s '
                    'airr_collections, a plain TSV of the source\'s own '
                    'columns otherwise.'),
    }
    for action, (helptext, description) in action_help.items():
        action_parser = actions.add_parser(action, help=helptext,
                                           description=description,
                                           formatter_class=CommonHelpFormatter)
        tables = action_parser.add_subparsers(dest='table', metavar='TABLE',
                                              required=True)

        entries = [(table, source.collection_help.get(table), table)
                   for table in source.collections]
        entries.append(('all', 'every table', None))

        for table_name, table_help, filter_table in entries:
            if filter_table is None:
                summary = 'every table this database offers'
            else:
                summary = 'the %s table (%s), optionally narrowed down with the filter flags below' % (table_name, table_help)
            leaf = tables.add_parser(
                table_name, help=table_help,
                description='%s %s.' % (action.capitalize(), summary),
                formatter_class=CommonHelpFormatter)
            if filter_table is not None:
                addFilterArgs(leaf, schema, name, filter_table)
            leaf.add_argument('--limit', type=int, default=None,
                              help='stop after this many rows'
                                   if filter_table is not None
                                   else 'stop after this many tables')
            if action == 'search':
                leaf.add_argument('-o', '--out', type=Path, default=None,
                                  help='write the result to a TSV file')
            else:
                leaf.add_argument('--outdir', type=Path, required=True,
                                  help='directory to write into: raw/ '
                                       'mirrors the upstream file(s), '
                                       'specificity/ holds the normalized '
                                       'TSV(s)')
                leaf.add_argument('--dry-run', action='store_true',
                                  help='report what would be fetched, then stop')
                leaf.add_argument('--no-resume', action='store_true',
                                  help='re-download in full rather than '
                                       'continuing a partly fetched file')
                leaf.add_argument('--strict-airr', action='store_true',
                                  help='for a table written as AIRR, drop '
                                       'columns the AIRR schema does not '
                                       'define; no effect on other tables')


def _addSpecificityGroupParser(commands):
    """
    Add the `specificity` subcommand: `list`, one per registered database, and
    a synthetic `all` that fans out to every registered database.
    """
    group = commands.add_parser(
        'specificity', help='search and download specificity databases',
        description='Search and download specificity databases (reference '
                    'tables of epitopes, assays and receptor sequences), as '
                    'opposed to repertoire sequencing sources.',
        formatter_class=CommonHelpFormatter)
    dbs = group.add_subparsers(dest='db', metavar='DB', required=True)

    dbs.add_parser(
        'list', help='list the specificity databases sourcerer knows about',
        description='List every specificity database sourcerer knows how to '
                    'fetch from, along with a one-line description and its '
                    'homepage.',
        formatter_class=CommonHelpFormatter)

    addBcrParser(dbs)

    for name, source in sorted(SPECIFICITY_REGISTRY.items()):
        _addSpecificityParser(dbs, name, source)

    all_group = dbs.add_parser(
        'all', help='every registered specificity database',
        description='Act on every registered specificity database at once.',
        formatter_class=CommonHelpFormatter)
    all_actions = all_group.add_subparsers(dest='action', metavar='ACTION',
                                           required=True)
    all_download = all_actions.add_parser(
        'download', help='download every table of every database',
        description='Download every table of every registered specificity '
                    'database. Filters cannot be applied here, since '
                    'different databases have unrelated fields; target one '
                    'database directly to filter it.',
        formatter_class=CommonHelpFormatter)
    all_download.add_argument('--outdir', type=Path, required=True,
                              help='directory to write into; one '
                                   'subdirectory per database')
    all_download.add_argument('--limit', type=int, default=None,
                              help='stop after this many tables, per database')
    all_download.add_argument('--dry-run', action='store_true',
                              help='report what would be fetched, then stop')
    all_download.add_argument('--strict-airr', action='store_true',
                              help='for a table written as AIRR, drop columns '
                                   'the AIRR schema does not define; no effect '
                                   'on other tables')


def makeClient(args):
    """Build the shared HTTP client."""
    return HttpClient()


def handleSources(args):
    """List registered sources."""
    for name, source in sorted(REGISTRY.items()):
        print('%-10s %s' % (name, source.description))
        print('%-10s %s' % ('', source.homepage))
        if source.license:
            print('%-10s license: %s' % ('', source.license))
        for paper in source.citation:
            print('%-10s cite: %s' % ('', paper))

    return 0


def handleSchemaShow(args):
    """Print a stored snapshot."""
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


def handleDownload(args):
    """Download and optionally convert matching data units."""
    client = makeClient(args)
    source = getSource(args.source, client)

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


def _resolveSpecificityTables(source, args):
    """
    Return the tables a specificity action should act on.

    Arguments:
      source (SourceBase): the specificity source.
      args (Namespace): parsed arguments, carrying `table` and `limit`.

    Returns:
      list: table names. `--limit` caps how many are used when `table` is
      'all'; for a single named table it has no effect, since that table
      alone is already the whole selection.
    """
    if args.table != 'all':
        return [args.table]

    tables = list(source.collections)

    return tables[:args.limit] if args.limit is not None else tables


def _writeSpecificityTsv(chunks, out):
    """
    Write normalized specificity records as a TSV, streaming chunk by chunk.

    Unlike `Convert.writeAirr`, this performs no AIRR schema validation:
    specificity columns are a source's own, not AIRR rearrangement fields.

    Arguments:
      chunks (iterable): DataFrames of normalized records.
      out (Path): output path.

    Returns:
      int: rows written.
    """
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)

    written = 0
    header_written = False
    with open(out, 'w', newline='') as handle:
        for frame in chunks:
            frame.to_csv(handle, sep='\t', index=False, header=not header_written)
            header_written = True
            written += len(frame)

    if not header_written:
        out.write_text('')

    return written


def handleSpecificityList(args):
    """List registered specificity databases."""
    for name, source in sorted(SPECIFICITY_REGISTRY.items()):
        print('%-10s %s' % (name, source.description))
        print('%-10s %s' % ('', source.homepage))
        if source.license:
            print('%-10s license: %s' % ('', source.license))
        for paper in source.citation:
            print('%-10s cite: %s' % ('', paper))

    return 0


def handleSpecificitySearch(args):
    """List what a table (or every table) contains, without downloading."""
    client = makeClient(args)
    source = getSpecificitySource(args.db, client)
    filters = collectFilters(args)
    tables = _resolveSpecificityTables(source, args)

    rows = []
    for table in tables:
        query = source.validateQuery(table, filters)
        for unit in source.searchUnits(query):
            rows.append({'unit_id': unit.unit_id, 'collection': unit.collection,
                        'url': unit.url})

    log.info('%d table(s)', len(rows))
    if args.out is not None:
        Catalog.saveCatalog(rows, args.out,
                            columns=('unit_id', 'collection', 'url'))
        log.info('wrote %s', args.out)
    else:
        for row in rows:
            print('%-16s %-24s %s' % (row['collection'], row['unit_id'], row['url']))

    return 0


def handleSpecificityDownload(args):
    """Download one or more tables: a raw mirror plus one normalized TSV each."""
    client = makeClient(args)
    source = getSpecificitySource(args.db, client)
    filters = collectFilters(args)
    tables = _resolveSpecificityTables(source, args)

    outdir = Path(args.outdir)
    raw_dir = outdir / 'raw'
    normalized_dir = outdir / 'specificity'

    records = []
    formats_written = set()
    for table in tables:
        query = source.validateQuery(table, filters)
        units = source.searchUnits(query)

        if args.dry_run:
            for unit in units:
                print('%-16s %s' % (table, unit.unit_id))
            continue

        for unit in units:
            result = source.fetchUnit(unit, raw_dir, resume=not args.no_resume)
            _, chunks, report = source.convertUnit(result.path, unit)

            dest = normalized_dir / ('%s.tsv' % table)
            if table in source.airr_collections:
                validation = Convert.writeAirr(chunks, dest,
                                               strict=args.strict_airr)
                Convert.writeValidationReport(validation, dest)
                log.info('%s: %d rows, %d invalid, %d rows in', dest.name,
                         validation['rows_checked'], validation['rows_invalid'],
                         report['rows_in'])
                output_format = 'airr'
            else:
                rows = _writeSpecificityTsv(chunks, dest)
                log.info('%s: %d rows in, %d rows written', dest.name,
                         report['rows_in'], rows)
                output_format = 'tsv'
            formats_written.add(output_format)

            records.append(Provenance.buildUnitRecord(
                unit, result, outdir, {output_format: dest}))

    if args.dry_run:
        log.info('dry run: nothing downloaded')
        return 0

    # A specificity download represents the database's current state, not one
    # more addition to a growing set the way OAS downloads are, so this always
    # replaces the metadata file rather than accumulating old runs into it --
    # see Specificity.Provenance for why that is a separate writer.
    writeSpecificityMetadata(
        outdir, args.db, args.table, filters, args.limit, records,
        schema=source.schema, license=source.license, citation=source.citation,
        formats=formats_written or ('tsv',))

    return 0


def handleSpecificityDownloadAll(args):
    """Download every table of every registered specificity database."""
    outdir = Path(args.outdir)

    for name in sorted(SPECIFICITY_REGISTRY):
        log.info('downloading specificity source %s', name)
        handleSpecificityDownload(Namespace(
            db=name, table='all', outdir=outdir / name, limit=args.limit,
            dry_run=args.dry_run, no_resume=False,
            strict_airr=args.strict_airr))

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

        if args.command in REGISTRY:
            # The action and collection levels are required subparsers, so
            # argparse has already rejected a commandline missing either.
            args.source = args.command
            if args.action == 'search':
                return handleSearch(args)
            if args.action == 'download':
                return handleDownload(args)

        if args.command == 'specificity':
            if args.db == 'list':
                return handleSpecificityList(args)
            if args.db == 'all':
                return handleSpecificityDownloadAll(args)
            if args.db == 'bcr':
                return handleBcr(args)
            # The action and table levels are required subparsers, so argparse
            # has already rejected a commandline missing either.
            if args.action == 'search':
                return handleSpecificitySearch(args)
            if args.action == 'download':
                return handleSpecificityDownload(args)

        if args.command == 'annotate':
            return handleAnnotate(args)

        parser.print_help(sys.stderr)
        return 1
    except SourcererError as error:
        log.error('%s', error)
        return 1


if __name__ == '__main__':
    sys.exit(main())
