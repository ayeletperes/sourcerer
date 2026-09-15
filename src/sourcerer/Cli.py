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
from datetime import UTC, datetime
from pathlib import Path

# Sourcerer imports
from sourcerer import Catalog, Contracts, Convert, Ncbi, Provenance, Reference
from sourcerer.Airrflow import buildSamplesheet, countUnresolvedSubjects
from sourcerer.Commandline import CommonHelpFormatter, setupLogging
from sourcerer.Exceptions import SourcererError
from sourcerer.Http import HttpClient
from sourcerer.Schema import PSEUDO_VALUES, loadSchema, saveSchema
from sourcerer.Sources import ALIASES, REGISTRY, canonicalName, getSource
from sourcerer.Sources.Oas import isNull
from sourcerer.Version import __date__, __version__

log = logging.getLogger('sourcerer')

#: Output formats the download and convert subcommands can produce.
FORMATS = ('raw', 'airr', 'fasta')

#: Pseudo-collection for germline sources: fetch every species sourcerer supports
#: into one reference_base. Not offered for OAS (its collections are paired and
#: unpaired, not species) nor for search (two species would merge two hit lists).
ALL_SPECIES = 'all'

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
        # A presence-only flag takes a fixed token set, so the usage line spells
        # it out the way argparse does for choices; the others take a value
        # from a vocabulary too long for the usage line, listed in the help.
        metavar = 'VALUE'
        if item.pseudo_values:
            summary = 'filter on whether %s is recorded' % item.name
            metavar = '{%s,%s}' % (item.wildcard, ','.join(sorted(PSEUDO_VALUES)))
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
                            metavar=metavar, default=None, help=summary)


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
    build.add_argument('--map', dest='map_file', type=Path, default=None,
                       metavar='MANIFEST', help='a manifest declaring the species and chain of files whose names do not say, one per line: <file> <species> <CHAIN> [aa]. It overrides the naming rule, so it can also correct a file the rule misreads')

    diff = actions.add_parser(
        'diff', help='compare two germline reference folders allele by allele',
        description='Compare two reference folders allele by allele and report '
                    'what is identical, added, removed or changed. Files are '
                    'matched by name in any layout, and sequences are compared '
                    'without gaps, so a gapped and an ungapped copy of the same '
                    'allele are not a false difference. Exits non-zero when the '
                    'two references differ.',
        formatter_class=CommonHelpFormatter)
    diff.add_argument('reference_a', type=Path,
                      help='the baseline reference folder')
    diff.add_argument('reference_b', type=Path,
                      help='the reference folder to compare against it')
    diff.add_argument('--species', nargs='+', choices=list(Reference.SPECIES),
                      default=None,
                      help='limit to these species; default is every species '
                           'found in either folder')
    diff.add_argument('--map', dest='map_file', type=Path, default=None,
                      metavar='MANIFEST', help='a manifest declaring the species and chain of files whose names do not say, one per line: <file> <species> <CHAIN> [aa]. It overrides the naming rule, so it can also correct a file the rule misreads')

    show = actions.add_parser(
        'show', help='report what a reference folder is and where it came from',
        description='Read a reference folder\'s provenance sidecars -- '
                    'IMGT.yaml, AIRRC.yaml and sourcerer_build.yaml -- and '
                    'report the release and sets it was built from, what was '
                    'built, and what the folder holds. Accepts a reference_base, '
                    'a directory containing one, or an igblast_base, which '
                    'carries copies of the same sidecars.',
        formatter_class=CommonHelpFormatter)
    show.add_argument('folder', type=Path,
                      help='the reference folder to describe')
    show.add_argument('--map', dest='map_file', type=Path, default=None,
                      metavar='MANIFEST', help='a manifest declaring the species and chain of files whose names do not say, one per line: <file> <species> <CHAIN> [aa]. It overrides the naming rule, so it can also correct a file the rule misreads')


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

    check = actions.add_parser(
        'check', help='classify drift between two snapshots',
        description='Compare the packaged snapshot against a stored one and '
                    'classify every difference by severity: additive (new '
                    'values or units), anomaly (internal inconsistency), '
                    'removed (something users may pin disappeared) or '
                    'structural (the shape sourcerer parses changed). The '
                    'exit status reflects the overall level when it reaches '
                    'the --fail-on threshold.',
        formatter_class=CommonHelpFormatter)
    check.add_argument('--source', required=True,
                       choices=sorted(REGISTRY) + sorted(ALIASES),
                       help='which source to check')
    check.add_argument('--against', default='git:HEAD',
                       help='what to compare the packaged snapshot to: '
                            'git:REV reads the snapshot committed at that '
                            'revision, anything else is a snapshot directory')
    check.add_argument('--report', type=Path, default=None,
                       help='write the findings as JSON to this file')
    check.add_argument('--markdown', type=Path, default=None,
                       help='write the findings as markdown to this file')
    check.add_argument('--fail-on', default='structural',
                       choices=['never', 'additive', 'anomaly', 'removed',
                                'structural'],
                       help='exit non-zero when the overall level is at or '
                            'above this; never always exits 0')
    check.add_argument('--no-probe', action='store_true',
                       help='skip the live URL rule probe, for offline use')


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
        names = list(source.collections)
        if action == 'download' and source.output == 'reference' and len(names) > 1:
            names.append(ALL_SPECIES)

        for collection in names:
            # Passing help is what makes argparse list the collection at all.
            collection_help = source.collection_help.get(
                collection,
                'every species sourcerer supports for this source (%s), into '
                'one reference_base. Not every species the source publishes'
                % ', '.join(source.collections))
            leaf = collections.add_parser(
                collection, help=collection_help,
                description='%s the %s collection (%s), optionally narrowed '
                            'down with the filter flags below.'
                            % (action.capitalize(), collection, collection_help),
                formatter_class=CommonHelpFormatter)
            addFilterArgs(leaf, schema, name,
                           source.collections[0] if collection == ALL_SPECIES
                           else collection)
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
                    leaf.add_argument('--from', dest='from_ref', type=Path,
                                      default=None, metavar='REFERENCE',
                                      help='re-download the versions pinned in a '
                                           'reference_base (its IMGT.yaml / '
                                           'AIRRC.yaml), or one of those files, '
                                           'instead of the latest')
                    leaf.add_argument('--compare', type=Path, default=None,
                                      metavar='REFERENCE',
                                      help='after building, compare the result '
                                           'allele by allele against this '
                                           'reference folder and report; a '
                                           'difference is a non-zero exit')
                    leaf.add_argument('--resolve-doi', action='store_true',
                                      help='resolve each OGRDB set\'s Zenodo DOI '
                                           'into AIRRC.yaml (scrapes the OGRDB '
                                           'web UI; ignored by imgt)')
                else:
                    leaf.add_argument('--format', action='append', dest='formats',
                                      choices=FORMATS,
                                      help='what to write, repeatable to write '
                                           'several; raw mirrors the source files '
                                           'untouched and is always written whether '
                                           'or not it is requested, since airr and '
                                           'fasta are converted from it, so omitting '
                                           'this writes airr alone')
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
    report is the one thing this command produces, always with both NCBI's raw
    sample name and a suggested subject, so there was nothing left for a flag
    to switch between. It carries every column the input samplesheet had, evidence
    columns appended, so it is a drop-in airrflow input rather than a
    side file -- see buildEvidenceRow and NCBI_EVIDENCE_COLUMNS.
    """
    verify = actions.add_parser(
        'verify', help='cross-reference unresolved subjects against NCBI',
        description='Write an evidence report with one row per samplesheet '
                    'row, every input column carried through unchanged. '
                    'Every row whose sample_name yields a run/sample '
                    'accession (SRR/ERR/DRR/GSM) is looked up against NCBI, '
                    'whether or not OAS itself recorded a subject_id -- an '
                    'OAS-recorded subject can still be wrong (a typo, a '
                    'short code reused across studies, or a pooled/hashed '
                    'run naming several donors under one value), and '
                    'NCBI\'s own record is independent evidence either way. '
                    'ncbi_sample_name carries NCBI\'s raw sample name -- never '
                    'guessed at further, since some studies\' names need '
                    'study-specific reading to turn into a subject id (see '
                    'the module docstring in Ncbi.py); ncbi_subject_suggested '
                    'carries the same value with the handful of generic '
                    'patterns (a trailing locus or visit suffix) stripped, '
                    'the ones safe to normalize regardless of study. '
                    'subject_check reports how that compares to the '
                    'samplesheet\'s own subject_id: \'unresolved\' when OAS '
                    'recorded none, \'pooled\' when either side names more '
                    'than one donor, \'agrees\' or \'differs\' otherwise, or '
                    '\'unverified\' when NCBI itself could not resolve the '
                    'accession. A pooled/multi-donor run gets an '
                    'AMBIGUOUS_POOLED marker naming the donor codes in both '
                    'ncbi_ columns instead of a guessed single '
                    'subject. Because every input column survives, the '
                    'report can be pointed at directly as airrflow --input '
                    'once subject_id is filled in or corrected for any row '
                    'that needs it.',
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
                             'rate from 3/s to 10/s; prefer setting the '
                             'NCBI_API_KEY environment variable instead, '
                             'since a key given here ends up in shell history')


#: Request delay for the IgBLAST support mirror. The default is cautious because
#: IMGT and OGRDB are small academic servers; NCBI's file server and GitHub's raw
#: host are bulk services, and the mirror alone makes 100+ requests, so the
#: default would spend most of a build asleep rather than transferring.
MIRROR_DELAY = 0.05


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
    changed_any = changed

    wanted = args.collection or list(source.collections)
    catalogs = {}
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

        path, changed = Catalog.saveCatalog(rows, path)
        log.info('%s %s (%d units)',
                 'wrote' if changed else 'unchanged, left alone:', path, len(rows))
        changed_any = changed_any or changed
        catalogs[collection] = rows

    for name, (path, changed) in sorted(
            source.harvestArtifacts(out, schema, catalogs).items()):
        log.info('%s %s', 'wrote' if changed else 'unchanged, left alone:', path)
        changed_any = changed_any or changed

    # The provenance record moves only when the snapshot did: a quiet refresh
    # leaves every tracked file alone, which is what keeps the scheduled
    # workflow from opening a pull request on a quiet month.
    if changed_any:
        stamp = datetime.now(UTC).strftime('%Y-%m-%dT%H:%M:%SZ')
        record = Contracts.saveProvenance(out, stamp, __version__)
        log.info('wrote %s', record)

    return 0


def handleSchemaCheck(args):
    """Classify drift between the packaged snapshot and a stored one."""
    from sourcerer import Drift

    args.source = canonicalName(args.source)

    old = Drift.loadSnapshot(args.source, args.against)
    new = Drift.loadSnapshotDir(args.source)
    client = None if args.no_probe else makeClient(args)

    findings = Drift.checkDrift(old, new, client=client)
    report = Drift.buildReport(args.source, args.against, findings)

    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(Contracts.serializeJson(report))
        log.info('wrote %s', args.report)
    if args.markdown is not None:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(Drift.renderMarkdown(report))
        log.info('wrote %s', args.markdown)

    level = report['overall_level']
    if not findings:
        log.info('no drift against %s', args.against)
    else:
        for finding in findings:
            where = (' [%s]' % finding.collection) if finding.collection else ''
            log.info('%-10s %s%s: %s', finding.level, finding.category, where,
                     finding.message)
        log.info('overall level: %s (%d finding(s))', level, len(findings))

    return Drift.exitCode(findings, args.fail_on)


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
    elif units:
        print(formatUnitTable(units))

    return 0


#: Metadata shown beside each hit on stdout: enough to tell units apart and
#: pick the ones to download without first saving a catalog with --out.
SEARCH_COLUMNS = ('Species', 'Disease', 'Subject', 'BSource')


def formatUnitTable(units, columns=SEARCH_COLUMNS):
    """
    Render data units as an aligned text table for stdout.

    Arguments:
      units (list): DataUnit objects.
      columns (tuple): metadata keys to show after the identifier and count.

    Returns:
      str: the table, header first, without a trailing newline.
    """
    header = ['unit_id', 'n_unique_sequences'] + list(columns)
    rows = [[x.unit_id, str(x.n_sequences or '')]
            + [str(x.metadata.get(c, '') or '') for c in columns]
            for x in units]
    widths = [max(len(row[i]) for row in [header] + rows)
              for i in range(len(header))]

    lines = []
    for row in [header] + rows:
        cells = [row[0].ljust(widths[0]), row[1].rjust(widths[1])]
        cells += [cell.ljust(width) for cell, width in zip(row[2:], widths[2:])]
        lines.append('  '.join(cells).rstrip())

    return '\n'.join(lines)

def loadMap(args):
    """
    Read the --map manifest, if one was given.

    Arguments:
      args (Namespace): parsed arguments.

    Returns:
      dict: the manifest, or None.
    """
    if getattr(args, 'map_file', None) is None:
        return None

    return Reference.loadReferenceMap(args.map_file)


def handleReferenceDiff(args):
    """Compare two reference folders and report; non-zero exit if they differ."""
    for folder in (args.reference_a, args.reference_b):
        if not folder.is_dir():
            raise SourcererError('no such reference folder: %s' % folder)

    diff = Reference.diffReference(args.reference_a, args.reference_b,
                                   species=args.species, mapping=loadMap(args))
    print(diff.summary())

    return 0 if diff.same else 1


def handleReferenceShow(args):
    """Report what a reference folder is and where it came from."""
    print(Reference.describeReference(args.folder, mapping=loadMap(args)))

    return 0


def handleReference(args):
    """Validate a reference folder and, unless --check, build its IgBLAST base."""
    if args.action == 'diff':
        return handleReferenceDiff(args)
    if args.action == 'show':
        return handleReferenceShow(args)

    if not args.folder.is_dir():
        raise SourcererError('no such reference folder: %s' % args.folder)

    plan = Reference.planReference(args.folder, species=args.species,
                                   mapping=loadMap(args))
    print(plan.summary())

    if not plan.ok:
        raise SourcererError('no databases can be built from %s; check the file '
                             'names against <species>_<CHAIN>.fasta' % args.folder)

    if args.check:
        return 0

    if args.out is None:
        raise SourcererError('--out is required to build; pass --check to only '
                             'validate the folder')

    report = Reference.buildFromPlan(plan, args.out, HttpClient(delay=MIRROR_DELAY))
    for path in Reference.writeBuildMetadata(
            args.out, report, args.folder, Provenance.timestamp()[:10],
            'sourcerer %s' % __version__):
        log.info('wrote %s', path)
    log.info('wrote %s', args.out)

    return 0


def applyPins(source, from_ref, species):
    """
    Pin a reference source to the versions recorded in a reference_base.

    Reads the IMGT.yaml / AIRRC.yaml at --from and applies whichever pins the
    source can act on, so an imgt source takes the IMGT release, an ogrdb source
    the set versions, and the blend both. A --from whose pins none of the source
    can use is a mistake worth stopping on rather than silently fetching latest.

    Only the species being downloaded is read. A reference_base can hold several,
    downloaded weeks apart from different releases, and pinning a mouse download
    to human's release would quietly build something other than what was asked
    for.

    Arguments:
      source (ReferenceSource): the source about to download.
      from_ref (Path): a reference_base or an IMGT.yaml/AIRRC.yaml file.
      species (str): the species being downloaded.
    """
    from sourcerer.Sources.Germline import ReferenceSource

    pins = Reference.loadReferencePins(from_ref)
    imgt, airrc = pins.get('imgt') or {}, pins.get('airrc') or {}
    can_release = type(source).pinRelease is not ReferenceSource.pinRelease

    applied = []
    entry = (imgt.get('species') or {}).get(species) or {}
    if entry.get('release') and can_release:
        source.pinRelease(entry['release'])
        applied.append('IMGT release %s' % entry['release'])

    sets = [item for item in airrc.get('sets') or []
            if item.get('species') == species]
    if sets and hasattr(source, 'pinSets'):
        source.pinSets(sets)
        applied.append('%d OGRDB set version(s)' % len(sets))

    if not applied:
        raise SourcererError('the reference at %s records no %s versions that %s '
                             'can re-download' % (from_ref, species, source.name))
    log.info('re-downloading pinned: %s', '; '.join(applied))


def handleReferenceDownload(args, source):
    """Download one species, or all, and build an airrflow reference_base."""
    species = (list(source.collections) if args.collection == ALL_SPECIES
               else [args.collection])
    if args.resolve_doi and hasattr(source, 'enableDoi'):
        source.enableDoi()

    outdir = Path(args.outdir)
    reference_dir = outdir / 'reference_base'
    provenance = []
    # Each species is pinned, fetched and built on its own; the provenance
    # sidecars merge across them, so several land in one reference_base.
    for name in species:
        if args.from_ref is not None:
            applyPins(source, args.from_ref, name)
        query = source.validateQuery(name, collectFilters(args))
        if args.limit is not None:
            query = type(query)(collection=query.collection,
                                filters=query.filters, limit=args.limit)
        units = source.searchUnits(query)
        log.info('%d germline files for %s %s', len(units), args.source, name)
        if args.dry_run:
            for unit in units:
                print('%-48s %s' % (unit.unit_id, unit.url))
            continue
        entries = []
        for unit in units:
            result = source.fetchUnit(unit, outdir / 'raw',
                                      resume=not args.no_resume)
            entries.append((unit, result.path))
            provenance.append(Provenance.buildUnitRecord(unit, result, outdir, {}))
        source.buildReference(entries, reference_dir).logSummary()
        for path in source.writeReferenceMetadata(reference_dir, units):
            log.info('wrote %s', path)

    if args.dry_run:
        log.info('dry run: nothing downloaded')
        return 0

    log.info('wrote %s', reference_dir)

    formats = ['reference']
    if args.igblast:
        igblast_out = args.igblast_out or (outdir / 'igblast_base')
        report = Reference.buildIgblastBase(reference_dir, igblast_out,
                                            HttpClient(delay=MIRROR_DELAY),
                                            species=species)
        Reference.writeBuildMetadata(igblast_out, report, reference_dir,
                                     Provenance.timestamp()[:10],
                                     'sourcerer %s' % __version__)
        log.info('wrote %s', igblast_out)
        formats.append('igblast')

    record = Provenance.writeDownloadMetadata(
        outdir, args.source, args.collection, collectFilters(args), args.limit,
        formats, provenance, schema=source.schema, license=source.license,
        citation=source.citation)
    log.info('wrote %s', record)

    if args.compare is not None:
        if not args.compare.is_dir():
            raise SourcererError('no such reference folder: %s' % args.compare)
        diff = Reference.diffReference(args.compare, reference_dir,
                                       species=species)
        print(diff.summary())
        if not diff.same:
            log.warning('the built reference differs from %s', args.compare)
            return 1

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
    # Defaults to airr, not raw: sourcerer exists to hand Immcantation
    # something it can use directly, and the raw mirror alone (the previous
    # default) needs a second `download` run before airrflow can read
    # anything. Nothing is lost by defaulting this way -- raw is always
    # written regardless, see the bucket comment below.
    formats = args.formats or ['airr']
    total = sum(x.n_sequences or 0 for x in units)

    log.info('%d data units, %s sequences, formats: %s',
             len(units), format(total, ','), ', '.join(formats))

    if args.dry_run:
        if units:
            print(formatUnitTable(units))
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
        # One pass over the unit feeds every requested writer. Converting once
        # per format read and normalized the whole (multi-GB) file once per
        # format, so airr+fasta cost twice what airr alone did.
        writers = {}
        if 'airr' in formats:
            dest = outdir / 'airr' / ('%s.tsv' % stem)
            writers['airr'] = Convert.AirrWriter(dest, strict=args.strict_airr)
        if 'fasta' in formats:
            dest = outdir / 'fasta' / ('%s.fasta' % stem)
            writers['fasta'] = Convert.FastaWriter(dest)

        if writers:
            _, chunks, report = source.convertUnit(result.path, unit)
            try:
                for frame in chunks:
                    for writer in writers.values():
                        writer.write(frame)
            finally:
                for writer in writers.values():
                    writer.close()

            loci[unit.unit_id] = report['loci']
            for fmt, writer in writers.items():
                written[fmt].append((unit, writer.out))
                outputs[fmt] = writer.out

            if 'airr' in writers:
                validation = writers['airr'].validation
                Convert.writeValidationReport(validation, writers['airr'].out)
                log.info('%s: %d rows, %d invalid, %d rows in',
                         writers['airr'].out.name, validation['rows_checked'],
                         validation['rows_invalid'], report['rows_in'])

        provenance.append(
            Provenance.buildUnitRecord(unit, result, outdir, outputs))

    # A samplesheet is a derived artifact of a data format, so one is written per
    # converted format rather than one ambiguous sheet naming a single file.
    #
    # Computed once, from the raw bucket rather than per format: every unit is
    # in it regardless of which formats were requested, and a unit's Subject
    # metadata does not depend on which format its samplesheet row ends up in.
    unresolved = countUnresolvedSubjects(written['raw'])
    for fmt in ('airr', 'fasta'):
        if written.get(fmt):
            sheet = outdir / ('samplesheet_airrflow_%s.tsv' % fmt)
            buildSamplesheet(written[fmt], sheet, args.collection, outdir,
                             loci=loci)
            log.info('wrote %s', sheet)
            if unresolved:
                # The likeliest silent-wrong-analysis outcome this command can
                # produce: every one of these rows gets the same OAS sentinel
                # ('no'/'None') as subject_id, so airrflow would treat them as
                # one subject unless the user runs verify first.
                log.warning(
                    "%d of %d units have no subject recorded in OAS; run "
                    "'sourcerer oas verify %s' to cross-reference them "
                    'against NCBI before using this samplesheet',
                    unresolved, len(units), sheet)

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
NCBI_EVIDENCE_COLUMNS = ('ncbi_sample_name', 'ncbi_subject_suggested', 'status',
                         'subject_check', 'accession', 'biosample_accession',
                         'biosample_url', 'pooled_codes')


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


def normalizeForCompare(text):
    """
    Reduce text to bare alphanumerics for a formatting-insensitive comparison.

    Arguments:
      text (str): the text to normalize.

    Returns:
      str: lowercased, with everything but letters and digits stripped.
    """
    return ''.join(c for c in str(text or '').lower() if c.isalnum())


def subjectCheck(subject_id, found):
    """
    Compare a samplesheet's own subject_id against NCBI's evidence for it.

    Run unconditionally, even when subject_id already looks real: OAS's own
    value can still be wrong -- a typo, a short code reused across studies
    (see buildEvidenceRow's caller), or a pooled/hashed run that names
    several donors under one value (Ncbi.poolCodes is applied to subject_id
    itself here, not only to NCBI's text, because OAS's own Subject field
    uses the same semicolon-separated donor lists).

    Arguments:
      subject_id (str): the samplesheet row's own subject_id.
      found (Ncbi.Evidence): the NCBI lookup for this row's accession, or
        None if no accession could be read from sample_name.

    Returns:
      str: 'unresolved' if OAS recorded no subject at all, 'pooled' if
        either side names more than one donor, 'unverified' if NCBI could
        not resolve the accession (so there is nothing to compare against),
        otherwise 'agrees' or 'differs'.
    """
    if isNull(subject_id):
        return 'unresolved'

    from sourcerer.Ncbi import poolCodes
    if poolCodes(subject_id) or (found is not None and found.status == 'pooled'):
        return 'pooled'

    if found is None or found.status != 'ok':
        return 'unverified'

    ncbi_text = normalizeForCompare(found.sample_name)
    oas_text = normalizeForCompare(subject_id)
    if oas_text and ncbi_text and (oas_text in ncbi_text or ncbi_text in oas_text):
        return 'agrees'

    return 'differs'


def buildEvidenceRow(row, found):
    """
    Build one row of the verify evidence report.

    Every column already on `row` passes through verbatim -- this only adds
    NCBI_EVIDENCE_COLUMNS, it never edits or drops what was already on the
    samplesheet. ncbi_sample_name and ncbi_subject_suggested always come from
    NCBI, never from the samplesheet's own subject_id: a real-looking
    subject_id is not proof it is correct, which is exactly what
    subject_check is for. A pooled/multi-donor run gets an AMBIGUOUS_POOLED
    marker in both columns, since there is no single subject to suggest.
    Anything else takes ncbi_sample_name from NCBI's raw sample name and
    ncbi_subject_suggested from the generic-heuristic strip of it (see
    Ncbi.suggestSubject) -- both always written, so joining this report back
    onto a samplesheet never needs a flag to decide which one it gets.

    Arguments:
      row (dict): the samplesheet row.
      found (Ncbi.Evidence): the NCBI lookup for this row's accession, or
        None if no accession could be read from sample_name.

    Returns:
      dict: row, plus NCBI_EVIDENCE_COLUMNS.
    """
    subject_id = row.get('subject_id', '')

    if found is None:
        evidence = {'ncbi_sample_name': '', 'ncbi_subject_suggested': '',
                   'status': 'no_accession', 'accession': '', 'biosample_accession': '',
                   'biosample_url': '', 'pooled_codes': ''}
    else:
        pooled_codes = ';'.join(found.pooled_codes)
        if found.status == 'pooled':
            ncbi_sample_name = ncbi_subject_suggested = 'AMBIGUOUS_POOLED:%s' % pooled_codes
        elif found.status == 'ok':
            ncbi_sample_name = found.sample_name
            ncbi_subject_suggested = found.suggested_subject
        else:
            ncbi_sample_name = ncbi_subject_suggested = ''
        evidence = {'ncbi_sample_name': ncbi_sample_name,
                   'ncbi_subject_suggested': ncbi_subject_suggested,
                   'status': found.status, 'accession': found.accession,
                   'biosample_accession': found.biosample_accession,
                   'biosample_url': found.url, 'pooled_codes': pooled_codes}

    evidence['subject_check'] = subjectCheck(subject_id, found)

    return {**row, **evidence}


def warnReusedSubjects(rows):
    """
    Log a warning when the same subject_id is used by more than one study.

    airrflow keys a subject on subject_id alone, so two different studies
    reusing the same short code (e.g. 'Donor-2') silently merges two
    unrelated people into one subject downstream. This is exactly the kind
    of collision a row-by-row read of subject_check is unlikely to catch,
    since each individual row looks unremarkable on its own.

    Arguments:
      rows (list): samplesheet rows; only those carrying both subject_id and
        study are considered.
    """
    studies_by_subject = {}
    for row in rows:
        subject = row.get('subject_id', '')
        study = row.get('study', '')
        if isNull(subject) or not study:
            continue
        studies_by_subject.setdefault(subject, set()).add(study)

    for subject, studies in sorted(studies_by_subject.items()):
        if len(studies) > 1:
            log.warning("subject_id '%s' is used by %d different studies: %s",
                       subject, len(studies), ', '.join(sorted(studies)))


def handleOasVerify(args):
    """Cross-reference every samplesheet row's subject against NCBI."""
    fields, rows = readSamplesheetRows(args.samplesheet)
    warnReusedSubjects(rows)

    # Every row with a readable accession is looked up, not only rows OAS
    # left unresolved: a subject_id OAS did record is still worth checking
    # against NCBI's own record (see subjectCheck), so there is no 'pending'
    # subset here to restrict the lookup to.
    accession_by_sample = {}
    for row in rows:
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
            if args.action == 'check':
                return handleSchemaCheck(args)
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
