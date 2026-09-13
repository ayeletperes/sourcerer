"""
`sourcerer annotate` -- the commandline front end for Annotate.py

Kept apart from sourcerer.Cli (which only imports addAnnotateParser/
handleAnnotate and wires them into the top level parser and dispatcher) so
that the annotation step -- argument definitions, help text, and the handler
itself -- lives alongside the rest of the specificity annotation code it
operates on.

Reads two AIRR rearrangement TSVs already on disk -- typically a user's own
receptor table as `--query` and a specificity database's downloaded table
(e.g. `sourcerer specificity iedb download bcr --format airr`, or its raw
`sourcerer specificity iedb download bcr` output) as `--reference` -- and
writes one antigen specificity annotation per query row. It does not itself
download or convert anything; that is steps one and two of the pipeline,
already handled by `sourcerer specificity <db> download`.
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import logging
from pathlib import Path

import pandas

# Sourcerer imports
from sourcerer.Commandline import CommonHelpFormatter
from sourcerer.Exceptions import AnnotationError
from sourcerer.Sources.Specificity.Annotate import (DEFAULT_ID_COL,
                                                    runExactAnnotation,
                                                    runFuzzyAnnotation)

log = logging.getLogger(__name__)


def addAnnotateParser(commands):
    """
    Add the `annotate` subcommand to the top level parser.

    Arguments:
      commands (argparse._SubParsersAction): the top level subparsers
        object, as built by `sourcerer.Cli.getArgParser`.
    """
    parser = commands.add_parser(
        'annotate', help='annotate a receptor table against a reference table',
        description='Annotate each row of a query AIRR rearrangement table '
                    'with antigen specificity information from a fixed '
                    'reference AIRR rearrangement table -- the third step '
                    'after downloading and converting a specificity '
                    'database, e.g. `sourcerer specificity iedb download '
                    'bcr --format airr`.',
        formatter_class=CommonHelpFormatter)

    parser.add_argument('--query', type=Path, required=True,
                        help='AIRR rearrangement TSV to annotate, one row '
                             'per chain/sequence')
    parser.add_argument('--reference', type=Path, required=True,
                        help='AIRR rearrangement TSV to annotate against, '
                             'e.g. a specificity database\'s downloaded '
                             'table')
    parser.add_argument('--compare-cols', action='append', required=True,
                        dest='compare_cols', metavar='COLUMN',
                        help='column to compare. At --max-mismatches 0 '
                             '(exact, CompAIRR-backed) exactly one column is '
                             'allowed; above 0 (fuzzy) this is repeatable, '
                             'and mismatches are pooled across every column '
                             'given, e.g. --compare-cols cdr3_heavy '
                             '--compare-cols cdr3_light')
    parser.add_argument('--max-mismatches', type=float, default=0,
                        help='0 (default) requires exact agreement, '
                             'CompAIRR-backed; above 0 allows that many '
                             'substitutions, pooled across --compare-cols, '
                             'pure Python')
    parser.add_argument('--fraction', action='store_true',
                        help='treat --max-mismatches as a fraction of the '
                             'query\'s length instead of a fixed count '
                             '(e.g. 0.15 allows up to 15%% mismatches); no '
                             'effect at --max-mismatches 0')
    parser.add_argument('--exact-cols', action='append', default=None,
                        dest='exact_cols', metavar='COLUMN',
                        help='fuzzy annotation only (--max-mismatches above '
                             '0): column that must agree exactly before '
                             '--compare-cols are even compared, e.g. '
                             'v_call; repeatable. For exact annotation, use '
                             '--vgene-col/--jgene-col instead')
    parser.add_argument('--vgene-col', default=None, dest='vgene_col',
                        metavar='COLUMN',
                        help='exact annotation only (--max-mismatches 0): '
                             'column to require exact V-gene agreement on')
    parser.add_argument('--jgene-col', default=None, dest='jgene_col',
                        metavar='COLUMN',
                        help='exact annotation only (--max-mismatches 0): '
                             'column to require exact J-gene agreement on')
    parser.add_argument('--compairr-bin', default=None, dest='compairr_bin',
                        metavar='PATH',
                        help='exact annotation only (--max-mismatches 0): '
                             'path to the compairr binary; without it, PATH '
                             'is searched')
    parser.add_argument('--threads', type=int, default=1,
                        help='exact annotation only (--max-mismatches 0): '
                             'compairr threads to use')
    parser.add_argument('--id-col', default=DEFAULT_ID_COL, dest='id_col',
                        metavar='COLUMN',
                        help='column in --reference to report as each '
                             'hit\'s identifier')
    parser.add_argument('--locus', default=None,
                        help='restrict both --query and --reference to '
                             'this locus (e.g. IGH) before annotating, for '
                             'AIRR tables holding more than one chain per '
                             'row; without it every row is annotated as is')
    parser.add_argument('-o', '--out', type=Path, required=True,
                        help='write the annotation table here')


def _readAirrTsv(path):
    """
    Read an AIRR rearrangement TSV as plain text columns.

    Arguments:
      path (Path): the file to read.

    Returns:
      pandas.DataFrame: every column as a string, blank cells as '' rather
      than NaN -- matching how Convert.writeAirr wrote them, and how
      Annotate.py's row comparisons expect a value to always be a string.
    """
    return pandas.read_csv(path, sep='\t', dtype=str, na_filter=False)


def _filterByLocus(frame, locus, label):
    """
    Restrict a table to one locus, if a locus filter was requested.

    Arguments:
      frame (pandas.DataFrame): the table to filter.
      locus (str): the locus to keep, or None to leave `frame` unchanged.
      label (str): 'query' or 'reference', for the log line only.

    Returns:
      pandas.DataFrame: the filtered table.
    """
    if locus is None:
        return frame

    filtered = frame[frame['locus'] == locus]
    log.info('%s: %d/%d rows are locus %s', label, len(filtered), len(frame), locus)

    return filtered


def handleAnnotate(args):
    """
    Read --query and --reference, run the requested annotation, and write
    the result.

    Arguments:
      args (argparse.Namespace): parsed arguments, as built by
        addAnnotateParser.

    Returns:
      int: process exit status.
    """
    queries = _filterByLocus(_readAirrTsv(args.query), args.locus, 'query')
    reference = _filterByLocus(_readAirrTsv(args.reference), args.locus, 'reference')
    queries = queries.reset_index(drop=True)
    reference = reference.reset_index(drop=True)

    if args.max_mismatches == 0:
        if len(args.compare_cols) != 1:
            raise AnnotationError(
                'exact annotation (CompAIRR-backed) compares exactly one '
                'column; got %d via --compare-cols (%s). Use '
                '--max-mismatches above 0 for multi-column pooled '
                'annotation instead.'
                % (len(args.compare_cols), ', '.join(args.compare_cols)))
        result = runExactAnnotation(queries, reference, args.compare_cols[0],
                                    id_col=args.id_col,
                                    vgene_col=args.vgene_col,
                                    jgene_col=args.jgene_col,
                                    compairr_bin=args.compairr_bin,
                                    threads=args.threads)
    elif args.fraction:
        result = runFuzzyAnnotation(queries, reference, args.compare_cols,
                                    args.max_mismatches, id_col=args.id_col,
                                    exact_cols=args.exact_cols,
                                    threshold_is_fraction=True)
    else:
        result = runFuzzyAnnotation(queries, reference, args.compare_cols,
                                    args.max_mismatches, id_col=args.id_col,
                                    exact_cols=args.exact_cols)

    hit = int(result['n_hits_total'].gt(0).sum())
    log.info('%d query rows, %d with at least one hit, %d unique (0 or 1 hit)',
             len(result), hit, int(result['is_unique_hit'].sum()))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.out, sep='\t', index=False)
    log.info('wrote %s', args.out)

    return 0
