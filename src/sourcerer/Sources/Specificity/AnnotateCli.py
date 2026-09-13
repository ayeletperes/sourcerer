"""
`sourcerer annotate` -- the commandline front end for Annotate.py

Kept apart from sourcerer.Cli (which only imports addAnnotateParser/
handleAnnotate and wires them into the top level parser and dispatcher) so
that the annotation step -- argument definitions, help text, and the handler
itself -- lives alongside the rest of the specificity annotation code it
operates on.

Reads two AIRR rearrangement TSVs already on disk -- typically a user's own
receptor table as `--query` and a specificity database's downloaded table
(e.g. `sourcerer specificity bcr`'s output) as `--reference` -- and writes
one antigen specificity annotation per query row. It does not itself
download or convert anything; that is handled by `sourcerer specificity
bcr`.
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
from sourcerer.Sources.Specificity.Annotate import (DEFAULT_ID_COL, HEAVY_LOCUS,
                                                    runExactAnnotation,
                                                    runFuzzyAnnotation,
                                                    runLevenshteinAnnotation,
                                                    runPairedExactAnnotation)

log = logging.getLogger(__name__)

#: --fuzzy's fixed gene-gating columns -- V-gene and J-gene must always agree
#: exactly; only --compare-cols (e.g. cdr3_aa) gets any substitution tolerance.
#: Matches runPairedExactAnnotation's own vgene_col/jgene_col defaults.
FUZZY_EXACT_COLS = ['v_call', 'j_call']


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
                    'reference AIRR rearrangement table, e.g. `sourcerer '
                    'specificity bcr`\'s output.',
        formatter_class=CommonHelpFormatter)

    parser.add_argument('--query', type=Path, required=True,
                        help='AIRR rearrangement TSV to annotate, one row '
                             'per chain/sequence')
    parser.add_argument('--reference', type=Path, required=True,
                        help='AIRR rearrangement TSV to annotate against, '
                             'e.g. a specificity database\'s downloaded '
                             'table')
    parser.add_argument('--compare-cols', action='append', default=None,
                        dest='compare_cols', metavar='COLUMN',
                        help='column to compare; required for --exact/'
                             '--fuzzy/--identity/--levenshtein, ignored '
                             'under --paired. Under --exact/--levenshtein '
                             'exactly one column is allowed; under --fuzzy/'
                             '--identity this is repeatable, and mismatches '
                             'are pooled across every column given (e.g. '
                             '--compare-cols cdr1_aa --compare-cols '
                             'cdr3_aa), on top of the fixed V/J-gene gating '
                             '-- both sides are still the heavy chain, none '
                             'of these pool across chains')

    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--exact', action='store_true',
                      help='require exact agreement on a single '
                           '--compare-cols column (the default): '
                           'CompAIRR-backed')
    mode.add_argument('--fuzzy', type=int, default=None, metavar='N',
                      help='exact V-gene and J-gene agreement (v_call, '
                           'j_call), plus up to N substitutions on '
                           '--compare-cols (pooled if more than one column); '
                           'pure Python, no CompAIRR. Restricted to the '
                           'heavy chain (locus IGH) to bound its cost, which '
                           'grows with candidate pool size unlike --exact\'s '
                           'hashed lookup; --locus must be IGH or left unset')
    mode.add_argument('--identity', type=float, default=None, metavar='PCT',
                      help='same as --fuzzy, but the substitution tolerance '
                           'is expressed as a percent identity over the '
                           'query\'s length instead of a fixed count, e.g. '
                           '85 allows up to 15%% of positions to differ -- '
                           'the right choice for --compare-cols columns '
                           'whose length varies row to row, where a fixed '
                           '--fuzzy N would be stricter on short sequences '
                           'and looser on long ones')
    mode.add_argument('--levenshtein', type=int, default=None, metavar='N',
                      help='exact V-gene and J-gene agreement (v_call, '
                           'j_call), plus up to N substitutions, insertions '
                           'and deletions combined on a single --compare-cols '
                           'column -- indel-tolerant, unlike --fuzzy/'
                           '--identity (substitutions only). Pure Python '
                           '(rapidfuzz), no CompAIRR; restricted to the '
                           'heavy chain same as --fuzzy; does not pool '
                           'across multiple --compare-cols and has no \'X\' '
                           'wildcard handling')
    mode.add_argument('--paired', action='store_true',
                      help='exact V-gene, exact J-gene and exact CDR3 '
                           'agreement, required on both the heavy and the '
                           'light chain of the same reference antibody -- '
                           'the strict definition of an exact match for a '
                           'paired receptor. Classifies each query cell_id '
                           'as both/heavy_only/light_only/conflicting/'
                           'neither rather than listing hits per chain row; '
                           '--compare-cols, --vgene-col, --jgene-col, '
                           '--id-col and --locus have no effect here, since '
                           'the columns and both loci (cdr3_aa, v_call, '
                           'j_call, cell_id, IGH/IGK/IGL) are fixed by that '
                           'definition')

    parser.add_argument('--vgene-col', default=None, dest='vgene_col',
                        metavar='COLUMN',
                        help='--exact only: column to require exact V-gene '
                             'agreement on')
    parser.add_argument('--jgene-col', default=None, dest='jgene_col',
                        metavar='COLUMN',
                        help='--exact only: column to require exact J-gene '
                             'agreement on')
    parser.add_argument('--compairr-bin', default=None, dest='compairr_bin',
                        metavar='PATH',
                        help='--exact/--paired only: path to the compairr '
                             'binary; without it, PATH is searched')
    parser.add_argument('--threads', type=int, default=1,
                        help='--exact/--paired only: compairr threads to use')
    parser.add_argument('--id-col', default=DEFAULT_ID_COL, dest='id_col',
                        metavar='COLUMN',
                        help='column in --reference to report as each '
                             'hit\'s identifier')
    parser.add_argument('--locus', default=None,
                        help='restrict both --query and --reference to '
                             'this locus (e.g. IGH) before annotating, for '
                             'AIRR tables holding more than one chain per '
                             'row; without it every row is annotated as is. '
                             '--fuzzy/--identity/--levenshtein always '
                             'restrict to IGH regardless')
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
    if args.paired:
        queries = _readAirrTsv(args.query)
        reference = _readAirrTsv(args.reference)
        result = runPairedExactAnnotation(queries, reference,
                                          compairr_bin=args.compairr_bin,
                                          threads=args.threads)
        log.info('%d query antibodies: %s', len(result),
                 ', '.join('%s=%d' % (label, count) for label, count
                          in result['classification'].value_counts().items()))
    else:
        fuzzy_family = (args.fuzzy is not None or args.identity is not None
                        or args.levenshtein is not None)
        if fuzzy_family:
            # Every pure-Python path (fuzzy, identity, levenshtein) has a cost
            # that grows with candidate pool size, unlike --exact's hashed
            # CompAIRR lookup -- restricting to the heavy chain keeps that
            # cost bounded rather than doubling it (or worse) by also
            # scanning the light chain.
            if args.locus not in (None, HEAVY_LOCUS):
                raise AnnotationError(
                    "--fuzzy/--identity/--levenshtein are restricted to the "
                    "heavy chain; --locus must be '%s' or left unset, got "
                    "'%s'" % (HEAVY_LOCUS, args.locus))
            locus = HEAVY_LOCUS
        else:
            locus = args.locus

        queries = _filterByLocus(_readAirrTsv(args.query), locus, 'query')
        reference = _filterByLocus(_readAirrTsv(args.reference), locus,
                                   'reference')
        queries = queries.reset_index(drop=True)
        reference = reference.reset_index(drop=True)

        if not args.compare_cols:
            raise AnnotationError('--compare-cols is required for --exact/'
                                  '--fuzzy/--identity/--levenshtein')

        if args.identity is not None:
            if not 0 < args.identity < 100:
                raise AnnotationError(
                    '--identity is a percent identity in (0, 100); got %s. '
                    'Use --exact for 100%% identity.' % args.identity)
            result = runFuzzyAnnotation(queries, reference, args.compare_cols,
                                        (100 - args.identity) / 100,
                                        id_col=args.id_col,
                                        exact_cols=FUZZY_EXACT_COLS,
                                        threshold_is_fraction=True)
        elif args.fuzzy is not None:
            result = runFuzzyAnnotation(queries, reference, args.compare_cols,
                                        args.fuzzy, id_col=args.id_col,
                                        exact_cols=FUZZY_EXACT_COLS)
        elif args.levenshtein is not None:
            if len(args.compare_cols) != 1:
                raise AnnotationError(
                    '--levenshtein compares exactly one column; got %d via '
                    '--compare-cols (%s).'
                    % (len(args.compare_cols), ', '.join(args.compare_cols)))
            result = runLevenshteinAnnotation(queries, reference,
                                              args.compare_cols[0],
                                              args.levenshtein,
                                              id_col=args.id_col,
                                              exact_cols=FUZZY_EXACT_COLS)
        else:
            if len(args.compare_cols) != 1:
                raise AnnotationError(
                    'exact annotation (CompAIRR-backed) compares exactly one '
                    'column; got %d via --compare-cols (%s). Use --fuzzy N '
                    'for multi-column pooled annotation instead.'
                    % (len(args.compare_cols), ', '.join(args.compare_cols)))
            result = runExactAnnotation(queries, reference, args.compare_cols[0],
                                        id_col=args.id_col,
                                        vgene_col=args.vgene_col,
                                        jgene_col=args.jgene_col,
                                        compairr_bin=args.compairr_bin,
                                        threads=args.threads)

        hit = int(result['n_hits_total'].gt(0).sum())
        log.info('%d query rows, %d with at least one hit, %d unique (0 or 1 hit)',
                 len(result), hit, int(result['is_unique_hit'].sum()))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.out, sep='\t', index=False)
    log.info('wrote %s', args.out)

    return 0
