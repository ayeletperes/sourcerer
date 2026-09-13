"""
`sourcerer specificity bcr` -- the no-verb shortcut for fetching IEDB BCR

Sugar for `sourcerer specificity iedb download bcr`: BCR is currently the
only table this project annotates against, so this drops both the source
name and the `download` verb -- naming the table is enough. Calling it
fetches IEDB's receptor bulk export and converts it to AIRR in the same
step; there is no separate command for "just download" or "just convert."
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
from argparse import Namespace
from pathlib import Path

# Sourcerer imports
from sourcerer.Commandline import CommonHelpFormatter

#: Written into when no --outdir is given.
DEFAULT_OUTDIR = Path('tmp')


def addBcrParser(dbs):
    """
    Add the `bcr` shortcut to the `specificity` subcommand.

    Arguments:
      dbs (argparse._SubParsersAction): the `specificity` subparsers object
        (its `db` level), as built by
        `sourcerer.Cli._addSpecificityGroupParser`.
    """
    parser = dbs.add_parser(
        'bcr', help='fetch + convert IEDB BCR receptor sequences to AIRR',
        description='Download IEDB\'s BCR (antibody) receptor table and '
                    'convert it to AIRR in one step -- shorthand for '
                    '`sourcerer specificity iedb download bcr`.',
        formatter_class=CommonHelpFormatter)
    parser.add_argument('--outdir', type=Path, default=DEFAULT_OUTDIR,
                        help='directory to write into')
    parser.add_argument('--strict-airr', action='store_true',
                        help='drop columns the AIRR schema does not define')


def handleBcr(args):
    """
    Fetch and convert IEDB's `bcr` table.

    Delegates to `sourcerer.Cli.handleSpecificityDownload` rather than
    duplicating it -- imported locally to avoid a circular import, since
    `sourcerer.Cli` is what wires this module in.

    Arguments:
      args (argparse.Namespace): parsed arguments, as built by addBcrParser.

    Returns:
      int: process exit status.
    """
    from sourcerer.Cli import handleSpecificityDownload

    return handleSpecificityDownload(Namespace(
        db='iedb', table='bcr', outdir=args.outdir, limit=None,
        dry_run=False, no_resume=False, strict_airr=args.strict_airr))
