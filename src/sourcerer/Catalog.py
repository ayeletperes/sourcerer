"""
Data unit catalogs

A catalog is the checked-in index of every data unit a collection contains. It
is written as a sorted TSV rather than kept in the source's native form so that
the monthly refresh produces a line wise diff a human can read: "four new units
in Smith_2026" rather than one changed line of a multi megabyte JSON document.

For OAS paired data the catalog is not merely a convenience. The source publishes
no machine readable index of it at all, so this file is the only way to search
paired data without the search form being up.
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import csv
import logging
from pathlib import Path

log = logging.getLogger(__name__)

#: Columns every catalog carries, in output order. Trailing columns may be empty:
#: run is absent from most paired filenames, and the enrichment columns are only
#: populated once a unit's detail page has been read.
CATALOG_COLUMNS = ('unit_id', 'collection', 'url', 'dir_segment', 'study', 'run',
                   'n_unique_sequences', 'Species', 'Isotype', 'Chain', 'Disease',
                   'Vaccine', 'Subject', 'Age', 'Longitudinal', 'BSource', 'BType',
                   'Author', 'detail_status', 'detail_attempted')

#: Value of detail_status meaning the unit's detail page has been read.
DETAIL_OK = 'ok'


def loadCatalog(path):
    """
    Read a catalog from disk.

    Arguments:
      path (Path): the TSV file.

    Returns:
      list: one dict per data unit.
    """
    path = Path(path)
    if not path.exists():
        return []

    with open(path, newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def saveCatalog(rows, path, columns=CATALOG_COLUMNS):
    """
    Write a catalog deterministically.

    Rows are sorted by identifier and columns are fixed, so re-harvesting an
    unchanged collection reproduces the file byte for byte. Without that, every
    scheduled refresh would look like a change and the review workflow would stop
    meaning anything.

    Arguments:
      rows (iterable): dicts of column to value.
      path (Path): where to write.
      columns (tuple): the column order.

    Returns:
      Path: the file written.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    ordered = sorted(rows, key=lambda x: x.get('unit_id', ''))
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter='\t',
                                extrasaction='ignore', lineterminator='\n')
        writer.writeheader()
        for row in ordered:
            writer.writerow({x: row.get(x, '') for x in columns})

    return path


def mergeEnrichment(existing, fresh):
    """
    Carry forward enrichment that a new harvest did not manage to fetch.

    A detail page that fails once must not blank values a previous run already
    obtained, and the unit must stay eligible for another attempt rather than
    being written off because it is no longer new.

    Arguments:
      existing (list): rows from the stored catalog.
      fresh (list): rows from the current harvest.

    Returns:
      list: fresh rows with previously known enrichment preserved.
    """
    previous = {x['unit_id']: x for x in existing if x.get('unit_id')}

    merged = []
    for row in fresh:
        old = previous.get(row.get('unit_id'))
        if old is not None:
            for column in ('BSource', 'BType', 'Author'):
                if not row.get(column) and old.get(column):
                    row[column] = old[column]
            if not row.get('detail_status') and old.get('detail_status'):
                row['detail_status'] = old['detail_status']
                row['detail_attempted'] = old.get('detail_attempted', '')
        merged.append(row)

    return merged


def needsDetail(row):
    """
    Test whether a unit still needs its detail page read.

    Selection is by recorded status rather than by novelty. Choosing only new
    units would strand any unit whose page failed once, because it will never be
    new again and so would keep its empty BSource and BType forever.

    Arguments:
      row (dict): a catalog row.

    Returns:
      bool: True if the detail page should be fetched.
    """
    return row.get('detail_status') != DETAIL_OK


def isDefined(value):
    """
    Test whether a catalog cell carries information.

    Arguments:
      value (str): the cell.

    Returns:
      bool: True if the value is neither empty nor a null sentinel.
    """
    return str(value or '').strip().lower() not in ('', 'no', 'none', 'na',
                                                    'unknown', 'undefined')


def filterCatalog(rows, filters, wildcard='*'):
    """
    Select catalog rows matching a set of field filters.

    Some fields filter on presence rather than on value: the source offers only
    'defined' and 'undefined' for them, so comparing literally would match
    nothing.

    Arguments:
      rows (iterable): catalog rows.
      filters (dict): field to value; the wildcard matches everything.
      wildcard (str): the value meaning "all".

    Returns:
      list: matching rows.
    """
    selected = []
    for row in rows:
        keep = True
        for name, value in filters.items():
            if value == wildcard:
                continue
            if name not in row:
                # A field the catalog does not carry cannot be filtered offline.
                keep = False
                break

            if value == 'defined':
                keep = isDefined(row[name])
            elif value == 'undefined':
                keep = not isDefined(row[name])
            else:
                keep = row[name] == value

            if not keep:
                break

        if keep:
            selected.append(row)

    return selected
