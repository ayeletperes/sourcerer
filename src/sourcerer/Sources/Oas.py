"""
Observed Antibody Space (OAS)

OAS exposes no documented API. Discovery works differently for its two
collections, and the difference is not cosmetic:

- unpaired has a complete catalog as a single JSON document, so every data unit
  and its metadata can be listed without touching the search form;
- paired has no catalog at all. The only way to enumerate it is to submit the
  search form and read the download commands out of the JavaScript in the reply.

Paths are treated as opaque throughout. Paired data currently lives under two
different directory layouts and several filename patterns, and for most units the
run accession does not appear in the filename at all, so anything that rebuilt a
path from parsed components would mishandle the majority of the collection.
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import csv
import gzip
import hashlib
import json
import logging
import re
from datetime import UTC
from pathlib import Path
from urllib.parse import urlparse

import pandas
from bs4 import BeautifulSoup

# Sourcerer imports
from sourcerer.Catalog import DETAIL_OK, filterCatalog, loadCatalog, needsDetail
from sourcerer.Convert import coerceAirrTypes
from sourcerer.Exceptions import OasParseError
from sourcerer.Sources.Base import DataUnit, SourceBase

log = logging.getLogger(__name__)

#: Endpoints.
HOST = 'https://opig.stats.ox.ac.uk'
PAIRED_FORM_URL = HOST + '/webapps/oas/oas_paired/'
UNPAIRED_FORM_URL = HOST + '/webapps/oas/oas_unpaired/'
CATALOG_URL = HOST + '/webapps/ngsdb/oas_metadata_map.json'
DOWNLOAD_BASE = HOST + '/webapps/ngsdb/'
DETAIL_URL = HOST + '/webapps/oas/dataunit_%s'

#: Catalog keys are absolute server paths; this prefix is what maps them to URLs.
CATALOG_KEY_PREFIX = '/vols/naga-datasets/oas/'

#: The collections OAS offers, paired first.
COLLECTIONS = ('paired', 'unpaired')

#: What each collection contains, for `--help`.
COLLECTION_HELP = {
    'paired': 'single cell runs, heavy and light chain paired per cell',
    'unpaired': 'bulk runs, one chain per sequence and no pairing',
}

#: The number of results is reported in prose, not in a machine readable field.
COUNT_REGEX = (r'yielded\s*<b>([\d,]+)</b>\s*filtered sequences from\s*'
               r'<b>([\d,]+)</b>\s*studies')

#: The download commands are embedded in a JavaScript array.
CSV_ARRAY_MARKER = 'var CSV = ['
CSV_ARRAY_REGEX = r'var CSV\s*=\s*\[(.*?)\]\.join'
WGET_REGEX = r'"wget ([^"]+)"'

#: OAS spells the same concept differently in the search results, the search form
#: and the data unit metadata. Everything downstream sees the form spelling.
FIELD_ALIASES = {
    'Organism': 'Species',
    'Individual': 'Subject',
    'DS Name': 'Study',
    '#Unique Sequences': 'Unique sequences',
}

#: Values OAS uses to mean "not recorded".
NULL_TOKENS = frozenset(['', 'no', 'No', 'none', 'None', 'NA', 'n/a',
                         'unknown', 'undefined'])

#: Isotype values that are not real isotypes. Bulk appears in unpaired metadata
#: and All in paired metadata; writing either into c_call would be false data.
NON_ISOTYPES = frozenset(['Bulk', 'All'])


def isNull(value):
    """
    Test whether a metadata value is one of the source's null sentinels.

    Arguments:
      value: the value to test.

    Returns:
      bool: True if the value carries no information.
    """
    if value is None:
        return True

    return str(value).strip() in NULL_TOKENS


def unescapeOption(text):
    """
    Undo the escaping OAS applies to option labels.

    Some vocabulary values contain commas, which the page escapes as ``\\,``.

    Arguments:
      text (str): the raw option text.

    Returns:
      str: the value as the form will accept it.
    """
    return text.replace('\\,', ',').strip()


def parseFormSchema(html, collection):
    """
    Extract the searchable fields and their vocabularies from a search form.

    The options are rendered server side, so the full controlled vocabulary is
    present in the HTML and no JavaScript needs to be executed.

    Arguments:
      html (str): the search form page.
      collection (str): 'paired' or 'unpaired', used only in error messages.

    Returns:
      list: dicts with keys name, wildcard, values and pseudo_values, in the
      order the form presents them.

    Raises:
      OasParseError: if the page contains no form or no select elements.
    """
    soup = BeautifulSoup(html, 'html.parser')
    form = soup.find('form')
    if form is None:
        raise OasParseError('no <form> found on the OAS %s search page; the page '
                            'layout has changed' % collection)

    selects = form.find_all('select')
    if not selects:
        raise OasParseError('no <select> elements in the OAS %s search form; the '
                            'vocabulary can no longer be harvested' % collection)

    fields = []
    for select in selects:
        name = select.get('name')
        if not name:
            continue

        options = [unescapeOption(x.get_text()) for x in select.find_all('option')]
        options = [x for x in options if x]
        if not options:
            raise OasParseError("the %s field '%s' has no options" %
                                (collection, name))

        wildcard = options[0]
        values = options[1:]
        # Age, Subject and Longitudinal on the paired form offer only presence
        # flags. Recording them as a vocabulary would make the interactive builder
        # offer nonsense and would make a comparison against unpaired look like a
        # mass deletion.
        pseudo = set(values) == {'defined', 'undefined'}

        fields.append({'name': name, 'wildcard': wildcard,
                       'values': [] if pseudo else values,
                       'pseudo_values': pseudo})

    if not fields:
        raise OasParseError('no named <select> elements in the OAS %s search form'
                            % collection)

    return fields


def parseSearchTotals(html):
    """
    Read the reported sequence and study counts from a search reply.

    Arguments:
      html (str): the search results page.

    Returns:
      dict: sequences and studies as integers.

    Raises:
      OasParseError: if the count sentence is absent, which means either the
        search failed or the page wording changed.
    """
    match = re.search(COUNT_REGEX, html)
    if match is None:
        raise OasParseError(
            'could not find the result count sentence in the OAS search reply; '
            'either the search returned nothing or the page wording changed')

    return {'sequences': int(match.group(1).replace(',', '')),
            'studies': int(match.group(2).replace(',', ''))}


def parseDownloadUrls(html):
    """
    Extract the data unit download URLs from a search reply.

    OAS builds a shell script client side and stores it as a JavaScript array;
    these are the same URLs its bulk_download.sh would contain.

    Arguments:
      html (str): the search results page.

    Returns:
      list: absolute download URLs in page order.

    Raises:
      OasParseError: if the array is missing or contains no commands.
    """
    match = re.search(CSV_ARRAY_REGEX, html, re.DOTALL)
    if match is None:
        raise OasParseError(
            'no "%s" array in the OAS search reply; the download script is no '
            'longer embedded the way sourcerer expects' % CSV_ARRAY_MARKER)

    urls = re.findall(WGET_REGEX, match.group(1))
    if not urls:
        raise OasParseError('the OAS download script contained no wget commands')

    return [x.strip() for x in urls]


def unitIdFromUrl(url):
    """
    Derive the opaque unit identifier from a download URL.

    The identifier is the path below the collection directory, taken verbatim.
    It is never split into study, run or filename components: paired data uses
    several directory layouts and filename patterns, and most paired filenames
    contain no run accession, so any structured interpretation would be wrong for
    the majority of the collection.

    Arguments:
      url (str): an absolute data unit URL.

    Returns:
      tuple: (collection, unit_id).

    Raises:
      OasParseError: if the URL sits under no known collection directory.
    """
    parts = urlparse(url).path.strip('/').split('/')
    for collection in COLLECTIONS:
        if collection in parts:
            index = parts.index(collection)
            unit_id = '/'.join(parts[index + 1:])
            if not unit_id:
                break
            return collection, unit_id

    raise OasParseError(
        "cannot place '%s' under a known OAS collection (%s); the download URL "
        'layout has changed' % (url, ', '.join(COLLECTIONS)))


def urlFromUnitId(collection, unit_id):
    """
    Build the download URL for a unit identifier.

    Arguments:
      collection (str): 'paired' or 'unpaired'.
      unit_id (str): the opaque identifier.

    Returns:
      str: the absolute download URL.
    """
    return '%s%s/%s' % (DOWNLOAD_BASE, collection, unit_id)


def urlFromCatalogKey(key):
    """
    Map an unpaired catalog key to its download URL.

    Catalog keys are absolute paths on the OAS file server; replacing the mount
    prefix with the web root yields the public URL.

    Arguments:
      key (str): a key from the unpaired catalog JSON.

    Returns:
      str: the absolute download URL.

    Raises:
      OasParseError: if the key does not carry the expected prefix.
    """
    if not key.startswith(CATALOG_KEY_PREFIX):
        raise OasParseError(
            "catalog key '%s' does not start with '%s'; the rule mapping catalog "
            'keys to download URLs has changed' % (key, CATALOG_KEY_PREFIX))

    return DOWNLOAD_BASE + key[len(CATALOG_KEY_PREFIX):]


def parseSearchTable(html, collection='paired'):
    """
    Read the per unit metadata table from a search reply.

    Rows are matched to units through the detail link rather than by row order,
    and columns are read by header name rather than by position. Positional
    access is how the predecessor tool worked, and a single inserted column
    upstream would have silently relabelled every field.

    Arguments:
      html (str): the search results page.
      collection (str): the collection being searched.

    Returns:
      list: dicts of canonical field name to value, each including unit_id.

    Raises:
      OasParseError: if no results table is present.
    """
    soup = BeautifulSoup(html, 'html.parser')

    table, headers = None, None
    for candidate in soup.find_all('table'):
        names = [x.get_text().strip() for x in candidate.find_all('th')]
        if names and 'Details' in names:
            table, headers = candidate, names
            break

    if table is None:
        raise OasParseError('no results table found in the OAS %s search reply'
                            % collection)

    rows = []
    for row in table.find_all('tr'):
        cells = row.find_all('td')
        if not cells:
            continue

        link = row.find('a', href=re.compile(r'unit='))
        if link is None:
            continue

        record = {'unit_id': link['href'].split('unit=', 1)[1]}
        for name, cell in zip(headers, cells):
            if name == 'Details':
                continue
            record[FIELD_ALIASES.get(name, name)] = cell.get_text().strip()

        rows.append(record)

    if not rows:
        raise OasParseError('the OAS %s results table contained no data unit rows'
                            % collection)

    return rows


def parseDetailPage(html):
    """
    Read the fields a data unit's detail page carries but the results table lacks.

    BSource and BType are searchable on the paired form and are needed for the
    airrflow samplesheet, but the paired results table does not include them.

    Arguments:
      html (str): a dataunit detail page.

    Returns:
      dict: canonical field name to value for whatever the page exposes.

    Raises:
      OasParseError: if the page exposes no label/value rows at all.
    """
    soup = BeautifulSoup(html, 'html.parser')

    found = {}
    for row in soup.find_all('tr'):
        cells = row.find_all(['td', 'th'])
        if len(cells) < 2:
            continue
        # A row of nothing but header cells is the table's own heading, not data.
        if all(x.name == 'th' for x in cells[:2]):
            continue

        label = cells[0].get_text().strip().rstrip(':')
        value = cells[1].get_text().strip()
        if label and value:
            found[FIELD_ALIASES.get(label, label)] = value

    # A page that yields nothing is a layout change, not a unit that happens to
    # record no metadata: every detail page carries at least its own identifiers.
    # Returning an empty dict here would let the caller mark the unit enriched and
    # leave BSource and BType permanently blank, with nothing anywhere saying why.
    if not found:
        raise OasParseError(
            'the OAS data unit detail page exposed no label/value rows; the page '
            'layout has changed')

    return found


# ---------------------------------------------------------------------------
# Reading and normalizing data units
# ---------------------------------------------------------------------------

#: Chain suffixes used by paired data units.
CHAINS = ('heavy', 'light')

#: Matches a paired column and splits it into stem and chain.
CHAIN_COLUMN = re.compile(r'^(?P<stem>.+)_(?P<chain>heavy|light)$')

#: Paired identifiers are 10x barcodes with a contig suffix, as in
#: AAACCTGAGTCAATAG-1_contig_2. The barcode names the cell, the contig names the
#: chain, so removing the contig leaves the cell.
CONTIG_SUFFIX = re.compile(r'_contig_\d+$')

#: Width of the zero padded row counter in generated identifiers. Fixed rather
#: than derived from the unit's row count so that identifiers do not depend on
#: knowing the total in advance, which would force a counting pass over a
#: multi gigabyte file before conversion could start.
ID_WIDTH = 9

#: Columns consumed during normalization and not carried into the output.
CONSUMED_COLUMNS = frozenset(['Isotype', 'Redundancy'])

#: OAS writes single letter locus codes. AIRR requires the full gene locus, so
#: these are only a fallback for when v_call is empty.
LOCUS_LETTERS = {'H': 'IGH', 'K': 'IGK', 'L': 'IGL'}

#: AIRR boolean spellings accepted on input.
TRUE_TOKENS = frozenset(['T', 'TRUE', 'TRUE.', '1', 'YES', 'Y'])
FALSE_TOKENS = frozenset(['F', 'FALSE', 'FALSE.', '0', 'NO', 'N'])

#: Fields sourcerer adds. The prefix guarantees they cannot collide with a
#: current or future AIRR field name.
PROVENANCE_FIELDS = ('sourcerer_source', 'sourcerer_collection',
                     'sourcerer_unit_id', 'sourcerer_original_sequence_id',
                     'sourcerer_row_hash')


def unitStem(unit_id):
    """
    Build a filesystem safe, globally unique prefix for a unit's identifiers.

    The whole opaque identifier is used rather than just the filename. Paired
    filenames repeat across studies, so a shorter prefix would produce colliding
    identifiers once more than one unit is converted.

    Arguments:
      unit_id (str): the opaque unit identifier.

    Returns:
      str: the identifier prefix.
    """
    stem = re.sub(r'\.csv\.gz$', '', unit_id)

    return re.sub(r'[^A-Za-z0-9]+', '_', stem).strip('_')


def readDataUnit(path, chunksize=50000):
    """
    Open an OAS data unit and return its metadata and record chunks.

    The first line is a single quoted CSV field, and a quoted field may legally
    contain embedded newlines, so it is consumed with a csv.reader rather than by
    reading one physical line. The same handle is then passed to pandas, which
    continues at the header. Mixing iteration and reads on a text handle is well
    defined in Python 3, so the handoff is safe.

    Arguments:
      path (Path): the data unit file.
      chunksize (int): rows per chunk.

    Returns:
      tuple: (metadata dict, iterator of DataFrames).

    Raises:
      OasParseError: if the metadata line is missing or is not JSON.
    """
    handle = gzip.open(path, 'rt', newline='')
    try:
        reader = csv.reader(handle)
        try:
            first = next(reader)
        except StopIteration:
            raise OasParseError('%s is empty' % path)

        if len(first) != 1:
            raise OasParseError(
                '%s does not start with a single metadata field; got %d fields. '
                'The data unit layout has changed.' % (path, len(first)))

        try:
            metadata = json.loads(first[0])
        except ValueError as error:
            raise OasParseError(
                'the first record of %s is not JSON metadata (%s). The data unit '
                'layout has changed.' % (path, error))

        frames = pandas.read_csv(handle, chunksize=chunksize, dtype=str,
                                 na_filter=False)
    except Exception:
        handle.close()
        raise

    def chunks():
        try:
            yield from frames
        finally:
            handle.close()

    return metadata, chunks()


def toAirrBool(value):
    """
    Convert an OAS boolean spelling to an AIRR boolean.

    Arguments:
      value: the raw value.

    Returns:
      str: 'T', 'F', or '' when the value is absent or unrecognized.
    """
    text = str(value).strip().upper()
    if text in TRUE_TOKENS:
        return 'T'
    if text in FALSE_TOKENS:
        return 'F'

    return ''


def deriveLocus(v_call, fallback=''):
    """
    Determine the AIRR locus for a rearrangement.

    Taken from v_call rather than from the file's own locus column, which holds
    single letters such as H, K and L. Those are not valid AIRR locus values, and
    the paired Chain metadata is coarser still: it cannot tell IGK from IGL.

    Arguments:
      v_call (str): the V gene call, possibly a comma separated list.
      fallback (str): the file's locus column, used only when v_call is empty.

    Returns:
      str: an AIRR locus such as IGH, or '' when it cannot be determined.
    """
    if v_call:
        gene = str(v_call).split(',')[0].strip().upper()
        if len(gene) >= 3 and gene[:2] in ('IG', 'TR'):
            return gene[:3]

    letter = str(fallback).strip().upper()

    return LOCUS_LETTERS.get(letter, '')


def isotypeToCall(value):
    """
    Map an OAS Isotype value to an AIRR c_call.

    'Bulk' and 'All' are sentinels meaning the library was not isotype resolved.
    Writing them into c_call, as a straight copy would, invents a constant region
    call that the experiment never measured.

    Arguments:
      value (str): the Isotype value.

    Returns:
      str: the c_call, or '' when the isotype is unknown or a sentinel.
    """
    text = '' if value is None else str(value).strip()
    if not text or text in NON_ISOTYPES or isNull(text):
        return ''

    return text


def splitChains(frame):
    """
    Split a wide paired frame into one frame per chain.

    Every column in a paired data unit is suffixed, and each stem appears for
    both chains, so the split is total: one input row becomes exactly two output
    rows. Column sets differ between the two paired directory layouts, so the
    stems are discovered per file rather than assumed.

    Arguments:
      frame (pandas.DataFrame): the wide chunk.

    Returns:
      dict: chain name to a frame whose columns are the bare stems.

    Raises:
      OasParseError: if the columns are not symmetric across the two chains.
    """
    mapping = {x: {} for x in CHAINS}
    unsuffixed = []
    for column in frame.columns:
        match = CHAIN_COLUMN.match(column)
        if match is None:
            unsuffixed.append(column)
            continue
        mapping[match.group('chain')][match.group('stem')] = column

    if unsuffixed:
        raise OasParseError(
            'paired data unit has columns with no chain suffix (%s); the pivot '
            'assumption no longer holds' % ', '.join(sorted(unsuffixed)[:5]))

    heavy, light = set(mapping['heavy']), set(mapping['light'])
    if heavy != light:
        raise OasParseError(
            'paired chain columns are not symmetric; heavy only: %s, light only: '
            '%s' % (sorted(heavy - light)[:5], sorted(light - heavy)[:5]))

    return {chain: frame[list(cols.values())].rename(
                columns={v: k for k, v in cols.items()})
            for chain, cols in mapping.items()}


def cellBarcode(sequence_id):
    """
    Reduce a paired sequence identifier to the cell it came from.

    Arguments:
      sequence_id (str): an OAS paired identifier.

    Returns:
      str: the barcode, or '' when there is no identifier to reduce.
    """
    text = '' if sequence_id is None else str(sequence_id).strip()
    if not text:
        return ''

    return CONTIG_SUFFIX.sub('', text)


def rowHash(row):
    """
    Build a short content hash for a rearrangement.

    Lets a re-downloaded unit be checked row for row even if the upstream row
    order changed, which the positional identifier alone cannot do.

    Arguments:
      row (pandas.Series): a normalized row.

    Returns:
      str: the first 12 hex characters of a SHA-256 digest.
    """
    key = '|'.join(str(row.get(x, '')) for x in
                   ('sequence', 'v_call', 'j_call', 'junction'))

    return hashlib.sha256(key.encode('utf-8')).hexdigest()[:12]


def newReport():
    """
    Create a fresh conversion report.

    Returns:
      dict: zeroed counters, accumulated across chunks.
    """
    return {'rows_in': 0, 'rows_out': 0, 'missing_v_call': 0,
            'empty_sequence': 0, 'unresolved_locus': 0, 'missing_c_call': 0,
            'missing_duplicate_count': 0, 'cell_barcode_mismatch': 0,
            'loci': set()}


def normalizeChunk(metadata, chunk, unit_id, collection, offset=0, report=None,
                   prefix_ids=False):
    """
    Convert one chunk of an OAS data unit into long form AIRR records.

    Where the source supplies identifiers they are kept; where it does not, they
    are derived from the chunk's global offset rather than from the position of a
    row within its chunk, so converting a unit in chunks produces exactly the
    same output as converting it whole.

    Arguments:
      metadata (dict): the data unit's metadata line.
      chunk (pandas.DataFrame): raw records, all columns as strings.
      unit_id (str): the opaque unit identifier.
      collection (str): 'paired' or 'unpaired'.
      offset (int): index of this chunk's first row within the whole unit.
      report (dict): counters to accumulate into, from newReport().
      prefix_ids (bool): namespace identifiers with the unit stem, for output
        that combines several units into one file.

    Returns:
      pandas.DataFrame: normalized records with AIRR field names.
    """
    if report is None:
        report = newReport()

    report['rows_in'] += len(chunk)
    stem = unitStem(unit_id)

    if collection == 'paired':
        frame = _pairChunk(chunk, stem, offset, report)
        if prefix_ids:
            # Only paired identifiers need this: they are the source's own and
            # are unique only within a unit. Unpaired identifiers are synthesized
            # with the stem already in them.
            for column in ('sequence_id', 'cell_id'):
                frame[column] = stem + '_' + frame[column].astype(str)
    else:
        frame = chunk.copy()
        frame['_row'] = range(offset, offset + len(frame))
        # Unpaired units carry no sequence_id at all, so there is nothing to
        # preserve; paired units do, and it is kept in _pairChunk.
        frame['_source_sequence_id'] = ''
        frame['sequence_id'] = ['%s_%0*d' % (stem, ID_WIDTH, x)
                                for x in frame['_row']]

    return _finishChunk(frame, metadata, unit_id, collection, report)


def _chainBarcodes(frame, count):
    """
    Read one chain's barcodes, tolerating a layout that has no identifiers.

    Arguments:
      frame (pandas.DataFrame): one chain's records.
      count (int): how many rows to return.

    Returns:
      list: one barcode per row, '' where there is none.
    """
    if 'sequence_id' not in frame.columns:
        return [''] * count

    return [cellBarcode(x) for x in frame['sequence_id']]


def _cellIds(chains, stem, offset, report):
    """
    Resolve one cell identifier per input row.

    Heavy and light of the same cell must end up on the same cell_id, so it is
    resolved once from the row rather than derived separately from each chain's
    own identifier. Deriving it twice would split a cell in two whenever the two
    columns disagreed, which is a silent failure: nothing downstream can tell a
    split cell from a cell that genuinely had one chain.

    Arguments:
      chains (dict): chain name to that chain's records.
      stem (str): identifier prefix, used only by the fallback.
      offset (int): global index of the chunk's first row.
      report (dict): counters to accumulate into.

    Returns:
      list: one cell identifier per input row.
    """
    count = len(chains['heavy'])
    heavy = _chainBarcodes(chains['heavy'], count)
    light = _chainBarcodes(chains['light'], count)

    cells = []
    for index in range(count):
        first, second = heavy[index], light[index]
        if first and second and first != second:
            report['cell_barcode_mismatch'] += 1
        # A row with no identifier at all still needs one, and the row index is
        # the only thing left that is stable across chunk sizes.
        cells.append(first or second
                     or '%s_cell_%0*d' % (stem, ID_WIDTH, offset + index))

    return cells


def _pairChunk(chunk, stem, offset, report):
    """
    Pivot a wide paired chunk into two rows per cell.

    Arguments:
      chunk (pandas.DataFrame): the wide chunk.
      stem (str): identifier prefix for this unit.
      offset (int): global index of the chunk's first row.
      report (dict): counters to accumulate into.

    Returns:
      pandas.DataFrame: long form records carrying cell_id and sequence_id.
    """
    chains = splitChains(chunk)
    cells = _cellIds(chains, stem, offset, report)

    parts = []
    for chain in CHAINS:
        part = chains[chain].copy()
        part['_row'] = range(offset, offset + len(part))
        part['_cell'] = cells
        part['_chain'] = chain
        parts.append(part)

    frame = pandas.concat(parts, ignore_index=True)
    # Cells stay together and heavy always precedes light, so the output order is
    # a deterministic function of the input row index and not of chunking.
    frame['_rank'] = frame['_chain'].map({'heavy': 0, 'light': 1})
    frame = frame.sort_values(['_row', '_rank'], kind='stable')
    frame = frame.reset_index(drop=True).drop(columns=['_rank'])

    if 'sequence_id' in frame.columns:
        original = frame['sequence_id'].fillna('').astype(str)
    else:
        original = pandas.Series([''] * len(frame), index=frame.index,
                                 dtype=object)

    # The source identifier is kept verbatim. It is the real 10x barcode and
    # contig, it is what joins a row back to the file it came from, and a row
    # counter carries neither property. It is unique within a unit, which is what
    # one output file per unit requires; see OasSource.prefix_ids for combining.
    frame['_source_sequence_id'] = original
    frame['cell_id'] = frame['_cell']
    frame['sequence_id'] = [o if o else '%s_%s' % (c, x)
                            for o, c, x in zip(original, frame['_cell'],
                                               frame['_chain'])]

    if len(frame) != 2 * len(chunk):
        raise OasParseError(
            'paired pivot produced %d rows from %d input rows; expected exactly '
            'two per row' % (len(frame), len(chunk)))

    return frame


def _finishChunk(frame, metadata, unit_id, collection, report):
    """
    Apply the field mappings shared by both collections.

    Arguments:
      frame (pandas.DataFrame): records after any pivot.
      metadata (dict): the data unit's metadata line.
      unit_id (str): the opaque unit identifier.
      collection (str): 'paired' or 'unpaired'.
      report (dict): counters to accumulate into.

    Returns:
      pandas.DataFrame: the normalized chunk.
    """
    frame = frame.copy()

    if '_source_sequence_id' not in frame.columns:
        frame['_source_sequence_id'] = ''

    # duplicate_count: only some layouts carry Redundancy.
    if 'Redundancy' in frame.columns:
        counts = pandas.to_numeric(frame['Redundancy'], errors='coerce')
        frame['duplicate_count'] = counts.fillna(1).astype(int)
    else:
        frame['duplicate_count'] = 1
        report['missing_duplicate_count'] += len(frame)

    # c_call: per chain Isotype where the layout has it, otherwise the unit level
    # Isotype. Sentinels never become a call.
    if 'Isotype' in frame.columns:
        frame['c_call'] = frame['Isotype'].map(isotypeToCall)
    else:
        frame['c_call'] = isotypeToCall(metadata.get('Isotype'))
    report['missing_c_call'] += int((frame['c_call'] == '').sum())

    # locus: always recomputed from v_call. The file's own locus column holds
    # single letters (H, K, L), which are not valid AIRR locus values.
    blank = pandas.Series([''] * len(frame), index=frame.index, dtype=object)
    calls = frame['v_call'].fillna('') if 'v_call' in frame.columns else blank
    letters = frame['locus'].fillna('') if 'locus' in frame.columns else blank
    frame['locus'] = [deriveLocus(v, f) for v, f in zip(calls, letters)]

    for column in ('stop_codon', 'vj_in_frame', 'productive', 'rev_comp',
                   'complete_vdj', 'v_frameshift'):
        if column in frame.columns:
            frame[column] = frame[column].map(toAirrBool)

    if 'v_call' in frame.columns:
        report['missing_v_call'] += int((frame['v_call'].fillna('') == '').sum())
    if 'sequence' in frame.columns:
        report['empty_sequence'] += int((frame['sequence'].fillna('') == '').sum())
    report['unresolved_locus'] += int((frame['locus'] == '').sum())
    # Collected here so the samplesheet can derive pcr_target_locus from what the
    # data actually contains rather than from an assumption about the source.
    report['loci'].update(x for x in frame['locus'].unique() if x)

    frame['repertoire_id'] = unit_id
    frame['sourcerer_source'] = 'oas'
    frame['sourcerer_collection'] = collection
    frame['sourcerer_unit_id'] = unit_id
    # Recorded only when sequence_id is not already the source's own value.
    # Repeating an identical value in a second column of every row is noise, not
    # provenance; a value here means the identifier was rewritten.
    source_ids = frame['_source_sequence_id'].astype(str)
    frame['sourcerer_original_sequence_id'] = source_ids.where(
        source_ids != frame['sequence_id'].astype(str), '')
    frame['sourcerer_row_hash'] = frame.apply(rowHash, axis=1)

    drop = [x for x in frame.columns
            if x in CONSUMED_COLUMNS or x.startswith('_')]
    frame = frame.drop(columns=drop)

    frame = coerceAirrTypes(frame)

    report['rows_out'] += len(frame)

    return frame


# ---------------------------------------------------------------------------
# The source
# ---------------------------------------------------------------------------

class OasSource(SourceBase):
    """
    The Observed Antibody Space source.
    """

    name = 'oas'
    description = 'Observed Antibody Space: cleaned, annotated antibody repertoires'
    homepage = 'https://opig.stats.ox.ac.uk/webapps/oas/'
    collections = COLLECTIONS
    collection_help = COLLECTION_HELP

    #: OAS distributes its data under CC BY 4.0, per its homepage; the two
    #: papers below are what it asks to be cited in exchange for that license.
    license = 'CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)'
    citation = (
        'Kovaltsuk A, Leem J, Kelm S, Snowden J, Deane CM, Krawczyk K. '
        'Observed Antibody Space: A Resource for Data Mining Next-Generation '
        'Sequencing of Antibody Repertoires. J Immunol. 2018;201(8):2502-2509. '
        'doi:10.4049/jimmunol.1800708',
        'Olsen TH, Boyles F, Deane CM. Observed Antibody Space: A diverse '
        'database of cleaned, annotated, and translated unpaired and paired '
        'antibody sequences. Protein Sci. 2022;31(1):141-146. '
        'doi:10.1002/pro.4205',
    )

    #: Namespace generated identifiers with the unit stem. Off by default: one
    #: output file per unit needs no prefix, and the source's own barcodes are
    #: more useful bare. Anything writing several units into one file must turn
    #: this on. 10x barcodes are drawn from a fixed whitelist and therefore recur
    #: in every unit, so combining units without a prefix silently merges
    #: unrelated cells rather than failing.
    prefix_ids = False

    def formUrl(self, collection):
        """
        Return the search form URL for a collection.

        Arguments:
          collection (str): 'paired' or 'unpaired'.

        Returns:
          str: the form URL.
        """
        return PAIRED_FORM_URL if collection == 'paired' else UNPAIRED_FORM_URL

    def harvestSchema(self):
        """
        Fetch both search forms and build a fresh snapshot.

        Returns:
          SourceSchema: the harvested snapshot.
        """
        from datetime import datetime

        from sourcerer.Schema import Collection, Field, SourceSchema
        from sourcerer.Version import __version__

        collections = {}
        for collection in self.collections:
            html = self.client.get(self.formUrl(collection)).text
            fields = tuple(
                Field(name=x['name'], values=tuple(x['values']),
                      wildcard=x['wildcard'], pseudo_values=x['pseudo_values'])
                for x in parseFormSchema(html, collection))
            collections[collection] = Collection(name=collection, fields=fields)

        return SourceSchema(
            source=self.name,
            harvested=datetime.now(UTC).strftime('%Y-%m-%dT%H:%M:%SZ'),
            harvested_by='sourcerer %s' % __version__,
            source_urls={'paired_form': PAIRED_FORM_URL,
                         'unpaired_form': UNPAIRED_FORM_URL,
                         'catalog': CATALOG_URL,
                         'download_base': DOWNLOAD_BASE},
            url_rules={'catalog_key_prefix': CATALOG_KEY_PREFIX,
                       'download_prefix': DOWNLOAD_BASE},
            parse_contracts={'count_regex': COUNT_REGEX,
                             'csv_array_marker': CSV_ARRAY_MARKER,
                             'detail_link_pattern': r'\.\./dataunit_(paired|unpaired)\?unit='},
            field_aliases=dict(FIELD_ALIASES),
            collections=collections)

    def submitSearch(self, collection, filters):
        """
        Submit the search form and return the reply.

        The form is multipart encoded and needs an explicit wildcard for every
        field; an empty string is not accepted as "all".

        Arguments:
          collection (str): which collection to search.
          filters (dict): resolved field to value pairs.

        Returns:
          str: the HTML reply.
        """
        url = self.formUrl(collection)
        payload = {k: (None, v) for k, v in filters.items()}
        response = self.client.post(url, files=payload, headers={'Referer': url})

        return response.text

    def harvestCatalog(self, collection, schema=None):
        """
        Build a catalog of every data unit in a collection.

        Paired data has no published index, so this submits an unfiltered search
        and reads the download script out of the reply. Unpaired has a JSON
        catalog and is handled separately.

        Arguments:
          collection (str): which collection to catalog.
          schema (SourceSchema): the schema to take wildcards from. Passed
            explicitly during a refresh, when the freshly harvested schema is
            newer than any packaged one and may be the only one that exists.

        Returns:
          list: catalog rows.
        """
        if collection == 'unpaired':
            return self._harvestUnpairedCatalog()

        schema = schema if schema is not None else self.schema
        wildcards = {x.name: x.wildcard
                     for x in schema.getCollection(collection).fields}
        html = self.submitSearch(collection, wildcards)

        urls = parseDownloadUrls(html)
        rows = {}
        for url in urls:
            found, unit_id = unitIdFromUrl(url)
            rows[unit_id] = {'unit_id': unit_id, 'collection': found, 'url': url,
                             'dir_segment': url.split('/')[-2],
                             'study': unit_id.split('/')[0]}

        for record in parseSearchTable(html, collection):
            row = rows.get(record['unit_id'])
            if row is None:
                continue
            row['n_unique_sequences'] = record.get('Unique sequences', '')
            for name in ('Species', 'Isotype', 'Chain', 'Disease', 'Vaccine',
                         'Subject', 'Age', 'Longitudinal'):
                if name in record:
                    row[name] = record[name]

        return list(rows.values())

    def _harvestUnpairedCatalog(self):
        """
        Build the unpaired catalog from the published JSON index.

        Returns:
          list: catalog rows.
        """
        payload = self.client.get(CATALOG_URL).json()

        rows = []
        for key, meta in payload.items():
            url = urlFromCatalogKey(key)
            collection, unit_id = unitIdFromUrl(url)
            row = {'unit_id': unit_id, 'collection': collection, 'url': url,
                   'dir_segment': url.split('/')[-2],
                   'study': unit_id.split('/')[0],
                   'run': meta.get('Run', ''),
                   'n_unique_sequences': meta.get('Unique sequences', ''),
                   'Author': meta.get('Author', '')}
            for name in ('Species', 'Isotype', 'Chain', 'Disease', 'Vaccine',
                         'Subject', 'Age', 'Longitudinal', 'BSource', 'BType'):
                row[name] = meta.get(name, '')
            # The JSON index carries everything the detail pages would add.
            row['detail_status'] = DETAIL_OK
            rows.append(row)

        return rows

    def enrichCatalog(self, rows, limit=None, force=False):
        """
        Fill in the fields only a unit's detail page carries.

        The paired results table has no BSource or BType, but the paired form
        filters on both and the samplesheet needs them. Units are selected by
        recorded status rather than by novelty, so a page that failed once is
        retried later instead of staying blank forever.

        Arguments:
          rows (list): catalog rows, modified in place.
          limit (int): stop after this many fetches.
          force (bool): re-read every unit's detail page, including those already
            recorded as read. Costs one request per unit, so it is for recovering
            from a page layout change rather than for routine use.

        Returns:
          int: how many units were successfully enriched.
        """
        from datetime import datetime

        pending = list(rows) if force else [x for x in rows if needsDetail(x)]
        if limit is not None:
            pending = pending[:limit]

        stamp = datetime.now(UTC).strftime('%Y-%m-%dT%H:%M:%SZ')
        enriched = 0
        for row in pending:
            url = '%s?unit=%s' % (DETAIL_URL % row['collection'], row['unit_id'])
            row['detail_attempted'] = stamp
            try:
                found = parseDetailPage(self.client.get(url).text)
            except Exception as error:
                # Never destructive: whatever a previous run learned stays.
                log.warning('detail page for %s failed (%s); will retry',
                            row['unit_id'], error)
                row['detail_status'] = 'failed'
                continue

            for name in ('BSource', 'BType', 'Author'):
                if found.get(name):
                    row[name] = found[name]
            row['detail_status'] = DETAIL_OK
            enriched += 1

        return enriched

    def catalogPath(self, collection):
        """
        Return the packaged catalog location for a collection.

        Arguments:
          collection (str): the collection.

        Returns:
          Path: the catalog file inside the installed package.
        """
        from importlib import resources

        anchor = resources.files('sourcerer').joinpath(
            'data/schemas', self.name, '%s_catalog.tsv' % collection)

        return Path(str(anchor))

    def searchUnits(self, query):
        """
        Resolve a query to data units using the packaged catalog.

        Arguments:
          query (Query): the validated request.

        Returns:
          list: DataUnit objects, ordered by identifier.
        """
        rows = loadCatalog(self.catalogPath(query.collection))
        if not rows:
            raise OasParseError(
                "no packaged catalog for OAS %s; run 'sourcerer schema refresh "
                "--source oas'" % query.collection)

        selected = filterCatalog(rows, query.filters)
        if query.limit is not None:
            selected = selected[:query.limit]

        units = []
        for row in selected:
            counts = row.get('n_unique_sequences') or ''
            units.append(DataUnit(
                unit_id=row['unit_id'], collection=row['collection'],
                url=row['url'], metadata=dict(row),
                n_sequences=int(counts) if counts.isdigit() else None))

        return units

    def readUnit(self, path, unit, chunksize=50000):
        """
        Open a downloaded unit.

        Arguments:
          path (Path): the downloaded file.
          unit (DataUnit): what it is.
          chunksize (int): rows per chunk.

        Returns:
          tuple: (metadata dict, iterator of raw record chunks).
        """
        return readDataUnit(path, chunksize=chunksize)

    def normalizeChunk(self, metadata, chunk, unit, offset, report):
        """
        Map one chunk of raw records to AIRR named records.

        Arguments:
          metadata (dict): the unit's metadata.
          chunk (pandas.DataFrame): raw records.
          unit (DataUnit): what they came from.
          offset (int): index of the chunk's first row within the whole unit.
          report (dict): counters to accumulate into.

        Returns:
          pandas.DataFrame: normalized records.
        """
        return normalizeChunk(metadata, chunk, unit.unit_id, unit.collection,
                              offset, report, prefix_ids=self.prefix_ids)
