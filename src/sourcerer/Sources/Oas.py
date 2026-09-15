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
import os
import re
from datetime import UTC
from pathlib import Path
from urllib.parse import urlparse

import pandas
from bs4 import BeautifulSoup

# Sourcerer imports
from sourcerer.Catalog import DETAIL_OK, filterCatalog, loadCatalog, needsDetail
from sourcerer.Commandline import CommonHelpFormatter
from sourcerer.Convert import coerceAirrTypes
from sourcerer.Exceptions import OasParseError, SourcererError
from sourcerer.Http import HttpClient
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

#: Every searchable field the snapshot may carry, mapped to what sourcerer does
#: with it: the AIRR or samplesheet column it feeds, or an explicit note that it
#: is a search filter only. A contract test asserts every field in the packaged
#: snapshot appears here, so when OAS adds a field the monthly refresh PR fails
#: CI with its name instead of silently ignoring it.
KNOWN_FIELDS = {
    'Species': 'samplesheet species',
    'Age': 'samplesheet age',
    'BSource': 'samplesheet tissue',
    'BType': 'samplesheet cell_subset',
    'Vaccine': 'samplesheet intervention',
    'Disease': 'samplesheet disease_diagnosis',
    'Subject': 'samplesheet subject_id',
    'Longitudinal': 'samplesheet longitudinal',
    'Chain': 'search filter only; locus is derived from v_call per row',
    'Isotype': 'c_call, from unit metadata when no per chain Isotype column',
    'Primer': 'search filter only; no output column',
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


def clean(value, default=''):
    """
    Normalize an OAS metadata value, mapping its null sentinels to a default.

    Arguments:
      value: the raw value.
      default (str): what to use when the value carries no information.

    Returns:
      str: the cleaned value.
    """
    if isNull(value):
        return default

    return str(value).strip()


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
        # Detail pages escape a comma inside a value (e.g. BType's
        # 'Plasmablasts\, Memory B cells and activated T cells') the same way
        # the search form does, so it must be undone the same way: unescaped,
        # a value like that can never match --btype's validated (unescaped)
        # filter, and the exact-match lookup in Catalog.filterCatalog silently
        # returns zero hits -- the failure mode the snapshot exists to prevent.
        value = unescapeOption(cells[1].get_text())
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
# Data contracts and catalog fingerprints
# ---------------------------------------------------------------------------

#: AIRR stems every data unit layout must carry for conversion to work. Declared
#: rather than derived so a layout that loses one shows up as a contract change.
REQUIRED_AIRR_STEMS = ('sequence', 'locus', 'v_call', 'j_call', 'junction',
                       'sequence_alignment', 'germline_alignment')

#: Columns whose presence or absence is what distinguishes the known layouts:
#: the 158 paired csv/ units genuinely lack Redundancy, c_region and Isotype,
#: and unpaired units carry no sequence_id at all. Recorded per probe unit as
#: required_columns and absent_columns so a layout gaining or losing one is a
#: visible contract change rather than a silent conversion difference.
CONTRACT_MARKERS = ('Redundancy', 'c_region', 'Isotype', 'sequence_id',
                    'duplicate_count')

#: How many complete CSV records a contract probe must decode before it stops
#: extending its byte range.
PROBE_RECORDS = 2

#: Cap on distinct values per key recorded in the catalog fingerprint. A key
#: like Run is unique per unit, and enumerating fifteen thousand values would
#: bloat the fingerprint without describing anything the drift check compares.
MAX_FINGERPRINT_VALUES = 500


def filenamePattern(unit_id):
    """
    Normalize a unit's filename to its layout pattern.

    Digit runs are collapsed so that every run accession maps to the same
    pattern. The result is descriptive bookkeeping for the snapshot, never
    parsed back: paths stay opaque, and a novel pattern is additive drift.

    Arguments:
      unit_id (str): the opaque unit identifier.

    Returns:
      str: the pattern, e.g. 'SRR<n>_paired' or '<n>_S<n>__<n>_Paired_All'.
    """
    name = unit_id.rsplit('/', 1)[-1]
    name = re.sub(r'\.csv\.gz$', '', name)

    return re.sub(r'\d+', '<n>', name)


def pathLayouts(rows):
    """
    Summarize the directory and filename layouts a catalog contains.

    Arguments:
      rows (list): catalog rows.

    Returns:
      dict: observed directory segments and filename patterns, each with unit
      counts.
    """
    segments = {}
    patterns = {}
    for row in rows:
        segment = row.get('dir_segment', '')
        segments[segment] = segments.get(segment, 0) + 1
        pattern = filenamePattern(row['unit_id'])
        patterns[pattern] = patterns.get(pattern, 0) + 1

    return {'observed_dir_segments': dict(sorted(segments.items())),
            'observed_filename_patterns': dict(sorted(patterns.items()))}


def probeComplete(raw):
    """
    Decide whether a ranged probe has fetched enough of a data unit.

    Enough means the metadata member has decoded completely (a second member has
    started) and the CSV member contains the header plus PROBE_RECORDS complete
    lines. A fixed byte window would stop sufficing the moment the metadata or
    header grew, so completeness is judged on structure rather than on size.

    Arguments:
      raw (bytes): the accumulated prefix of the remote file.

    Returns:
      bool: True once the prefix contains what parseProbeFacts needs.
    """
    from sourcerer.Gzip import splitMembers

    members = splitMembers(raw)
    if len(members) < 2:
        return False

    text = members[1].decode('utf-8', errors='replace')

    return text.count('\n') >= 1 + PROBE_RECORDS


def jsonTypeName(value):
    """
    Name a JSON value's type for the contract record.

    Arguments:
      value: a value decoded from JSON.

    Returns:
      str: one of 'null', 'bool', 'int', 'float', 'str', 'list', 'dict'.
    """
    if value is None:
        return 'null'
    # bool subclasses int, so it has to be tested first.
    for kind, name in ((bool, 'bool'), (int, 'int'), (float, 'float'),
                       (str, 'str'), (list, 'list'), (dict, 'dict')):
        if isinstance(value, kind):
            return name

    return type(value).__name__


def parseProbeFacts(raw, unit_id, collection):
    """
    Extract the file format facts from a probed data unit prefix.

    Arguments:
      raw (bytes): the leading bytes of the remote file.
      unit_id (str): which unit was probed.
      collection (str): 'paired' or 'unpaired'.

    Returns:
      dict: the contract facts for this unit's layout.

    Raises:
      OasParseError: if the prefix does not have the expected structure.
    """
    import io

    from sourcerer.Gzip import splitMembers

    members = splitMembers(raw)
    if len(members) < 2:
        raise OasParseError(
            'probe of %s decoded %d gzip member(s) where 2 were expected '
            '(metadata, then CSV); the data unit framing has changed'
            % (unit_id, len(members)))

    meta_text = members[0].decode('utf-8', errors='replace')
    record = next(csv.reader(io.StringIO(meta_text)))
    if len(record) != 1:
        raise OasParseError(
            'probe of %s: the metadata member holds %d CSV fields where 1 was '
            'expected' % (unit_id, len(record)))

    try:
        metadata = json.loads(record[0])
    except ValueError as error:
        raise OasParseError(
            'probe of %s: the metadata member is not JSON (%s)'
            % (unit_id, error))

    columns = next(csv.reader(io.StringIO(
        members[1].decode('utf-8', errors='replace'))))

    facts = {
        'unit_id': unit_id,
        'gzip_members': len(members),
        'metadata_keys': sorted(metadata),
        'metadata_value_types': {k: jsonTypeName(v) for k, v in metadata.items()},
        'n_columns': len(columns),
        'columns': sorted(columns),
        'required_airr_stems': list(REQUIRED_AIRR_STEMS),
    }

    if collection == 'paired':
        stems = {chain: set() for chain in CHAINS}
        unsuffixed = []
        for column in columns:
            match = CHAIN_COLUMN.match(column)
            if match is not None:
                stems[match['chain']].add(match['stem'])
            else:
                unsuffixed.append(column)
        shared = stems['heavy'] & stems['light']
        facts['suffix_pairing'] = {
            'suffixes': ['_%s' % x for x in CHAINS],
            'unsuffixed_columns': sorted(unsuffixed),
            'paired_stems': len(shared),
            'unmatched_stems': sorted(stems['heavy'] ^ stems['light'])}
        present = shared
    else:
        present = set(columns)

    facts['required_columns'] = sorted(x for x in CONTRACT_MARKERS
                                       if x in present)
    facts['absent_columns'] = sorted(x for x in CONTRACT_MARKERS
                                     if x not in present)

    return facts


def chooseProbeUnits(rows, pinned):
    """
    Select one probe unit per path layout, honoring existing pins.

    A pinned unit still present in the catalog is kept: re-choosing every run
    would churn the snapshot diff and untie the recorded contract from the file
    it was measured on. A fresh unit is chosen only for a layout with no live
    pin, preferring the smallest by sequence count so the probe stays cheap; the
    replacement then appears in the snapshot diff, where a reviewer sees it.

    Arguments:
      rows (list): catalog rows.
      pinned (iterable): unit_ids pinned by the existing contracts.

    Returns:
      dict: dir_segment to the chosen catalog row, in sorted segment order.
    """
    def size(row):
        counts = str(row.get('n_unique_sequences') or '')
        return (0, int(counts)) if counts.isdigit() else (1, 0)

    chosen = {}
    pins = set(pinned)
    for row in rows:
        segment = row.get('dir_segment', '')
        if row['unit_id'] in pins and segment not in chosen:
            chosen[segment] = row

    groups = {}
    for row in rows:
        groups.setdefault(row.get('dir_segment', ''), []).append(row)
    for segment, group in groups.items():
        if segment not in chosen:
            chosen[segment] = min(group, key=lambda x: (size(x), x['unit_id']))

    return dict(sorted(chosen.items()))


def buildFingerprint(content, headers, payload):
    """
    Condense the raw unpaired catalog into the facts the drift check compares.

    The unpaired catalog is a single ~7 MB JSON document, one entry per data
    unit. Diffing it whole every month would be slow to fetch and unreadable
    to review, so this reduces it once, at harvest time, to the handful of
    aggregate facts Drift.py actually looks at, written to
    catalog_fingerprint.json and compared against the version committed at
    the last drift check.

    Returned fields, and what each is for:

    - ``collection``: the collection the fingerprinted rows belong to
      ('unpaired' today, since that is the only machine readable index OAS
      publishes), or None if the payload named more than one -- see the
      comment on that key below. Drift.py labels its findings from this
      rather than assuming a name.
    - ``sha256`` / ``etag`` / ``last_modified``: cheap validators for "did the
      document change at all", checked before the more expensive comparisons
      below run.
    - ``n_units``: the catalog's total row count; growth or shrinkage between
      two fingerprints is reported directly from this.
    - ``key_counts``: how many units carry each metadata key. A key whose
      count is neither 0 nor n_units is present on a strict subset of units --
      by definition partial -- which is the anomaly the drift check reports
      (the real world case: a stray 'Organism' key on exactly 1 of 15,631
      unpaired units).
    - ``value_types``: the distinct JSON types seen for each key's values
      (see jsonTypeName), so a key silently changing shape (e.g. a count
      written as a string instead of a number) is caught even though the key
      itself did not change.
    - ``value_counts``: how many units carry each distinct value of each key,
      capped per key at MAX_FINGERPRINT_VALUES entries -- past that (an
      accession-like key with one value per unit, say) the key's entry is
      None instead, since enumerating unique-per-unit values would make the
      fingerprint scale with the catalog rather than stay a summary.
    - ``study_index``: unit count per study, both a human readable summary and
      what a drift finding's "N new units in <study>" grouping is built from.
    - ``sample_unit``: one verbatim catalog entry (the first key in sorted
      order), kept as a live example of the document's actual shape and as a
      lightweight parser fixture.

    Arguments:
      content (bytes): the catalog document as served.
      headers: the response headers.
      payload (dict): the parsed catalog.

    Returns:
      dict: the fingerprint, with the fields described above.
    """
    key_counts = {}
    value_types = {}
    value_counts = {}
    studies = {}
    collections = set()
    for key, meta in payload.items():
        collection, unit_id = unitIdFromUrl(urlFromCatalogKey(key))
        collections.add(collection)
        study = unit_id.split('/')[0]
        studies[study] = studies.get(study, 0) + 1
        for name, value in meta.items():
            key_counts[name] = key_counts.get(name, 0) + 1
            value_types.setdefault(name, set()).add(jsonTypeName(value))
            counts = value_counts.setdefault(name, {})
            if counts is not None:
                counts[str(value)] = counts.get(str(value), 0) + 1
                if len(counts) > MAX_FINGERPRINT_VALUES:
                    value_counts[name] = None

    return {
        # The Drift findings this fingerprint feeds label themselves from this
        # field rather than assuming a fixed collection name. Today this is
        # always 'unpaired' -- the unpaired catalog is the only one OAS
        # publishes as a machine readable index -- but a mixed document, were
        # one ever fed in, is a fact worth recording rather than guessing a
        # label from, so it is left unset (None) rather than picking one side.
        'collection': next(iter(collections)) if len(collections) == 1 else None,
        'sha256': hashlib.sha256(content).hexdigest(),
        'n_units': len(payload),
        'etag': headers.get('ETag', ''),
        'last_modified': headers.get('Last-Modified', ''),
        'key_counts': dict(sorted(key_counts.items())),
        'value_types': {k: sorted(v) for k, v in sorted(value_types.items())},
        # None marks a key with more distinct values than the fingerprint
        # enumerates, such as per unit accessions.
        'value_counts': {k: (dict(sorted(v.items())) if v is not None else None)
                         for k, v in sorted(value_counts.items())},
        'study_index': dict(sorted(studies.items())),
        'sample_unit': {'key': min(payload), 'metadata': payload[min(payload)]},
    }


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
# The airrflow samplesheet row
# ---------------------------------------------------------------------------

def countUnresolvedSubjects(entries):
    """
    Count data units whose Subject metadata is one of OAS's own null sentinels.

    A unit like this writes that raw sentinel ('no', 'None', ...) into the
    samplesheet's subject_id column rather than a real identifier -- see
    samplesheetRow's own comment on why the value is kept raw rather than
    collapsed to a placeholder. `sourcerer oas verify` is what turns those
    into real evidence from NCBI; this count is what a caller uses to decide
    whether it is worth telling the user to run it.

    Arguments:
      entries (list): (DataUnit, Path) pairs, the same shape
        Airrflow.buildSamplesheet takes; only the unit's metadata is read
        here.

    Returns:
      int: how many units carry no recorded subject.
    """
    return sum(1 for unit, _ in entries
              if isNull((unit.metadata or {}).get('Subject')))


def samplesheetRow(unit):
    """
    Map an OAS data unit's metadata to airrflow samplesheet columns.

    Arguments:
      unit (DataUnit): the converted unit.

    Returns:
      dict: the samplesheet columns OAS metadata can fill in.
    """
    metadata = unit.metadata or {}

    # Subject is passed through raw rather than via clean(): a value like OAS's
    # own "no" carries real information (subject identity was not recorded) and
    # must not be collapsed and then replaced by the study name, which would
    # falsely tell airrflow that every otherwise-unidentified unit in the study
    # is the same subject. The study is used only when OAS supplies no value
    # for Subject at all.
    raw_subject = metadata.get('Subject')
    subject = (str(raw_subject).strip() if raw_subject not in (None, '') else '')
    subject = subject or clean(metadata.get('study'))

    return {
        'subject_id': subject.replace(' ', '_'),
        'species': clean(metadata.get('Species'), 'human').lower(),
        'tissue': clean(metadata.get('BSource'), 'unknown'),
        # Not derivable from OAS, but airrflow requires the column to be
        # populated and asks for NA when it is unknown. NA is a placeholder
        # here too, so a hand-edited value still survives a later merge.
        'sex': 'NA',
        'age': clean(metadata.get('Age'), 'NA'),
        'biomaterial_provider': clean(metadata.get('Author'),
                                      clean(metadata.get('study'))),
        # Driven by the collection rather than hardcoded: only paired data is
        # single cell, and the R implementation assumed TRUE because it only
        # ever handled paired.
        'single_cell': 'TRUE' if unit.collection == 'paired' else 'FALSE',
        'disease_diagnosis': clean(metadata.get('Disease')),
        'intervention': clean(metadata.get('Vaccine')),
        # Like Age, OAS records this as a presence flag ("no" when the
        # study carries no longitudinal design), so the same null-token
        # collapse to 'NA' applies, unlike Subject's raw pass-through.
        'longitudinal': clean(metadata.get('Longitudinal'), 'NA'),
        'cell_subset': clean(metadata.get('BType')),
        'study': clean(metadata.get('study')) or unit.study,
    }


# ---------------------------------------------------------------------------
# oas verify: cross-referencing unresolved subjects against NCBI
# ---------------------------------------------------------------------------

#: The columns `sourcerer oas verify` adds, appended after whatever columns
#: the input samplesheet already had (never inserted among them) -- so the
#: report is the input samplesheet plus evidence, and can be used in its
#: place as airrflow input, rather than a separate file missing the columns
#: airrflow actually reads (`filename`, `species`, `pcr_target_locus`, ...).
NCBI_EVIDENCE_COLUMNS = ('ncbi_sample_name', 'ncbi_subject_suggested', 'status',
                         'subject_check', 'accession', 'biosample_accession',
                         'biosample_url', 'pooled_codes')


def _addVerifyAction(actions):
    """
    Add oas verify: cross-reference a samplesheet's unresolved subjects
    against NCBI.

    Not a search/download action: it takes no COLLECTION or filter flags,
    since it reads a samplesheet already on disk rather than querying OAS.
    Every row with subject_id 'no' or 'None' (OAS's own null sentinels; see
    isNull above) names an SRA run or GEO sample accession in its
    sample_name column, and that accession's NCBI BioSample record usually
    names the sample plainly enough for a human to read off the subject.

    Deliberately two flags, not the larger surface an earlier version of this
    command had (--apply, --use-suggestion, --evidence-out, --limit): the
    report is the one thing this command produces, always with both NCBI's raw
    sample name and a suggested subject, so there was nothing left for a flag
    to switch between. It carries every column the input samplesheet had, evidence
    columns appended, so it is a drop-in airrflow input rather than a
    side file -- see buildEvidenceRow and NCBI_EVIDENCE_COLUMNS.

    Arguments:
      actions: the oas subcommand's action subparsers (search and download
        are already on it).
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
                    'NCBI columns instead of a guessed single '
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
    from sourcerer.Ncbi import (
        DEFAULT_DELAY,
        KEYED_DELAY,
        accessionFromText,
        gatherEvidence,
    )

    fields, rows = readSamplesheetRows(args.samplesheet)
    warnReusedSubjects(rows)

    # Every row with a readable accession is looked up, not only rows OAS
    # left unresolved: a subject_id OAS did record is still worth checking
    # against NCBI's own record (see subjectCheck), so there is no 'pending'
    # subset here to restrict the lookup to.
    accession_by_sample = {}
    for row in rows:
        accession = accessionFromText(row.get('sample_name', ''))
        if accession is not None:
            accession_by_sample[row['sample_id']] = accession

    accessions = set(accession_by_sample.values())
    api_key = args.ncbi_api_key or os.environ.get('NCBI_API_KEY')
    delay = KEYED_DELAY if api_key else DEFAULT_DELAY
    client = HttpClient(delay=delay)
    evidence = (gatherEvidence(client, accessions, api_key=api_key)
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

    #: Fingerprint of the raw unpaired catalog, set as a side effect of
    #: harvesting it, so that a refresh fingerprints the same bytes it indexed.
    catalog_fingerprint = None

    catalog_columns = ('unit_id', 'collection', 'url', 'dir_segment', 'study',
                       'run', 'n_unique_sequences', 'Species', 'Isotype',
                       'Chain', 'Disease', 'Vaccine', 'Subject', 'Age',
                       'Longitudinal', 'BSource', 'BType', 'Author',
                       'detail_status', 'detail_attempted')
    search_columns = ('Species', 'Disease', 'Subject', 'BSource')
    enrichment_columns = ('BSource', 'BType', 'Author')

    def newReport(self):
        """See SourceBase.newReport."""
        return newReport()

    def samplesheetRow(self, unit):
        """See SourceBase.samplesheetRow."""
        return samplesheetRow(unit)

    def countUnresolvedSubjects(self, entries):
        """See SourceBase.countUnresolvedSubjects."""
        return countUnresolvedSubjects(entries)

    @classmethod
    def addActions(cls, actions):
        """Add oas verify. See SourceBase.addActions."""
        _addVerifyAction(actions)

        return {'verify': handleOasVerify}

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

        Also fingerprints the raw document while it is in hand, so that the
        drift check can compare the catalog's shape without re-fetching it.

        Returns:
          list: catalog rows.
        """
        response = self.client.get(CATALOG_URL)
        payload = json.loads(response.content)
        self.catalog_fingerprint = buildFingerprint(response.content,
                                                    response.headers, payload)

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

            for name in self.enrichment_columns:
                if found.get(name):
                    row[name] = found[name]
            row['detail_status'] = DETAIL_OK
            enriched += 1

        return enriched

    def harvestContracts(self, catalogs, existing=None):
        """
        Probe pinned data units and record the downloaded file format.

        One unit per path layout is fetched by progressive byte ranges, decoded
        far enough to see the metadata member and the CSV header, and reduced to
        the facts conversion depends on. Collections not harvested this run keep
        their existing entries, so a partial refresh cannot silently drop a
        contract.

        Arguments:
          catalogs (dict): collection name to catalog rows.
          existing (dict): the stored contracts, for probe unit pins.

        Returns:
          dict: the data contracts payload.

        Raises:
          ProbeIncompleteError: if a probe hit its byte cap. This is a harvest
            failure, not drift; a slow or truncated response must not read as a
            format change.
          OasParseError: if a probed prefix does not have the expected shape.
        """
        from datetime import datetime

        from sourcerer.Contracts import CONTRACTS_VERSION
        from sourcerer.Version import __version__

        collections = dict((existing or {}).get('collections') or {})
        for collection, rows in sorted(catalogs.items()):
            pinned = [x['unit_id'] for x in
                      (collections.get(collection) or {}).get('probe_units', [])]
            probes = []
            for segment, row in chooseProbeUnits(rows, pinned).items():
                log.info('probing %s %s unit %s', self.name, collection,
                         row['unit_id'])
                raw = self.client.readRanges(row['url'], probeComplete)
                facts = parseProbeFacts(raw, row['unit_id'], collection)
                facts['dir_segment'] = segment
                probes.append(facts)
            collections[collection] = {'path_layouts': pathLayouts(rows),
                                       'probe_units': probes}

        return {'schema_version': CONTRACTS_VERSION,
                'harvested': datetime.now(UTC).strftime('%Y-%m-%dT%H:%M:%SZ'),
                'harvested_by': 'sourcerer %s' % __version__,
                'collections': collections}

    def harvestArtifacts(self, out, schema, catalogs):
        """
        Write the OAS specific snapshot artifacts.

        Arguments:
          out (Path): the snapshot directory being written.
          schema (SourceSchema): the freshly harvested schema.
          catalogs (dict): collection name to the catalog rows harvested this
            run.

        Returns:
          dict: artifact name to (path, changed).
        """
        from sourcerer.Contracts import (
            loadContracts,
            saveContracts,
            saveFingerprint,
        )

        # Pins come from the directory being written when it already holds
        # contracts, otherwise from the packaged snapshot, so a refresh into a
        # fresh --out directory still keeps the committed pins.
        existing = loadContracts(self.name, path=out) or loadContracts(self.name)

        written = {}
        contracts = self.harvestContracts(catalogs, existing=existing)
        written['data_contracts'] = saveContracts(contracts, out)

        if self.catalog_fingerprint is not None:
            written['catalog_fingerprint'] = saveFingerprint(
                self.catalog_fingerprint, out)

        return written

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
