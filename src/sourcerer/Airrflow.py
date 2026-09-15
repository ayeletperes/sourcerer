"""
nf-core/airrflow samplesheets

The column set and the hygiene rules follow the mapping already in use in
Rmd/airrflow.Rmd, with two corrections noted at the point they are applied.

A samplesheet is a derived artifact of a data format, not a format in its own
right: its filename column has to name exactly one file per sample. Requesting
both AIRR and FASTA therefore produces two samplesheets, each internally
consistent, rather than one that has to choose.

A samplesheet describes every unit converted into its output directory, not just
the units of the run that happened to write it last. Downloading in several
passes is the normal way to assemble a dataset, so writing it fresh each time
would silently orphan the files of every earlier pass. Rewriting is therefore a
merge: see mergeSamplesheet for what that guarantees.
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import csv
import logging
from pathlib import Path

# Sourcerer imports
from sourcerer.Exceptions import SourcererError
from sourcerer.Sources.Oas import isNull

log = logging.getLogger(__name__)

#: Columns airrflow expects, in order. sample_name is added by sourcerer so that
#: the original identifier survives the rewrite of sample_id.
SAMPLESHEET_COLUMNS = ('sample_id', 'filename', 'subject_id', 'species',
                       'pcr_target_locus', 'tissue', 'sex', 'age',
                       'biomaterial_provider', 'single_cell', 'intervention',
                       'disease_diagnosis', 'cell_subset', 'study',
                       'sample_name')

#: Loci that map to each airrflow pcr_target_locus value.
IG_LOCI = frozenset(['IGH', 'IGK', 'IGL'])
TR_LOCI = frozenset(['TRA', 'TRB', 'TRG', 'TRD'])

#: Values this module writes when it has nothing to say. They are stand ins for a
#: missing value rather than data, so a merge treats them as absent and lets a
#: real value replace them. Without this a field with a non-empty default, such
#: as tissue, could never be improved by a later run.
PLACEHOLDER_VALUES = frozenset(['', 'unknown', 'NA'])


def isPlaceholder(value):
    """
    Report whether a samplesheet value carries no information.

    Arguments:
      value (str): the value to test.

    Returns:
      bool: True if the value is missing or a stand in for missing.
    """
    return value is None or str(value).strip() in PLACEHOLDER_VALUES


def targetLocus(loci):
    """
    Collapse observed loci to the receptor class airrflow asks for.

    This collapse belongs here and nowhere else. Doing it to the rearrangement
    file's own locus column, as the R implementation does, produces a file that
    is not valid AIRR, because locus must name the actual gene locus.

    Arguments:
      loci (iterable): observed locus values.

    Returns:
      str: 'IG', 'TR', or ''.

    Raises:
      ValueError: if a single unit mixes immunoglobulin and T cell receptor loci.
    """
    found = {x for x in loci if x}
    ig, tr = found & IG_LOCI, found & TR_LOCI

    if ig and tr:
        raise ValueError('data unit mixes IG (%s) and TR (%s) loci'
                         % (sorted(ig), sorted(tr)))
    if ig:
        return 'IG'
    if tr:
        return 'TR'

    return ''


def clean(value, default=''):
    """
    Normalize a metadata value, mapping the source's null sentinels to a default.

    Arguments:
      value: the raw value.
      default (str): what to use when the value carries no information.

    Returns:
      str: the cleaned value.
    """
    if isNull(value):
        return default

    return str(value).strip()


def loadSamplesheet(path):
    """
    Read a samplesheet sourcerer wrote earlier.

    A file whose header is not the one written here is not something this code
    may rewrite, so it raises rather than merging into it. Overwriting a
    hand-built samplesheet because it happened to occupy the expected filename
    would destroy work that cannot be regenerated.

    Arguments:
      path (Path): the samplesheet to read.

    Returns:
      list: rows in file order, empty if the file does not exist.

    Raises:
      SourcererError: if the file exists but was not written by sourcerer.
    """
    path = Path(path)
    if not path.exists():
        return []

    with open(path, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = reader.fieldnames or []
        if list(fields) != list(SAMPLESHEET_COLUMNS):
            raise SourcererError(
                '%s does not look like a sourcerer samplesheet (expected columns '
                '%s), refusing to overwrite it'
                % (path, ', '.join(SAMPLESHEET_COLUMNS)))

        return [dict(row) for row in reader]


def mergeSamplesheet(existing, fresh):
    """
    Merge freshly converted units into the rows already in a samplesheet.

    Three guarantees, all aimed at making a repeated download additive rather
    than destructive:

    - A unit already described keeps its sample_id, so identifiers a user has
      already referenced downstream stay valid when the sheet grows.
    - New units are appended in resolved order with the next free sample_id,
      so the numbering never renumbers what came before.
    - For a unit seen again, a value already in the file wins over the newly
      derived one wherever it carries information. Re-running fills in blanks and
      replaces placeholders, for instance once detail page enrichment supplies a
      tissue, but never overwrites a field edited by hand. sex is not derivable
      from the source at all, so this is the only thing protecting it.

    Arguments:
      existing (list): rows already in the file, in order.
      fresh (list): rows for this run's units, without a sample_id.

    Returns:
      list: the merged rows, in file order.
    """
    merged = list(existing)
    seen = {row.get('sample_name'): index for index, row in enumerate(merged)}
    used = {row.get('sample_id') for row in merged}
    counter = 0

    for row in fresh:
        index = seen.get(row['sample_name'])
        if index is not None:
            kept = merged[index]
            merged[index] = {
                column: (row.get(column, '') if isPlaceholder(kept.get(column))
                         else kept[column])
                for column in SAMPLESHEET_COLUMNS}
            merged[index]['sample_id'] = kept.get('sample_id', '')
            continue

        # Skip identifiers already in the file rather than assuming the existing
        # rows are a contiguous ssr_1..ssr_n run: they may have been edited, and
        # a duplicate sample_id would make airrflow merge two samples.
        while True:
            counter += 1
            candidate = 'ssr_%d' % counter
            if candidate not in used:
                break

        used.add(candidate)
        row['sample_id'] = candidate
        seen[row['sample_name']] = len(merged)
        merged.append(row)

    for row in merged:
        # Deferred until now because the fallback names the sample_id, which is
        # only known once the merge has assigned it.
        if not row.get('subject_id'):
            row['subject_id'] = '%s_subj' % row['sample_id']

    return merged


def buildSamplesheet(entries, out, collection, root=None, loci=None):
    """
    Write an airrflow samplesheet describing converted data units.

    One row per data unit, since one unit is one repertoire is one sample.

    An existing samplesheet at `out` is merged into rather than replaced, so
    building a dataset over several downloads accumulates. See mergeSamplesheet.

    Arguments:
      entries (list): (DataUnit, Path) pairs naming the converted output.
      out (Path): where to write the samplesheet.
      collection (str): the collection the units came from.
      root (Path): if given, filenames are written relative to it.
      loci (dict): unit_id to the loci observed in its converted output, used to
        derive pcr_target_locus.

    Returns:
      Path: the file written.
    """
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    loci = loci or {}

    rows = []
    for unit, path in entries:
        metadata = unit.metadata or {}
        filename = Path(path)
        if root is not None:
            try:
                filename = filename.relative_to(Path(root))
            except ValueError:
                pass

        # sample_id is left for the merge to assign, since it depends on what the
        # samplesheet already contains. The real identifier is preserved in
        # sample_name, which is what the merge keys on.
        subject = clean(metadata.get('Subject')) or clean(metadata.get('study'))

        rows.append({
            'sample_id': '',
            'filename': str(filename),
            'subject_id': subject.replace(' ', '_'),
            'species': clean(metadata.get('Species'), 'human').lower(),
            'pcr_target_locus': targetLocus(loci.get(unit.unit_id, [])),
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
            'single_cell': 'TRUE' if collection == 'paired' else 'FALSE',
            'disease_diagnosis': clean(metadata.get('Disease')),
            'intervention': clean(metadata.get('Vaccine')),
            'cell_subset': clean(metadata.get('BType')),
            'study': clean(metadata.get('study')) or unit.study,
            'sample_name': unit.unit_id,
        })

    existing = loadSamplesheet(out)
    merged = mergeSamplesheet(existing, rows)
    if existing:
        log.info('%s: %d units already described, %d now in total',
                 out.name, len(existing), len(merged))

    with open(out, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(SAMPLESHEET_COLUMNS),
                                delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(merged)

    return out
