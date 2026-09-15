"""
Download provenance

What was fetched, from where, when, and what it hashed to. The raw mirror is the
one artifact a user cannot regenerate from anything else in the output directory,
so it is the one that most needs a record of its own origin.

This is deliberately not a samplesheet. A samplesheet is an airrflow input and
names one file per sample; raw source files are not valid airrflow input in
either mode, so describing them in that shape would invite feeding them to a
pipeline that cannot read them. This file answers a different question: what is
in this directory and where did it come from.
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import logging
import sys
from datetime import UTC, datetime
from pathlib import Path

import yaml

# Sourcerer imports
from sourcerer.Exceptions import SourcererError
from sourcerer.Version import __version__

log = logging.getLogger(__name__)

#: Name of the record written at the root of a download directory. It contains
#: neither 'samplesheet' nor 'airrflow' because it is neither: nothing consumes
#: it as a pipeline input.
DOWNLOAD_METADATA = 'download_metadata.yml'

#: Format version of this file, so a later reader can tell what it is looking at.
METADATA_VERSION = 1


def timestamp():
    """
    Returns:
      str: the current time as an ISO 8601 UTC string.
    """
    return datetime.now(UTC).strftime('%Y-%m-%dT%H:%M:%SZ')


def relativize(path, root):
    """
    Express a path relative to the output root where possible.

    Recorded paths are relative so that moving or renaming the download
    directory does not invalidate the record.

    Arguments:
      path (Path): the path to express.
      root (Path): the download root.

    Returns:
      str: the relative path, or the original if it lies outside root.
    """
    try:
        return str(Path(path).relative_to(Path(root)))
    except ValueError:
        return str(path)


def loadMetadata(path):
    """
    Read a metadata file written by an earlier run.

    A file that is not one of ours is not something this code may rewrite, so it
    raises rather than merging into it.

    Arguments:
      path (Path): the file to read.

    Returns:
      dict: the parsed record, or a fresh skeleton if the file does not exist.

    Raises:
      SourcererError: if the file exists but was not written by sourcerer.
    """
    path = Path(path)
    if not path.exists():
        return {'sourcerer_metadata_version': METADATA_VERSION,
                'runs': [], 'units': []}

    with open(path) as handle:
        record = yaml.safe_load(handle)

    if not isinstance(record, dict) or 'sourcerer_metadata_version' not in record:
        raise SourcererError('%s was not written by sourcerer, refusing to '
                             'overwrite it' % path)

    record.setdefault('runs', [])
    record.setdefault('units', [])

    return record


def mergeUnits(existing, fresh):
    """
    Merge this run's units into those already recorded.

    Keyed on unit_id, because assembling a dataset over several downloads is
    normal and rewriting the file from one run's units alone would drop every
    earlier one. A unit fetched again replaces its entry, since the newer digest
    and timestamp describe what is actually on disk now, but its recorded outputs
    accumulate: converting to FASTA today does not unrecord the AIRR file written
    yesterday, which is still there.

    Arguments:
      existing (list): unit records already in the file.
      fresh (list): unit records from this run.

    Returns:
      list: the merged records, earlier units keeping their position.
    """
    merged = list(existing)
    seen = {x.get('unit_id'): i for i, x in enumerate(merged)}

    for unit in fresh:
        index = seen.get(unit['unit_id'])
        if index is None:
            seen[unit['unit_id']] = len(merged)
            merged.append(unit)
            continue

        outputs = dict(merged[index].get('outputs') or {})
        outputs.update(unit.get('outputs') or {})
        merged[index] = {**unit, 'outputs': outputs}

    return merged


def buildUnitRecord(unit, result, root, outputs=None):
    """
    Describe one downloaded unit.

    Arguments:
      unit (DataUnit): what was fetched.
      result (DownloadResult): the outcome of fetching it.
      root (Path): the download root, for relative paths.
      outputs (dict): format name to written path, for converted formats.

    Returns:
      dict: the record.
    """
    return {
        'unit_id': unit.unit_id,
        'collection': unit.collection,
        'url': unit.url,
        'raw': relativize(result.path, root),
        'sha256': result.sha256,
        'size_bytes': result.size_bytes,
        'n_sequences': unit.n_sequences,
        'downloaded': timestamp(),
        'outputs': {k: relativize(v, root) for k, v in (outputs or {}).items()},
    }


def commandLine():
    """
    Reconstruct the invocation under the tool's own name.

    argv[0] is whatever launched the process, which may be an absolute path to
    a module inside the installed package. Substituting the entry point name
    keeps the recorded command both readable and runnable.

    Returns:
      str: the command line.
    """
    return ' '.join(['sourcerer'] + sys.argv[1:])


def writeDownloadMetadata(out, source, collection, filters, limit, formats,
                          units, schema=None, license=None, citation=None):
    """
    Write or update the provenance record for a download directory.

    Arguments:
      out (Path): the download root. The file is written at its top level.
      source (str): the source name.
      collection (str): the collection downloaded.
      filters (dict): the filters the user supplied, so the download can be
        repeated. Not the resolved query, whose unconstrained fields are filled
        in with a source specific wildcard that carries no information here.
      limit (int): the unit cap, or None.
      formats (list): formats written this run.
      units (list): unit records from buildUnitRecord.
      schema (SourceSchema): the snapshot the query resolved against, if known.
      license (str): the source's data license, if known. Recorded so a reader
        of this directory alone, with no access to sourcerer's own docs, still
        knows the terms the data was obtained under.
      citation (tuple): the source's requested citation(s), if known, for the
        same reason.

    Returns:
      Path: the file written.
    """
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / DOWNLOAD_METADATA

    record = loadMetadata(path)

    run = {
        'finished': timestamp(),
        # The invocation, so a reader can repeat or amend the download without
        # reconstructing the filters from the unit list.
        'command': commandLine(),
        'collection': collection,
        'filters': dict(filters or {}),
        'limit': limit,
        'formats': list(formats),
        'units': len(units),
    }
    if schema is not None:
        run['schema_harvested'] = schema.harvested
        run['schema_harvested_by'] = schema.harvested_by

    # Rebuilt in a fixed order rather than updated in place, so the header keys
    # stay at the top of the file however the loaded record was ordered.
    merged = {
        'sourcerer_metadata_version': METADATA_VERSION,
        'source': source,
        'generated_by': 'sourcerer %s' % __version__,
    }
    if license is not None:
        merged['data_license'] = license
    if citation:
        merged['data_citation'] = list(citation)
    merged['runs'] = list(record.get('runs') or []) + [run]
    merged['units'] = mergeUnits(record.get('units') or [], units)

    with open(path, 'w') as handle:
        yaml.safe_dump(merged, handle, sort_keys=False, default_flow_style=False)

    return path
