"""
Downloaded-file contracts and catalog fingerprints

The schema snapshot records what a source can be searched on; it says nothing
about the inside of a downloaded file, which is the riskiest contract in the
package: the wide to long pivot, the gzip member framing and every AIRR field
mapping depend on it, and the source could restructure a file without touching a
single form field. The data contracts file pins that format, per path layout,
from a few pinned probe units. The catalog fingerprint condenses a source's raw
catalog into the facts the drift check compares, including per key unit counts,
which is what surfaces a key present on only some units.

Everything here is written deterministically and only when its content actually
changed, because a file that changes on every run would make the scheduled
refresh open a pull request every month whether or not anything drifted.
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import json
from importlib import resources
from pathlib import Path

import yaml

# Sourcerer imports
from sourcerer.Exceptions import SchemaError

#: File names within a source's snapshot directory.
CONTRACTS_FILE = 'data_contracts.yaml'
FINGERPRINT_FILE = 'catalog_fingerprint.json'
PROVENANCE_FILE = 'provenance.json'

#: Contract format version this code understands.
CONTRACTS_VERSION = 1

#: Keys that record when a harvest ran rather than what it found, excluded when
#: deciding whether a file actually changed. ETag and Last-Modified are server
#: bookkeeping that can flip on a redeploy with identical content; treating that
#: as a change would open a pull request that reviews nothing.
VOLATILE_KEYS = ('harvested', 'harvested_by', 'etag', 'last_modified')


def _stable(payload):
    """Drop the volatile keys, for change comparison."""
    return {k: v for k, v in (payload or {}).items() if k not in VOLATILE_KEYS}


def writeIfChanged(path, text, previous_stable=None, stable=None):
    """
    Write a file only when its content meaningfully changed.

    Arguments:
      path (Path): where to write.
      text (str): the serialized content.
      previous_stable: the existing file's content with volatile keys removed,
        or None when the file does not exist.
      stable: the new content with volatile keys removed.

    Returns:
      bool: True if the file was written.
    """
    path = Path(path)
    if path.exists() and previous_stable == stable:
        return False

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)

    return True


def serializeYaml(payload):
    """
    Render a payload as deterministic YAML.

    Arguments:
      payload (dict): plain data.

    Returns:
      str: the serialized payload.
    """
    return yaml.safe_dump(payload, sort_keys=True, default_flow_style=False,
                          width=88)


def serializeJson(payload):
    """
    Render a payload as deterministic JSON.

    Arguments:
      payload (dict): plain data.

    Returns:
      str: the serialized payload.
    """
    return json.dumps(payload, sort_keys=True, indent=2) + '\n'


def _read(source, name, path=None):
    """
    Read one snapshot file, from the packaged data by default.

    Arguments:
      source (str): the source name.
      name (str): the file name within the snapshot directory.
      path (Path): a directory to read instead of the packaged snapshot.

    Returns:
      str: the file's text, or None if it does not exist.
    """
    if path is not None:
        handle = Path(path) / name
        return handle.read_text() if handle.exists() else None

    anchor = resources.files('sourcerer').joinpath('data/schemas', source, name)

    return anchor.read_text() if anchor.is_file() else None


def loadContracts(source, path=None):
    """
    Load a source's data contracts.

    Arguments:
      source (str): the source name.
      path (Path): a directory to read instead of the packaged snapshot.

    Returns:
      dict: the contracts, or None if none are stored.

    Raises:
      SchemaError: if the stored contracts are malformed or too new.
    """
    text = _read(source, CONTRACTS_FILE, path=path)
    if text is None:
        return None

    payload = yaml.safe_load(text)
    if not isinstance(payload, dict):
        raise SchemaError('%s for %s is not a mapping' % (CONTRACTS_FILE, source))

    version = payload.get('schema_version')
    if version is None or version > CONTRACTS_VERSION:
        raise SchemaError(
            '%s for %s declares schema_version %s but this sourcerer understands '
            'at most %s' % (CONTRACTS_FILE, source, version, CONTRACTS_VERSION))

    return payload


def saveContracts(payload, path):
    """
    Write data contracts, but only when their content actually changed.

    Arguments:
      payload (dict): the contracts.
      path (Path): the snapshot directory.

    Returns:
      tuple: (Path, changed) where changed is False if the file was left alone.
    """
    handle = Path(path) / CONTRACTS_FILE
    previous = None
    if handle.exists():
        previous = _stable(yaml.safe_load(handle.read_text()))

    changed = writeIfChanged(handle, serializeYaml(payload),
                             previous_stable=previous, stable=_stable(payload))

    return handle, changed


def loadFingerprint(source, path=None):
    """
    Load a source's catalog fingerprint.

    Arguments:
      source (str): the source name.
      path (Path): a directory to read instead of the packaged snapshot.

    Returns:
      dict: the fingerprint, or None if none is stored.
    """
    text = _read(source, FINGERPRINT_FILE, path=path)

    return json.loads(text) if text is not None else None


def saveFingerprint(payload, path):
    """
    Write a catalog fingerprint, but only when its content actually changed.

    Arguments:
      payload (dict): the fingerprint.
      path (Path): the snapshot directory.

    Returns:
      tuple: (Path, changed) where changed is False if the file was left alone.
    """
    handle = Path(path) / FINGERPRINT_FILE
    previous = None
    if handle.exists():
        previous = _stable(json.loads(handle.read_text()))

    changed = writeIfChanged(handle, serializeJson(payload),
                             previous_stable=previous, stable=_stable(payload))

    return handle, changed


def saveProvenance(path, timestamp, version):
    """
    Record that the snapshot changed.

    Deliberately carries no "last checked" stamp: a tracked file touched on every
    run would dirty the working tree monthly and defeat the no-drift, no-pull-
    request policy. The record of a quiet check lives in the workflow's job
    summary instead. This is therefore only called after a harvest that changed
    something.

    Arguments:
      path (Path): the snapshot directory.
      timestamp (str): ISO 8601 UTC time of the change.
      version (str): the sourcerer version that produced it.

    Returns:
      Path: the file written.
    """
    handle = Path(path) / PROVENANCE_FILE
    handle.parent.mkdir(parents=True, exist_ok=True)
    handle.write_text(serializeJson({'last_changed': timestamp,
                                     'sourcerer_version': version}))

    return handle
