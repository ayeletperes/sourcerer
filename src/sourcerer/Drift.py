"""
Schema drift detection

Compares two snapshots of a source -- the stored one and a freshly harvested
one -- and classifies every difference. A refresh normally produces many
findings at once, so severity is an ordered threshold rather than a category:
each finding carries its own level, the run's overall level is the maximum
among them, and that alone determines the exit code.

The levels, in order:

- ``additive``: something new appeared upstream (a vocabulary value, a study,
  data units, a path layout). Existing scripts keep working.
- ``anomaly``: something is internally inconsistent (a catalog key present on
  only some units, a form value no unit carries). Worth a look, breaks nothing.
- ``removed``: something users may be pinning disappeared (a value, a unit, a
  pinned probe unit). Their scripts break; sourcerer itself keeps working.
- ``structural``: the shape sourcerer parses changed (a field added or removed,
  a parse contract broken, a data unit format deviation, an unresolvable URL).
  Code changes may be needed, and the contract tests on the refresh PR say so.
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import io
import json
import logging
import subprocess
from dataclasses import dataclass
from dataclasses import field as dcField
from datetime import UTC, datetime
from pathlib import Path

import yaml

# Sourcerer imports
from sourcerer.Catalog import loadCatalog
from sourcerer.Contracts import (
    CONTRACTS_FILE,
    FINGERPRINT_FILE,
    loadContracts,
    loadFingerprint,
)
from sourcerer.Exceptions import SourcererError
from sourcerer.Schema import fromDict, loadSchema
from sourcerer.Version import __version__

log = logging.getLogger(__name__)

#: Severity levels in ascending order. The exit codes leave room between
#: levels for any future refinement without renumbering.
LEVELS = ('none', 'additive', 'anomaly', 'removed', 'structural')
LEVEL_CODES = {'none': 0, 'additive': 10, 'anomaly': 20, 'removed': 25,
               'structural': 30}

#: How many units the URL rule probe samples per collection.
PROBE_SAMPLE = 3

#: How many items a finding message enumerates before eliding the rest.
LIST_CAP = 8


@dataclass(frozen=True)
class Finding:
    """
    One classified difference between two snapshots.

    Arguments:
      level (str): one of LEVELS[1:].
      category (str): what kind of thing changed, e.g. 'field', 'unit'.
      message (str): human readable description.
      collection (str): which collection it concerns, or None for source wide.
    """
    level: str
    category: str
    message: str
    collection: str = None


@dataclass
class Snapshot:
    """
    Everything one version of a source's snapshot contains.

    Arguments:
      schema (SourceSchema): the search schema, or None if absent.
      catalogs (dict): collection name to catalog rows.
      contracts (dict): the data contracts, or None if absent.
      fingerprint (dict): the catalog fingerprint, or None if absent.
    """
    schema: object = None
    catalogs: dict = dcField(default_factory=dict)
    contracts: dict = None
    fingerprint: dict = None


def listSome(values):
    """
    Render a list for a finding message, eliding a long tail.

    Arguments:
      values (iterable): the items.

    Returns:
      str: up to LIST_CAP items, with a count of the remainder.
    """
    values = sorted(str(x) for x in values)
    shown = ', '.join(values[:LIST_CAP])
    if len(values) > LIST_CAP:
        shown += ', ... %d more' % (len(values) - LIST_CAP)

    return shown


# ---------------------------------------------------------------------------
# Loading snapshots
# ---------------------------------------------------------------------------

def packagedDir(source):
    """
    Locate the packaged snapshot directory for a source.

    Arguments:
      source (str): the source name.

    Returns:
      Path: the snapshot directory.
    """
    from importlib import resources

    return Path(str(resources.files('sourcerer').joinpath('data/schemas',
                                                          source)))


def loadSnapshotDir(source, path=None):
    """
    Load a snapshot from a directory, the packaged one by default.

    Arguments:
      source (str): the source name.
      path (Path): the directory, or None for the packaged snapshot.

    Returns:
      Snapshot: whatever the directory contains.
    """
    schema = loadSchema(source, path=path)

    catalogs = {}
    root = Path(path) if path is not None else packagedDir(source)
    for collection in schema.collection_names:
        handle = root / ('%s_catalog.tsv' % collection)
        if handle.exists():
            catalogs[collection] = loadCatalog(handle)

    return Snapshot(schema=schema, catalogs=catalogs,
                    contracts=loadContracts(source, path=path),
                    fingerprint=loadFingerprint(source, path=path))


def gitShow(root, rev, relpath):
    """
    Read one file from a git revision.

    Arguments:
      root (Path): the repository working tree.
      rev (str): the revision.
      relpath (str): repository relative path.

    Returns:
      str: the file's text, or None if it does not exist at that revision.
    """
    result = subprocess.run(
        ['git', '-C', str(root), 'show', '%s:%s' % (rev, relpath)],
        capture_output=True, text=True)

    return result.stdout if result.returncode == 0 else None


def loadSnapshotGit(source, rev):
    """
    Load a snapshot from a git revision of the working checkout.

    Requires the packaged snapshot directory to live inside a git working tree,
    which is the case for a development checkout or an editable install; it is
    how the scheduled workflow compares its fresh harvest against the last
    committed state.

    Arguments:
      source (str): the source name.
      rev (str): the git revision, e.g. 'HEAD'.

    Returns:
      Snapshot: whatever that revision contains.

    Raises:
      SourcererError: if the snapshot directory is not inside a git repository
        or the revision has no snapshot.
    """
    packaged = packagedDir(source)
    probe = subprocess.run(
        ['git', '-C', str(packaged), 'rev-parse', '--show-toplevel'],
        capture_output=True, text=True)
    if probe.returncode != 0:
        raise SourcererError(
            "cannot compare against 'git:%s': the packaged snapshot at %s is "
            'not inside a git checkout; pass a snapshot directory instead'
            % (rev, packaged))

    root = Path(probe.stdout.strip())
    prefix = packaged.resolve().relative_to(root.resolve()).as_posix()

    text = gitShow(root, rev, '%s/schema.yaml' % prefix)
    if text is None:
        raise SourcererError("no schema.yaml for '%s' at git revision '%s'"
                             % (source, rev))
    schema = fromDict(yaml.safe_load(text))

    catalogs = {}
    for collection in schema.collection_names:
        text = gitShow(root, rev, '%s/%s_catalog.tsv' % (prefix, collection))
        if text is not None:
            import csv
            catalogs[collection] = list(csv.DictReader(io.StringIO(text),
                                                       delimiter='\t'))

    contracts = gitShow(root, rev, '%s/%s' % (prefix, CONTRACTS_FILE))
    fingerprint = gitShow(root, rev, '%s/%s' % (prefix, FINGERPRINT_FILE))

    return Snapshot(
        schema=schema, catalogs=catalogs,
        contracts=yaml.safe_load(contracts) if contracts is not None else None,
        fingerprint=json.loads(fingerprint) if fingerprint is not None else None)


def loadSnapshot(source, location):
    """
    Load a snapshot from wherever a --against argument points.

    Arguments:
      source (str): the source name.
      location (str): 'git:REV' or a directory path.

    Returns:
      Snapshot: the loaded snapshot.
    """
    if location.startswith('git:'):
        return loadSnapshotGit(source, location[len('git:'):])

    return loadSnapshotDir(source, path=Path(location))


# ---------------------------------------------------------------------------
# Comparisons
# ---------------------------------------------------------------------------

def compareFields(collection_name, old, new):
    """
    Compare one collection's field list and vocabularies.

    Arguments:
      collection_name (str): which collection.
      old (Collection): the stored version.
      new (Collection): the fresh version.

    Returns:
      list: findings.
    """
    findings = []
    old_names = set(old.field_names)
    new_names = set(new.field_names)

    for name in sorted(new_names - old_names):
        findings.append(Finding(
            'structural', 'field', "field '%s' appeared" % name,
            collection_name))
    for name in sorted(old_names - new_names):
        findings.append(Finding(
            'structural', 'field', "field '%s' disappeared" % name,
            collection_name))

    for name in sorted(old_names & new_names):
        before, after = old.getField(name), new.getField(name)

        if before.wildcard != after.wildcard:
            findings.append(Finding(
                'structural', 'field',
                "field '%s' changed its wildcard from '%s' to '%s'"
                % (name, before.wildcard, after.wildcard), collection_name))
        if before.pseudo_values != after.pseudo_values:
            findings.append(Finding(
                'structural', 'field',
                "field '%s' %s a presence-only field"
                % (name, 'became' if after.pseudo_values else 'is no longer'),
                collection_name))

        added = set(after.values) - set(before.values)
        gone = set(before.values) - set(after.values)
        if added:
            findings.append(Finding(
                'additive', 'value',
                "field '%s' gained %d value(s): %s"
                % (name, len(added), listSome(added)), collection_name))
        if gone:
            findings.append(Finding(
                'removed', 'value',
                "field '%s' lost %d value(s): %s"
                % (name, len(gone), listSome(gone)), collection_name))

    return findings


def compareSchemas(old, new):
    """
    Compare two search schemas.

    Arguments:
      old (SourceSchema): the stored version.
      new (SourceSchema): the fresh version.

    Returns:
      list: findings.
    """
    findings = []

    for name, label in (('parse_contracts', 'parse contract'),
                        ('url_rules', 'URL rule'),
                        ('field_aliases', 'field alias'),
                        ('source_urls', 'source URL')):
        before, after = getattr(old, name), getattr(new, name)
        for key in sorted(set(before) | set(after)):
            if before.get(key) != after.get(key):
                findings.append(Finding(
                    'structural', 'contract',
                    "%s '%s' changed from %r to %r"
                    % (label, key, before.get(key), after.get(key))))

    old_names = set(old.collection_names)
    new_names = set(new.collection_names)
    for name in sorted(new_names - old_names):
        findings.append(Finding('structural', 'collection',
                                "collection '%s' appeared" % name))
    for name in sorted(old_names - new_names):
        findings.append(Finding('structural', 'collection',
                                "collection '%s' disappeared" % name))

    for name in sorted(old_names & new_names):
        findings += compareFields(name, old.getCollection(name),
                                  new.getCollection(name))

    return findings


def compareCatalogs(collection, old_rows, new_rows):
    """
    Compare two versions of one collection's catalog.

    Arguments:
      collection (str): which collection.
      old_rows (list): the stored catalog.
      new_rows (list): the fresh catalog.

    Returns:
      list: findings.
    """
    findings = []
    before = {x['unit_id']: x for x in old_rows}
    after = {x['unit_id']: x for x in new_rows}

    added = set(after) - set(before)
    gone = set(before) - set(after)

    def byStudy(unit_ids, source):
        groups = {}
        for unit_id in unit_ids:
            study = source[unit_id].get('study', '')
            groups.setdefault(study, []).append(unit_id)
        return sorted(groups.items())

    for study, units in byStudy(added, after):
        findings.append(Finding(
            'additive', 'unit',
            '%d new unit(s) in %s' % (len(units), study or '(no study)'),
            collection))
    for study, units in byStudy(gone, before):
        findings.append(Finding(
            'removed', 'unit',
            '%d unit(s) disappeared from %s: %s'
            % (len(units), study or '(no study)', listSome(units)), collection))

    return findings


def compareProbeUnits(collection, old_units, new_units):
    """
    Compare the pinned probe units of one collection's data contracts.

    Any fact deviation is structural: the probe facts are exactly what the
    conversion code assumes about the file format.

    Arguments:
      collection (str): which collection.
      old_units (list): probe unit entries from the stored contracts.
      new_units (list): probe unit entries from the fresh contracts.

    Returns:
      list: findings.
    """
    findings = []
    before = {x.get('dir_segment', ''): x for x in old_units}
    after = {x.get('dir_segment', ''): x for x in new_units}

    for segment in sorted(set(before) - set(after)):
        findings.append(Finding(
            'removed', 'probe-unit',
            "pinned probe unit '%s' has no counterpart: the '%s' layout "
            'disappeared' % (before[segment].get('unit_id'), segment),
            collection))

    for segment in sorted(set(before) & set(after)):
        old_unit, new_unit = before[segment], after[segment]
        if old_unit.get('unit_id') != new_unit.get('unit_id'):
            findings.append(Finding(
                'removed', 'probe-unit',
                "pinned probe unit '%s' disappeared; '%s' now probes the "
                "'%s' layout" % (old_unit.get('unit_id'),
                                 new_unit.get('unit_id'), segment), collection))
            # The two units may legitimately differ in row level facts, so the
            # fact comparison only applies to the same unit.
            continue

        for key in sorted(set(old_unit) | set(new_unit)):
            if key in ('unit_id', 'dir_segment'):
                continue
            if old_unit.get(key) != new_unit.get(key):
                findings.append(Finding(
                    'structural', 'data-contract',
                    "probe unit '%s': %s changed from %s to %s"
                    % (new_unit.get('unit_id'), key,
                       json.dumps(old_unit.get(key), sort_keys=True),
                       json.dumps(new_unit.get(key), sort_keys=True)),
                    collection))

    return findings


def compareContracts(old, new):
    """
    Compare two data contracts documents.

    Arguments:
      old (dict): the stored contracts.
      new (dict): the fresh contracts.

    Returns:
      list: findings.
    """
    findings = []
    before = (old or {}).get('collections') or {}
    after = (new or {}).get('collections') or {}

    for collection in sorted(set(before) - set(after)):
        findings.append(Finding(
            'structural', 'data-contract',
            'the data contracts no longer cover this collection', collection))

    for collection in sorted(set(before) & set(after)):
        old_layouts = before[collection].get('path_layouts') or {}
        new_layouts = after[collection].get('path_layouts') or {}
        for key, label in (('observed_dir_segments', 'directory layout'),
                           ('observed_filename_patterns', 'filename pattern')):
            olds = set(old_layouts.get(key) or {})
            news = set(new_layouts.get(key) or {})
            # Paths are opaque, so a novel layout costs nothing: additive.
            for item in sorted(news - olds):
                findings.append(Finding(
                    'additive', 'layout',
                    "new %s '%s'" % (label, item), collection))
            for item in sorted(olds - news):
                findings.append(Finding(
                    'removed', 'layout',
                    "%s '%s' disappeared" % (label, item), collection))

        findings += compareProbeUnits(
            collection, before[collection].get('probe_units') or [],
            after[collection].get('probe_units') or [])

    return findings


def compareFingerprints(old, new):
    """
    Compare two catalog fingerprints.

    Arguments:
      old (dict): the stored fingerprint.
      new (dict): the fresh fingerprint.

    Returns:
      list: findings.
    """
    findings = []

    # The fingerprint names its own collection (buildFingerprint stamps it
    # from the catalog rows actually fingerprinted) rather than this module
    # assuming one; the fresh side wins when both are known, since it
    # describes the catalog the findings below are actually about.
    collection = (new or {}).get('collection') or (old or {}).get('collection')

    before_keys = set((old or {}).get('key_counts') or {})
    after_keys = set((new or {}).get('key_counts') or {})
    for key in sorted(after_keys - before_keys):
        count = new['key_counts'][key]
        # A key on every unit is a new column of the catalog; a key on a few is
        # more likely upstream inconsistency, and the anomaly check reports it.
        level = 'structural' if count == new.get('n_units') else 'additive'
        findings.append(Finding(
            level, 'catalog-key',
            "catalog key '%s' appeared on %d of %d units"
            % (key, count, new.get('n_units', 0)), collection))
    for key in sorted(before_keys - after_keys):
        findings.append(Finding(
            'structural', 'catalog-key',
            "catalog key '%s' disappeared" % key, collection))

    before_types = (old or {}).get('value_types') or {}
    after_types = (new or {}).get('value_types') or {}
    for key in sorted(before_types.keys() & after_types.keys()):
        if before_types[key] != after_types[key]:
            findings.append(Finding(
                'structural', 'value-type',
                "catalog key '%s' changed type from %s to %s"
                % (key, before_types[key], after_types[key]), collection))

    old_units = (old or {}).get('n_units')
    new_units = (new or {}).get('n_units')
    if old_units is not None and new_units is not None:
        if new_units > old_units:
            findings.append(Finding(
                'additive', 'catalog',
                'catalog grew from %d to %d units' % (old_units, new_units),
                collection))
        elif new_units < old_units:
            findings.append(Finding(
                'removed', 'catalog',
                'catalog shrank from %d to %d units' % (old_units, new_units),
                collection))

    return findings


# ---------------------------------------------------------------------------
# Checks on the new snapshot alone
# ---------------------------------------------------------------------------

def findAnomalies(snapshot):
    """
    Report internal inconsistencies of a snapshot.

    Arguments:
      snapshot (Snapshot): the snapshot to inspect.

    Returns:
      list: findings.
    """
    findings = []

    fingerprint = snapshot.fingerprint or {}
    collection = fingerprint.get('collection')
    n_units = fingerprint.get('n_units', 0)
    for key, count in sorted((fingerprint.get('key_counts') or {}).items()):
        if 0 < count < n_units:
            findings.append(Finding(
                'anomaly', 'partial-key',
                "catalog key '%s' is present on %d of %d units"
                % (key, count, n_units), collection))

    if snapshot.schema is None:
        return findings

    for collection_name in snapshot.schema.collection_names:
        rows = snapshot.catalogs.get(collection_name)
        if not rows:
            continue
        columns = rows[0].keys()
        for item in snapshot.schema.getCollection(collection_name).fields:
            if item.pseudo_values or item.name not in columns:
                continue
            observed = {x.get(item.name, '') for x in rows}
            unseen = set(item.values) - observed
            if unseen:
                findings.append(Finding(
                    'anomaly', 'unseen-value',
                    "form value(s) for '%s' matched by no cataloged unit: %s"
                    % (item.name, listSome(unseen)), collection_name))

    return findings


def checkPathSafety(snapshot):
    """
    Assert the only two path facts the code depends on.

    Paths are opaque, so a novel layout is fine -- but a unit whose relative
    path cannot be preserved locally (absolute, escaping upward, or colliding
    with another unit) breaks the mirror, and that is structural.

    Arguments:
      snapshot (Snapshot): the snapshot to inspect.

    Returns:
      list: findings.
    """
    from pathlib import PurePosixPath

    findings = []
    for collection, rows in sorted(snapshot.catalogs.items()):
        seen = {}
        for row in rows:
            unit_id = row.get('unit_id', '')
            parts = PurePosixPath(unit_id).parts
            if PurePosixPath(unit_id).is_absolute() or '..' in parts:
                findings.append(Finding(
                    'structural', 'path',
                    "unit_id '%s' cannot be preserved as a local relative path"
                    % unit_id, collection))
            normalized = str(PurePosixPath(unit_id))
            if normalized in seen and seen[normalized] != unit_id:
                findings.append(Finding(
                    'structural', 'path',
                    "unit_ids '%s' and '%s' collide after normalization"
                    % (seen[normalized], unit_id), collection))
            seen.setdefault(normalized, unit_id)

    return findings


def probeUrls(snapshot, client, sample=PROBE_SAMPLE):
    """
    Verify the URL rewrite rule against the live server.

    No vocabulary diff can catch a change to the catalog-key-to-URL mapping, so
    a few units per collection are probed with a HEAD, falling back to a one
    byte ranged GET. Nothing whole is downloaded.

    Arguments:
      snapshot (Snapshot): the snapshot whose catalogs supply the URLs.
      client (HttpClient): the client to probe with.
      sample (int): units probed per collection.

    Returns:
      list: findings.
    """
    findings = []
    for collection, rows in sorted(snapshot.catalogs.items()):
        chosen = sorted(rows, key=lambda x: x.get('unit_id', ''))[:sample]
        for row in chosen:
            if not client.probeAlive(row['url']):
                findings.append(Finding(
                    'structural', 'url-rule',
                    "the URL derived from unit_id '%s' does not resolve: %s"
                    % (row.get('unit_id'), row.get('url')), collection))

    return findings


# ---------------------------------------------------------------------------
# Orchestration and reporting
# ---------------------------------------------------------------------------

def checkDrift(old, new, client=None):
    """
    Compare two snapshots and classify every difference.

    Arguments:
      old (Snapshot): the stored snapshot.
      new (Snapshot): the fresh snapshot.
      client (HttpClient): probe the URL rule live when given; skipped when
        None, which is what offline runs pass.

    Returns:
      list: findings, in a stable order.
    """
    findings = []
    findings += compareSchemas(old.schema, new.schema)

    for collection in sorted(set(old.catalogs) | set(new.catalogs)):
        findings += compareCatalogs(collection,
                                    old.catalogs.get(collection, []),
                                    new.catalogs.get(collection, []))

    if old.contracts is not None and new.contracts is not None:
        findings += compareContracts(old.contracts, new.contracts)
    if old.fingerprint is not None and new.fingerprint is not None:
        findings += compareFingerprints(old.fingerprint, new.fingerprint)

    findings += findAnomalies(new)
    findings += checkPathSafety(new)

    if client is not None:
        findings += probeUrls(new, client)

    return findings


def overallLevel(findings):
    """
    Reduce findings to the run's overall level.

    Arguments:
      findings (list): the findings.

    Returns:
      str: the maximum severity present, or 'none'.
    """
    worst = 'none'
    for finding in findings:
        if LEVELS.index(finding.level) > LEVELS.index(worst):
            worst = finding.level

    return worst


def exitCode(findings, fail_on):
    """
    Map findings and a threshold to a process exit status.

    Arguments:
      findings (list): the findings.
      fail_on (str): 'never' or a level name; the run fails when its overall
        level is at or above this.

    Returns:
      int: 0, or the overall level's code.
    """
    level = overallLevel(findings)
    if fail_on == 'never' or level == 'none':
        return 0
    if LEVELS.index(level) >= LEVELS.index(fail_on):
        return LEVEL_CODES[level]

    return 0


def buildReport(source, against, findings):
    """
    Assemble the machine readable drift report.

    Arguments:
      source (str): the source checked.
      against (str): what the snapshot was compared to.
      findings (list): the findings.

    Returns:
      dict: the report.
    """
    return {
        'source': source,
        'against': against,
        'checked': datetime.now(UTC).strftime('%Y-%m-%dT%H:%M:%SZ'),
        'checked_by': 'sourcerer %s' % __version__,
        'overall_level': overallLevel(findings),
        'findings': [{'level': x.level, 'category': x.category,
                      'collection': x.collection, 'message': x.message}
                     for x in findings],
    }


def renderMarkdown(report):
    """
    Render a drift report as markdown, for the job summary and the PR body.

    Arguments:
      report (dict): the report from buildReport.

    Returns:
      str: the markdown.
    """
    lines = ['# Schema drift: %s' % report['source'], '']
    lines.append('Compared against `%s` at %s.'
                 % (report['against'], report['checked']))
    lines.append('')

    findings = report['findings']
    if not findings:
        lines.append('No drift detected.')
        return '\n'.join(lines) + '\n'

    lines.append('**Overall level: %s** (%d finding(s))'
                 % (report['overall_level'], len(findings)))

    for level in reversed(LEVELS[1:]):
        subset = [x for x in findings if x['level'] == level]
        if not subset:
            continue
        lines += ['', '## %s (%d)' % (level, len(subset)), '']
        for finding in subset:
            where = (' `[%s]`' % finding['collection']
                     if finding['collection'] else '')
            lines.append('- **%s**%s %s'
                         % (finding['category'], where, finding['message']))

    return '\n'.join(lines) + '\n'
