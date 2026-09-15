"""
Stored schema snapshots

A snapshot is the checked-in record of what a remote source looked like at a
point in time: which fields it can be searched on, what values those fields
accept, and the structural invariants the parsing code relies on.

Everything user facing is generated from this, so no field list is hardcoded
anywhere in the package. That is deliberate: the tool this replaces carried a
literal list of paired search fields that had silently gone out of date, and
nothing detected it.
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import difflib
from dataclasses import dataclass
from dataclasses import field as dcField
from importlib import resources
from pathlib import Path

import yaml

# Sourcerer imports
from sourcerer.Exceptions import SchemaError

#: Snapshot format version this code understands. A snapshot declaring a higher
#: version is refused rather than misread.
SCHEMA_VERSION = 1

#: Values that mean "this field is recorded but has no vocabulary", used by OAS
#: paired search for Age, Subject and Longitudinal.
PSEUDO_VALUES = frozenset(['defined', 'undefined'])


@dataclass(frozen=True)
class Field:
    """
    One searchable field and its controlled vocabulary.

    Arguments:
      name (str): the field name as the remote source spells it.
      values (tuple): allowed values, excluding the wildcard.
      wildcard (str): the value meaning "all".
      pseudo_values (bool): True when values are presence flags rather than a
        real vocabulary, in which case offering them as choices would mislead.
    """
    name: str
    values: tuple = ()
    wildcard: str = '*'
    pseudo_values: bool = False

    @property
    def flag(self):
        """str: the commandline flag generated for this field."""
        return '--%s' % self.name.lower().replace(' ', '-').replace('_', '-')

    def accepts(self, value):
        """
        Test whether a value is valid for this field.

        Arguments:
          value (str): the candidate value.

        Returns:
          bool: True if the value may be sent.
        """
        if value == self.wildcard:
            return True

        # A presence-only field has no vocabulary of its own; what it accepts is
        # the presence tokens, which is also what the rejection message advises.
        if self.pseudo_values:
            return value in PSEUDO_VALUES

        return value in self.values


@dataclass(frozen=True)
class Collection:
    """
    One searchable collection within a source, such as OAS paired or unpaired.

    Arguments:
      name (str): the collection name.
      fields (tuple): Field objects in the order the remote form presents them.
      reported_totals (dict): counts the source reported at harvest time.
    """
    name: str
    fields: tuple = ()
    reported_totals: dict = dcField(default_factory=dict)

    def getField(self, name):
        """
        Look up a field by name, case insensitively.

        Arguments:
          name (str): the field name.

        Returns:
          Field: the matching field, or None.
        """
        for item in self.fields:
            if item.name.lower() == name.lower():
                return item

        return None

    @property
    def field_names(self):
        """tuple: the field names in document order."""
        return tuple(x.name for x in self.fields)


@dataclass(frozen=True)
class SourceSchema:
    """
    A complete snapshot of one remote source.

    Arguments:
      source (str): the source name, e.g. 'oas'.
      schema_version (int): snapshot format version.
      harvested (str): ISO 8601 UTC timestamp of the harvest.
      harvested_by (str): the tool version that produced it.
      source_urls (dict): named endpoints.
      url_rules (dict): how catalog keys map to download URLs.
      parse_contracts (dict): structural invariants asserted by the drift check.
      field_aliases (dict): remote synonyms mapped to canonical field names.
      collections (dict): name to Collection.
    """
    source: str
    schema_version: int = SCHEMA_VERSION
    harvested: str = None
    harvested_by: str = None
    source_urls: dict = dcField(default_factory=dict)
    url_rules: dict = dcField(default_factory=dict)
    parse_contracts: dict = dcField(default_factory=dict)
    field_aliases: dict = dcField(default_factory=dict)
    collections: dict = dcField(default_factory=dict)

    @property
    def collection_names(self):
        """tuple: collection names in sorted order."""
        return tuple(sorted(self.collections))

    def getCollection(self, name):
        """
        Look up a collection by name.

        Arguments:
          name (str): the collection name.

        Returns:
          Collection: the matching collection.

        Raises:
          SchemaError: if the collection is not in the snapshot.
        """
        if name not in self.collections:
            raise SchemaError("unknown collection '%s' for source '%s'; known: %s"
                              % (name, self.source,
                                 ', '.join(self.collection_names)))

        return self.collections[name]

    def canonicalField(self, name):
        """
        Resolve a remote synonym to the canonical field name.

        The same concept is spelled differently in different parts of a source.
        OAS search results say Organism and Individual where the form and the
        data unit metadata say Species and Subject.

        Arguments:
          name (str): a field name as it appears upstream.

        Returns:
          str: the canonical name, or the input unchanged if it is not an alias.
        """
        return self.field_aliases.get(name, name)

    def validateFilters(self, collection_name, filters):
        """
        Check filters against the snapshot and fill in wildcards.

        Unknown fields and unknown values are errors, not warnings. Sending an
        unrecognized value upstream returns zero results rather than an error,
        which reads as "no data matched" instead of "you asked for something that
        does not exist" -- the failure mode this validation exists to prevent.

        Arguments:
          collection_name (str): which collection is being searched.
          filters (dict): user supplied field to value pairs.

        Returns:
          dict: every field in the collection, defaulted to its wildcard.

        Raises:
          SchemaError: if a field or value is not in the snapshot.
        """
        collection = self.getCollection(collection_name)
        resolved = {x.name: x.wildcard for x in collection.fields}

        for name, value in filters.items():
            target = collection.getField(self.canonicalField(name))
            if target is None:
                raise SchemaError(
                    "unknown filter field '%s' for %s %s; available fields: %s"
                    % (name, self.source, collection_name,
                       ', '.join(collection.field_names)))

            if not target.accepts(value):
                raise SchemaError(self._badValueMessage(collection_name, target,
                                                        value))

            resolved[target.name] = value

        return resolved

    def _badValueMessage(self, collection_name, target, value):
        """
        Build an actionable message for a rejected filter value.

        Arguments:
          collection_name (str): the collection being searched.
          target (Field): the field the value was offered for.
          value (str): the rejected value.

        Returns:
          str: the error message, including close matches where any exist.
        """
        if target.pseudo_values:
            return ("'%s' is not valid for %s in %s; this field only filters on "
                    "presence, so use one of: %s, %s"
                    % (value, target.name, collection_name, target.wildcard,
                       ', '.join(sorted(PSEUDO_VALUES))))

        message = ("'%s' is not an available value for %s in %s (%d values known)"
                   % (value, target.name, collection_name, len(target.values)))

        close = difflib.get_close_matches(value, target.values, n=3, cutoff=0.6)
        if close:
            message += '; did you mean: %s' % ', '.join(close)
        message += ("; run 'sourcerer schema show --source %s --collection %s "
                    "--field %s' to list them" % (self.source, collection_name,
                                                  target.name))

        return message


def fromDict(payload):
    """
    Build a SourceSchema from parsed YAML.

    Arguments:
      payload (dict): the deserialized snapshot.

    Returns:
      SourceSchema: the snapshot.

    Raises:
      SchemaError: if the snapshot is malformed or too new to understand.
    """
    version = payload.get('schema_version')
    if version is None:
        raise SchemaError('snapshot is missing schema_version')

    if version > SCHEMA_VERSION:
        raise SchemaError(
            'snapshot declares schema_version %s but this sourcerer understands '
            'at most %s; upgrade sourcerer rather than reading it partially'
            % (version, SCHEMA_VERSION))

    collections = {}
    for name, body in (payload.get('collections') or {}).items():
        fields = tuple(
            Field(name=x['name'],
                  values=tuple(x.get('values') or ()),
                  wildcard=x.get('wildcard', '*'),
                  pseudo_values=bool(x.get('pseudo_values', False)))
            for x in (body.get('fields') or []))
        collections[name] = Collection(
            name=name, fields=fields,
            reported_totals=body.get('reported_totals') or {})

    return SourceSchema(
        source=payload['source'],
        schema_version=version,
        harvested=payload.get('harvested'),
        harvested_by=payload.get('harvested_by'),
        source_urls=payload.get('source_urls') or {},
        url_rules=payload.get('url_rules') or {},
        parse_contracts=payload.get('parse_contracts') or {},
        field_aliases=payload.get('field_aliases') or {},
        collections=collections)


def toDict(schema):
    """
    Serialize a SourceSchema to plain data for YAML output.

    Ordering is fixed and keys are sorted on dump so that re-harvesting an
    unchanged source produces a byte identical file. Without that, every
    scheduled refresh would appear to be a change.

    Arguments:
      schema (SourceSchema): the snapshot.

    Returns:
      dict: plain data ready for yaml.safe_dump.
    """
    collections = {}
    for name in sorted(schema.collections):
        collection = schema.collections[name]
        collections[name] = {
            'reported_totals': dict(collection.reported_totals),
            'fields': [{'name': x.name,
                        'wildcard': x.wildcard,
                        'pseudo_values': x.pseudo_values,
                        'values': list(x.values)}
                       for x in collection.fields]}

    return {'source': schema.source,
            'schema_version': schema.schema_version,
            'harvested': schema.harvested,
            'harvested_by': schema.harvested_by,
            'source_urls': dict(schema.source_urls),
            'url_rules': dict(schema.url_rules),
            'parse_contracts': dict(schema.parse_contracts),
            'field_aliases': dict(schema.field_aliases),
            'collections': collections}


def loadSchema(source, path=None):
    """
    Load a snapshot, from the packaged data by default.

    Arguments:
      source (str): the source name, e.g. 'oas'.
      path (Path): a directory to read instead of the packaged snapshot.

    Returns:
      SourceSchema: the snapshot.

    Raises:
      SchemaError: if the snapshot is missing or unreadable.
    """
    if path is not None:
        handle = Path(path) / 'schema.yaml'
        if not handle.exists():
            raise SchemaError('no schema.yaml under %s' % path)
        text = handle.read_text()
    else:
        # importlib.resources rather than a path relative to __file__, so that the
        # snapshot is read from the installed wheel and a packaging mistake fails
        # here instead of silently reading the source tree.
        anchor = resources.files('sourcerer').joinpath('data/schemas', source,
                                                       'schema.yaml')
        if not anchor.is_file():
            raise SchemaError(
                "no packaged schema for source '%s'; run 'sourcerer schema "
                "refresh --source %s'" % (source, source))
        text = anchor.read_text()

    payload = yaml.safe_load(text)
    if not isinstance(payload, dict):
        raise SchemaError('schema.yaml for %s is not a mapping' % source)

    return fromDict(payload)


#: Fields that record when a harvest ran rather than what it found. They are
#: excluded when deciding whether a snapshot actually changed.
HARVEST_STAMPS = ('harvested', 'harvested_by')


def serializeSchema(schema):
    """
    Render a snapshot as YAML.

    Arguments:
      schema (SourceSchema): the snapshot.

    Returns:
      str: the serialized snapshot.
    """
    return yaml.safe_dump(toDict(schema), sort_keys=True,
                          default_flow_style=False, width=88)


def sameContent(left, right):
    """
    Compare two serialized snapshots ignoring the harvest timestamp.

    Arguments:
      left (str): one serialized snapshot.
      right (str): the other.

    Returns:
      bool: True if they describe the same source state.
    """
    def strip(text):
        payload = yaml.safe_load(text) or {}
        return {k: v for k, v in payload.items() if k not in HARVEST_STAMPS}

    return strip(left) == strip(right)


def saveSchema(schema, path):
    """
    Write a snapshot, but only when its content actually changed.

    A harvest that finds nothing new must leave the working tree untouched.
    Stamping a fresh timestamp on every run would make the scheduled refresh
    modify a tracked file every month, and the automation that opens a pull
    request whenever the tree is dirty would then open one every month with no
    change in it. Keeping the previous timestamp is what makes "no drift, no
    pull request" true by construction rather than by an extra condition.

    Arguments:
      schema (SourceSchema): the snapshot to write.
      path (Path): the directory to write schema.yaml into.

    Returns:
      tuple: (Path, changed) where changed is False if the file was left alone.
    """
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    handle = path / 'schema.yaml'

    payload = serializeSchema(schema)
    if handle.exists():
        existing = handle.read_text()
        if sameContent(existing, payload):
            return handle, False

    handle.write_text(payload)

    return handle, True
