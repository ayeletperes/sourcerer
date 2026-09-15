"""
Unit tests for schema snapshots
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import tempfile
import unittest
from pathlib import Path

# Sourcerer imports
from sourcerer.Exceptions import SchemaError
from sourcerer.Schema import (
    SCHEMA_VERSION,
    Collection,
    Field,
    SourceSchema,
    fromDict,
    loadSchema,
    saveSchema,
    toDict,
)


def makeSchema(harvested='2026-08-04T00:00:00Z'):
    """Build a small two collection snapshot with disjoint vocabularies."""
    return SourceSchema(
        source='demo', harvested=harvested, harvested_by='test',
        field_aliases={'Organism': 'Species'},
        collections={
            'paired': Collection(name='paired', fields=(
                Field(name='Species', values=('human', 'rat_SD')),
                Field(name='Age', pseudo_values=True))),
            'unpaired': Collection(name='unpaired', fields=(
                Field(name='Species', values=('rabbit', 'camel')),
                Field(name='Primer', values=('p1',))))})


class TestValidation(unittest.TestCase):
    """
    Tests for filter validation
    """

    def setUp(self):
        self.schema = makeSchema()

    def test_defaults_to_wildcards(self):
        resolved = self.schema.validateFilters('paired', {})
        self.assertEqual(resolved, {'Species': '*', 'Age': '*'})

    def test_accepts_a_known_value(self):
        resolved = self.schema.validateFilters('paired', {'Species': 'human'})
        self.assertEqual(resolved['Species'], 'human')

    def test_rejects_an_unknown_value_with_suggestions(self):
        """
        An unknown value is an error, not a warning.

        Sending it upstream returns zero results, which reads as "no data
        matched" rather than "that value does not exist".
        """
        with self.assertRaises(SchemaError) as raised:
            self.schema.validateFilters('paired', {'Species': 'humn'})

        self.assertIn('human', str(raised.exception))

    def test_rejects_an_unknown_field(self):
        with self.assertRaises(SchemaError):
            self.schema.validateFilters('paired', {'Nonesuch': 'x'})

    def test_aliases_resolve(self):
        resolved = self.schema.validateFilters('paired', {'Organism': 'human'})
        self.assertEqual(resolved['Species'], 'human')

    def test_presence_only_fields_reject_real_values(self):
        with self.assertRaises(SchemaError):
            self.schema.validateFilters('paired', {'Age': '35'})

        self.assertEqual(
            self.schema.validateFilters('paired', {'Age': 'defined'})['Age'],
            'defined')

    def test_collections_never_share_a_vocabulary(self):
        """
        Validating one collection must never consult another's values.

        The predecessor tool validated paired searches against an index that
        contains only unpaired units, so paired-only values such as rat_SD were
        reported as unknown and then sent anyway.
        """
        self.assertEqual(
            self.schema.validateFilters('paired', {'Species': 'rat_SD'})['Species'],
            'rat_SD')

        with self.assertRaises(SchemaError):
            self.schema.validateFilters('unpaired', {'Species': 'rat_SD'})

        with self.assertRaises(SchemaError):
            self.schema.validateFilters('paired', {'Species': 'camel'})

    def test_unknown_collection_raises(self):
        with self.assertRaises(SchemaError):
            self.schema.validateFilters('nonesuch', {})


class TestRoundTrip(unittest.TestCase):
    """
    Tests for serialization
    """

    def test_round_trip(self):
        schema = makeSchema()
        again = fromDict(toDict(schema))

        self.assertEqual(again.collection_names, schema.collection_names)
        self.assertEqual(again.getCollection('unpaired').field_names,
                         ('Species', 'Primer'))

    def test_a_newer_snapshot_is_refused(self):
        """
        A snapshot from a future version is refused rather than partly read.

        Reading only the parts this version understands would silently drop
        fields and make searches quietly narrower than the user asked for.
        """
        payload = toDict(makeSchema())
        payload['schema_version'] = SCHEMA_VERSION + 1

        with self.assertRaises(SchemaError):
            fromDict(payload)


class TestQuietWrite(unittest.TestCase):
    """
    Tests that an unchanged harvest leaves the working tree alone
    """

    def test_unchanged_content_is_not_rewritten(self):
        """
        Re-harvesting an unchanged source must not modify the file.

        The scheduled refresh opens a pull request whenever the tree is dirty, so
        stamping a new timestamp every month would produce twelve empty pull
        requests a year and train everyone to ignore them.
        """
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            written, changed = saveSchema(makeSchema(), path)
            self.assertTrue(changed)
            before = written.read_text()

            # Same source state, harvested later.
            _, changed = saveSchema(makeSchema(harvested='2026-09-01T00:00:00Z'),
                                    path)

            self.assertFalse(changed)
            self.assertEqual(written.read_text(), before)

    def test_real_change_is_written(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            saveSchema(makeSchema(), path)

            changed_schema = makeSchema()
            changed_schema.collections['paired'] = Collection(
                name='paired', fields=(Field(name='Species',
                                             values=('human', 'rat_SD', 'newt')),))
            _, changed = saveSchema(changed_schema, path)

            self.assertTrue(changed)


class TestPackagedSnapshot(unittest.TestCase):
    """
    Tests for the snapshot that ships with the package
    """

    def test_oas_snapshot_loads_from_the_installed_package(self):
        """
        Loading goes through importlib.resources, so a packaging mistake fails
        here rather than on a user's machine.
        """
        schema = loadSchema('oas')

        self.assertEqual(schema.source, 'oas')
        self.assertEqual(set(schema.collection_names), {'paired', 'unpaired'})

    def test_the_two_oas_collections_differ_as_the_live_forms_do(self):
        schema = loadSchema('oas')
        paired = set(schema.getCollection('paired').field_names)
        unpaired = set(schema.getCollection('unpaired').field_names)

        self.assertNotIn('Isotype', paired)
        self.assertIn('Isotype', unpaired)
        self.assertIn('Primer', unpaired)
        self.assertTrue({'BSource', 'BType', 'Subject'} <= paired)


if __name__ == '__main__':
    unittest.main()
