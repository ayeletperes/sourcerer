"""
Unit tests for packaging metadata
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import os
import re
import tomllib
import unittest

# Sourcerer imports
import sourcerer

test_path = os.path.dirname(os.path.realpath(__file__))
repo_path = os.path.dirname(test_path)


def parseRequirement(line):
    """
    Split a requirement line into its distribution name and specifier.

    Arguments:
      line (str): a single requirements.txt or PEP 621 dependency entry.

    Returns:
      tuple: (name, specifier) with the name lowercased and dashes normalized.
    """
    match = re.match(r'^\s*([A-Za-z0-9._-]+)\s*(.*)$', line.strip())
    name, spec = match.group(1), match.group(2).strip()

    return name.lower().replace('_', '-'), spec


class TestDependencies(unittest.TestCase):
    """
    Tests keeping the two dependency lists in agreement
    """

    def setUp(self):
        with open(os.path.join(repo_path, 'pyproject.toml'), 'rb') as handle:
            self.pyproject = tomllib.load(handle)

        with open(os.path.join(repo_path, 'requirements.txt')) as handle:
            lines = [x for x in handle if x.strip() and not x.startswith('#')]

        self.requirements = dict(parseRequirement(x) for x in lines)
        self.dependencies = dict(parseRequirement(x)
                                 for x in self.pyproject['project']['dependencies'])

    def test_lists_agree(self):
        """
        requirements.txt matches [project] dependencies.

        Only pyproject.toml installs anything for a wheel user; requirements.txt
        exists for the CI minimum-dependency job. If they drift, CI tests a
        different dependency set than users actually get.
        """
        self.assertEqual(self.dependencies, self.requirements)

    def test_airr_floor_is_two(self):
        """
        The airr floor stays at 2.0.

        The streaming validation report calls RearrangementSchema.validate_header
        and validate_row, whose behaviour was verified against airr 2.0.0. Lowering
        this floor would silently disable validation reporting.
        """
        self.assertEqual(self.dependencies['airr'], '>=2.0')


class TestVersion(unittest.TestCase):
    """
    Tests for version metadata
    """

    def test_version_is_exposed(self):
        """The package exposes a PEP 440 style version."""
        self.assertRegex(sourcerer.__version__, r'^\d+\.\d+\.\d+')

    def test_hatch_reads_version_file(self):
        """Hatchling sources the version from Version.py, not a duplicate literal."""
        with open(os.path.join(repo_path, 'pyproject.toml'), 'rb') as handle:
            pyproject = tomllib.load(handle)

        self.assertEqual(pyproject['tool']['hatch']['version']['path'],
                         'src/sourcerer/Version.py')
        self.assertIn('version', pyproject['project']['dynamic'])


if __name__ == '__main__':
    unittest.main()
