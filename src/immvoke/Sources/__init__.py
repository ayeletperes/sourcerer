"""
Source registry

A plain dictionary. Entry point based plugin discovery is a small change to make
later, once a second source exists to justify it.
"""

# Info
__author__ = 'Susanna Marquez'

# Immvoke imports
from immvoke.Sources.AirrcImgt import AirrcImgtSource
from immvoke.Sources.Imgt import ImgtSource
from immvoke.Sources.Oas import OasSource
from immvoke.Sources.Ogrdb import OgrdbSource

#: Every source immvoke knows about, by canonical commandline name.
REGISTRY = {source.name: source
            for source in (OasSource, ImgtSource, OgrdbSource, AirrcImgtSource)}

#: Alternative names that resolve to a canonical source, e.g. 'airrc' -> 'ogrdb'.
ALIASES = {alias: source.name
           for source in REGISTRY.values()
           for alias in source.aliases}


def canonicalName(name):
    """
    Resolve an alias to the canonical source name, or return it unchanged.

    Arguments:
      name (str): a source name or alias.

    Returns:
      str: the canonical source name.
    """
    return ALIASES.get(name, name)


def getSource(name, client, schema=None):
    """
    Instantiate a source by name or alias.

    Arguments:
      name (str): the source name or alias.
      client (HttpClient): the shared HTTP client.
      schema (SourceSchema): a preloaded snapshot, or None to load on demand.

    Returns:
      SourceBase: the source.

    Raises:
      KeyError: if the name is not a known source or alias.
    """
    name = canonicalName(name)
    if name not in REGISTRY:
        raise KeyError("unknown source '%s'; known sources: %s"
                       % (name, ', '.join(sorted(REGISTRY))))

    return REGISTRY[name](client, schema=schema)
