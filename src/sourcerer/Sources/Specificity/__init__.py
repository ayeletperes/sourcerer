"""
Specificity source registry

A second, separate registry from sourcerer.Sources: specificity databases
(IEDB, CoV-AbDab, ASD, and others to follow) are SourceBase subclasses too,
but are kept apart from repertoire sources so the commandline can mount them
under their own `sourcerer specificity <db>` tree instead of `sourcerer
<source>`.
"""

# Info
__author__ = 'Pramod Shinde'

# Sourcerer imports
from sourcerer.Sources.Specificity.Asd import AsdSource
from sourcerer.Sources.Specificity.CovAbDab import CovAbDabSource
from sourcerer.Sources.Specificity.Iedb import IedbSource

#: Every specificity source sourcerer knows about, by commandline name.
REGISTRY = {IedbSource.name: IedbSource, CovAbDabSource.name: CovAbDabSource,
           AsdSource.name: AsdSource}


def getSpecificitySource(name, client, schema=None):
    """
    Instantiate a specificity source by name.

    Arguments:
      name (str): the source name.
      client (HttpClient): the shared HTTP client.
      schema (SourceSchema): a preloaded snapshot, or None to load on demand.

    Returns:
      SourceBase: the source.

    Raises:
      KeyError: if the name is not registered.
    """
    if name not in REGISTRY:
        raise KeyError("unknown specificity source '%s'; known: %s"
                       % (name, ', '.join(sorted(REGISTRY))))

    return REGISTRY[name](client, schema=schema)
