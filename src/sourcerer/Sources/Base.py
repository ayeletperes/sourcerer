"""
The source interface

Five methods separate a database from everything else in the package. Below
normalizeChunk nothing knows where the data came from, so a second database
inherits the AIRR writer, the FASTA writer, the samplesheet builder, provenance
and the interactive builder without changes.
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath

log = logging.getLogger(__name__)


@dataclass(frozen=True)
class DataUnit:
    """
    One downloadable file and whatever the source says about it.

    Arguments:
      unit_id (str): the source's own opaque identifier. Never parsed into parts.
      collection (str): which collection it belongs to.
      url (str): absolute download URL.
      metadata (dict): source native metadata, unmapped.
      size_bytes (int): file size if known.
      n_sequences (int): sequence count if the source reports one.
    """
    unit_id: str
    collection: str
    url: str
    metadata: dict = field(default_factory=dict)
    size_bytes: int = None
    n_sequences: int = None

    @property
    def relpath(self):
        """PurePosixPath: mirror location, preserving the upstream path exactly."""
        return PurePosixPath(self.collection) / self.unit_id

    @property
    def study(self):
        """
        str: the leading path component, or '' when there is not one.

        Only ever used for grouping and display. Nothing depends on it being
        meaningful, because for some layouts it is the only interpretable part of
        the path and for others it is not present at all.
        """
        parts = PurePosixPath(self.unit_id).parts

        return parts[0] if len(parts) > 1 else ''


@dataclass(frozen=True)
class Query:
    """
    A resolved search request.

    Arguments:
      collection (str): which collection to search.
      filters (dict): field to value, already validated against the snapshot.
      limit (int): stop after this many units, or None for all.
    """
    collection: str
    filters: dict = field(default_factory=dict)
    limit: int = None


@dataclass(frozen=True)
class DownloadResult:
    """
    The outcome of fetching one unit.

    Arguments:
      unit (DataUnit): what was fetched.
      path (Path): where it landed.
      sha256 (str): digest of the completed file, read back from disk.
      size_bytes (int): size on disk.
      resumed (bool): whether an interrupted transfer was continued.
      skipped (bool): whether it was already present.
    """
    unit: DataUnit
    path: Path
    sha256: str
    size_bytes: int
    resumed: bool = False
    skipped: bool = False


class SourceBase(ABC):
    """
    Base class for every remote source.
    """

    #: Short name used on the commandline and as the schema directory name.
    name = None
    #: Alternative commandline names for the same source, e.g. ('airrc',) for
    #: OGRDB. They share the source's subcommand, flags and schema; the canonical
    #: ``name`` is what schema and provenance are keyed on.
    aliases = ()
    #: One line description for `sourcerer sources list`.
    description = ''
    #: Where a human can read about the source.
    homepage = ''
    #: Collections this source offers, in the order they should be presented.
    collections = ()
    #: Collection name to one line description, shown in `--help`.
    collection_help = {}
    #: License the source distributes its data under, e.g. 'CC BY 4.0'.
    #: Shown in `sourcerer sources list` and recorded in download provenance,
    #: since it is what governs how downloaded data may be reused.
    license = ''
    #: How to cite this source, one string per paper, oldest first. Shown in
    #: `sourcerer sources list` and recorded in download provenance, so that
    #: the record of what was downloaded travels with a reminder of how to
    #: give the source credit for it.
    citation = ()
    #: What the source produces, and therefore which output path `download`
    #: drives. 'dataset' sources are repertoires: they convert to AIRR/FASTA and
    #: write an airrflow samplesheet. 'reference' sources are germline sets: they
    #: build an airrflow germline reference_base instead, and never touch the
    #: rearrangement conversion path. See sourcerer.Reference.ReferenceSource.
    output = 'dataset'

    def __init__(self, client, schema=None):
        """
        Arguments:
          client (HttpClient): the shared HTTP client.
          schema (SourceSchema): the loaded snapshot, or None to load on demand.
        """
        self.client = client
        self._schema = schema

    @property
    def schema(self):
        """SourceSchema: the stored snapshot, loaded lazily."""
        if self._schema is None:
            from sourcerer.Schema import loadSchema
            self._schema = loadSchema(self.name)

        return self._schema

    @abstractmethod
    def harvestSchema(self):
        """
        Contact the live source and build a fresh snapshot.

        Returns:
          SourceSchema: the newly harvested snapshot.
        """

    @abstractmethod
    def searchUnits(self, query):
        """
        Resolve a query to concrete data units.

        Arguments:
          query (Query): the validated request.

        Returns:
          list: DataUnit objects.
        """

    @abstractmethod
    def readUnit(self, path, unit):
        """
        Open a downloaded unit.

        Arguments:
          path (Path): the downloaded file.
          unit (DataUnit): what it is.

        Returns:
          tuple: (metadata dict, iterator of raw record chunks).
        """

    @abstractmethod
    def normalizeChunk(self, metadata, chunk, unit, offset, report):
        """
        Map one chunk of raw records to AIRR named records.

        Arguments:
          metadata (dict): the unit's metadata.
          chunk: raw records.
          unit (DataUnit): what they came from.
          offset (int): index of the chunk's first row within the whole unit.
          report (dict): counters to accumulate into.

        Returns:
          pandas.DataFrame: normalized records.
        """

    def validateQuery(self, collection, filters):
        """
        Validate filters against the stored snapshot.

        Arguments:
          collection (str): the collection being searched.
          filters (dict): user supplied filters.

        Returns:
          Query: the resolved query.
        """
        resolved = self.schema.validateFilters(collection, filters)

        return Query(collection=collection, filters=resolved)

    def resolveOutputPath(self, unit, outdir):
        """
        Build the mirror location for a unit.

        The upstream relative path is preserved verbatim rather than rebuilt, so
        layouts the code has never seen still land somewhere sensible.

        Arguments:
          unit (DataUnit): the unit.
          outdir (Path): the mirror root.

        Returns:
          Path: where the unit belongs.
        """
        return Path(outdir) / unit.relpath

    def fetchUnit(self, unit, outdir, resume=True, progress=True):
        """
        Download one unit into the mirror.

        Arguments:
          unit (DataUnit): what to fetch.
          outdir (Path): the mirror root.
          resume (bool): continue an interrupted transfer if possible.
          progress (bool): show a progress bar.

        Returns:
          DownloadResult: the outcome.
        """
        dest = self.resolveOutputPath(unit, outdir)
        outcome = self.client.fetch(unit.url, dest, resume=resume,
                                    progress=progress)

        return DownloadResult(unit=unit, path=outcome.path,
                              sha256=outcome.sha256,
                              size_bytes=outcome.size_bytes,
                              resumed=outcome.resumed, skipped=outcome.skipped)

    def convertUnit(self, path, unit, chunksize=50000):
        """
        Read and normalize a downloaded unit, one chunk at a time.

        Arguments:
          path (Path): the downloaded file.
          unit (DataUnit): what it is.
          chunksize (int): rows per chunk.

        Returns:
          tuple: (metadata, generator of normalized chunks, report dict). The
          report is filled in as the generator is consumed.
        """
        from sourcerer.Sources.Oas import newReport

        metadata, chunks = self.readUnit(path, unit)
        report = newReport()

        def normalized():
            offset = 0
            for chunk in chunks:
                yield self.normalizeChunk(metadata, chunk, unit, offset, report)
                offset += len(chunk)

        return metadata, normalized(), report
