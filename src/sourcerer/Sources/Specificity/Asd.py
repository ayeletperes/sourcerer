"""
ASD (Antigen-Specific antibody Database) - Ab-Ag pairs as a Delta Lake

ASD (https://naturalantibody.com/agab/) publishes its whole collection as a
Delta Lake (Parquet partitions plus a _delta_log transaction log) staged in a
Google Drive folder behind a Colab notebook; there is no plain HTTP file to
GET. fetchUnit therefore reaches for gdown rather than the shared
HttpClient, and what it downloads is one folder rather than one file.

A table is not an index into many separately downloadable files: the Delta
Lake itself *is* the data unit, so there is exactly one. ASD has no HTML
search form, so harvestSchema returns a hand curated snapshot rather than
one built from a live page; the one filterable field (dataset) is the
controlled vocabulary observed directly in the Delta Lake, since the Drive
folder offers no per-dataset export to filter server side.
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import hashlib
import json
import logging
from datetime import UTC
from pathlib import Path

import pandas

# Sourcerer imports
from sourcerer.Exceptions import AsdFetchError, AsdParseError
from sourcerer.Http import hashFile
from sourcerer.Sources.Base import DataUnit, DownloadResult, SourceBase

log = logging.getLogger(__name__)

#: Endpoints. There is no per-file download URL; the whole Delta Lake is one
#: Google Drive folder, reached through the Colab notebook linked from HOME_URL.
HOME_URL = 'https://naturalantibody.com/agab/'
DRIVE_FOLDER_URL = ('https://drive.google.com/drive/folders/'
                    '1HY2GOVj-HR8t6Jmhe5MIJ6xmBCOxC45q')

#: The one table this source offers; the whole Delta Lake is the data unit.
COLLECTION = 'asd'

#: dataset has no server-side filter on the Drive folder download, so it is
#: applied client-side after fetching; this is its known controlled
#: vocabulary, observed directly in the Delta Lake.
_DATASETS = (
    'aae', 'aatp', 'ab-bind', 'abbd', 'abdesign', 'alphaseq', 'biomap',
    'buzz', 'covid-19', 'dlgo', 'flab_hie2022', 'flab_koenig2017',
    'flab_rosace2023', 'flab_shanehsazzadeh2023', 'flab_warszawski2019',
    'genbank', 'hiv', 'literature', 'met', 'osh', 'patents', 'rmna',
    'skempiv2', 'structures-antibodies', 'structures-nanobodies',
)

#: Columns whose values are nested (Arrow struct/list), which pandas reads
#: back as Python dicts/lists. A TSV cell that is str()'d from a dict is not
#: valid JSON and cannot be round tripped, so these are serialized instead.
_NESTED_COLUMNS = ('metadata',)


def _hashDirectory(path):
    """
    Build a content digest for a downloaded Delta Lake directory.

    HttpClient.fetch hashes a single completed file; a Delta Lake is a
    directory of them, so this combines a hash of every file's relative path
    and content into one digest, computed in a fixed (sorted) order so the
    result is stable regardless of the order gdown wrote files in.

    Arguments:
      path (Path): the Delta Lake directory.

    Returns:
      str: hexadecimal SHA-256 digest.
    """
    digest = hashlib.sha256()
    for item in sorted(Path(path).rglob('*')):
        if item.is_file():
            digest.update(str(item.relative_to(path)).encode('utf-8'))
            digest.update(hashFile(item).encode('utf-8'))

    return digest.hexdigest()


def _directorySize(path):
    """
    Total size in bytes of every file under a Delta Lake directory.

    Arguments:
      path (Path): the Delta Lake directory.

    Returns:
      int: total size in bytes.
    """
    return sum(item.stat().st_size for item in Path(path).rglob('*')
              if item.is_file())


def _rowHash(row):
    """
    Build a short content hash for an Ab-Ag pair row.

    Arguments:
      row (pandas.Series): a row of upstream columns, before provenance
        columns are added.

    Returns:
      str: the first 12 hex characters of a SHA-256 digest.
    """
    key = '|'.join(str(x) for x in row.tolist())

    return hashlib.sha256(key.encode('utf-8')).hexdigest()[:12]


class AsdSource(SourceBase):
    """
    The ASD specificity source.
    """

    name = 'asd'
    description = 'ASD: antigen-specific antibody and nanobody Ab-Ag pairs'
    homepage = HOME_URL
    collections = (COLLECTION,)
    collection_help = {
        COLLECTION: ('the full Ab-Ag pair Delta Lake (structures, affinity '
                     'assays, literature and patent mined pairs, …), '
                     'optionally narrowed to one dataset'),
    }

    #: TODO: confirm license wording against https://naturalantibody.com/agab/
    #: before this is relied on for redistribution terms.
    license = ('CC BY-NC 4.0; see the LICENSE.txt distributed alongside the '
               'Delta Lake in the Drive folder')
    citation = (
        'naturalantibody.com. ASD: Antigen-Specific antibody Database. '
        'https://naturalantibody.com/agab/',
    )

    def harvestSchema(self):
        """
        Build the hand curated snapshot.

        There is no live form to scrape, so this returns the same field
        definitions packaged in schema.yaml, timestamped now. `sourcerer
        schema refresh --source asd` exists mainly so the timestamp and any
        future edits to this method stay reproducible the same way a real
        harvest would.

        Returns:
          SourceSchema: the snapshot.
        """
        from datetime import datetime

        from sourcerer.Schema import Collection, Field, SourceSchema
        from sourcerer.Version import __version__

        collections = {
            COLLECTION: Collection(name=COLLECTION, fields=(
                Field(name='dataset', values=_DATASETS),
            )),
        }

        return SourceSchema(
            source=self.name,
            harvested=datetime.now(UTC).strftime('%Y-%m-%dT%H:%M:%SZ'),
            harvested_by=('sourcerer %s (hand curated, not scraped: ASD '
                          'exposes no HTML form)' % __version__),
            source_urls={'home': HOME_URL, 'drive_folder': DRIVE_FOLDER_URL},
            collections=collections)

    def searchUnits(self, query):
        """
        Resolve a query to the one DataUnit the Delta Lake represents.

        Arguments:
          query (Query): the validated request.

        Returns:
          list: a single DataUnit.
        """
        return [DataUnit(unit_id=COLLECTION, collection=query.collection,
                         url=DRIVE_FOLDER_URL, metadata=dict(query.filters))]

    def fetchUnit(self, unit, outdir, resume=True, progress=True):
        """
        Download the Delta Lake folder from Google Drive via gdown.

        Unlike HttpClient.fetch, gdown has no resumable-range support, so
        `resume` only controls whether an already complete looking directory
        (one that already has *.parquet files) is trusted and left alone,
        not whether a partial transfer is continued.

        Arguments:
          unit (DataUnit): the asd table.
          outdir (Path): the mirror root.
          resume (bool): reuse a directory that already looks complete.
          progress (bool): show gdown's own progress output.

        Returns:
          DownloadResult: the outcome.

        Raises:
          AsdFetchError: if gdown is not installed, the folder could not be
            fetched, or the result contains no *.parquet files.
        """
        dest = self.resolveOutputPath(unit, outdir)

        if resume and dest.exists() and any(dest.rglob('*.parquet')):
            log.info('%s already present; skipping', dest)
        else:
            self._downloadDeltaLake(unit.url, dest, progress=progress)

        return DownloadResult(unit=unit, path=dest,
                              sha256=_hashDirectory(dest),
                              size_bytes=_directorySize(dest))

    @staticmethod
    def _downloadDeltaLake(url, dest, progress=True):
        """
        Download one Google Drive folder into dest with gdown.

        Arguments:
          url (str): the Drive folder URL.
          dest (Path): where to write it.
          progress (bool): show gdown's own progress output.

        Raises:
          AsdFetchError: if gdown is not installed, the download failed, or
            the folder contained no *.parquet files.
        """
        try:
            import gdown
        except ImportError as error:
            raise AsdFetchError(
                "gdown is required to download ASD's Delta Lake from Google "
                "Drive; install it with 'pip install gdown'") from error

        dest.mkdir(parents=True, exist_ok=True)
        try:
            gdown.download_folder(url, output=str(dest), quiet=not progress,
                                  use_cookies=False)
        except Exception as error:
            raise AsdFetchError(
                'failed to download the ASD Delta Lake folder from %s: %s'
                % (url, error)) from error

        if not any(dest.rglob('*.parquet')):
            raise AsdFetchError(
                'gdown reported success but no *.parquet files were found '
                'under %s; the ASD Drive folder layout has changed' % dest)

    def readUnit(self, path, unit, chunksize=None):
        """
        Open the downloaded Delta Lake, one partition file per chunk.

        Arguments:
          path (Path): the downloaded Delta Lake directory.
          unit (DataUnit): what it is.
          chunksize (int): unused; each Parquet partition is already a
            manageably sized chunk. Kept for interface parity with SourceBase.

        Returns:
          tuple: ({}, iterator of raw record chunks). ASD carries no
          per-unit metadata line.

        Raises:
          AsdParseError: if the directory contains no *.parquet partitions.
        """
        partitions = sorted(Path(path).glob('*.parquet'))
        if not partitions:
            raise AsdParseError(
                'no *.parquet files under %s; the ASD Delta Lake layout has '
                'changed' % path)

        return {}, (pandas.read_parquet(x) for x in partitions)

    def normalizeChunk(self, metadata, chunk, unit, offset, report):
        """
        Pass a chunk's columns through unchanged, applying provenance and any
        dataset filter the fetch step could not apply server-side.

        Arguments:
          metadata (dict): unused; ASD carries no per-unit metadata.
          chunk (pandas.DataFrame): raw records, ASD's own columns.
          unit (DataUnit): what they came from.
          offset (int): unused; kept for interface parity with SourceBase.
          report (dict): counters to accumulate into.

        Returns:
          pandas.DataFrame: the chunk with provenance columns added.
        """
        report['rows_in'] += len(chunk)

        frame = self._applyClientFilter(chunk, unit)
        frame = frame.copy()

        for column in _NESTED_COLUMNS:
            if column in frame.columns:
                frame[column] = frame[column].apply(
                    lambda x: json.dumps(x, default=str) if x is not None else '')

        # apply(..., axis=1) on an empty frame cannot infer a row shape, so it
        # is skipped rather than trusted to return an empty result.
        if len(frame):
            row_hash = frame.apply(_rowHash, axis=1)
        else:
            row_hash = pandas.Series([], dtype=str, index=frame.index)

        frame['sourcerer_source'] = self.name
        frame['sourcerer_collection'] = unit.collection
        frame['sourcerer_unit_id'] = unit.unit_id
        frame['sourcerer_row_hash'] = row_hash

        report['rows_out'] += len(frame)

        return frame

    def _applyClientFilter(self, frame, unit):
        """
        Apply dataset filtering that the Drive folder download cannot do itself.

        Arguments:
          frame (pandas.DataFrame): raw records for one chunk.
          unit (DataUnit): the table, carrying the resolved filters.

        Returns:
          pandas.DataFrame: the filtered chunk, or the chunk unchanged when no
          filter applies.
        """
        wanted = (unit.metadata or {}).get('dataset')
        if not wanted or wanted == '*' or 'dataset' not in frame.columns:
            return frame

        return frame[frame['dataset'] == wanted]
