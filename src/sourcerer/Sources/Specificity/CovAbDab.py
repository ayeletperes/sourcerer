"""
CoV-AbDab - coronavirus-binding antibody and nanobody entries

CoV-AbDab (https://opig.stats.ox.ac.uk/webapps/covabdab/) publishes its whole
collection as a single bulk CSV, so the one table this source offers *is*
that CSV rather than an index into many separately downloadable files. The
CSV is not served from a fixed endpoint, though: its filename is stamped
with the release date (e.g. CoV-AbDab_080224.csv), so the current download
URL has to be read off the homepage, and that costs an HTTP request every
time the table is resolved.

There is no search form to build a controlled vocabulary from, so
harvestSchema returns a hand curated snapshot rather than one scraped from a
live page.
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import hashlib
import logging
from datetime import UTC

import pandas
from bs4 import BeautifulSoup

# Sourcerer imports
from sourcerer.Exceptions import CovAbDabParseError
from sourcerer.Sources.Base import DataUnit, SourceBase

log = logging.getLogger(__name__)

#: Endpoints.
HOST = 'https://opig.stats.ox.ac.uk'
HOME_URL = HOST + '/webapps/covabdab/'

#: The one table this source offers; the whole CSV is the data unit.
COLLECTION = 'antibodies'


def findCsvUrl(client):
    """
    Read the current bulk CSV download link off the CoV-AbDab homepage.

    Arguments:
      client (HttpClient): the shared HTTP client.

    Returns:
      str: the absolute CSV URL.

    Raises:
      CovAbDabParseError: if no download link can be found, meaning the page
        layout has changed.
    """
    html = client.get(HOME_URL).text
    soup = BeautifulSoup(html, 'html.parser')

    for link in soup.find_all('a', href=True):
        href = link['href']
        if 'download' in href.lower() and href.lower().endswith('.csv'):
            return href if href.startswith('http') else HOST + href

    raise CovAbDabParseError(
        'no CSV download link found on %s; the homepage layout has changed'
        % HOME_URL)


def _rowHash(row):
    """
    Build a short content hash for an antibody row.

    Arguments:
      row (pandas.Series): a row of upstream columns, before provenance
        columns are added.

    Returns:
      str: the first 12 hex characters of a SHA-256 digest.
    """
    key = '|'.join(str(x) for x in row.tolist())

    return hashlib.sha256(key.encode('utf-8')).hexdigest()[:12]


class CovAbDabSource(SourceBase):
    """
    The CoV-AbDab specificity source.
    """

    name = 'covabdab'
    description = ('CoV-AbDab: coronavirus-binding antibody and nanobody '
                   'sequences')
    homepage = HOME_URL
    collections = (COLLECTION,)
    collection_help = {
        COLLECTION: ('antibody and nanobody entries against SARS-CoV-2 and '
                     'related coronaviruses, from the bulk CSV export'),
    }

    license = 'CC BY-NC 4.0 (https://creativecommons.org/licenses/by-nc/4.0/)'
    citation = (
        'Raybould MIJ, Kovaltsuk A, Marks C, Deane CM. CoV-AbDab: the '
        'Coronavirus Antibody Database. Bioinformatics. 2021;37(5):734-735. '
        'doi:10.1093/bioinformatics/btaa739',
    )

    def harvestSchema(self):
        """
        Build the hand curated snapshot.

        There is no live form to scrape for a controlled vocabulary, so this
        returns the same field definitions packaged in schema.yaml, timestamped
        now. `sourcerer schema refresh --source covabdab` exists mainly so the
        timestamp and any future edits to this method stay reproducible the
        same way a real harvest would.

        Returns:
          SourceSchema: the snapshot.
        """
        from datetime import datetime

        from sourcerer.Schema import Collection, SourceSchema
        from sourcerer.Version import __version__

        collections = {COLLECTION: Collection(name=COLLECTION)}

        return SourceSchema(
            source=self.name,
            harvested=datetime.now(UTC).strftime('%Y-%m-%dT%H:%M:%SZ'),
            harvested_by=('sourcerer %s (hand curated, not scraped: '
                          'CoV-AbDab exposes no HTML form)' % __version__),
            source_urls={'home': HOME_URL},
            collections=collections)

    def searchUnits(self, query):
        """
        Resolve a query to the one DataUnit the current CSV represents.

        The CSV's filename and URL change with every release, so this reads
        the current download link off the homepage rather than caching one.

        Arguments:
          query (Query): the validated request.

        Returns:
          list: a single DataUnit.
        """
        url = findCsvUrl(self.client)
        unit_id = url.rsplit('/', 1)[-1]

        return [DataUnit(unit_id=unit_id, collection=query.collection, url=url)]

    def readUnit(self, path, unit, chunksize=50000):
        """
        Open the downloaded CSV.

        Arguments:
          path (Path): the downloaded file.
          unit (DataUnit): what it is.
          chunksize (int): rows per chunk.

        Returns:
          tuple: ({}, iterator of raw record chunks). CoV-AbDab carries no
          per-unit metadata line.
        """
        # utf-8-sig: the export is saved with a leading BOM.
        chunks = pandas.read_csv(path, chunksize=chunksize, dtype=str,
                                 na_filter=False, encoding='utf-8-sig')

        return {}, chunks

    def normalizeChunk(self, metadata, chunk, unit, offset, report):
        """
        Pass a chunk's columns through unchanged, applying provenance.

        Arguments:
          metadata (dict): unused; CoV-AbDab carries no per-unit metadata.
          chunk (pandas.DataFrame): raw records, CoV-AbDab's own columns.
          unit (DataUnit): what they came from.
          offset (int): unused; kept for interface parity with SourceBase.
          report (dict): counters to accumulate into.

        Returns:
          pandas.DataFrame: the chunk with provenance columns added.
        """
        report['rows_in'] += len(chunk)

        # apply(..., axis=1) on an empty frame cannot infer a row shape, so it
        # is skipped rather than trusted to return an empty result.
        if len(chunk):
            row_hash = chunk.apply(_rowHash, axis=1)
        else:
            row_hash = pandas.Series([], dtype=str, index=chunk.index)

        frame = chunk.copy()
        frame['sourcerer_source'] = self.name
        frame['sourcerer_collection'] = unit.collection
        frame['sourcerer_unit_id'] = unit.unit_id
        frame['sourcerer_row_hash'] = row_hash

        report['rows_out'] += len(frame)

        return frame
