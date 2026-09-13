"""
IEDB - B-cell/antibody receptor sequences and assay annotations

IEDB (https://www.iedb.org/) publishes receptor sequences as a single bulk
ZIP export and assay data through a PostgREST API
(https://query-api.iedb.org/api/v1). A table is not an index into many
separately downloadable files: the table itself *is* the data unit, so each
table produces exactly one DataUnit.

IEDB has no HTML search form to scrape, so harvestSchema returns a hand
curated snapshot rather than one built from a live page.

The receptor bulk CSVs (bcr, tcr) are written with two header rows: a row of
category labels (Receptor, Chain 1, Chain 2, ...) over a row of field names,
and one file row describes a whole receptor -- both chains at once, side by
side, rather than one row per sequence. normalizeChunk splits each receptor
row into one AIRR rearrangement record per chain it actually carries, linked
by cell_id, which is the shape `sourcerer specificity iedb download bcr`
needs to write a real rearrangement TSV instead of a source-column dump.
bcell and bcr_to_bcell carry no sequences -- one is assay results, the other
a join table -- so they keep their own column names.
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import hashlib
import json
import logging
import zipfile
from datetime import UTC
from pathlib import Path

import pandas

# Sourcerer imports
from sourcerer.Sources.Specificity.Paginate import pageByRange
from sourcerer.Convert import coerceAirrTypes
from sourcerer.Exceptions import IedbParseError
from sourcerer.Http import hashFile
from sourcerer.Sources.Base import DataUnit, DownloadResult, SourceBase

log = logging.getLogger(__name__)

#: Endpoints.
API_BASE = 'https://query-api.iedb.org/api/v1'
BULK_URL = 'https://www.iedb.org/downloader.php?file_name=doc/receptor_full_v3.zip'

#: Rows requested per page from the PostgREST API.
API_PAGE_SIZE = 2000

#: table -> (unit_id, url). bcr and tcr share one bulk ZIP; bcell and
#: bcr_to_bcell are each their own API endpoint.
_TABLE_SOURCE = {
    'bcr': ('bcr_full_v3.csv', BULK_URL),
    'tcr': ('tcr_full_v3.csv', BULK_URL),
    'bcell': ('bcell_search.json', API_BASE + '/bcell_search'),
    'bcr_to_bcell': ('bcr_to_bcell.json', API_BASE + '/bcr_to_bcell'),
}

#: Tables served from the shared receptor bulk ZIP rather than the API.
_BULK_TABLES = frozenset(['bcr', 'tcr'])

#: qualitative_measure has no server-side filter on bcell_search, so it is
#: applied client-side after fetching; this is its known controlled
#: vocabulary.
_QUALITATIVE_MEASURES = ('Positive', 'Positive-High', 'Positive-Intermediate',
                         'Positive-Low', 'Negative')

#: The two chain slots a receptor CSV row carries, in file order, numbered
#: for building each chain's sequence_id and sort rank.
_CHAIN_LABELS = (('Chain 1', 1), ('Chain 2', 2))

#: A chain's own "Type" value, as IEDB spells it, mapped to an AIRR locus.
#: 'light' (kappa vs lambda undetermined), 'IgNAR' (shark heavy-chain-only),
#: 'construct' and 'scFv' (engineered, not a single natural locus), and a
#: blank (no chain in this slot) all resolve to '' rather than a guess.
_LOCUS = {
    'heavy': 'IGH', 'kappa_light': 'IGK', 'lambda_light': 'IGL',
    'alpha': 'TRA', 'beta': 'TRB', 'gamma': 'TRG', 'delta': 'TRD',
}

#: Receptor-level columns, shared by both of a receptor's chain records,
#: mapped from their (category, field name) header pair to a plain name.
_SHARED_COLUMNS = {
    ('Receptor', 'Group IRI'): 'receptor_group_id',
    ('Receptor', 'IEDB Receptor ID'): 'iedb_receptor_id',
    ('Receptor', 'Reference Name'): 'reference_name',
    ('Receptor', 'Type'): 'receptor_type',
    ('Reference', 'IEDB IRI'): 'reference_iri',
    ('Epitope', 'IEDB IRI'): 'epitope_iri',
    ('Epitope', 'Name'): 'epitope_name',
    ('Epitope', 'Source Molecule'): 'epitope_source_molecule',
    ('Epitope', 'Source Organism'): 'epitope_source_organism',
    ('Assay', 'Type'): 'assay_type',
    ('Assay', 'IEDB IDs'): 'assay_iedb_ids',
    ('Assay', 'MHC Allele Names'): 'assay_mhc_allele_names',
}

#: Per-chain fields IEDB records twice, once curated (expert reviewed) and
#: once calculated (automated); the curated value is kept where present, the
#: calculated one otherwise. Maps the AIRR field name to (curated column,
#: calculated column).
_CHAIN_PREFERRED = {
    'v_call': ('Curated V Gene', 'Calculated V Gene'),
    'd_call': ('Curated D Gene', 'Calculated D Gene'),
    'j_call': ('Curated J Gene', 'Calculated J Gene'),
    'cdr1_aa': ('CDR1 Curated', 'CDR1 Calculated'),
    'cdr1_start': ('CDR1 Start Curated', 'CDR1 Start Calculated'),
    'cdr1_end': ('CDR1 End Curated', 'CDR1 End Calculated'),
    'cdr2_aa': ('CDR2 Curated', 'CDR2 Calculated'),
    'cdr2_start': ('CDR2 Start Curated', 'CDR2 Start Calculated'),
    'cdr2_end': ('CDR2 End Curated', 'CDR2 End Calculated'),
    'cdr3_aa': ('CDR3 Curated', 'CDR3 Calculated'),
    'cdr3_start': ('CDR3 Start Curated', 'CDR3 Start Calculated'),
    'cdr3_end': ('CDR3 End Curated', 'CDR3 End Calculated'),
}

#: Per-chain fields with one source column each, renamed straight to AIRR.
_CHAIN_DIRECT = {
    'sequence': 'Nucleotide Sequence',
    'sequence_aa': 'Protein Sequence',
    # IEDB's "Junction Calculated" runs Cys-to-Trp/Phe, the AIRR junction
    # convention; there is no separate nucleotide junction column.
    'junction_aa': 'Junction Calculated',
}

#: Per-chain fields kept verbatim under a sourcerer-owned name: no core AIRR
#: field represents them, but they are real information, not noise.
_CHAIN_EXTRA = {
    'iedb_chain_type': 'Type',
    'protein_iri': 'Protein IRI',
    'v_domain_calculated_aa': 'V Domain Calculated',
}

#: Final column order for a receptor's AIRR rearrangement records: AIRR named
#: fields first, then the sourcerer-owned extras, then the shared receptor,
#: epitope and assay context every chain of the same receptor repeats.
_AIRR_COLUMNS = (
    'sequence_id', 'cell_id', 'locus', 'sequence', 'sequence_aa',
    'v_call', 'd_call', 'j_call',
    'cdr1_aa', 'cdr1_start', 'cdr1_end',
    'cdr2_aa', 'cdr2_start', 'cdr2_end',
    'cdr3_aa', 'cdr3_start', 'cdr3_end',
    'junction_aa', 'junction_aa_length',
    'iedb_chain_type', 'protein_iri', 'v_domain_calculated_aa',
    'receptor_group_id', 'iedb_receptor_id', 'receptor_type',
    'reference_name', 'reference_iri',
    'epitope_iri', 'epitope_name', 'epitope_source_molecule',
    'epitope_source_organism',
    'assay_type', 'assay_iedb_ids', 'assay_mhc_allele_names',
)


def _readReceptorCsv(path, chunksize):
    """
    Open a receptor bulk CSV, honoring its two header rows.

    Arguments:
      path (Path): the downloaded bcr_full_v3.csv or tcr_full_v3.csv.
      chunksize (int): rows per chunk.

    Returns:
      iterator: DataFrames with a (category, field) MultiIndex for columns.
    """
    return pandas.read_csv(path, header=[0, 1], chunksize=chunksize, dtype=str,
                           na_filter=False)


def _preferCurated(chain, curated, calculated):
    """
    Pick a chain field's curated value, falling back to the calculated one.

    Arguments:
      chain (pandas.DataFrame): one chain's columns, bare field names.
      curated (str): the curated column name.
      calculated (str): the calculated column name.

    Returns:
      pandas.Series: the curated value where present, the calculated one
      elsewhere.
    """
    return chain[curated].where(chain[curated] != '', chain[calculated])


def _chainFrame(chunk, chain_label, chain_number, shared):
    """
    Build one chain's AIRR rearrangement records out of a receptor chunk.

    Arguments:
      chunk (pandas.DataFrame): raw receptor rows, (category, field) columns.
      chain_label (str): 'Chain 1' or 'Chain 2'.
      chain_number (int): 1 or 2, used in sequence_id and the sort rank.
      shared (pandas.DataFrame): the receptor-level columns, same index as
        chunk, already renamed by _SHARED_COLUMNS.

    Returns:
      pandas.DataFrame: one row per input row, including rows with no chain in
      this slot; the caller drops those using the '_present' column.
    """
    chain = chunk[chain_label]

    frame = pandas.DataFrame(index=chunk.index)
    for airr_field, (curated, calculated) in _CHAIN_PREFERRED.items():
        frame[airr_field] = _preferCurated(chain, curated, calculated)
    for airr_field, column in _CHAIN_DIRECT.items():
        frame[airr_field] = chain[column]
    for name, column in _CHAIN_EXTRA.items():
        frame[name] = chain[column]

    frame['locus'] = frame['iedb_chain_type'].map(_LOCUS).fillna('')

    has_junction = frame['junction_aa'] != ''
    frame['junction_aa_length'] = ''
    frame.loc[has_junction, 'junction_aa_length'] = \
        frame.loc[has_junction, 'junction_aa'].str.len().astype(str)

    frame['sequence_id'] = shared['iedb_receptor_id'] + ('_%d' % chain_number)
    frame['cell_id'] = shared['iedb_receptor_id']
    for name in _SHARED_COLUMNS.values():
        frame[name] = shared[name]

    # A chain slot with nothing in it at all (single-chain receptors are
    # common: a heavy-only submission, a TCR with only its beta sequenced)
    # must not become a fabricated empty rearrangement record.
    frame['_present'] = ((chain['Type'] != '') | (chain['Nucleotide Sequence'] != '')
                         | (chain['Protein Sequence'] != ''))
    # A stable sort key: chunk.index is the file's global row number, so this
    # keeps a receptor's chains adjacent and chain 1 before chain 2, the same
    # order regardless of how the file was chunked.
    frame['_row'] = chunk.index
    frame['_chain'] = chain_number

    return frame


def _receptorChunkToAirr(chunk, report):
    """
    Convert one chunk of receptor rows into AIRR rearrangement records.

    Arguments:
      chunk (pandas.DataFrame): raw receptor rows, (category, field) columns.
      report (dict): counters to accumulate into.

    Returns:
      pandas.DataFrame: one row per chain actually present, sorted so a
      receptor's chains stay adjacent with chain 1 first.
    """
    report['rows_in'] += len(chunk)

    shared = pandas.DataFrame(index=chunk.index)
    for key, name in _SHARED_COLUMNS.items():
        shared[name] = chunk[key]

    chains = [_chainFrame(chunk, label, number, shared)
             for label, number in _CHAIN_LABELS]

    frame = pandas.concat(chains, ignore_index=True)
    frame = frame[frame['_present']]
    frame = frame.sort_values(['_row', '_chain'], kind='stable')
    frame = frame[list(_AIRR_COLUMNS)].reset_index(drop=True)
    frame = coerceAirrTypes(frame)

    report['rows_out'] += len(frame)

    return frame


def _rowHash(row):
    """
    Build a short content hash for an annotation row.

    Arguments:
      row (pandas.Series): a row of upstream columns, before provenance
        columns are added.

    Returns:
      str: the first 12 hex characters of a SHA-256 digest.
    """
    key = '|'.join(str(x) for x in row.tolist())

    return hashlib.sha256(key.encode('utf-8')).hexdigest()[:12]


def _findZipMember(archive, filename):
    """
    Locate a member of the receptor bulk ZIP by its filename.

    Arguments:
      archive (zipfile.ZipFile): the opened archive.
      filename (str): the filename to find, ignoring any directory prefix.

    Returns:
      str: the matching member name.

    Raises:
      IedbParseError: if no member has that filename.
    """
    for member in archive.namelist():
        if Path(member).name == filename:
            return member

    raise IedbParseError(
        "'%s' not found in %s; the receptor bulk export layout has changed"
        % (filename, BULK_URL))


class IedbSource(SourceBase):
    """
    The IEDB specificity source.
    """

    name = 'iedb'
    description = ('IEDB: curated B-cell/antibody receptor sequences and '
                   'assay annotations')
    homepage = 'https://www.iedb.org/'
    collections = ('bcr', 'tcr', 'bcell', 'bcr_to_bcell')
    collection_help = {
        'bcr': 'BCR (antibody) receptor sequences, from the receptor bulk export',
        'tcr': 'TCR receptor sequences, from the same bulk export, kept for reference',
        'bcell': 'B-cell assay records: antigen, epitope and qualitative outcome',
        'bcr_to_bcell': 'join table linking BCR receptor groups to B-cell assay records',
    }
    airr_collections = _BULK_TABLES

    #: TODO: confirm exact license wording against https://www.iedb.org/about
    #: before this is relied on for redistribution terms.
    license = 'Freely available; see https://www.iedb.org/about'
    citation = (
        'Vita R, Mahajan S, Overton JA, Dhanda SK, Martini S, Cantrell JR, '
        'Wheeler DK, Sette A, Peters B. The Immune Epitope Database (IEDB): '
        '2018 update. Nucleic Acids Res. 2019;47(D1):D339-D343. '
        'doi:10.1093/nar/gky1006',
    )

    def harvestSchema(self):
        """
        Build the hand curated snapshot.

        There is no live form to scrape, so this returns the same field
        definitions packaged in schema.yaml, timestamped now. `sourcerer
        schema refresh --source iedb` exists mainly so the timestamp and any
        future edits to this method stay reproducible the same way a real
        harvest would.

        Returns:
          SourceSchema: the snapshot.
        """
        from datetime import datetime

        from sourcerer.Schema import Collection, Field, SourceSchema
        from sourcerer.Version import __version__

        collections = {
            'bcr': Collection(name='bcr'),
            'tcr': Collection(name='tcr'),
            'bcell': Collection(name='bcell', fields=(
                Field(name='qualitative_measure',
                      values=_QUALITATIVE_MEASURES),
            )),
            'bcr_to_bcell': Collection(name='bcr_to_bcell'),
        }

        return SourceSchema(
            source=self.name,
            harvested=datetime.now(UTC).strftime('%Y-%m-%dT%H:%M:%SZ'),
            harvested_by=('sourcerer %s (hand curated, not scraped: IEDB '
                          'exposes no HTML form)' % __version__),
            source_urls={'receptor_bulk': BULK_URL, 'api': API_BASE},
            collections=collections)

    def searchUnits(self, query):
        """
        Resolve a query to the one DataUnit its table represents.

        A table is not an index of many files; the table itself is the unit,
        so there is always exactly one.

        Arguments:
          query (Query): the validated request.

        Returns:
          list: a single DataUnit.
        """
        unit_id, url = _TABLE_SOURCE[query.collection]

        return [DataUnit(unit_id=unit_id, collection=query.collection, url=url,
                         metadata=dict(query.filters))]

    def fetchUnit(self, unit, outdir, resume=True, progress=True):
        """
        Download one table.

        Arguments:
          unit (DataUnit): what to fetch.
          outdir (Path): the mirror root.
          resume (bool): continue an interrupted transfer if possible.
          progress (bool): show a progress bar.

        Returns:
          DownloadResult: the outcome.
        """
        if unit.collection in _BULK_TABLES:
            return self._fetchFromBulkZip(unit, outdir, resume=resume,
                                          progress=progress)

        return self._fetchFromApi(unit, outdir)

    def _fetchFromBulkZip(self, unit, outdir, resume=True, progress=True):
        """
        Fetch the shared receptor bulk ZIP and extract one table from it.

        The ZIP is cached at a fixed path under outdir: a rerun with the file
        already present reuses it (HttpClient.fetch skips a download whose
        destination already exists), so fetching bcr and then tcr in the same
        output directory downloads the ZIP only once.

        Arguments:
          unit (DataUnit): the bcr or tcr table.
          outdir (Path): the mirror root.
          resume (bool): continue an interrupted ZIP transfer if possible.
          progress (bool): show a progress bar for the ZIP transfer.

        Returns:
          DownloadResult: the outcome, describing the extracted table file.
        """
        zip_dest = Path(outdir) / 'receptor_full_v3.zip'
        self.client.fetch(unit.url, zip_dest, resume=resume, progress=progress)

        dest = self.resolveOutputPath(unit, outdir)
        dest.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(zip_dest) as archive:
            member = _findZipMember(archive, unit.unit_id)
            with archive.open(member) as source, open(dest, 'wb') as target:
                target.write(source.read())

        return DownloadResult(unit=unit, path=dest, sha256=hashFile(dest),
                              size_bytes=dest.stat().st_size)

    def _fetchFromApi(self, unit, outdir):
        """
        Page the PostgREST API for one table and write it as JSON.

        Arguments:
          unit (DataUnit): the bcell or bcr_to_bcell table.
          outdir (Path): the mirror root.

        Returns:
          DownloadResult: the outcome.
        """
        records = []
        for batch in pageByRange(self.client, unit.url, page_size=API_PAGE_SIZE):
            records.extend(batch)

        dest = self.resolveOutputPath(unit, outdir)
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(json.dumps(records))

        return DownloadResult(unit=unit, path=dest, sha256=hashFile(dest),
                              size_bytes=dest.stat().st_size)

    def readUnit(self, path, unit, chunksize=50000):
        """
        Open a downloaded table.

        Arguments:
          path (Path): the downloaded file.
          unit (DataUnit): what it is.
          chunksize (int): rows per chunk, for the CSV tables. The JSON
            tables are modest enough to read as a single chunk.

        Returns:
          tuple: ({}, iterator of raw record chunks). IEDB tables carry no
          per-unit metadata line. bcr and tcr chunks carry a (category,
          field) MultiIndex for columns rather than flat column names, since
          the source CSV has two header rows.
        """
        if unit.collection in _BULK_TABLES:
            return {}, _readReceptorCsv(path, chunksize)

        with open(path) as handle:
            records = json.load(handle)

        return {}, iter([pandas.DataFrame.from_records(records)])

    def normalizeChunk(self, metadata, chunk, unit, offset, report):
        """
        Normalize a chunk, applying provenance and any filter the fetch step
        could not apply server-side.

        bcr and tcr are mapped to AIRR rearrangement records, one per chain;
        bcell and bcr_to_bcell pass their own columns through unchanged, since
        they carry no sequences to name AIRR fields after.

        Arguments:
          metadata (dict): unused; IEDB tables carry no per-unit metadata.
          chunk (pandas.DataFrame): raw records, IEDB's own columns.
          unit (DataUnit): what they came from.
          offset (int): unused; kept for interface parity with SourceBase.
          report (dict): counters to accumulate into.

        Returns:
          pandas.DataFrame: the normalized chunk, with provenance columns
          added.
        """
        if unit.collection in _BULK_TABLES:
            frame = _receptorChunkToAirr(chunk, report)
        else:
            frame = self._normalizeAssayChunk(chunk, unit, report)

        # apply(..., axis=1) on an empty frame cannot infer a row shape, so it
        # is skipped rather than trusted to return an empty result.
        if len(frame):
            row_hash = frame.apply(_rowHash, axis=1)
        else:
            row_hash = pandas.Series([], dtype=str, index=frame.index)

        frame = frame.copy()
        frame['sourcerer_source'] = self.name
        frame['sourcerer_collection'] = unit.collection
        frame['sourcerer_unit_id'] = unit.unit_id
        frame['sourcerer_row_hash'] = row_hash

        return frame

    def _normalizeAssayChunk(self, chunk, unit, report):
        """
        Pass an assay or join table chunk's columns through unchanged,
        applying any filter the fetch step could not apply server-side.

        Arguments:
          chunk (pandas.DataFrame): raw records, IEDB's own columns.
          unit (DataUnit): what they came from.
          report (dict): counters to accumulate into.

        Returns:
          pandas.DataFrame: the filtered chunk.
        """
        report['rows_in'] += len(chunk)

        frame = self._applyClientFilter(chunk, unit)

        report['rows_out'] += len(frame)

        return frame

    def _applyClientFilter(self, frame, unit):
        """
        Apply qualitative_measure filtering that bcell_search cannot do itself.

        Arguments:
          frame (pandas.DataFrame): raw records for one chunk.
          unit (DataUnit): the table, carrying the resolved filters.

        Returns:
          pandas.DataFrame: the filtered chunk, or the chunk unchanged for
          tables and filters this does not apply to.
        """
        wanted = (unit.metadata or {}).get('qualitative_measure')
        if not wanted or wanted == '*' or 'qualitative_measure' not in frame.columns:
            return frame

        return frame[frame['qualitative_measure'] == wanted]
