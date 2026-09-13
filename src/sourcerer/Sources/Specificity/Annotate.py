"""
Antigen specificity annotation: query-vs-reference receptor comparison

For each row of a query AIRR table, find every row of a fixed reference AIRR
table it agrees with -- one annotation per query row (`query_row`,
`n_hits_total`, `is_unique_hit`, `hit_ids`). This is the third pipeline
stage: download a specificity database (`sourcerer specificity iedb download
bcr`), which already lands as a real AIRR rearrangement TSV, then annotate a
user's own AIRR table against it as the reference.

Three ways to decide a hit today:

    runExactAnnotation        a single column, byte-identical, CompAIRR
                              accelerated (hash-based lookup, optionally
                              gated by exact V/J-gene agreement)
    runFuzzyAnnotation        one or more columns, up to `max_mismatches`
                              substitutions allowed (pooled across columns
                              if more than one), pure Python
    runLevenshteinAnnotation  a single column, up to `max_edits`
                              substitutions/insertions/deletions combined,
                              pure Python (rapidfuzz)

runPairedExactAnnotation builds on runExactAnnotation for the case where a
query row is one chain of an antibody rather than a self-contained unit: it
runs the exact criterion once per chain (heavy, light) and reconciles the
two per antibody, classifying each as agreeing with the same reference
antibody on both chains, one chain only, or neither.

runExactAnnotation requires the compairr binary
(https://github.com/uio-bmi/compairr, Rognes et al. 2022, Bioinformatics),
built from source -- no external deps beyond a C++11 compiler and make:

    git clone https://github.com/uio-bmi/compairr.git
    cd compairr && make

then either put the resulting `compairr` binary on PATH, or pass its path as
`compairr_bin=`. compairr compares sequences literally and has no concept of
the "X" wildcard runFuzzyAnnotation gives unresolved (e.g. PDB-derived)
residues, and hard-errors on an empty/missing sequence value rather than
treating it as a length-0 field -- see runExactAnnotation's docstring for how
those rows are handled instead.

runLevenshteinAnnotation has neither the wildcard handling nor the
length-bucketed index the other two share: Levenshtein distance lets a query
and a candidate differ in length, so length can't prefilter candidates, and a
wildcard has no clean meaning once an alignment can shift position. It also
takes only a single column, unlike runFuzzyAnnotation -- see its section
header for why edit distance does not pool across columns the way Hamming
distance does.

Substitution-matrix-weighted (BLOSUM) annotation is intentionally not
included in this first pass -- see the sibling project this was ported from
(`analysis/matching/matching.py` in specificity-annotations) if it is needed
later.
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas
from rapidfuzz.distance import Levenshtein

# Sourcerer imports
from sourcerer.Exceptions import AnnotationError

#: Reference column read as each hit's identifier when the caller does not
#: say otherwise -- the AIRR rearrangement primary key.
DEFAULT_ID_COL = 'sequence_id'
#: A position where either side has this character never counts as a
#: mismatch. Only meaningful to runFuzzyAnnotation and the pure-Python
#: fallback runExactAnnotation uses for rows compairr cannot handle --
#: compairr itself has no wildcard concept.
WILDCARD = 'X'

#: runPairedExactAnnotation's fixed column names -- one antibody's heavy and
#: light chain are two rows of the same AIRR table, linked by `cell_id` and
#: told apart by `locus`; there is nothing to make configurable here beyond
#: what the AIRR schema already fixes these fields to mean.
CELL_ID_COL = 'cell_id'
LOCUS_COL = 'locus'
HEAVY_LOCUS = 'IGH'
LIGHT_LOCI = ('IGK', 'IGL')


def _valueDistance(a, b):
    """
    Count of substitution mismatches between `a` and `b`.

    Arguments:
      a (str): the first value.
      b (str): the second value.

    Returns:
      int: mismatch count, or None if the lengths differ -- this only models
      substitutions, not insertions/deletions. A WILDCARD on either side of a
      position is never a mismatch.
    """
    if len(a) != len(b):
        return None

    return sum(ca != cb and ca != WILDCARD and cb != WILDCARD
              for ca, cb in zip(a, b))


def _rowDistance(a, b):
    """
    Total substitution mismatches summed across every comparison column.

    Arguments:
      a (tuple): one row's values, one per comparison column.
      b (tuple): the row it is being compared to, same shape.

    Returns:
      int: the summed distance, or None if any column pair differs in length.
    """
    total = 0
    for x, y in zip(a, b):
        distance = _valueDistance(x, y)
        if distance is None:
            return None
        total += distance

    return total


def _mismatchPositions(a, b):
    """
    Every substitution mismatch between `a` and `b`.

    Arguments:
      a (tuple): one row's values, one per comparison column.
      b (tuple): the row it is being compared to, same shape.

    Returns:
      list: (column index, 0-based position) pairs, same wildcard rule as
      _valueDistance.
    """
    return [(col, pos)
           for col, (x, y) in enumerate(zip(a, b))
           for pos, (cx, cy) in enumerate(zip(x, y))
           if cx != cy and cx != WILDCARD and cy != WILDCARD]


def _formatPositions(positions, n_cols):
    """
    Render mismatch positions as a comma-separated string.

    Arguments:
      positions (list): (column index, position) pairs, from
        _mismatchPositions.
      n_cols (int): how many comparison columns are in play.

    Returns:
      str: bare positions ("3,7") for the common single-column case,
      column-prefixed ("0:3,1:2") when there is more than one -- a bare
      position is ambiguous once more than one column is being compared.
    """
    if n_cols == 1:
        return ','.join(str(pos) for _, pos in positions)

    return ','.join('%d:%d' % (col, pos) for col, pos in positions)


def _rowValues(frame, columns):
    """
    Read one or more columns as a list of per-row tuples.

    Arguments:
      frame (pandas.DataFrame): the table to read from.
      columns (list): column names, in the order values should be compared.

    Returns:
      list: one tuple per row, values in `columns` order, blank standing in
      for a missing value.
    """
    return list(frame[columns].fillna('').itertuples(index=False, name=None))


def _bucketKey(row_values, exact_values):
    """
    Build the index bucket a row belongs to.

    Arguments:
      row_values (tuple): the row's comparison column values.
      exact_values (tuple): the row's exact column values, empty if none.

    Returns:
      tuple: exact column values followed by the length of each comparison
      column value -- substitutions can't change length, so two rows can
      only agree when every comparison column is the same length.
    """
    return exact_values + tuple(len(value) for value in row_values)


def buildAnnotationIndex(reference, compare_cols, id_col=DEFAULT_ID_COL,
                         exact_cols=None):
    """
    Bucket a reference table for fast lookup.

    Arguments:
      reference (pandas.DataFrame): the fixed table being annotated against.
      compare_cols (list): column names to compare.
      id_col (str): column holding each reference row's identifier.
      exact_cols (list): column names that must agree exactly before
        `compare_cols` are even compared, e.g. `['v_call']` to only ever
        compare CDR3s within the same V-gene. None for no such gate.

    Returns:
      dict: bucket key (see _bucketKey) to a list of (row_values, id) pairs.
    """
    values = _rowValues(reference, compare_cols)
    exact_values = (_rowValues(reference, exact_cols) if exact_cols
                    else [()] * len(reference))
    ids = reference[id_col].tolist()

    index = {}
    for row_values, row_exact, row_id in zip(values, exact_values, ids):
        index.setdefault(_bucketKey(row_values, row_exact), []).append(
            (row_values, row_id))

    return index


def runFuzzyAnnotation(queries, reference, compare_cols, max_mismatches,
                       id_col=DEFAULT_ID_COL, exact_cols=None,
                       threshold_is_fraction=False):
    """
    Annotate each query row against `reference` allowing up to
    `max_mismatches` substitutions across `compare_cols` (pooled if more
    than one) instead of requiring exact agreement.

    `exact_cols` matters beyond speed: an unconstrained CDR3-only fuzzy
    annotation is prone to coincidental hits (especially on short CDR3s,
    where 1-2 substitutions is a large fraction of the sequence) that don't
    share a V-gene with the query and are unlikely to be a real related
    receptor.

    Arguments:
      queries (pandas.DataFrame): the input rows to annotate.
      reference (pandas.DataFrame): the fixed table being annotated against.
      compare_cols (list): column names to compare.
      max_mismatches (float): the mismatch tolerance -- an absolute count, or
        a fraction of the query's length when `threshold_is_fraction` is True.
      id_col (str): column in `reference` holding each row's identifier.
      exact_cols (list): column names that must agree exactly first, e.g.
        `['v_call']`. None for no such gate.
      threshold_is_fraction (bool): treat `max_mismatches` as a fraction of
        the query's summed comparison column length (e.g. 0.15 allows up to
        15% mismatches) instead of a fixed count. Resolved to an absolute
        count once per query row, since every candidate in its bucket
        already shares its exact length.

    Returns:
      pandas.DataFrame: one row per query row, columns `query_row`,
      `n_hits_total`, `min_mismatches`, `is_unique_hit`, `hit_ids`,
      `mismatch_positions` (per id in `hit_ids`, same order; column
      prefixed when `compare_cols` has more than one column).
    """
    index = buildAnnotationIndex(reference, compare_cols, id_col, exact_cols)
    values = _rowValues(queries, compare_cols)
    exact_values = (_rowValues(queries, exact_cols) if exact_cols
                    else [()] * len(values))

    rows = []
    for query_row, (row_values, row_exact) in enumerate(zip(values, exact_values)):
        candidates = index.get(_bucketKey(row_values, row_exact), [])
        max_distance = (int(max_mismatches * sum(len(v) for v in row_values))
                       if threshold_is_fraction else max_mismatches)

        best_by_id = {}
        for cand_values, row_id in candidates:
            distance = _rowDistance(row_values, cand_values)
            if distance is not None and distance <= max_distance:
                if row_id not in best_by_id or distance < best_by_id[row_id][0]:
                    best_by_id[row_id] = (distance, cand_values)

        hit_ids = sorted(best_by_id, key=lambda x: (best_by_id[x][0], str(x)))
        mismatch_strs = [
            _formatPositions(_mismatchPositions(row_values, best_by_id[x][1]),
                            len(compare_cols))
            for x in hit_ids
        ]
        rows.append({
            'query_row': query_row,
            'n_hits_total': len(hit_ids),
            # float('nan'), not None: keeps the column a uniform float64
            # rather than a mixed-type object column when some query rows
            # have no hit at all.
            'min_mismatches': float(min(d for d, _ in best_by_id.values()))
                              if best_by_id else float('nan'),
            'is_unique_hit': len(hit_ids) <= 1,
            'hit_ids': ';'.join(str(x) for x in hit_ids),
            'mismatch_positions': ';'.join(mismatch_strs),
        })

    return pandas.DataFrame(rows)


# ── CompAIRR-accelerated exact annotation (single column) ──────────────────

def _findCompairr(compairr_bin):
    """
    Resolve the compairr binary to invoke.

    Arguments:
      compairr_bin (str): an explicit path, or None to search PATH.

    Returns:
      str: the resolved path.

    Raises:
      AnnotationError: if compairr is not on PATH and no path was given.
    """
    path = compairr_bin or shutil.which('compairr')
    if path is None:
        raise AnnotationError(
            "compairr binary not found. Build it from source and either put "
            "it on PATH or pass compairr_bin=<path>:\n"
            "  git clone https://github.com/uio-bmi/compairr.git\n"
            "  cd compairr && make")

    return path


def _writeAirrTsv(frame, id_col, seq_col, repertoire_id, path, vgene_col,
                  jgene_col=None):
    """
    Write one side of a compairr comparison as a minimal AIRR TSV.

    Arguments:
      frame (pandas.DataFrame): the rows to write.
      id_col (str): column to write as `sequence_id`.
      seq_col (str): column to write as `junction_aa`.
      repertoire_id (str): constant `repertoire_id` value -- compairr
        compares one sequence field per repertoire; both sides here are
        exactly one repertoire each.
      path (Path): where to write the TSV.
      vgene_col (str): column to write as `v_call`, or None.
      jgene_col (str): column to write as `j_call`, or None.
    """
    out = pandas.DataFrame({
        'repertoire_id': repertoire_id,
        'sequence_id': frame[id_col].astype(str),
        'junction_aa': frame[seq_col],
        'duplicate_count': 1,
    })
    if vgene_col is not None or jgene_col is not None:
        # compairr requires v_call *and* j_call unless -g is given; a shared
        # 'na' constant on whichever side has no column keeps that side
        # always agreeing, so only the column(s) actually given constrain
        # hits (vgene_col alone, jgene_col alone, or both).
        out['v_call'] = frame[vgene_col] if vgene_col is not None else 'na'
        out['j_call'] = frame[jgene_col] if jgene_col is not None else 'na'
    out.to_csv(path, sep='\t', index=False)


def _runCompairrPairs(query_path, reference_path, use_genes, compairr_bin,
                      threads=1):
    """
    Run compairr in exact (`-d 0`), cross-comparison (`-x`) mode.

    Arguments:
      query_path (Path): the query side, written by _writeAirrTsv.
      reference_path (Path): the reference side, written by _writeAirrTsv.
      use_genes (bool): whether v_call/j_call were written and should gate
        agreement; False passes compairr's `-g` (ignore genes).
      compairr_bin (str): the resolved binary path.
      threads (int): compairr's own `-t`, 1-256.

    Returns:
      pandas.DataFrame: compairr's pairs output (`-p`): sequence_id_1
      (query), sequence_id_2 (reference), distance (always 0 here).
    """
    with tempfile.TemporaryDirectory() as tmp:
        out_path = Path(tmp) / 'existence.tsv'
        pairs_path = Path(tmp) / 'pairs.tsv'
        cmd = [
            compairr_bin, '-x', str(query_path), str(reference_path),
            '-d', '0',
            '-f',  # ignore duplicate_count -- no real abundances are tracked
            '--distance',
            '-o', str(out_path),
            '-p', str(pairs_path),
            '-t', str(threads),
        ]
        if not use_genes:
            cmd.append('-g')
        subprocess.run(cmd, check=True, capture_output=True, text=True)

        return pandas.read_csv(pairs_path, sep='\t')


def _pairsToHits(pairs):
    """
    Fold compairr's pairs output into a hit set per query row.

    Arguments:
      pairs (pandas.DataFrame): _runCompairrPairs's output.

    Returns:
      dict: query_row (int) to set of reference ids (str).
    """
    hits = {}
    for row in pairs.itertuples(index=False):
        query_row = int(row.sequence_id_1)
        hits.setdefault(query_row, set()).add(str(row.sequence_id_2))

    return hits


def _fallbackExactAnnotation(queries_subset, reference_subset, seq_col, id_col,
                             exact_cols):
    """
    Pure-Python exact annotation for the rows compairr cannot handle: an
    'X' (unresolved residue) or blank `seq_col` value on either side.

    Arguments:
      queries_subset (pandas.DataFrame): query rows to check, keeping their
        original position (query_row) as the DataFrame index.
      reference_subset (pandas.DataFrame): reference rows to check against.
      seq_col (str): the column being compared.
      id_col (str): column in `reference_subset` holding each row's
        identifier.
      exact_cols (list): vgene_col/jgene_col to gate on, or None.

    Returns:
      dict: query_row (int) to set of reference ids (str), omitting query
      rows with no hit.
    """
    if len(queries_subset) == 0 or len(reference_subset) == 0:
        return {}

    index = buildAnnotationIndex(reference_subset, [seq_col], id_col, exact_cols)
    values = _rowValues(queries_subset, [seq_col])
    exact_values = (_rowValues(queries_subset, exact_cols) if exact_cols
                    else [()] * len(values))

    hits = {}
    for query_row, row_values, row_exact in zip(queries_subset.index, values,
                                                 exact_values):
        candidates = index.get(_bucketKey(row_values, row_exact), [])
        matched = {str(row_id) for cand_values, row_id in candidates
                  if _rowDistance(row_values, cand_values) == 0}
        if matched:
            hits[int(query_row)] = matched

    return hits


def runExactAnnotation(queries, reference, seq_col, id_col=DEFAULT_ID_COL,
                       vgene_col=None, jgene_col=None, compairr_bin=None,
                       threads=1):
    """
    Annotate each query row against `reference` on exact identity in
    `seq_col`, optionally gated by exact `vgene_col`/`jgene_col` agreement
    (either, both, or neither) -- CompAIRR accelerated: hash-based lookup
    instead of a per-bucket scan, built for datasets far larger than a
    pure-Python bucketed comparison handles comfortably.

    'X' (e.g. an unresolved PDB-derived residue) and blank values have no
    meaning to compairr, which errors on them outright, so any query or
    reference row containing either is routed around compairr entirely and
    resolved with the same wildcard-aware pure-Python comparison
    runFuzzyAnnotation uses (module docstring): once for such query rows
    against the *whole* reference set, and once for clean query rows against
    the reference's own X/blank rows (compairr's bulk path below already
    covers clean-vs-clean and can't see X-containing reference rows at all).
    Results from all three sources are merged per query row.

    Arguments:
      queries (pandas.DataFrame): the input rows to annotate.
      reference (pandas.DataFrame): the fixed table being annotated against.
      seq_col (str): the single column to compare.
      id_col (str): column in `reference` holding each row's identifier.
      vgene_col (str): column requiring exact agreement, e.g. a V-gene call.
      jgene_col (str): column requiring exact agreement, e.g. a J-gene call.
      compairr_bin (str): path to the compairr binary, or None to search
        PATH.
      threads (int): compairr's own `-t` (1-256); only applies to the
        clean-vs-clean bulk path, the Python fallback is single-threaded
        regardless.

    Returns:
      pandas.DataFrame: one row per query row, columns `query_row`,
      `n_hits_total`, `is_unique_hit`, `hit_ids` (semicolon joined, sorted).

    Raises:
      AnnotationError: if the compairr binary cannot be found.
    """
    compairr_bin = _findCompairr(compairr_bin)
    exact_cols = [c for c in (vgene_col, jgene_col) if c is not None] or None
    use_genes = vgene_col is not None or jgene_col is not None

    queries = queries.reset_index(drop=True)
    queries = queries.assign(**{seq_col: queries[seq_col].fillna('').astype(str)})
    reference = reference.assign(**{seq_col: reference[seq_col].fillna('').astype(str)})
    n = len(queries)

    query_has_x = (queries[seq_col].str.contains(WILDCARD, regex=False)
                  | (queries[seq_col] == ''))
    reference_has_x = (reference[seq_col].str.contains(WILDCARD, regex=False)
                       | (reference[seq_col] == ''))

    clean_queries = queries[~query_has_x]
    x_queries = queries[query_has_x]
    clean_reference = reference[~reference_has_x]
    x_reference = reference[reference_has_x]

    hits = {}

    def _merge(new):
        for query_row, ids in new.items():
            hits.setdefault(query_row, set()).update(ids)

    # Clean query rows vs clean reference rows: bulk path, via compairr.
    if len(clean_queries) and len(clean_reference):
        with tempfile.TemporaryDirectory() as tmp:
            query_path = Path(tmp) / 'query.tsv'
            reference_path = Path(tmp) / 'reference.tsv'
            _writeAirrTsv(clean_queries.assign(__query_row__=clean_queries.index),
                         '__query_row__', seq_col, 'query', query_path,
                         vgene_col, jgene_col)
            _writeAirrTsv(clean_reference, id_col, seq_col, 'reference',
                         reference_path, vgene_col, jgene_col)
            pairs = _runCompairrPairs(query_path, reference_path, use_genes,
                                      compairr_bin, threads)
        _merge(_pairsToHits(pairs))

    # X-containing or blank query rows vs the *whole* reference set.
    _merge(_fallbackExactAnnotation(x_queries, reference, seq_col, id_col, exact_cols))

    # Clean query rows vs X-containing/blank reference rows only.
    _merge(_fallbackExactAnnotation(clean_queries, x_reference, seq_col, id_col,
                                    exact_cols))

    rows = []
    for query_row in range(n):
        hit_ids = sorted(hits.get(query_row, set()))
        rows.append({
            'query_row': query_row,
            'n_hits_total': len(hit_ids),
            'is_unique_hit': len(hit_ids) <= 1,
            'hit_ids': ';'.join(hit_ids),
        })

    return pandas.DataFrame(rows)


# ── paired heavy+light exact annotation ─────────────────────────────────────

def _hitsByCellId(query_chain, chain_result, cell_id_col):
    """
    Fold one chain's runExactAnnotation output into a hit set per antibody.

    Arguments:
      query_chain (pandas.DataFrame): the query rows runExactAnnotation was
        called on for this chain, in the same order (its index reset to
        0..n-1, matching `chain_result`'s `query_row`).
      chain_result (pandas.DataFrame): runExactAnnotation's output for
        `query_chain`.
      cell_id_col (str): column in `query_chain` identifying the antibody a
        chain row belongs to.

    Returns:
      dict: query antibody id to the set of reference antibody ids (as read
      off `id_col` when runExactAnnotation was called -- expected to be
      `cell_id_col` on the reference side too) whose corresponding chain
      agreed with it.
    """
    cell_ids = query_chain[cell_id_col].reset_index(drop=True)
    hits = {}
    for row in chain_result.itertuples(index=False):
        if not row.hit_ids:
            continue
        cell_id = cell_ids.iloc[row.query_row]
        hits.setdefault(cell_id, set()).update(row.hit_ids.split(';'))

    return hits


def runPairedExactAnnotation(queries, reference, cdr3_col='cdr3_aa',
                             vgene_col='v_call', jgene_col='j_call',
                             compairr_bin=None, threads=1):
    """
    Classify each query antibody by whether the *same* reference antibody
    agrees with it on the heavy chain, the light chain, both, or neither.

    An antibody is two rows of an AIRR table -- a heavy chain and, where
    present, a light chain -- sharing one `cell_id` and told apart by
    `locus` (`IGH` for heavy, `IGK`/`IGL` for light). Agreement on a chain
    is runExactAnnotation's criterion (exact V-gene, exact J-gene, exact
    CDR3), run once for the heavy rows and once for the light rows; a
    reference antibody counts toward "both" only if its `cell_id` appears in
    *both* chains' hit sets for the query antibody, not merely in each
    independently -- a query whose heavy chain happens to match one
    reference antibody and whose light chain happens to match an unrelated
    one is `conflicting`, not `both`.

    Arguments:
      queries (pandas.DataFrame): query AIRR rows, both chains, with `locus`
        and `cell_id` columns.
      reference (pandas.DataFrame): reference AIRR rows, same shape.
      cdr3_col (str): the CDR3 amino acid column.
      vgene_col (str): the V-gene call column.
      jgene_col (str): the J-gene call column.
      compairr_bin (str): path to the compairr binary, or None to search
        PATH.
      threads (int): compairr's own `-t`; see runExactAnnotation.

    Returns:
      pandas.DataFrame: one row per query `cell_id`: `classification` (one
      of `both`, `heavy_only`, `light_only`, `conflicting`, `neither`),
      `heavy_hit_ids`, `light_hit_ids` (every reference antibody agreeing on
      that chain alone, semicolon joined, sorted) and `both_hit_ids` (their
      intersection -- non-empty only when `classification` is `both`).
    """
    query_heavy = queries[queries[LOCUS_COL] == HEAVY_LOCUS]
    query_light = queries[queries[LOCUS_COL].isin(LIGHT_LOCI)]
    reference_heavy = reference[reference[LOCUS_COL] == HEAVY_LOCUS]
    reference_light = reference[reference[LOCUS_COL].isin(LIGHT_LOCI)]

    query_heavy = query_heavy.reset_index(drop=True)
    query_light = query_light.reset_index(drop=True)
    reference_heavy = reference_heavy.reset_index(drop=True)
    reference_light = reference_light.reset_index(drop=True)

    heavy_result = runExactAnnotation(query_heavy, reference_heavy, cdr3_col,
                                      id_col=CELL_ID_COL, vgene_col=vgene_col,
                                      jgene_col=jgene_col,
                                      compairr_bin=compairr_bin, threads=threads)
    light_result = runExactAnnotation(query_light, reference_light, cdr3_col,
                                      id_col=CELL_ID_COL, vgene_col=vgene_col,
                                      jgene_col=jgene_col,
                                      compairr_bin=compairr_bin, threads=threads)

    heavy_hits = _hitsByCellId(query_heavy, heavy_result, CELL_ID_COL)
    light_hits = _hitsByCellId(query_light, light_result, CELL_ID_COL)

    rows = []
    for cell_id in sorted(set(queries[CELL_ID_COL])):
        heavy = heavy_hits.get(cell_id, set())
        light = light_hits.get(cell_id, set())
        both = heavy & light

        if both:
            classification = 'both'
        elif heavy and light:
            classification = 'conflicting'
        elif heavy:
            classification = 'heavy_only'
        elif light:
            classification = 'light_only'
        else:
            classification = 'neither'

        rows.append({
            'cell_id': cell_id,
            'classification': classification,
            'heavy_hit_ids': ';'.join(sorted(heavy)),
            'light_hit_ids': ';'.join(sorted(light)),
            'both_hit_ids': ';'.join(sorted(both)),
        })

    return pandas.DataFrame(rows)


# ── indel-tolerant annotation (Levenshtein) ─────────────────────────────────
#
# Unlike runFuzzyAnnotation, there is no length-bucketed index and no
# 'X'-wildcard handling here: Levenshtein allows a query and a candidate to
# differ in length, so a per-column length key can't prefilter candidates, and
# a wildcard has no clean meaning once an alignment can shift positions.
# Candidates are only pre-filtered by `exact_cols` (e.g. V-gene, J-gene), then
# compared to every query brute-force within that block. Also unlike
# runFuzzyAnnotation, this only ever takes a single seq_col: edit distance
# summed across multiple columns independently aligned doesn't correspond to
# any one alignment of the concatenated sequence, so it does not generalize
# the way Hamming distance pools cleanly across columns.

def runLevenshteinAnnotation(queries, reference, seq_col, max_edits,
                             id_col=DEFAULT_ID_COL, exact_cols=None):
    """
    Annotate each query row against `reference` allowing up to `max_edits`
    substitutions, insertions and deletions combined in `seq_col` -- the
    indel-tolerant counterpart to runFuzzyAnnotation's substitution-only
    criterion, via rapidfuzz's Levenshtein distance.

    Arguments:
      queries (pandas.DataFrame): the input rows to annotate.
      reference (pandas.DataFrame): the fixed table being annotated against.
      seq_col (str): the single column to compare.
      max_edits (int): the edit-distance tolerance (substitutions +
        insertions + deletions combined).
      id_col (str): column in `reference` holding each row's identifier.
      exact_cols (list): column names that must agree exactly first, e.g.
        `['v_call', 'j_call']`. None for no such gate.

    Returns:
      pandas.DataFrame: one row per query row, columns `query_row`,
      `n_hits_total`, `min_mismatches` (edit distance here, not substitution
      count), `is_unique_hit`, `hit_ids`.
    """
    reference_exact = (_rowValues(reference, exact_cols) if exact_cols
                       else [()] * len(reference))
    reference_seqs = reference[seq_col].fillna('').astype(str).tolist()
    reference_ids = reference[id_col].tolist()

    buckets = {}
    for seq, exact_values, row_id in zip(reference_seqs, reference_exact,
                                         reference_ids):
        buckets.setdefault(exact_values, []).append((seq, row_id))

    query_exact = (_rowValues(queries, exact_cols) if exact_cols
                  else [()] * len(queries))
    query_seqs = queries[seq_col].fillna('').astype(str).tolist()

    rows = []
    for query_row, (query_seq, row_exact) in enumerate(zip(query_seqs, query_exact)):
        best_by_id = {}
        for cand_seq, row_id in buckets.get(row_exact, []):
            distance = Levenshtein.distance(query_seq, cand_seq,
                                            score_cutoff=max_edits)
            if distance <= max_edits:
                if row_id not in best_by_id or distance < best_by_id[row_id]:
                    best_by_id[row_id] = distance

        hit_ids = sorted(best_by_id, key=lambda x: (best_by_id[x], str(x)))
        rows.append({
            'query_row': query_row,
            'n_hits_total': len(hit_ids),
            'min_mismatches': float(min(best_by_id.values()))
                              if best_by_id else float('nan'),
            'is_unique_hit': len(hit_ids) <= 1,
            'hit_ids': ';'.join(str(x) for x in hit_ids),
        })

    return pandas.DataFrame(rows)
