# Specificity annotations (work in progress)

Here we provide specificity annotation against public databases, currently
IEDB: download the reference table, convert it to AIRR, then match your own
receptors against it.

## IEDB: download and convert

[IEDB](https://www.iedb.org/) (`Iedb.py`) publishes BCR (antibody) receptor
sequences as a single table, `bcr`, inside a shared bulk ZIP export.

`bcr` downloads as a schema-validated AIRR rearrangement TSV, with a
`.validation.txt` report next to it. One row per chain, linked by `cell_id`,
so a paired receptor becomes two rows; where IEDB records both a curated
(expert reviewed) and a calculated value, the curated one wins.

```bash
# download and convert
sourcerer specificity bcr --outdir tmp

# drop columns the AIRR schema does not define
sourcerer specificity bcr --outdir tmp --strict-airr
```

## Annotating against it

`sourcerer annotate` takes a query AIRR table and a reference one -- the
`tmp/specificity/bcr.tsv` written above -- and reports, for every query row,
which reference rows it matches. One flag picks the criterion:

| Flag | Criterion | `--compare-cols` |
| --- | --- | --- |
| `--exact` (default) | byte-identical, CompAIRR-backed, optional V/J-gene gate | one |
| `--fuzzy N` | exact V+J gene, up to N substitutions | one or more, pooled |
| `--identity PCT` | same as `--fuzzy`, tolerance as percent identity rather than a count | one or more, pooled |
| `--levenshtein N` | exact V+J gene, up to N substitutions/insertions/deletions | exactly one |
| `--paired` | exact V-gene, J-gene and CDR3 on *both* chains of the same reference antibody | fixed (`cdr3_aa`) |

The three fuzzy modes are pure Python, always gate on exact V-gene and
J-gene agreement, and are always restricted to the heavy chain (`--locus`
may only be unset or `IGH`) -- their cost grows with the candidate pool,
unlike `--exact`'s hashed CompAIRR lookup. `--levenshtein` takes a single
`--compare-cols` column, since edit distance doesn't pool across columns the
way Hamming distance does.

`--exact` runs through [CompAIRR](https://github.com/uio-bmi/compairr),
built from source and put on `PATH` (or pointed at with `--compairr-bin`).
CompAIRR has no `X` wildcard and errors on a blank sequence, so those rows
are compared in Python instead, where `X` matches anything.

Heavy and light are separate rows, so `--exact` matches them as two calls
via `--locus`; the fuzzy modes filter to `IGH` themselves.

```bash
# VDJ exact match, heavy chain
sourcerer annotate --query mine.tsv --reference tmp/specificity/bcr.tsv \
    --compare-cols sequence_aa --locus IGH -o annotated_heavy.tsv

# same for the light chain (IGK for kappa, IGL for lambda)
sourcerer annotate --query mine.tsv --reference tmp/specificity/bcr.tsv \
    --compare-cols sequence_aa --locus IGK -o annotated_light.tsv

# V+J gene match, CDR3 85% amino acid identity (Hamming)
sourcerer annotate --query mine.tsv --reference tmp/specificity/bcr.tsv \
    --compare-cols cdr3_aa --identity 85 -o annotated.tsv

# up to 2 substitutions on CDR3
sourcerer annotate --query mine.tsv --reference tmp/specificity/bcr.tsv \
    --compare-cols cdr3_aa --fuzzy 2 -o annotated.tsv

# pooled across CDR1 and CDR3 (one shared tolerance, not one each)
sourcerer annotate --query mine.tsv --reference tmp/specificity/bcr.tsv \
    --compare-cols cdr1_aa --compare-cols cdr3_aa --fuzzy 3 -o annotated.tsv

# up to 2 edits on CDR3, tolerating a length difference --fuzzy cannot
sourcerer annotate --query mine.tsv --reference tmp/specificity/bcr.tsv \
    --compare-cols cdr3_aa --levenshtein 2 -o annotated.tsv
```

Every query row comes back with `n_hits_total`, `is_unique_hit` and
`hit_ids`. The fuzzy modes add `min_mismatches` (edit distance for
`--levenshtein`, substitution count otherwise); `--fuzzy`/`--identity` also
add `mismatch_positions`, which `--levenshtein` can't produce since an indel
shifts every position after it. BLOSUM is not implemented.

### Paired heavy+light exact match

`--paired` is stricter than the per-row modes: an antibody matches only if
the *same* reference antibody agrees on V-gene, J-gene and CDR3 for both its
heavy and its light chain -- a heavy hit on one reference antibody and an
unrelated light hit on another doesn't count. `--compare-cols` and `--locus`
don't apply; the columns are fixed by that definition.

```bash
sourcerer annotate --query mine.tsv --reference tmp/specificity/bcr.tsv \
    --paired -o annotated.tsv
```

One row per query antibody (`cell_id`), classified as `both`, `heavy_only`,
`light_only`, `conflicting` (heavy and light each hit, but different
reference antibodies), or `neither`, alongside `heavy_hit_ids`,
`light_hit_ids` and `both_hit_ids` (their intersection).

Reference command: `sourcerer specificity bcr --help`,
`sourcerer annotate --help`.
