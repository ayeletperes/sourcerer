# Specificity annotations (work in progress)

Here we provide specificity annotation against public databases, currently
IEDB: download the reference table, convert it to AIRR, then match your own
receptors against it.

## IEDB: download and convert

[IEDB](https://www.iedb.org/) (`Iedb.py`) publishes BCR (antibody) receptor
sequences as a single table, `bcr`, inside a shared bulk ZIP export.

`bcr` downloads as a  schema-validated AIRR rearrangement TSV (plus a
`.validation.txt` report next to it).
The bulk export writes one file row per receptor, both chains side by side
under a two-row header (a category row over a field name row); each row
becomes one rearrangement record per chain it actually carries (one for a
single-chain submission, two for a paired one), linked by `cell_id`. Curated
values (expert reviewed) are preferred over calculated ones where IEDB
records both.

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
| `--identity PCT` | same as `--fuzzy`, tolerance expressed as percent identity instead of a count | one or more, pooled |
| `--levenshtein N` | exact V+J gene, up to N substitutions/insertions/deletions | exactly one |
| `--paired` | exact V-gene, J-gene and CDR3 on *both* chains of the same reference antibody | fixed (`cdr3_aa`) |

`--fuzzy`, `--identity` and `--levenshtein` are pure Python, always gate on
exact V-gene and J-gene agreement, and are always restricted to the heavy
chain -- `--locus` may only be left unset or set to `IGH` under any of them.
None of that is optional, since their cost grows with the candidate pool
(unlike `--exact`'s hashed CompAIRR lookup) and limiting to one chain keeps
it bounded. `--levenshtein` is the odd one out among the three: it tolerates
insertions/deletions, not just substitutions, and for that reason only ever
takes a single `--compare-cols` column -- edit distance doesn't pool across
columns the way Hamming distance does.

`--exact` runs through [CompAIRR](https://github.com/uio-bmi/compairr),
which has to be built from source and put on `PATH` (or pointed at with
`--compairr-bin`). CompAIRR has no concept of the `X` wildcard for an
unresolved residue and errors on a blank sequence, so such rows are compared
in Python instead, where `X` matches anything.

`bcr.tsv` is one row per chain, not per antibody: `--exact` matches heavy
and light as two separate calls via `--locus`; the pure-Python modes filter
to `IGH` on their own.

```bash
# VDJ exact match (heavy chain only -- "VDJ" implies the D segment, which
# only the heavy chain rearranges; light chain is V+J alone)
sourcerer annotate --query mine.tsv --reference tmp/specificity/bcr.tsv \
    --compare-cols sequence_aa --locus IGH -o annotated_heavy.tsv

# exact, full light chain (IGK for kappa, IGL for lambda)
sourcerer annotate --query mine.tsv --reference tmp/specificity/bcr.tsv \
    --compare-cols sequence_aa --locus IGK -o annotated_light.tsv

# V+J gene match, CDR3 85% amino acid identity (Hamming distance) -- heavy
# chain only, gene gating and locus restriction both automatic under --fuzzy
sourcerer annotate --query mine.tsv --reference tmp/specificity/bcr.tsv \
    --compare-cols cdr3_aa --identity 85 -o annotated.tsv

# up to 2 substitutions on CDR3 -- --locus IGH is implicit and could be left out
sourcerer annotate --query mine.tsv --reference tmp/specificity/bcr.tsv \
    --compare-cols cdr3_aa --fuzzy 2 -o annotated.tsv

# substitutions pooled across CDR1 and CDR3 together (one shared tolerance,
# not one each)
sourcerer annotate --query mine.tsv --reference tmp/specificity/bcr.tsv \
    --compare-cols cdr1_aa --compare-cols cdr3_aa --fuzzy 3 -o annotated.tsv

# up to 2 edits (substitutions/insertions/deletions) on CDR3 -- tolerates a
# query and reference CDR3 of different lengths, which --fuzzy cannot at all
sourcerer annotate --query mine.tsv --reference tmp/specificity/bcr.tsv \
    --compare-cols cdr3_aa --levenshtein 2 -o annotated.tsv
```

Every query row comes back with `n_hits_total`, `is_unique_hit` and
`hit_ids`. The three fuzzy modes add `min_mismatches` (edit distance for
`--levenshtein`, substitution count otherwise); `--fuzzy`/`--identity` also
add `mismatch_positions`, which `--levenshtein` can't produce since an indel
shifts every position after it. BLOSUM is not implemented.

### Paired heavy+light exact match

`--exact`/`--fuzzy` treat each chain row on its own. `--paired` is stricter:
an antibody counts as a match only if the *same* reference antibody agrees
on V-gene, J-gene and CDR3 for both its heavy and its light chain -- a heavy
hit on one reference antibody and an unrelated light hit on another doesn't
count. `--compare-cols` and `--locus` don't apply here; the columns
(`cdr3_aa`, `v_call`, `j_call`, `cell_id`, `locus`) are fixed by that
definition.

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
