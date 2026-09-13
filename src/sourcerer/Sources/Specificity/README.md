# Specificity databases (work in progress)

Specificity databases are reference tables  (epitopes, assays, receptor
sequences) rather than repertoire sequencing runs, so they are grouped under
their own `sourcerer specificity` subcommand instead of being one more
top-level source. Each registered database exposes `search` and `download`,
each taking a table as a further subcommand, plus a synthetic `all` table
that acts on every table the database offers.

## IEDB

[IEDB](https://www.iedb.org/) (`src/sourcerer/Sources/Specificity/Iedb.py`)
offers four tables:

| table            | contents                                                       | source                           |
| ---------------- | -------------------------------------------------------------- | -------------------------------- |
| `bcr`          | BCR (antibody) receptor sequences                              | shared bulk ZIP export           |
| `tcr`          | TCR receptor sequences, kept for reference                     | shared bulk ZIP export           |
| `bcell`        | B-cell assay records: antigen, epitope, qualitative outcome    | PostgREST API (`bcell_search`) |
| `bcr_to_bcell` | join table linking BCR receptor groups to B-cell assay records | PostgREST API (`bcr_to_bcell`) |

`bcr` and `tcr` download as a real, schema-validated AIRR rearrangement TSV
(plus a `.validation.txt` report next to it), not a dump of IEDB's own column
names. The bulk export writes one file row per receptor, both chains side by
side under a two-row header (a category row over a field name row); each row
becomes one rearrangement record per chain it actually carries (one for a
single-chain submission, two for a paired one), linked by `cell_id`. Curated
values (expert reviewed) are preferred over calculated ones where IEDB
records both. `bcell` and `bcr_to_bcell` are assay results and a join table,
not sequences, so they keep their own column names as a plain TSV.

```bash
# download one table
sourcerer specificity iedb download bcr --outdir tmp

# narrow a table down with its filter flags (see --help for what a table supports)
sourcerer specificity iedb download bcell --qualitative-measure Positive --outdir tmp

# drop columns the AIRR schema does not define, for bcr/tcr
sourcerer specificity iedb download bcr --outdir tmp --strict-airr

# download every iedb table
sourcerer specificity iedb download all --outdir tmp
```

## CoV-AbDab

[CoV-AbDab](https://opig.stats.ox.ac.uk/webapps/covabdab/)
(`src/sourcerer/Sources/Specificity/CovAbDab.py`) publishes its whole
collection as a single bulk CSV, so it offers one table:

| table         | contents                                                          | source          |
| ------------- | ------------------------------------------------------------------ | --------------- |
| `antibodies` | antibody and nanobody entries against SARS-CoV-2 and related coronaviruses | bulk CSV export |

The CSV's filename and download URL change with every release, so `search`
and `download` both read the current link off the CoV-AbDab homepage before
resolving the table.

```bash
sourcerer specificity covabdab download antibodies --outdir tmp
sourcerer specificity covabdab search antibodies -o hits.tsv
```

## ASD

[ASD](https://naturalantibody.com/agab/) (Antigen-Specific antibody
Database, `src/sourcerer/Sources/Specificity/Asd.py`) publishes its whole
collection as a Delta Lake (Parquet partitions plus a `_delta_log`
transaction log), staged in a Google Drive folder rather than served over
plain HTTP. Downloading it goes through `gdown` instead of the shared HTTP
client. It offers one table:

| table | contents                                                                    | source            |
| ----- | ---------------------------------------------------------------------------- | ----------------- |
| `asd` | Ab-Ag pairs (structures, affinity assays, literature and patent mined pairs) | Google Drive folder (Delta Lake) |

`asd` can be narrowed with `--dataset`, whose known values (e.g. `hiv`,
`covid-19`, `structures-antibodies`) are recorded in `schema.yaml`; the whole
Delta Lake is still fetched, since Drive offers no server-side per-dataset
export, and the filter is applied while normalizing the downloaded records.

```bash
sourcerer specificity asd download asd --outdir tmp
sourcerer specificity asd download asd --dataset hiv --outdir tmp
```

## Every database at once

```bash
# list every registered specificity database
sourcerer specificity list

# download every table of every registered specificity database
sourcerer specificity all download --outdir tmp

# see what would be fetched without downloading anything
sourcerer specificity iedb download bcr --outdir tmp --dry-run
```

Reference command: `sourcerer specificity <db> download <table> --help`.

## Antigen specificity annotation

Once a specificity database's receptor table has been downloaded as AIRR
(`sourcerer specificity iedb download bcr --format airr`), `sourcerer
annotate` (`Annotate.py`, `AnnotateCli.py`) is the third pipeline step:
annotate a user's own AIRR rearrangement table (`--query`) against that
reference table (`--reference`), one hit list per query row. It works on any
two AIRR tables, not only IEDB's, but IEDB is the only source this has been
wired up against end to end so far.

Two ways to decide a hit, picked automatically from `--max-mismatches` alone:

- **`--max-mismatches 0` (exact, the default).** A single `--compare-cols`
  column, byte-identical, [CompAIRR](https://github.com/uio-bmi/compairr)
  accelerated (hash-based lookup, not a per-bucket scan) -- optionally gated
  by exact `--vgene-col`/`--jgene-col` agreement. **Requires the compairr
  binary**, built from source and put on PATH or pointed at with
  `--compairr-bin`:

  ```bash
  git clone https://github.com/uio-bmi/compairr.git
  cd compairr && make
  ```

  compairr has no concept of the `X` wildcard runFuzzyAnnotation gives
  unresolved (e.g. PDB-derived) residues and errors on an empty sequence
  value outright, so any query or reference row containing either is routed
  around it and resolved with the same wildcard-aware pure-Python
  comparison instead -- nothing is silently dropped, it just does not get
  compairr's speed.

- **`--max-mismatches` above 0 (fuzzy).** One or more `--compare-cols`
  columns, pooled if more than one, up to that many substitutions, pure
  Python -- optionally gated by exact `--exact-cols` agreement (unrelated to
  `--vgene-col`/`--jgene-col`, which only apply to the exact, compairr-backed
  path). `--fraction` treats `--max-mismatches` as a fraction of the query's
  length instead of a fixed count.

Indel-tolerant (Levenshtein) and substitution-matrix-weighted (BLOSUM)
annotation are not included yet.

```bash
# exact annotation on CDR3, compairr-backed
sourcerer annotate --query mine.tsv --reference iedb_bcr.tsv \
    --compare-cols cdr3_aa -o annotated.tsv

# same, gated to exact V-gene agreement, using 8 compairr threads
sourcerer annotate --query mine.tsv --reference iedb_bcr.tsv \
    --compare-cols cdr3_aa --vgene-col v_call --threads 8 -o annotated.tsv

# up to 2 substitutions, gated to the same V-gene, only the heavy chain rows
sourcerer annotate --query mine.tsv --reference iedb_bcr.tsv \
    --compare-cols cdr3_aa --max-mismatches 2 --exact-cols v_call \
    --locus IGH -o annotated.tsv

# mismatch budget as a fraction of the query's length instead of a fixed count
sourcerer annotate --query mine.tsv --reference iedb_bcr.tsv \
    --compare-cols cdr3_aa --max-mismatches 0.15 --fraction -o annotated.tsv
```
