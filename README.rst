sourcerer
================================================================================

.. image:: https://img.shields.io/badge/lifecycle-experimental-orange.svg
    :target: https://www.tidyverse.org/lifecycle/#experimental
    :alt: Lifecycle: experimental

**Experimental.** This is early-stage, unreleased software. Breaking changes --
including changes to the package name -- can and will happen without notice.
Do not depend on it for production pipelines yet.

``sourcerer`` downloads data from online immune repertoire databases and formats
it for use with the Immcantation_ framework and `nf-core/airrflow`_. Each
external source is a module, and sources come in two kinds:

- *dataset* sources such as OAS_ (Observed Antibody Space) download sequencing
  data and write an airrflow samplesheet;
- *germline reference* sources -- IMGT_, OGRDB_, and an ``airrc-imgt`` blend of
  the two -- download germline sets and build the ``reference_base`` and IgBLAST
  databases airrflow consumes. ``sourcerer reference`` can also validate and
  build those databases from a reference folder you already have.

.. _Immcantation: https://immcantation.readthedocs.io
.. _nf-core/airrflow: https://nf-co.re/airrflow
.. _OAS: https://opig.stats.ox.ac.uk/webapps/oas/
.. _IMGT: https://www.imgt.org/genedb/
.. _OGRDB: https://ogrdb.airr-community.org/

Why
--------------------------------------------------------------------------------

Remote databases change. They add organisms, diseases and studies, they rename
search fields, and occasionally they change the format of the files themselves.
Code that hardcodes those details breaks quietly, producing empty result sets or
mislabelled columns rather than errors.

``sourcerer`` keeps a checked-in snapshot of each remote source's search schema
*and* of its downloaded file format. The command line is generated from that
snapshot, so no field list is hardcoded, and a scheduled job re-harvests the
remote schema and opens a pull request whenever it drifts. Breakage shows up as
a reviewable diff and a failing test, not as a silently wrong download.

Usage
--------------------------------------------------------------------------------

Datasets (OAS), producing an airrflow samplesheet:

.. code-block:: bash

    sourcerer --version
    sourcerer oas download paired --species human --limit 1 --outdir tmp    # To convert to fasta, rerun (hashes any file already on disk)
    sourcerer oas download paired --species human --limit 3 --outdir tmp --format fasta
    cd tmp
    nextflow run nf-core/airrflow -r 5.1.0 \
        -profile docker \
        --mode assembled \
        --input samplesheet_airrflow_fasta.tsv \
        --outdir airrflow_out -c ../airrflow.config \
        --clonal_threshold 0.2 -resume

Germline references, producing the ``reference_base`` and IgBLAST databases:

.. code-block:: bash

    # IMGT germline for a species, and (with --igblast) the IgBLAST databases
    sourcerer imgt download human --outdir ref --igblast

    # the AIRR-C sets blended with IMGT (immunoglobulin from OGRDB, TR and the
    # remaining constants from IMGT) -- the airrflow airrc-imgt reference
    sourcerer airrc-imgt download human --outdir ref --igblast

    # validate a germline folder someone provided and build its databases;
    # --check validates only, without makeblastdb
    sourcerer reference build ref/reference_base --out igblast_base --check

nf-core/airrflow uses the result either way: point ``--reference_fasta`` and
``--reference_igblast`` at ``ref/reference_base`` and ``igblast_base`` with
``--fetch_germlines none``.


License
--------------------------------------------------------------------------------

This work is licensed under the `GNU Affero General Public License 3 (AGPL-3)
<https://www.gnu.org/licenses/agpl-3.0.en.html>`_.
