Release Notes
================================================================================

Version 0.1.0: 2026.08.04
-------------------------------------------------------------------------------

Initial release.

General:

+ Added the ``sourcerer`` commandline tool, with one module per external
  immune-repertoire database and OAS_ (Observed Antibody Space) as the first
  source.
+ Added ``sourcerer oas search`` and ``sourcerer oas download``, for both the
  ``paired`` and ``unpaired`` OAS collections, with commandline filter flags
  generated from a checked-in snapshot of the OAS search form rather than
  hardcoded.
+ Added conversion of downloaded data to an AIRR rearrangement TSV, with a
  streaming validation report, and to FASTA.
+ Added ``nf-core/airrflow`` samplesheet generation, one samplesheet per
  converted format, that merges across repeated download runs into the same
  output directory instead of overwriting.
+ Added ``sourcerer schema show`` and ``sourcerer schema refresh``, to inspect
  and re-harvest the stored snapshot of a source's search fields and data unit
  catalog.
+ Added a download provenance record (what was fetched, from where, when, and
  its hash) written alongside every download.

Germline references:

+ Added germline reference sources IMGT_ (``sourcerer imgt``) and OGRDB_
  (``sourcerer ogrdb``, also reachable as ``sourcerer airrc``), and an
  ``airrc-imgt`` blend that takes immunoglobulin V, D and J from OGRDB's AIRR-C
  sets and the T-cell receptor and remaining constants from IMGT.
+ Added ``sourcerer <source> download <species>``, which writes the germline
  ``reference_base`` in the `nf-core/airrflow`_ layout, and ``--igblast`` to also
  build the IgBLAST databases (``makeblastdb`` plus the NCBI internal_data and
  optional_file trees).
+ Added ``sourcerer reference build``, to validate a germline reference folder
  -- in the ``reference_base`` layout or a flat folder of FASTAs -- and build
  its IgBLAST databases, with ``--check`` to validate without building.
+ nf-core/airrflow can fetch germlines through ``sourcerer`` for its ``imgt``
  and ``airrc-imgt`` database types, or consume a ``sourcerer``-built reference
  passed to ``--reference_fasta`` / ``--reference_igblast``.

.. _OAS: https://opig.stats.ox.ac.uk/webapps/oas/
.. _IMGT: https://www.imgt.org/genedb/
.. _OGRDB: https://ogrdb.airr-community.org/
.. _nf-core/airrflow: https://nf-co.re/airrflow
