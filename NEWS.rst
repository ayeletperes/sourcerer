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
  optional_file trees). ``download all`` fetches every species the source covers
  into one ``reference_base``, with one ``IMGT.yaml`` and one ``AIRRC.yaml``
  describing the whole of it. The OAS source has no equivalent: its collections
  are paired and unpaired, which have different fields and different outputs.
+ Added ``sourcerer reference build``, to validate a germline reference folder
  -- in the ``reference_base`` layout or a flat folder of FASTAs -- and build
  its IgBLAST databases, with ``--check`` to validate without building. The
  build writes a ``sourcerer_build.yaml`` recording what it produced and carries
  any ``IMGT.yaml`` / ``AIRRC.yaml`` from the source reference into the
  ``igblast_base``, so a custom reference keeps its provenance.
+ Added ``--map``, a manifest naming the species and chain of reference files
  whose own names do not say -- an OGRDB set downloaded as ``IGH_VDJ_V.fasta``,
  for one. It overrides the naming rule rather than guessing, since a chain
  inferred wrongly does not fail, it files the alleles under the wrong locus.
  ``.fa`` and ``.fna`` are read alongside ``.fasta``.
+ Added ``sourcerer reference show``, to report what a reference folder is and
  where it came from: the release and sets behind it, what was built, and what
  it holds, read back from the provenance sidecars.
+ An allele name longer than the 50 characters ``makeblastdb -parse_seqids``
  accepts is now shortened rather than failing the build, keeping the head of
  the name and appending a digest of the whole of it. The mapping back to the
  original is written to ``shortened_alleles.tsv`` in the ``igblast_base``, since
  the shortened name is what IgBLAST reports into a ``v_call``. Only
  VDJbase-style novel allele names reach the limit; nothing IMGT or OGRDB
  publishes is close.
+ A build now reports J alleles that the mirrored NCBI auxiliary file does not
  name, and records them in ``sourcerer_build.yaml``. IgBLAST looks a J germline
  up in that file by name, so an allele missing from it gets no CDR3 and no
  productivity call, silently. This is what OGRDB's mouse sets hit: they name
  their J alleles ``IGKJ0-4JXG*00`` and NCBI's ``mouse_gl.aux`` lists none of
  them. The recommendation is to build an auxiliary file from the reference and
  pass it to ``igblastn -auxiliary_data``; sourcerer reports which alleles need
  rows rather than building the file itself.
+ A download now records its provenance inside ``reference_base``: ``IMGT.yaml``
  with the GENE-DB release, and ``AIRRC.yaml`` with each OGRDB set's version and
  release date, and its Zenodo DOI under ``--resolve-doi``. A download fetches
  one species but a ``reference_base`` holds as many as were downloaded into it,
  so both sidecars merge rather than replace: downloading mouse into a folder
  that already holds human keeps human. The IMGT release is recorded per
  species, since the two can be fetched weeks apart from different builds, and
  ``--from`` pins only the species being downloaded.
+ Added ``download --from <reference>``, to re-download the exact versions a
  ``reference_base`` was built from: OGRDB sets through the versioned API, and an
  IMGT release from the genedb-releases_ archive (the nearest release, with a
  warning, when the exact one is not archived, since IMGT serves only the
  current build). A substituted release is recorded as one: ``IMGT.yaml`` keeps
  the release asked for beside the release used, so a later ``--from`` cannot
  quietly pin the neighbour as though it were the original.
+ Added ``download --compare`` and ``sourcerer reference diff``, to compare two
  reference folders allele by allele -- identical, added, removed or changed --
  so a re-downloaded release can be checked against the original.
+ nf-core/airrflow can fetch germlines through ``sourcerer`` for its ``imgt``
  and ``airrc-imgt`` database types, or consume a ``sourcerer``-built reference
  passed to ``--reference_fasta`` / ``--reference_igblast``.

.. _OAS: https://opig.stats.ox.ac.uk/webapps/oas/
.. _genedb-releases: https://github.com/JamieHeather/genedb-releases
.. _IMGT: https://www.imgt.org/genedb/
.. _OGRDB: https://ogrdb.airr-community.org/
.. _nf-core/airrflow: https://nf-co.re/airrflow
