.. _UsageReference:

sourcerer reference
================================================================================

Validate a folder of germline FASTAs and build the IgBLAST databases from it,
for a reference someone supplies rather than one sourcerer downloaded. Files are
recognised by name in any directory layout --
``[<prefix>_][aa_]<species>_<CHAIN>.fasta``, for example ``human_IGHV.fasta`` or
``imgt_human_IGHV.fasta`` -- so a nested ``reference_base`` and a flat folder both
work, and ``.fa`` and ``.fna`` are read alongside ``.fasta``. ``--check``
validates and reports what would build without building anything, and needs no
``makeblastdb``.

A reference assembled by hand often does not follow that naming -- an OGRDB set
downloaded as ``IGH_VDJ_V.fasta`` says nothing a filename rule could read.
``--map`` takes a manifest that says it outright, rather than guessing, because
a chain inferred wrongly does not fail: it builds a database with the alleles
filed under the wrong locus. The manifest is tab- or whitespace-separated, with
``#`` comments ignored::

    #file                  species  chain  [aa]
    IGH_VDJ_V.fasta        human    IGHV
    C57BL-6_IGH_V.fasta    mouse    IGHV
    my_translated.fasta    human    IGHV   aa

Files are matched on their path relative to the reference folder, or failing
that on their basename. The manifest overrides the naming rule, so it can also
correct a file the rule misreads, and it applies to ``diff`` and ``show`` too.

An allele name longer than the 50 characters ``makeblastdb -parse_seqids``
accepts is shortened rather than failing the build: the head of the name is kept
and a digest of the whole name replaces the tail. Because the shortened name is
what IgBLAST reports into a ``v_call``, the mapping back is written to
``shortened_alleles.tsv`` beside the databases. In practice only VDJbase-style
novel allele names get near the limit.

A build also reports any J allele the mirrored NCBI auxiliary file does not
name. IgBLAST looks a J germline up in that file by name, so an allele it does
not list gets no CDR3 and no productivity call without any error. The names are
recorded under ``aux_not_covered`` in ``sourcerer_build.yaml``.

Where that happens, build an auxiliary file from the reference itself and pass it
to ``igblastn`` with ``-auxiliary_data``. Sourcerer does not build one -- that
belongs with the pipeline running IgBLAST -- but it names the alleles that need
rows, so the file can be built for exactly those.

``sourcerer reference show`` reports what a folder is and where it came from,
reading the ``IMGT.yaml``, ``AIRRC.yaml`` and ``sourcerer_build.yaml`` sidecars
a download or build leaves behind.

.. autoprogram:: sourcerer.Cli:getArgParser()
   :prog: sourcerer
   :start_command: reference
