.. _UsageReference:

sourcerer reference
================================================================================

Validate a folder of germline FASTAs and build the IgBLAST databases from it,
for a reference someone supplies rather than one sourcerer downloaded. Files are
recognised by name in any directory layout --
``[<prefix>_][aa_]<species>_<CHAIN>.fasta``, for example ``human_IGHV.fasta`` or
``imgt_human_IGHV.fasta`` -- so a nested ``reference_base`` and a flat folder both
work. ``--check`` validates and reports what would build without building
anything, and needs no ``makeblastdb``.

.. autoprogram:: sourcerer.Cli:getArgParser()
   :prog: sourcerer
   :start_command: reference
