.. _UsageImgt:

sourcerer imgt
================================================================================

`IMGT/GENE-DB <https://www.imgt.org/genedb/>`__: germline V, D, J and C
reference sequences. Offers ``human`` and ``mouse`` collections, each narrowed
with the ``--locus`` and ``--segment`` filters below. ``download`` writes an
airrflow ``reference_base``; ``--igblast`` additionally builds the IgBLAST
databases, which needs ``makeblastdb`` on the path.

.. autoprogram:: sourcerer.Cli:getArgParser()
   :prog: sourcerer
   :start_command: imgt
