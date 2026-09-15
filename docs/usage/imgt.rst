.. _UsageImgt:

sourcerer imgt
================================================================================

`IMGT/GENE-DB <https://www.imgt.org/genedb/>`__: germline V, D, J and C
reference sequences. Offers ``human`` and ``mouse`` collections, each narrowed
with the ``--locus`` and ``--segment`` filters below. ``download`` writes an
airrflow ``reference_base``; ``--igblast`` additionally builds the IgBLAST
databases, which needs ``makeblastdb`` on the path.

Downloading ``all`` instead of a single species fetches every species
**sourcerer supports** into one ``reference_base``, described by one
``IMGT.yaml`` and one ``AIRRC.yaml``. That is not every species the source
publishes: sourcerer covers human and mouse, while OGRDB also carries rhesus
macaque, deer mouse and rainbow trout, and IMGT many more.

.. autoprogram:: sourcerer.Cli:getArgParser()
   :prog: sourcerer
   :start_command: imgt
