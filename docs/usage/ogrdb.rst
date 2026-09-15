.. _UsageOgrdb:

sourcerer ogrdb
================================================================================

`OGRDB <https://ogrdb.airr-community.org/>`__: AIRR Community curated
immunoglobulin germline sets. Offers ``human`` and ``mouse`` collections,
narrowed with the ``--locus`` filter below. ``download`` writes an airrflow
``reference_base``; ``--igblast`` additionally builds the IgBLAST databases,
which needs ``makeblastdb`` on the path.

OGRDB is the AIRR Community database, so ``ogrdb`` also answers to the alias
``airrc`` (``sourcerer airrc download ...``). For a reference that additionally
fills in the T-cell receptor and the remaining constants from IMGT, use the
``airrc-imgt`` source instead.

.. autoprogram:: sourcerer.Cli:getArgParser()
   :prog: sourcerer
   :start_command: ogrdb
