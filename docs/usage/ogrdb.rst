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

Downloading ``all`` instead of a single species fetches every species
**sourcerer supports** into one ``reference_base``, described by one
``IMGT.yaml`` and one ``AIRRC.yaml``. That is not every species the source
publishes: sourcerer covers human and mouse, while OGRDB also carries rhesus
macaque, deer mouse and rainbow trout, and IMGT many more.

.. autoprogram:: sourcerer.Cli:getArgParser()
   :prog: sourcerer
   :start_command: ogrdb
