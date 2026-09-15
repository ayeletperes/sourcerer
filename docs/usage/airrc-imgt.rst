.. _UsageAirrcImgt:

sourcerer airrc-imgt
================================================================================

The AIRR-C germline sets blended with IMGT: immunoglobulin V, D and J from
OGRDB, and everything OGRDB does not cover -- all of the T-cell receptor, and
the immunoglobulin constants without a published set -- from IMGT. Offers
``human`` and ``mouse`` collections. ``download`` writes an airrflow
``reference_base`` mixing ``airrc_`` and ``imgt_`` files; ``--igblast``
additionally builds the IgBLAST databases.

Downloading ``all`` instead of a single species fetches every species
**sourcerer supports** into one ``reference_base``, described by one
``IMGT.yaml`` and one ``AIRRC.yaml``. That is not every species the source
publishes: sourcerer covers human and mouse, while OGRDB also carries rhesus
macaque, deer mouse and rainbow trout, and IMGT many more.

.. autoprogram:: sourcerer.Cli:getArgParser()
   :prog: sourcerer
   :start_command: airrc-imgt
