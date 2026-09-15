.. _UsageAirrcImgt:

sourcerer airrc-imgt
================================================================================

The AIRR-C germline sets blended with IMGT: immunoglobulin V, D and J from
OGRDB, and everything OGRDB does not cover -- all of the T-cell receptor, and
the immunoglobulin constants without a published set -- from IMGT. Offers
``human`` and ``mouse`` collections. ``download`` writes an airrflow
``reference_base`` mixing ``airrc_`` and ``imgt_`` files; ``--igblast``
additionally builds the IgBLAST databases.

.. autoprogram:: sourcerer.Cli:getArgParser()
   :prog: sourcerer
   :start_command: airrc-imgt
