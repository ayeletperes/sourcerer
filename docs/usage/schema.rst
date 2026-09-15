.. _UsageSchema:

sourcerer schema
================================================================================

Inspects the schema snapshot checked into the package, or re-harvests it
from a remote source when the upstream API changes (new organisms, fields,
values, ...). See :ref:`API` for what a snapshot contains.

.. autoprogram:: sourcerer.Cli:getArgParser()
   :prog: sourcerer
   :start_command: schema
