.. _UsageSchema:

immvoke schema
================================================================================

Inspects the schema snapshot checked into the package, or re-harvests it
from a remote source when the upstream API changes (new organisms, fields,
values, ...). See :ref:`API` for what a snapshot contains.

.. autoprogram:: immvoke.Cli:getArgParser()
   :prog: immvoke
   :start_command: schema
