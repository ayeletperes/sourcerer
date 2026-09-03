.. _UsageOas:

sourcerer oas
================================================================================

`Observed Antibody Space <https://opig.stats.ox.ac.uk/webapps/oas/>`__:
cleaned, annotated antibody repertoires. Offers ``paired`` and ``unpaired``
collections, each searchable and downloadable with the filter flags below.
``verify`` cross-references a downloaded samplesheet's unresolved subjects
(OAS's own null sentinels, ``no`` and ``None``) against NCBI, from the run or
sample accession embedded in each row's ``sample_name``.

.. autoprogram:: sourcerer.Cli:getArgParser()
   :prog: sourcerer
   :start_command: oas
