.. _Usage:

Commandline Usage
================================================================================

``sourcerer`` is a single command with a subcommand tree: one subcommand per
external source (``oas``, ``imgt``, ``ogrdb`` and ``airrc-imgt`` today), a
``schema`` subcommand
to inspect and re-harvest the stored snapshot each source is built from, a
``reference`` subcommand to validate a germline reference folder and build its
IgBLAST databases, and a ``sources`` subcommand that lists what is registered.

Sources come in two kinds. A repertoire source such as ``oas`` downloads
sequencing data and writes an airrflow samplesheet; a germline reference source
such as ``imgt`` and ``ogrdb`` downloads germline sets and writes an airrflow
``reference_base``, optionally building the IgBLAST databases with ``--igblast``.
The ``reference`` subcommand does that same build for a reference folder supplied
by hand, in either the ``reference_base`` layout or a flat folder of FASTAs.

Every source subcommand exposes the same two actions, ``search`` and
``download``, each taking a collection (for example ``paired`` and ``unpaired``
for OAS, or a species for the germline sources) as a further subcommand. The
filter flags under a collection — ``--species``, ``--locus``, and so on — are
not hardcoded: they are generated at parser-construction time from the checked-in
schema snapshot described in :ref:`API`, which is also why they appear below
exactly as they would in ``--help`` on the machine building these docs.

Commands are documented one page per top level subcommand, mirroring how the
commandline itself groups them: :doc:`sources` and :doc:`schema` apply
across every source, and each remaining page below documents one source
module's own ``search``/``download`` tree. Adding a new source means adding
one page here alongside it.

.. toctree::
   :maxdepth: 2

   sources
   schema
   reference
   oas
   imgt
   ogrdb
   airrc-imgt
