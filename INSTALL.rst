Installation
================================================================================

``sourcerer`` is not yet published on PyPI (the name is already taken by an
unrelated package). Until that is resolved, install directly from GitHub or
from a local checkout.

The current development build can be installed using pip and git::

    > pip3 install git+https://github.com/immcantation/sourcerer@master --user

If you currently have a development version installed, then you will likely
need to add the arguments ``--upgrade --no-deps --force-reinstall`` to the
pip3 command.

To install from a local checkout instead::

    > git clone https://github.com/immcantation/sourcerer
    > cd sourcerer
    > pip3 install . --user

For development, install in editable mode with the ``dev`` extra, which adds
``ruff``::

    > pip3 install -e ".[dev]"

Requirements
--------------------------------------------------------------------------------

+  `Python 3.11 <https://python.org>`__
+  `requests 2.28 <https://requests.readthedocs.io>`__
+  `beautifulsoup4 4.11 <https://www.crummy.com/software/BeautifulSoup>`__
+  `PyYAML 6.0 <https://pyyaml.org>`__
+  `pandas 2.2.3 <https://pandas.pydata.org>`__
+  `airr 2.0 <https://airr-standards.readthedocs.io>`__
+  `tqdm 4.64 <https://tqdm.github.io>`__

All of the above are installed automatically by pip; there is nothing to
install by hand.

Optional
--------------------------------------------------------------------------------

``sourcerer`` itself has no dependency on Nextflow or Docker. They are only
needed to run the ``nf-core/airrflow`` pipeline on the samplesheets
``sourcerer`` writes:

+  `Nextflow <https://www.nextflow.io>`__
+  `Docker <https://www.docker.com>`__ or another Nextflow-supported container
   engine
