Contact
--------------------------------------------------------------------------------

If you need help or have any questions, please contact the `Immcantation Group <mailto:immcantation@googlegroups.com>`__.

If you have discovered a bug or have a feature request, you can open an issue using the
`issue tracker <http://github.com/immcantation/sourcerer/issues>`__.

To receive alerts about Immcantation releases, news, events, and tutorials, join the `Immcantation News <https://groups.google.com/g/immcantation-news>`__ Google Group. `Membership settings <https://groups.google.com/g/immcantation-news/membership>`__ can be adjusted to change the frequency of email updates.


Citing downloaded data
--------------------------------------------------------------------------------

``sourcerer``'s own license (below) covers the tool, not the data it downloads.
Each remote source distributes its data under its own license and asks to be
cited in its own way; downloading through ``sourcerer`` does not change either
obligation. Run ``sourcerer sources list`` to print the license and citation
for every source, and check ``data_license`` / ``data_citation`` in the
``download_metadata.yml`` written alongside a download for the record tied to
that specific dataset.

OAS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Data from `OAS <https://opig.stats.ox.ac.uk/webapps/oas/>`__ is distributed
under a `CC BY 4.0
<https://creativecommons.org/licenses/by/4.0/>`__ license. In exchange, OAS
asks that both of the following be cited:

- Kovaltsuk A, Leem J, Kelm S, Snowden J, Deane CM, Krawczyk K. Observed
  Antibody Space: A Resource for Data Mining Next-Generation Sequencing of
  Antibody Repertoires. *J Immunol*. 2018;201(8):2502-2509.
  doi:`10.4049/jimmunol.1800708 <https://doi.org/10.4049/jimmunol.1800708>`__
- Olsen TH, Boyles F, Deane CM. Observed Antibody Space: A diverse database of
  cleaned, annotated, and translated unpaired and paired antibody sequences.
  *Protein Sci*. 2022;31(1):141-146.
  doi:`10.1002/pro.4205 <https://doi.org/10.1002/pro.4205>`__

IMGT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Data from `IMGT <https://www.imgt.org/genedb/>`__ is governed by the `IMGT terms
of use <https://www.imgt.org/about/termsofuse.php>`__: it is free for academic
research on condition that IMGT is cited. IMGT asks that the following be cited:

- Lefranc MP, Giudicelli V, Duroux P, et al. IMGT, the international
  ImMunoGeneTics information system 25 years on. *Nucleic Acids Res*.
  2015;43(Database issue):D413-D422.
  doi:`10.1093/nar/gku1056 <https://doi.org/10.1093/nar/gku1056>`__

The ``airrc-imgt`` blend uses IMGT data too, so its downloads carry this
obligation as well as the OGRDB one below.

OGRDB
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Data from `OGRDB <https://ogrdb.airr-community.org/>`__ is distributed under a
`CC BY 4.0 <https://creativecommons.org/licenses/by/4.0/>`__ license. In
exchange, OGRDB asks that the following be cited:

- Lees WD, Busse CE, Corcoran M, et al. OGRDB: a reference database of inferred
  immune receptor genes. *Nucleic Acids Res*. 2020;48(D1):D964-D970.
  doi:`10.1093/nar/gkz822 <https://doi.org/10.1093/nar/gkz822>`__


License
--------------------------------------------------------------------------------

This work is licensed under the
`GNU Affero General Public License Version 3 (AGPL-3) <https://www.gnu.org/licenses/agpl-3.0.en.html>`__.
This covers the ``sourcerer`` codebase; data downloaded through it carries the
source's own license, see `Citing downloaded data`_ above.
