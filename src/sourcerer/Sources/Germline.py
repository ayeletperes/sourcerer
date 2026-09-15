"""
Germline reference sources

A germline reference is a set of allele sequences. IMGT and OGRDB provide them,
and ReferenceSource turns a download into the airrflow reference tree through
buildReference, reusing SourceBase's download machinery: searchUnits, fetchUnit,
the schema and the HTTP client. The output shaping itself lives in
sourcerer.Reference.

SourceBase also declares readUnit and normalizeChunk, the conversion step for
sources that emit AIRR records; a germline source has nothing to convert, so
those two are filled in here with a clear error rather than left for each source
to repeat.
"""

# Info
__author__ = 'Ayelet Peres'

# Sourcerer imports
from sourcerer.Exceptions import SourcererError
from sourcerer.Reference import referenceFastaPath, writeFastaText
from sourcerer.Sources.Base import SourceBase


class ReferenceSource(SourceBase):
    """
    Base class for germline reference sources such as IMGT and OGRDB.
    """

    output = 'reference'

    #: The tag written into reference_base filenames, e.g. 'imgt' gives
    #: imgt_human_IGHV.fasta. Subclasses set this.
    prefix = None

    def readUnit(self, path, unit):
        """Unused: a germline source has nothing to convert."""
        raise SourcererError('%s builds a germline reference; there is nothing '
                             'to convert' % self.name)

    def normalizeChunk(self, metadata, chunk, unit, offset, report):
        """Unused: a germline source has nothing to convert."""
        raise SourcererError('%s builds a germline reference; there is nothing '
                             'to convert' % self.name)

    def buildReference(self, entries, reference_dir):
        """
        Turn downloaded germline files into an airrflow reference_base.

        Arguments:
          entries (list): (DataUnit, Path) pairs, one per fetched file.
          reference_dir (Path): the reference_base root to write into.

        Returns:
          ReferenceReport: the files written.
        """
        raise NotImplementedError

    def writeChain(self, reference_dir, species, kind, chain, records):
        """
        Write one chain's FASTA into the reference tree.

        Arguments:
          reference_dir (Path): the reference_base root.
          species (str): the species.
          kind (str): the subdirectory (KIND_VDJ, KIND_CONSTANT or KIND_AA).
          chain (str): the chain, e.g. 'IGHV'.
          records (iterable): (header, sequence) tuples to write verbatim.

        Returns:
          Path: the file written.
        """
        path = referenceFastaPath(reference_dir, self.prefix, species, kind, chain)
        writeFastaText(path, records)

        return path
