"""
Output writers

These are source agnostic: they take normalized, AIRR named records and know
nothing about where the data came from. A second database gets them for free.
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import collections
import logging
from pathlib import Path

log = logging.getLogger(__name__)

#: Identifier written when a record has no sequence_id. Only the identifier is
#: positional, so this is the one place a literal is still needed.
FASTA_NULL = 'NA'

#: Annotations carried in the FASTA header, in order, as AIRR column names. The
#: pRESTO convention is a sequence identifier followed by pipe separated
#: key=value pairs, and the keys are kept identical to the AIRR column names so
#: that a header and the rearrangement TSV never disagree about what a field is
#: called.
#:
#: Only cell_id is carried. Anything listed here is re-parsed by Change-O's
#: MakeDb.py and written over the value IgBLAST assigned, so carrying a field the
#: aligner derives for itself, such as c_call or locus, replaces an alignment
#: based call with the source's. cell_id is the only annotation IgBLAST cannot
#: recover from the sequence, and is what makes single cell data usable in
#: airrflow's assembled mode.
FASTA_ANNOTATIONS = ('cell_id',)


def coerceAirrTypes(frame):
    """
    Cast columns to the types the AIRR schema declares for them.

    OAS writes whole numbers in float form, so an alignment coordinate arrives as
    '1.0' where AIRR requires '1'. Every integer field in every row fails
    validation without this, which would bury any real problem under tens of
    thousands of spurious errors.

    Everything is held as text: these files are streamed straight to TSV, and
    round tripping through a numeric dtype would reintroduce the same '1.0' on
    output as well as turning empty cells into NaN.

    Arguments:
      frame (pandas.DataFrame): normalized records, columns as strings.

    Returns:
      pandas.DataFrame: the frame with numeric columns reformatted in place.
    """
    from airr.schema import RearrangementSchema

    for column in frame.columns:
        spec = RearrangementSchema.properties.get(column)
        if spec is None or spec.get('type') != 'integer':
            continue
        frame[column] = frame[column].map(_asIntegerText)

    return frame


def _asIntegerText(value):
    """
    Render a value as an integer string, or empty if it is not a whole number.

    Arguments:
      value: the raw value.

    Returns:
      str: the integer as text, or ''.
    """
    text = '' if value is None else str(value).strip()
    if not text:
        return ''

    try:
        number = float(text)
    except ValueError:
        return text

    if number != int(number):
        return text

    return str(int(number))


def newValidation():
    """
    Create a fresh validation report.

    Returns:
      dict: zeroed counters and an empty error tally.
    """
    return {'header_valid': None, 'header_error': None, 'rows_checked': 0,
            'rows_invalid': 0, 'errors': collections.Counter()}


def summarizeValidation(validation):
    """
    Render a validation report as text.

    Arguments:
      validation (dict): a report from newValidation().

    Returns:
      str: a human readable summary.
    """
    lines = ['header_valid: %s' % validation['header_valid']]
    if validation['header_error']:
        lines.append('header_error: %s' % validation['header_error'])
    lines.append('rows_checked: %d' % validation['rows_checked'])
    lines.append('rows_invalid: %d' % validation['rows_invalid'])
    for message, count in validation['errors'].most_common():
        lines.append('  %6d  %s' % (count, message))

    return '\n'.join(lines) + '\n'


class AirrWriter:
    """
    Stream normalized records into an AIRR rearrangement TSV, one chunk at a time.

    A push-style writer rather than a function over an iterable, so that one
    pass over a converted unit can feed this and a FastaWriter at once: the
    chunks come from a generator that reads and normalizes a multi-GB file,
    and pulling it once per output format meant doing that work once per
    format too.

    Validation is performed row by row as the records stream past, using the AIRR
    schema primitives, and is reported rather than enforced. There is no
    validating writer in the airr package: RearrangementWriter takes no validate
    argument, so a report has to be accumulated here.

    The writer is created with the airr default base=1. Change-O passes base=0
    because its internal model is zero based; OAS coordinates are already one
    based, so copying that would shift every start and end coordinate by one.

    Arguments:
      out (Path): output path.
      strict (bool): if True, drop columns the AIRR schema does not define.
      validation (dict): a report to accumulate into, from newValidation().

    Attributes:
      validation (dict): the validation report, complete once closed.
    """

    def __init__(self, out, strict=False, validation=None):
        self.out = Path(out)
        self.strict = strict
        self.validation = newValidation() if validation is None else validation
        self._writer = None
        self._handle = None

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        self.close()

    def write(self, frame):
        """
        Write one chunk, opening the file on the first.

        Arguments:
          frame (pandas.DataFrame): normalized records.
        """
        import airr
        from airr.schema import RearrangementSchema, ValidationError

        if self._writer is None:
            fields = list(frame.columns)
            if self.strict:
                fields = [x for x in fields if x in RearrangementSchema.properties]
            self.out.parent.mkdir(parents=True, exist_ok=True)
            self._handle = open(self.out, 'w')
            self._writer = airr.io.RearrangementWriter(self._handle, fields=fields,
                                                       base=1)
            try:
                RearrangementSchema.validate_header(self._writer.fields)
                self.validation['header_valid'] = True
            except ValidationError as error:
                self.validation['header_valid'] = False
                self.validation['header_error'] = str(error)

        for record in frame.to_dict('records'):
            self.validation['rows_checked'] += 1
            try:
                RearrangementSchema.validate_row(record)
            except ValidationError as error:
                self.validation['rows_invalid'] += 1
                self.validation['errors'][str(error)] += 1
            self._writer.write(record)

    def close(self):
        """Flush and close; with no chunks written, still leave an empty file."""
        if self._writer is not None:
            self._writer.close()
        elif self._handle is not None:
            self._handle.close()
        else:
            # No chunks at all: still produce a file so downstream steps do
            # not have to special case its absence.
            self.out.parent.mkdir(parents=True, exist_ok=True)
            self.out.write_text('')
        self._writer, self._handle = None, None


def writeAirr(chunks, out, strict=False, validation=None):
    """
    Write normalized records as an AIRR rearrangement TSV.

    A convenience over AirrWriter for a caller with a single output; see the
    class for the details.

    Arguments:
      chunks (iterable): DataFrames of normalized records.
      out (Path): output path.
      strict (bool): if True, drop columns the AIRR schema does not define.
      validation (dict): a report to accumulate into, from newValidation().

    Returns:
      dict: the validation report.
    """
    with AirrWriter(out, strict=strict, validation=validation) as writer:
        for frame in chunks:
            writer.write(frame)

    return writer.validation


def writeValidationReport(validation, out):
    """
    Write a validation summary next to its rearrangement file.

    Arguments:
      validation (dict): the report.
      out (Path): the rearrangement file the report describes.

    Returns:
      Path: the report file.
    """
    report = Path(str(out) + '.validation.txt')
    report.write_text(summarizeValidation(validation))

    return report


class FastaWriter:
    """
    Stream normalized records into a FASTA file with pRESTO style headers.

    Push-style for the same reason as AirrWriter: so both formats can be fed
    from one pass over a converted unit.

    Headers carry the cell, so that the pairing of heavy and light chains
    survives a format that has no other place to put it. See FASTA_ANNOTATIONS
    for why nothing else is carried.

    An annotation with no value is left out rather than written as a placeholder.
    pRESTO parses the header into a dictionary keyed by annotation name, so an
    absent key reads as absent, whereas a placeholder is indistinguishable from
    data and ends up in the rearrangement table as a literal 'NA'. Bulk records
    have no cell_id at all, and must not appear to have one.

    Arguments:
      out (Path): output path.
      sequence_field (str): which column holds the sequence.
      annotations (tuple): AIRR column names to carry as key=value pairs.

    Attributes:
      written (int): the number of records written so far.
    """

    def __init__(self, out, sequence_field='sequence',
                 annotations=FASTA_ANNOTATIONS):
        self.out = Path(out)
        self.sequence_field = sequence_field
        self.annotations = annotations
        self.written = 0
        self.out.parent.mkdir(parents=True, exist_ok=True)
        # Opened up front, unlike AirrWriter, because the FASTA header needs
        # no column list from the first chunk; an empty unit still gets a file.
        self._handle = open(self.out, 'w')

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        self.close()

    def write(self, frame):
        """
        Write one chunk.

        Arguments:
          frame (pandas.DataFrame): normalized records.
        """
        for record in frame.to_dict('records'):
            sequence = str(record.get(self.sequence_field, '') or '')
            if not sequence:
                continue

            fields = ['%s=%s' % (x, record[x]) for x in self.annotations
                      if record.get(x)]
            header = '|'.join([record.get('sequence_id') or FASTA_NULL] + fields)
            self._handle.write('>%s\n%s\n' % (header, sequence))
            self.written += 1

    def close(self):
        """Close the file."""
        if self._handle is not None:
            self._handle.close()
            self._handle = None


def writeFasta(chunks, out, sequence_field='sequence',
               annotations=FASTA_ANNOTATIONS):
    """
    Write normalized records as FASTA with pRESTO style headers.

    A convenience over FastaWriter for a caller with a single output; see the
    class for the details.

    Arguments:
      chunks (iterable): DataFrames of normalized records.
      out (Path): output path.
      sequence_field (str): which column holds the sequence.
      annotations (tuple): AIRR column names to carry as key=value pairs.

    Returns:
      int: the number of records written.
    """
    with FastaWriter(out, sequence_field=sequence_field,
                     annotations=annotations) as writer:
        for frame in chunks:
            writer.write(frame)

    return writer.written
