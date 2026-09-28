"""
Incremental gzip decoding

OAS data units are multi member gzip streams: the first member holds the JSON
metadata line and the second holds the CSV. Python's gzip.open joins members
transparently, so ordinary reading is unaffected, but anything decoding raw bytes
must handle it explicitly. A single decompressobj stops at the end of the first
member and reports success, which yields the metadata line and an apparently
empty file rather than an error.
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import zlib

#: gzip wrapped deflate; see zlib.decompressobj.
GZIP_WBITS = 31


def decompressPrefix(data):
    """
    Decode as much of a possibly truncated multi member gzip stream as possible.

    Used on byte ranges fetched from the head of a remote file, where the final
    member is expected to be cut off mid stream. A truncated tail is normal here
    and is not an error.

    Arguments:
      data (bytes): the leading bytes of a gzip stream.

    Returns:
      bytes: everything that could be decoded.
    """
    out = bytearray()
    rest = data

    while rest:
        decoder = zlib.decompressobj(GZIP_WBITS)
        try:
            out += decoder.decompress(rest)
        except zlib.error:
            # A member that cannot be started at all; keep what we already have.
            break

        if not decoder.eof:
            # The stream ran out inside this member, which is the expected end
            # state for a ranged prefix.
            break

        rest = decoder.unused_data

    return bytes(out)


def countCompleteLines(data, encoding='utf-8'):
    """
    Count newline terminated lines in decoded bytes.

    Arguments:
      data (bytes): decoded text.
      encoding (str): text encoding.

    Returns:
      int: the number of complete lines.
    """
    return data.decode(encoding, errors='replace').count('\n')


def hasCompleteLines(count, encoding='utf-8'):
    """
    Build a predicate for HttpClient.readRanges.

    Arguments:
      count (int): how many complete lines are needed.
      encoding (str): text encoding.

    Returns:
      callable: given accumulated raw bytes, returns True once the decoded
      prefix contains at least count complete lines.
    """
    def predicate(raw):
        return countCompleteLines(decompressPrefix(raw), encoding) >= count

    return predicate
