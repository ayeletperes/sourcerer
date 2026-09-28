"""
Typed exceptions

Parsing and probing failures are always raised, never swallowed into an empty
result. A scraper that returns {} on failure produces confidently wrong output;
these exceptions exist so that every failure names what was expected and where.
"""

# Info
__author__ = 'Susanna Marquez'


class ImmvokeError(Exception):
    """Base class for all immvoke errors."""
    pass


class HttpError(ImmvokeError):
    """A request failed after exhausting retries."""
    pass


class ProbeIncompleteError(ImmvokeError):
    """
    A progressive range probe hit its byte cap without decoding what it needed.

    This is a harvest failure, not schema drift. Conflating the two would make a
    slow or truncated response look like an upstream format change.
    """
    pass


class ParseError(ImmvokeError):
    """Remote content did not match the structure the code expects."""
    pass


class OasParseError(ParseError):
    """OAS content did not match the expected structure."""
    pass


class ImgtParseError(ParseError):
    """IMGT content did not match the expected structure."""
    pass


class OgrdbParseError(ParseError):
    """OGRDB content did not match the expected structure."""
    pass


class SchemaError(ImmvokeError):
    """A stored schema snapshot is missing, malformed or too new to understand."""
    pass


class ConversionError(ImmvokeError):
    """A data unit could not be converted."""
    pass
