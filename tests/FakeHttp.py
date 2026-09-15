"""
A scriptable stand-in for requests.Session

Injecting a fake session rather than a fake HttpClient means the tests exercise
the real retry, resume and integrity logic instead of a parallel implementation
of it. Fault injection lives here so that each test can describe the failure it
cares about in one line.
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import requests
from requests.structures import CaseInsensitiveDict


class FakeResponse:
    """
    A minimal stand-in for requests.Response.
    """

    def __init__(self, status_code=200, body=b'', headers=None):
        """
        Arguments:
          status_code (int): HTTP status to report.
          body (bytes): the response body.
          headers (dict): response headers.
        """
        self.status_code = status_code
        self.content = body
        self.headers = CaseInsensitiveDict(headers or {})
        self.closed = False

    def iter_content(self, chunk_size=1024):
        """Yield the body in chunk_size pieces."""
        for start in range(0, len(self.content), chunk_size):
            yield self.content[start:start + chunk_size]

    def close(self):
        """Mark the response closed."""
        self.closed = True


class FakeSession:
    """
    A requests.Session replacement driven by a handler callable.
    """

    def __init__(self, handler):
        """
        Arguments:
          handler (callable): called as handler(method, url, headers, call_index)
            and returns a FakeResponse, or raises to simulate a transport error.
        """
        self.headers = CaseInsensitiveDict()
        self.handler = handler
        self.calls = []

    def request(self, method, url, stream=False, timeout=None, headers=None,
                **kwargs):
        """Record the call and delegate to the handler."""
        merged = CaseInsensitiveDict(self.headers)
        merged.update(headers or {})
        self.calls.append({'method': method, 'url': url, 'headers': merged,
                           'stream': stream, 'kwargs': kwargs})

        return self.handler(method, url, merged, len(self.calls) - 1)


def sequenceHandler(responses):
    """
    Build a handler that returns each response in turn.

    Arguments:
      responses (list): FakeResponse objects or exception instances to raise.

    Returns:
      callable: a handler suitable for FakeSession.
    """
    def handler(method, url, headers, index):
        item = responses[min(index, len(responses) - 1)]
        if isinstance(item, Exception):
            raise item
        return item

    return handler


def rangeHandler(body, etag='"v1"', support_range=True, serve=None):
    """
    Build a handler that serves a body, honoring Range requests.

    Arguments:
      body (bytes): the full object.
      etag (str): the ETag to report.
      support_range (bool): if False, always answer 200 with the whole body,
        which is how a server that ignores Range behaves.
      serve (bytes): an alternative body to serve, used to simulate the object
        changing underneath a partial download.

    Returns:
      callable: a handler suitable for FakeSession.
    """
    def handler(method, url, headers, index):
        payload = body if serve is None else serve
        common = {'ETag': etag, 'Last-Modified': 'Mon, 04 Aug 2026 00:00:00 GMT'}

        range_header = headers.get('Range')
        if range_header and support_range:
            start = int(range_header.split('=')[1].split('-')[0])
            end_text = range_header.split('-')[1]
            end = int(end_text) if end_text else len(payload) - 1
            end = min(end, len(payload) - 1)
            chunk = payload[start:end + 1]
            common['Content-Range'] = 'bytes %d-%d/%d' % (start, end, len(payload))
            common['Content-Length'] = str(len(chunk))
            return FakeResponse(206, chunk, common)

        common['Content-Length'] = str(len(payload))
        return FakeResponse(200, payload, common)

    return handler


class Boom(requests.ConnectionError):
    """A transport level failure."""
    pass
