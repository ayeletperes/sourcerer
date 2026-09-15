"""
HTTP client

Every network call in sourcerer goes through HttpClient. That single seam is what
lets the test suite run with no network at all, and it is where politeness,
retries, timeouts and download integrity are enforced once rather than at each
call site.
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import hashlib
import json
import logging
import os
import random
import re
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import requests
from tqdm import tqdm

# Sourcerer imports
from sourcerer.Exceptions import HttpError, ProbeIncompleteError
from sourcerer.Version import __version__

log = logging.getLogger(__name__)

#: Status codes worth retrying. Everything else is a definitive answer.
RETRY_STATUS = frozenset([429, 500, 502, 503, 504])

#: Identify honestly. Spoofing a browser user agent is what gets an academic
#: host to block you, and it is bad citizenship toward a group providing free
#: data. See OASTools.py:65 for the anti-pattern this replaces.
USER_AGENT = ('sourcerer/%s (+https://github.com/immcantation/sourcerer; '
              'immcantation@googlegroups.com)' % __version__)

#: Read chunk size for streamed bodies.
CHUNK_BYTES = 1 << 16


@dataclass(frozen=True)
class Validators:
    """
    Server-supplied identity of a remote object at a point in time.

    Arguments:
      etag (str): the ETag header, or None.
      last_modified (str): the Last-Modified header, or None.
      size_bytes (int): the total size in bytes, or None if the server did not say.
    """
    etag: str = None
    last_modified: str = None
    size_bytes: int = None

    @property
    def tag(self):
        """str: the value to send in If-Range, preferring the strong validator."""
        return self.etag or self.last_modified

    def matches(self, other):
        """
        Test whether two validator sets describe the same remote object.

        A missing validator on either side is treated as a mismatch: resuming
        against an object we cannot identify risks concatenating two different
        bodies, which is far worse than re-fetching.

        Arguments:
          other (Validators): the validators to compare against.

        Returns:
          bool: True if the objects are provably the same.
        """
        if other is None or self.tag is None or other.tag is None:
            return False

        return self.tag == other.tag and self.size_bytes == other.size_bytes


@dataclass(frozen=True)
class FetchOutcome:
    """
    The result of fetching one remote object to disk.

    Arguments:
      path (Path): the completed file.
      sha256 (str): digest computed by reading the finished file from disk.
      size_bytes (int): size of the completed file.
      validators (Validators): server validators recorded at download time.
      resumed (bool): whether an interrupted transfer was continued.
      skipped (bool): whether the file was already present and complete.
    """
    path: Path
    sha256: str
    size_bytes: int
    validators: Validators
    resumed: bool = False
    skipped: bool = False


def hashFile(path, algorithm='sha256'):
    """
    Compute a digest by reading a completed file from disk.

    Always hash the finished file rather than the bytes streamed in this session.
    After a resumed download the streamed bytes are only the tail, so a digest
    accumulated while streaming would describe a fragment while appearing to
    describe the whole file.

    Arguments:
      path (Path): file to hash.
      algorithm (str): hashlib algorithm name.

    Returns:
      str: the hexadecimal digest.
    """
    digest = hashlib.new(algorithm)
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(CHUNK_BYTES), b''):
            digest.update(block)

    return digest.hexdigest()


def parseContentRangeTotal(value):
    """
    Extract the total object size from a Content-Range header.

    Arguments:
      value (str): the header value, e.g. 'bytes 100-199/1234'.

    Returns:
      int: the total size, or None if absent or unknown ('*').
    """
    if not value:
        return None

    match = re.search(r'/(\d+)\s*$', value)

    return int(match.group(1)) if match else None


class HttpClient:
    """
    A polite, retrying HTTP client with resumable downloads.
    """

    def __init__(self, user_agent=USER_AGENT, delay=0.5, max_retries=3,
                 backoff=1.0, connect_timeout=10, read_timeout=120,
                 body_timeout=300, session=None):
        """
        Arguments:
          user_agent (str): value of the User-Agent header.
          delay (float): minimum seconds between requests to the same host.
          max_retries (int): attempts after the first before giving up.
          backoff (float): base seconds for exponential backoff.
          connect_timeout (float): seconds to wait for a connection.
          read_timeout (float): seconds to wait for page-sized responses.
          body_timeout (float): seconds to wait for streamed body responses.
          session (requests.Session): an existing session, mainly for testing.
        """
        self.delay = delay
        self.max_retries = max_retries
        self.backoff = backoff
        self.connect_timeout = connect_timeout
        self.read_timeout = read_timeout
        self.body_timeout = body_timeout

        self.session = session if session is not None else requests.Session()
        self.session.headers.update({'User-Agent': user_agent})

        self._last_request = 0.0

    def _sleepForPoliteness(self):
        """Space requests out by at least self.delay seconds."""
        if self.delay <= 0:
            return

        elapsed = time.monotonic() - self._last_request
        if elapsed < self.delay:
            time.sleep(self.delay - elapsed)

    def _retryDelay(self, attempt, response):
        """
        Compute how long to wait before the next attempt.

        Honors Retry-After when the server supplies it, otherwise uses exponential
        backoff with jitter so that concurrent clients do not resynchronize.

        Arguments:
          attempt (int): zero-based attempt number that just failed.
          response (requests.Response): the failed response, or None.

        Returns:
          float: seconds to wait.
        """
        if response is not None:
            retry_after = response.headers.get('Retry-After')
            if retry_after and retry_after.strip().isdigit():
                return float(retry_after.strip())

        return self.backoff * (2 ** attempt) + random.uniform(0, self.backoff)

    def request(self, method, url, stream=False, **kwargs):
        """
        Issue a request, retrying transient failures.

        Arguments:
          method (str): HTTP method.
          url (str): absolute URL.
          stream (bool): if True do not preload the body.
          kwargs: passed through to requests.

        Returns:
          requests.Response: the successful response.

        Raises:
          HttpError: if every attempt failed.
        """
        timeout = kwargs.pop('timeout', None)
        if timeout is None:
            read = self.body_timeout if stream else self.read_timeout
            timeout = (self.connect_timeout, read)

        last_error = None
        for attempt in range(self.max_retries + 1):
            self._sleepForPoliteness()
            response = None
            try:
                response = self.session.request(method, url, stream=stream,
                                                timeout=timeout, **kwargs)
            except requests.RequestException as error:
                last_error = error
            finally:
                self._last_request = time.monotonic()

            if response is not None and response.status_code not in RETRY_STATUS:
                return response

            if response is not None:
                last_error = HttpError('%s %s returned HTTP %d'
                                       % (method, url, response.status_code))

            if attempt < self.max_retries:
                wait = self._retryDelay(attempt, response)
                log.warning('%s %s failed (%s); retrying in %.1fs',
                            method, url, last_error, wait)
                time.sleep(wait)

        raise HttpError('%s %s failed after %d attempts: %s'
                        % (method, url, self.max_retries + 1, last_error))

    def get(self, url, **kwargs):
        """Issue a GET request. Returns a requests.Response."""
        return self.request('GET', url, **kwargs)

    def post(self, url, **kwargs):
        """Issue a POST request. Returns a requests.Response."""
        return self.request('POST', url, **kwargs)

    def head(self, url, **kwargs):
        """Issue a HEAD request. Returns a requests.Response."""
        return self.request('HEAD', url, **kwargs)

    def getRange(self, url, start, end=None, **kwargs):
        """
        Request a byte range.

        Arguments:
          url (str): absolute URL.
          start (int): first byte offset, inclusive.
          end (int): last byte offset inclusive, or None for open-ended.

        Returns:
          requests.Response: the response, typically 206 but 200 if the server
          ignored the range.
        """
        headers = dict(kwargs.pop('headers', {}))
        headers['Range'] = 'bytes=%d-%s' % (start, '' if end is None else end)

        return self.request('GET', url, headers=headers, **kwargs)

    def probeAlive(self, url):
        """
        Check that a URL resolves without downloading it.

        Tries HEAD first as the cheapest option, then falls back to a one byte
        ranged GET. Some hosts handle HEAD unreliably, so a HEAD failure alone is
        not evidence that the URL is dead.

        Arguments:
          url (str): absolute URL.

        Returns:
          bool: True if the object appears to exist.
        """
        try:
            response = self.head(url)
            if response.status_code < 400:
                return True
        except HttpError:
            log.debug('HEAD unsupported or failed for %s; falling back to range', url)

        try:
            response = self.getRange(url, 0, 0, stream=True)
        except HttpError:
            return False

        try:
            return response.status_code in (200, 206)
        finally:
            response.close()

    def readRanges(self, url, is_complete, initial=1 << 18, cap=1 << 23):
        """
        Fetch increasing contiguous byte ranges until enough data has arrived.

        Used to inspect the head of a large gzipped file without downloading it.
        A fixed window would stop being large enough the moment the remote file's
        header grew, turning benign upstream growth into a false structural
        difference, so the window extends on demand instead.

        Arguments:
          url (str): absolute URL.
          is_complete (callable): given the bytes accumulated so far, returns True
            when no more data is needed.
          initial (int): size of the first range request.
          cap (int): maximum total bytes to fetch before giving up.

        Returns:
          bytes: the accumulated prefix of the remote object.

        Raises:
          ProbeIncompleteError: if the cap was reached without is_complete
            returning True, or the server stopped supplying data early.
        """
        buffer = b''
        want = initial

        while len(buffer) < cap:
            end = min(len(buffer) + want, cap) - 1
            response = self.getRange(url, len(buffer), end, stream=True)

            if response.status_code not in (200, 206):
                response.close()
                raise ProbeIncompleteError(
                    'range probe of %s returned HTTP %d' % (url, response.status_code))

            body = response.content
            served_whole = response.status_code == 200
            response.close()

            if served_whole:
                # Server ignored Range and sent everything; nothing left to ask for.
                buffer = body
                if is_complete(buffer):
                    return buffer
                raise ProbeIncompleteError(
                    'whole body of %s did not contain the expected structure' % url)

            if not body:
                raise ProbeIncompleteError(
                    'range probe of %s returned no data at offset %d'
                    % (url, len(buffer)))

            buffer += body
            if is_complete(buffer):
                return buffer

            want = min(want * 2, cap)

        raise ProbeIncompleteError(
            'range probe of %s reached the %d byte cap without finding the '
            'expected structure' % (url, cap))

    def fetch(self, url, dest, resume=True, progress=True, expected_sha256=None):
        """
        Download a URL to a path, resuming safely and verifying from disk.

        The remote source publishes no checksums, so integrity here means
        self-consistency: the digest is computed by re-reading the finished file,
        and a partial transfer is only continued when the server proves the object
        is unchanged. If it cannot prove that, the partial file is discarded and
        the download restarts, because a silently concatenated old-plus-new file
        is far worse than spending the bandwidth again.

        Arguments:
          url (str): absolute URL.
          dest (Path): final output path. Parent directories are created.
          resume (bool): whether to continue an interrupted transfer.
          progress (bool): whether to show a progress bar.
          expected_sha256 (str): if given and dest already matches, skip the
            download.

        Returns:
          FetchOutcome: what happened, including the digest and validators.
        """
        dest = Path(dest)
        dest.parent.mkdir(parents=True, exist_ok=True)
        temp = dest.with_name(dest.name + '.tmp')
        sidecar = dest.with_name(dest.name + '.tmp.json')

        if dest.exists():
            digest = hashFile(dest)
            if expected_sha256 is None or digest == expected_sha256:
                log.info('%s already present; skipping', dest.name)
                return FetchOutcome(path=dest, sha256=digest,
                                    size_bytes=dest.stat().st_size,
                                    validators=Validators(), skipped=True)
            log.warning('%s exists but does not match the expected digest; '
                        're-downloading', dest.name)
            dest.unlink()

        offset, known = self._resumeState(temp, sidecar) if resume else (0, None)

        response, offset, resumed = self._openStream(url, offset, known)
        try:
            validators = self._responseValidators(response, offset)

            self._writeSidecar(sidecar, url, validators)
            total = validators.size_bytes

            mode = 'ab' if resumed else 'wb'
            # A progress bar redrawn into a pipe or a log file emits one line per
            # update, which buries everything else. Only show it on a terminal.
            show = progress and sys.stderr.isatty()
            with open(temp, mode) as handle:
                bar = tqdm(total=total, initial=offset, unit='B', unit_scale=True,
                           desc=dest.name, disable=not show, leave=False)
                with bar:
                    for block in response.iter_content(chunk_size=CHUNK_BYTES):
                        if block:
                            handle.write(block)
                            bar.update(len(block))
        finally:
            response.close()

        os.replace(temp, dest)
        sidecar.unlink(missing_ok=True)

        digest = hashFile(dest)
        size = dest.stat().st_size

        if total is not None and size != total:
            log.warning('%s finished at %d bytes but the server reported %d',
                        dest.name, size, total)

        return FetchOutcome(path=dest, sha256=digest, size_bytes=size,
                            validators=validators, resumed=resumed)

    def _resumeState(self, temp, sidecar):
        """
        Recover how far a previous attempt got, and what it was downloading.

        Arguments:
          temp (Path): the partial file.
          sidecar (Path): the JSON record written alongside it.

        Returns:
          tuple: (offset, Validators) where offset is 0 when resuming is not
          possible.
        """
        if not temp.exists() or not sidecar.exists():
            return 0, None

        try:
            with open(sidecar) as handle:
                stored = json.load(handle)
            known = Validators(**stored['validators'])
        except (OSError, ValueError, KeyError, TypeError):
            log.warning('unreadable resume record %s; starting over', sidecar.name)
            return 0, None

        return temp.stat().st_size, known

    def _openStream(self, url, offset, known):
        """
        Open the response to stream from, resuming only when it is provably safe.

        A 206 whose body starts mid-object cannot be written from byte zero, so
        when a partial response turns out not to match the file on disk this
        re-issues the request without a Range header rather than reusing the
        response it already holds.

        Arguments:
          url (str): absolute URL.
          offset (int): bytes already on disk.
          known (Validators): validators recorded when the partial file was made.

        Returns:
          tuple: (response, offset, resumed). offset is 0 when starting over.

        Raises:
          HttpError: if the response status is not usable.
        """
        headers = {}
        if offset and known is not None and known.tag:
            headers['Range'] = 'bytes=%d-' % offset
            headers['If-Range'] = known.tag

        response = self.request('GET', url, stream=True, headers=headers)

        if response.status_code == 206:
            total = parseContentRangeTotal(response.headers.get('Content-Range'))
            if known is not None and known.size_bytes is not None \
                    and total is not None and total != known.size_bytes:
                response.close()
                log.warning('%s changed size since the partial download; '
                            'restarting from the beginning', url)
                return self.request('GET', url, stream=True), 0, False
            return response, offset, True

        if response.status_code == 200:
            # Either we did not ask for a range, or If-Range failed because the
            # object changed. Both mean start from scratch, and a 200 body always
            # begins at byte zero so it is safe to use directly.
            if offset:
                log.warning('%s no longer matches the partial download; '
                            'restarting from the beginning', url)
            return response, 0, False

        response.close()
        raise HttpError('GET %s returned HTTP %d' % (url, response.status_code))

    def _responseValidators(self, response, offset):
        """
        Build validators describing the whole object, not just this response.

        Arguments:
          response (requests.Response): the opened response.
          offset (int): bytes already on disk.

        Returns:
          Validators: etag, last modified and total size where known.
        """
        total = parseContentRangeTotal(response.headers.get('Content-Range'))
        if total is None:
            length = response.headers.get('Content-Length')
            if length is not None and length.isdigit():
                total = int(length) + offset

        return Validators(etag=response.headers.get('ETag'),
                          last_modified=response.headers.get('Last-Modified'),
                          size_bytes=total)

    def _writeSidecar(self, sidecar, url, validators):
        """
        Record what the partial file is, so a later run can resume it safely.

        Arguments:
          sidecar (Path): where to write.
          url (str): the URL being downloaded.
          validators (Validators): the server validators for the object.
        """
        payload = {'url': url, 'validators': asdict(validators)}
        with open(sidecar, 'w') as handle:
            json.dump(payload, handle, sort_keys=True, indent=2)
