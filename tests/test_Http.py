"""
Unit tests for the HTTP client
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

# Sourcerer imports
from sourcerer.Exceptions import HttpError, ProbeIncompleteError
from sourcerer.Http import HttpClient, Validators, hashFile, parseContentRangeTotal
from tests.FakeHttp import Boom, FakeResponse, FakeSession, rangeHandler, sequenceHandler

URL = 'https://example.org/unit.csv.gz'
BODY = bytes(range(256)) * 40  # 10,240 bytes, non-repeating enough to catch splices


def makeClient(handler):
    """Build an HttpClient with no politeness delay over a scripted session."""
    return HttpClient(delay=0, backoff=0, session=FakeSession(handler))


class TestHelpers(unittest.TestCase):
    """
    Tests for the module level helpers
    """

    def test_parse_content_range_total(self):
        self.assertEqual(parseContentRangeTotal('bytes 0-99/1234'), 1234)
        self.assertIsNone(parseContentRangeTotal('bytes 0-99/*'))
        self.assertIsNone(parseContentRangeTotal(None))

    def test_validators_require_a_tag_to_match(self):
        """
        Validators without a tag never match.

        Resuming against an object we cannot identify risks concatenating two
        different bodies, so an unidentifiable object must not compare equal.
        """
        anonymous = Validators(size_bytes=10)
        self.assertFalse(anonymous.matches(Validators(size_bytes=10)))
        tagged = Validators(etag='"v1"', size_bytes=10)
        self.assertTrue(tagged.matches(Validators(etag='"v1"', size_bytes=10)))
        self.assertFalse(tagged.matches(Validators(etag='"v2"', size_bytes=10)))


class TestRetry(unittest.TestCase):
    """
    Tests for transient failure handling
    """

    def setUp(self):
        patcher = mock.patch('sourcerer.Http.time.sleep')
        self.sleep = patcher.start()
        self.addCleanup(patcher.stop)

    def test_retries_then_succeeds(self):
        """A 500 followed by a 200 is retried, not surfaced."""
        client = makeClient(sequenceHandler([
            FakeResponse(500), FakeResponse(500), FakeResponse(200, b'ok')]))
        response = client.get(URL)

        self.assertEqual(response.status_code, 200)
        self.assertEqual(len(client.session.calls), 3)

    def test_retries_connection_errors(self):
        """Transport failures are retried like retryable statuses."""
        client = makeClient(sequenceHandler([Boom('reset'), FakeResponse(200, b'ok')]))

        self.assertEqual(client.get(URL).status_code, 200)

    def test_gives_up_and_raises(self):
        """Exhausting retries raises HttpError rather than returning a bad response."""
        client = makeClient(sequenceHandler([FakeResponse(503)]))

        with self.assertRaises(HttpError):
            client.get(URL)

        self.assertEqual(len(client.session.calls), 4)

    def test_does_not_retry_client_errors(self):
        """A 404 is a definitive answer and is returned immediately."""
        client = makeClient(sequenceHandler([FakeResponse(404)]))

        self.assertEqual(client.get(URL).status_code, 404)
        self.assertEqual(len(client.session.calls), 1)

    def test_honors_retry_after(self):
        """Retry-After overrides the computed backoff."""
        client = makeClient(sequenceHandler([
            FakeResponse(429, headers={'Retry-After': '7'}), FakeResponse(200)]))
        client.get(URL)

        self.assertIn(7.0, [call.args[0] for call in self.sleep.call_args_list])


class TestFetch(unittest.TestCase):
    """
    Tests for downloading, resuming and verifying
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dest = Path(self.tmp.name) / 'unit.csv.gz'
        self.digest = hashlib.sha256(BODY).hexdigest()

        patcher = mock.patch('sourcerer.Http.time.sleep')
        patcher.start()
        self.addCleanup(patcher.stop)

    def stagePartial(self, nbytes, etag='"v1"', size=len(BODY)):
        """Write a partial download plus its resume record, as an interruption would."""
        temp = self.dest.with_name(self.dest.name + '.tmp')
        temp.write_bytes(BODY[:nbytes])
        sidecar = self.dest.with_name(self.dest.name + '.tmp.json')
        sidecar.write_text(json.dumps({
            'url': URL,
            'validators': {'etag': etag, 'last_modified': None, 'size_bytes': size}}))

    def test_downloads_and_hashes_from_disk(self):
        client = makeClient(rangeHandler(BODY))
        outcome = client.fetch(URL, self.dest, progress=False)

        self.assertEqual(self.dest.read_bytes(), BODY)
        self.assertEqual(outcome.sha256, self.digest)
        self.assertEqual(outcome.size_bytes, len(BODY))
        self.assertFalse(outcome.resumed)

    def test_cleans_up_temporary_files(self):
        client = makeClient(rangeHandler(BODY))
        client.fetch(URL, self.dest, progress=False)

        siblings = sorted(p.name for p in self.dest.parent.iterdir())
        self.assertEqual(siblings, ['unit.csv.gz'])

    def test_resumes_and_matches_a_single_shot_download(self):
        """
        A resumed transfer produces the same bytes and digest as a whole one.

        This is the property that matters: hashing only the streamed tail would
        pass a naive test but produce a digest describing a fragment.
        """
        self.stagePartial(4096)
        client = makeClient(rangeHandler(BODY))
        outcome = client.fetch(URL, self.dest, progress=False)

        self.assertTrue(outcome.resumed)
        self.assertEqual(self.dest.read_bytes(), BODY)
        self.assertEqual(outcome.sha256, self.digest)

    def test_restarts_when_server_ignores_range(self):
        """
        A 200 answer to a Range request means start over, not append.

        Appending a whole body to a partial file is the silent corruption this
        guards against.
        """
        self.stagePartial(4096)
        client = makeClient(rangeHandler(BODY, support_range=False))
        outcome = client.fetch(URL, self.dest, progress=False)

        self.assertFalse(outcome.resumed)
        self.assertEqual(self.dest.read_bytes(), BODY)
        self.assertEqual(outcome.sha256, self.digest)

    def test_restarts_when_the_object_changed_size(self):
        """
        A 206 whose total disagrees with the partial file triggers a clean restart.

        The partial response body starts mid-object, so it must not be reused;
        the client has to re-request from byte zero.
        """
        replacement = BODY + b'appended'
        self.stagePartial(4096, size=len(BODY))
        client = makeClient(rangeHandler(replacement, serve=replacement))
        outcome = client.fetch(URL, self.dest, progress=False)

        self.assertFalse(outcome.resumed)
        self.assertEqual(self.dest.read_bytes(), replacement)
        self.assertEqual(outcome.sha256, hashlib.sha256(replacement).hexdigest())

    def test_skips_an_existing_file(self):
        self.dest.write_bytes(BODY)
        client = makeClient(rangeHandler(BODY))
        outcome = client.fetch(URL, self.dest, progress=False)

        self.assertTrue(outcome.skipped)
        self.assertEqual(outcome.sha256, self.digest)
        self.assertEqual(len(client.session.calls), 0)

    def test_redownloads_when_the_digest_disagrees(self):
        self.dest.write_bytes(b'stale')
        client = makeClient(rangeHandler(BODY))
        outcome = client.fetch(URL, self.dest, progress=False,
                               expected_sha256=self.digest)

        self.assertFalse(outcome.skipped)
        self.assertEqual(self.dest.read_bytes(), BODY)

    def test_hash_file_matches_hashlib(self):
        self.dest.write_bytes(BODY)
        self.assertEqual(hashFile(self.dest), self.digest)


class TestProbes(unittest.TestCase):
    """
    Tests for liveness and ranged head probes
    """

    def setUp(self):
        patcher = mock.patch('sourcerer.Http.time.sleep')
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_probe_alive_uses_head_when_it_works(self):
        client = makeClient(sequenceHandler([FakeResponse(200)]))

        self.assertTrue(client.probeAlive(URL))
        self.assertEqual(client.session.calls[0]['method'], 'HEAD')

    def test_probe_alive_falls_back_to_range(self):
        """
        A host that mishandles HEAD is not a dead URL.

        HEAD is tried first because it is cheap, but a failure only means the
        method is unreliable, so the probe falls back to a one byte GET.
        """
        def handler(method, url, headers, index):
            if method == 'HEAD':
                raise Boom('no HEAD here')
            return FakeResponse(206, b'x', {'Content-Range': 'bytes 0-0/100'})

        client = makeClient(handler)

        self.assertTrue(client.probeAlive(URL))
        self.assertEqual(client.session.calls[-1]['method'], 'GET')

    def test_read_ranges_extends_until_complete(self):
        """The window grows on demand rather than assuming a fixed prefix size."""
        client = makeClient(rangeHandler(BODY))
        prefix = client.readRanges(URL, lambda buf: len(buf) >= 3000, initial=1024)

        self.assertGreaterEqual(len(prefix), 3000)
        self.assertEqual(prefix, BODY[:len(prefix)])
        self.assertGreater(len(client.session.calls), 1)

    def test_read_ranges_reports_incomplete_rather_than_drift(self):
        """
        Hitting the cap raises ProbeIncompleteError.

        This must stay distinct from a structural difference: a truncated probe is
        a harvest failure, and treating it as drift would open a pull request
        claiming the remote format changed when it did not.
        """
        client = makeClient(rangeHandler(BODY))

        with self.assertRaises(ProbeIncompleteError):
            client.readRanges(URL, lambda buf: False, initial=512, cap=2048)


if __name__ == '__main__':
    unittest.main()
