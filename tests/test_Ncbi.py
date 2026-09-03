"""
Unit tests for the NCBI cross-reference
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import unittest

from sourcerer.Http import HttpClient
from sourcerer.Ncbi import (
    Evidence,
    accessionFromText,
    fetchBiosamples,
    gatherEvidence,
    poolCodes,
    resolveAccessions,
    suggestSubject,
)
from tests.FakeHttp import FakeResponse, FakeSession


def esearchXml(*ids):
    """An esearch response naming the given UIDs, in order."""
    id_list = ''.join('<Id>%s</Id>' % uid for uid in ids)
    return ('<?xml version="1.0"?><eSearchResult><IdList>%s</IdList>'
           '</eSearchResult>' % id_list).encode()


def docsumXml(uid, title, biosample, run_acc):
    """One esummary DocSum, shaped like a real SRA docsum's escaped ExpXml."""
    exp_xml = ('&lt;Summary&gt;&lt;Title&gt;%s&lt;/Title&gt;&lt;/Summary&gt;'
              '&lt;Biosample&gt;%s&lt;/Biosample&gt;' % (title, biosample))
    runs_xml = '&lt;Run acc="%s" total_spots="1"/&gt;' % run_acc
    return ('<DocSum><Id>%s</Id>'
           '<Item Name="ExpXml" Type="String">%s</Item>'
           '<Item Name="Runs" Type="String">%s</Item></DocSum>'
           % (uid, exp_xml, runs_xml))


def esummaryXml(*docsums):
    """An esummary response wrapping the given DocSum fragments."""
    return ('<?xml version="1.0"?><eSummaryResult>%s</eSummaryResult>'
           % ''.join(docsums)).encode()


def biosampleXml(accession, sample_name=None, title=''):
    """One BioSample fragment, as efetch --rettype full --retmode xml returns."""
    id_xml = '<Id db="BioSample">%s</Id>' % accession
    if sample_name is not None:
        id_xml += '<Id db_label="Sample name">%s</Id>' % sample_name
    return ('<BioSample accession="%s"><Ids>%s</Ids>'
           '<Description><Title>%s</Title></Description></BioSample>'
           % (accession, id_xml, title))


def efetchXml(*biosamples):
    """A biosample efetch response wrapping the given BioSample fragments."""
    return ('<?xml version="1.0"?><BioSampleSet>%s</BioSampleSet>'
           % ''.join(biosamples)).encode()


def routedHandler(routes):
    """
    Build a FakeSession handler that answers by URL substring.

    Arguments:
      routes (dict): substring (e.g. 'esearch') to FakeResponse.

    Returns:
      callable: a handler suitable for FakeSession.
    """
    def handler(method, url, headers, index):
        for substring, response in routes.items():
            if substring in url:
                return response
        raise AssertionError('no route matched %s' % url)

    return handler


def makeClient(routes):
    """Build an HttpClient wired to a FakeSession answering by URL substring."""
    return HttpClient(delay=0, backoff=0, session=FakeSession(routedHandler(routes)))


class TestAccessionFromText(unittest.TestCase):
    """
    Tests for pulling a run/sample accession out of OAS's unit id
    """

    def test_finds_an_sra_run_accession(self):
        """An SRR accession embedded in an OAS unit id is found."""
        text = 'Corinaldesi_2024/csv_paired/SRR25557617_1_Paired_All.csv.gz'
        self.assertEqual(accessionFromText(text), 'SRR25557617')

    def test_finds_a_geo_sample_accession(self):
        """A GSM accession embedded in an OAS unit id is found."""
        text = 'McIntire_2024/csv_paired/GSM6504685_1_Paired_All.csv.gz'
        self.assertEqual(accessionFromText(text), 'GSM6504685')

    def test_accession_immediately_followed_by_an_underscore_still_matches(self):
        """
        A trailing \\b would never fire here: \\d and _ are both word
        characters, so an accession run straight into '_1_Paired_All' (the
        normal OAS shape) needs no boundary after the digits to stop at them.
        """
        self.assertEqual(accessionFromText('SRR123_1_Paired_All'), 'SRR123')

    def test_returns_none_when_no_accession_is_present(self):
        """Text with no recognizable accession returns None, not a false match."""
        self.assertIsNone(accessionFromText('not-an-accession-XSRR123'))


class TestPoolCodes(unittest.TestCase):
    """
    Tests for detecting a multi-donor pooled/hashed sample
    """

    def test_reads_off_codes_from_a_hashed_sample(self):
        """A hashed sample naming several donors reports all of them."""
        text = 'Hashed scBCR sample (FA007, FA048)'
        self.assertEqual(poolCodes(text), ('FA007', 'FA048'))

    def test_ordinary_sample_is_not_pooled(self):
        """Text with no pooling language is never treated as pooled."""
        self.assertEqual(poolCodes('BL-110_VDJ'), ())

    def test_a_single_parenthesized_code_is_not_pooled(self):
        """One code in parentheses names a single donor, not a pool."""
        self.assertEqual(poolCodes('Hashed scBCR sample (FA007)'), ())


class TestSuggestSubject(unittest.TestCase):
    """
    Tests for the generic, study-agnostic subject suggestion
    """

    def test_strips_a_trailing_locus_token(self):
        """A trailing assay/locus marker is generic enough to strip."""
        self.assertEqual(suggestSubject('BL-110_VDJ'), 'BL-110')

    def test_strips_a_comma_separated_locus_token(self):
        """The separator before the locus token may be a comma, not just _."""
        self.assertEqual(suggestSubject('ED8, VDJ'), 'ED8')

    def test_strips_a_trailing_visit_marker(self):
        """A trailing _V<n> visit marker is generic enough to strip."""
        self.assertEqual(suggestSubject('TT04_subj6_V3'), 'TT04_subj6')

    def test_leaves_a_bare_trailing_digit_alone(self):
        """
        A bare trailing number is exactly the kind of study-specific detail
        this function must never guess at: 'Donor_5' might be the fifth donor,
        not a locus/visit marker, and only the specific study's convention
        could tell those apart.
        """
        self.assertEqual(suggestSubject('Donor_5'), 'Donor_5')

    def test_leaves_unrecognized_text_unchanged(self):
        """Text matching no generic pattern is returned as-is."""
        self.assertEqual(suggestSubject('P05_FNA_d0_1_Y1'), 'P05_FNA_d0_1_Y1')


class TestResolveAccessions(unittest.TestCase):
    """
    Tests for the SRA run/sample -> BioSample lookup
    """

    def test_resolves_an_sra_native_run(self):
        """A run submitted straight to SRA resolves via its own Run acc."""
        client = makeClient({
            'esearch': FakeResponse(200, esearchXml('1')),
            'esummary': FakeResponse(200, esummaryXml(docsumXml(
                '1', 'BL-110_VDJ', 'SAMN36877371', 'SRR25557617'))),
        })

        resolved = resolveAccessions(client, ['SRR25557617'])

        self.assertEqual(resolved['SRR25557617'], ('SAMN36877371', 'BL-110_VDJ'))

    def test_resolves_a_geo_mediated_sample_by_its_gsm_token(self):
        """
        A GEO-submitted sample's docsum names no matching Run acc (the run
        accession differs from the GSM the caller searched for), so
        attribution falls back to the GSM token embedded in the docsum's own
        title.
        """
        client = makeClient({
            'esearch': FakeResponse(200, esearchXml('1')),
            'esummary': FakeResponse(200, esummaryXml(docsumXml(
                '1', 'GSM7688951: BL-110_VDJ', 'SAMN36877371', 'SRR25557617'))),
        })

        resolved = resolveAccessions(client, ['GSM7688951'])

        self.assertEqual(resolved['GSM7688951'], ('SAMN36877371', 'BL-110_VDJ'))

    def test_does_not_attribute_a_docsum_to_an_accession_outside_the_batch(self):
        """
        A docsum is matched by its own content, not by request order: one
        that names neither a batch accession's Run acc nor its GSM token
        resolves nothing, even though NCBI still returned it.
        """
        client = makeClient({
            'esearch': FakeResponse(200, esearchXml('1')),
            'esummary': FakeResponse(200, esummaryXml(docsumXml(
                '1', 'unrelated sample', 'SAMN00000000', 'SRR00000000'))),
        })

        resolved = resolveAccessions(client, ['SRR25557617'])

        self.assertEqual(resolved, {})

    def test_accession_with_no_sra_hit_is_absent_from_the_result(self):
        """An accession esearch finds nothing for is simply not in the result."""
        client = makeClient({'esearch': FakeResponse(200, esearchXml())})

        resolved = resolveAccessions(client, ['SRR99999999'])

        self.assertEqual(resolved, {})


class TestFetchBiosamples(unittest.TestCase):
    """
    Tests for reading a BioSample record's sample name
    """

    def test_prefers_the_sample_name_identifier(self):
        """The dedicated 'Sample name' Id wins over the record's own Title."""
        client = makeClient({'efetch': FakeResponse(200, efetchXml(biosampleXml(
            'SAMN15865971', sample_name='32_CSF_uns_5pIGSEQ_2',
            title='Human sample from Homo sapiens')))})

        names = fetchBiosamples(client, ['SAMN15865971'])

        self.assertEqual(names['SAMN15865971'], '32_CSF_uns_5pIGSEQ_2')

    def test_falls_back_to_title_when_no_sample_name_is_set(self):
        """A record with no Sample name Id falls back to its Title."""
        client = makeClient({'efetch': FakeResponse(200, efetchXml(biosampleXml(
            'SAMN36877371', title='BL-110_VDJ')))})

        names = fetchBiosamples(client, ['SAMN36877371'])

        self.assertEqual(names['SAMN36877371'], 'BL-110_VDJ')


class TestGatherEvidence(unittest.TestCase):
    """
    Tests for the combined accession -> Evidence pipeline
    """

    def test_ok_run_carries_biosample_and_url(self):
        """A cleanly resolved run reports its BioSample, name and a check link."""
        client = makeClient({
            'esearch': FakeResponse(200, esearchXml('1')),
            'esummary': FakeResponse(200, esummaryXml(docsumXml(
                '1', 'GSM7688951: BL-110_VDJ', 'SAMN36877371', 'SRR25557617'))),
            'efetch': FakeResponse(200, efetchXml(
                biosampleXml('SAMN36877371', title='BL-110_VDJ'))),
        })

        evidence = gatherEvidence(client, ['SRR25557617'])['SRR25557617']

        self.assertEqual(evidence, Evidence(
            accession='SRR25557617', status='ok',
            biosample_accession='SAMN36877371', sample_name='BL-110_VDJ',
            url='https://www.ncbi.nlm.nih.gov/biosample/SAMN36877371'))
        self.assertEqual(evidence.suggested_subject, 'BL-110')

    def test_pooled_run_is_flagged_rather_than_resolved(self):
        """A multi-donor hashed run is marked pooled, not guessed at."""
        client = makeClient({
            'esearch': FakeResponse(200, esearchXml('1')),
            'esummary': FakeResponse(200, esummaryXml(docsumXml(
                '1', 'GSM6339422: Hashed scBCR sample (FA007, FA048)',
                'SAMN29756173', 'SRR20210850'))),
            'efetch': FakeResponse(200, efetchXml(biosampleXml(
                'SAMN29756173', title='Hashed scBCR sample (FA007, FA048)'))),
        })

        evidence = gatherEvidence(client, ['SRR20210850'])['SRR20210850']

        self.assertEqual(evidence.status, 'pooled')
        self.assertEqual(evidence.pooled_codes, ('FA007', 'FA048'))
        self.assertEqual(evidence.suggested_subject, '')

    def test_unresolved_accession_still_gets_an_entry(self):
        """
        An accession absent from NCBI still appears in the result, so callers
        never have to guard a missing key.
        """
        client = makeClient({'esearch': FakeResponse(200, esearchXml())})

        evidence = gatherEvidence(client, ['SRR99999999'])['SRR99999999']

        self.assertEqual(evidence.status, 'not_found')
        self.assertEqual(evidence.biosample_accession, '')
        self.assertEqual(evidence.url, '')

    def test_falls_back_to_the_sra_description_with_no_biosample_text_at_all(self):
        """
        A record with neither a Sample name Id nor any Title text still gets
        a usable sample_name, from the SRA experiment's own description --
        the same text a person would see on the run's trace page without
        ever following the BioSample link at all. Note this is a narrower
        fallback than "the Title happens to be generic": fetchBiosamples
        already prefers Sample name over Title (see its own tests), so any
        BioSample record contributing *some* text, however boilerplate, wins
        over this description.
        """
        client = makeClient({
            'esearch': FakeResponse(200, esearchXml('1')),
            'esummary': FakeResponse(200, esummaryXml(docsumXml(
                '1', 'IG-Seq of CSF unsorted cells from sample 32, 5prime replicate',
                'SAMN15865971', 'SRR12483427'))),
            'efetch': FakeResponse(200, efetchXml(biosampleXml(
                'SAMN15865971', title=''))),
        })

        evidence = gatherEvidence(client, ['SRR12483427'])['SRR12483427']

        self.assertEqual(evidence.status, 'ok')
        self.assertEqual(evidence.sample_name,
                         'IG-Seq of CSF unsorted cells from sample 32, 5prime replicate')


if __name__ == '__main__':
    unittest.main()
