"""
NCBI cross-reference

Some OAS studies do not record a Subject in their own metadata (`sourcerer oas
download` passes that through raw, as the OAS null sentinel 'no', rather than
guessing — see Airrflow.buildSamplesheet). Every such run still has an SRA run
accession or GEO sample accession embedded in its OAS unit id, and that
accession's BioSample record on NCBI usually names the sample plainly, e.g.
BioSample SAMN36877371 for run SRR25557617 gives the sample name 'BL-110_VDJ'.

This module is the deterministic half of closing that gap: given a batch of
accessions, look up each one's BioSample and return its raw sample name plus a
link a person can open to check the record themselves. It does not decide what
part of that raw text is the subject's identity — 'BL-110_VDJ' strips to
BL-110 in one study while 'TT04_subj6_V3' strips to TT04_subj6 in another, and
telling those apart needs to know the specific study's naming convention, not
just read the string. suggestSubject applies the handful of patterns generic
enough to be safe across studies; everything else is left for
`sourcerer oas verify`'s evidence report to surface for a human to decide.
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import logging
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from html import unescape
from urllib.parse import urlencode

log = logging.getLogger(__name__)

#: E-utilities base URL. See https://www.ncbi.nlm.nih.gov/books/NBK25501/
EUTILS_BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'

#: Where a person can check a BioSample record by hand, the same page
#: `sourcerer oas verify`'s evidence links point at.
BIOSAMPLE_URL = 'https://www.ncbi.nlm.nih.gov/biosample/%s'

#: Run/sample accession formats OAS unit ids embed: SRA runs (SRR/ERR/DRR) and
#: GEO samples (GSM), the latter for studies submitted through GEO rather than
#: directly to SRA. No trailing \b: OAS unit ids run the accession straight
#: into a '_1_Paired_All...' suffix, and \d is a word character, so a
#: boundary would never fire between the last digit and that underscore.
#: The digits themselves are exactly what \d+ stops matching on, so leaving
#: it off costs nothing.
ACCESSION_RE = re.compile(r'\b(SRR\d+|ERR\d+|DRR\d+|GSM\d+)')

#: NCBI's politeness guideline without a key is 3 requests/second; an API key
#: (see https://www.ncbi.nlm.nih.gov/account/settings/) raises that to 10.
#: Not used inside this module: HttpClient already enforces a minimum delay
#: between requests at the one seam every call goes through (see Http.py's
#: module docstring), so the caller passes one of these to HttpClient's own
#: delay argument when building the client rather than this module pacing
#: requests a second time on top of it.
DEFAULT_DELAY = 0.34
KEYED_DELAY = 0.11

#: How many accessions/ids to fold into one esearch OR-query or one
#: esummary/efetch batch. Comfortably under any URL length limit while still
#: cutting hundreds of runs down to a handful of requests.
CHUNK_SIZE = 40

#: Sample descriptions naming more than one donor in a single 10x cell-hashing
#: or multiplexed run. A BioSample record like this genuinely does not
#: identify one subject; suggestSubject and gatherEvidence both refuse to
#: guess which one instead of silently picking a donor.
POOLED_RE = re.compile(r'\b(hash(?:ed)?|pool(?:ed)?|multiplex(?:ed)?)\b',
                       re.IGNORECASE)

#: Trailing tokens generic enough, across studies, to be safe to strip when
#: suggesting a subject id: assay/locus markers and visit/replicate numbers.
#: Nothing study-specific (e.g. a particular cohort's prefix convention)
#: belongs here — see the module docstring for why.
#: The 5'/5-prime alternatives require an explicit GEX/prime marker rather
#: than accepting a bare trailing '5' on its own — a plain digit suffix (e.g.
#: 'Donor_5', a replicate or cohort number) is exactly the kind of study
#: specific subject-identifying detail this function must never strip.
_LOCUS_SUFFIX_RE = re.compile(
    r'''[_,-]\s*(?:
        VDJ | VJ | TCR | BCR | IGH | IGK | IGL | GEX
        | 5[\'′](?:[- ]?(?:GEX|prime))?
        | 5[- ]?(?:GEX|prime)
    )\s*$''',
    re.IGNORECASE | re.VERBOSE)
_VISIT_SUFFIX_RE = re.compile(r'[_-][Vv]\d+$')


@dataclass(frozen=True)
class Evidence:
    """
    What NCBI says about one run/sample accession.

    Arguments:
      accession (str): the SRR/ERR/DRR/GSM accession looked up.
      status (str): 'ok' if a BioSample sample name was found, 'pooled' if it
        names more than one donor, 'not_found' if SRA has no record of the
        accession at all.
      biosample_accession (str): the SAMN accession, or '' if status is
        'not_found'.
      sample_name (str): the raw text NCBI associates with the sample: the
        BioSample's own 'Sample name' identifier when the submitter set one,
        its Title otherwise, or (rarely) the SRA experiment's own title as a
        last resort. Never normalized — see suggestSubject for that.
      url (str): a BioSample page a person can open to check this by hand, or
        '' if status is 'not_found'.
      pooled_codes (tuple): the donor codes named in sample_name, if status is
        'pooled'; empty otherwise.
    """
    accession: str
    status: str = 'not_found'
    biosample_accession: str = ''
    sample_name: str = ''
    url: str = ''
    pooled_codes: tuple = ()

    @property
    def suggested_subject(self):
        """str: suggestSubject(sample_name), or '' if there is nothing to suggest."""
        if self.status != 'ok':
            return ''
        return suggestSubject(self.sample_name)


def accessionFromText(text):
    """
    Pull the first SRR/ERR/DRR/GSM accession out of free text.

    Meant for OAS's unit id (e.g.
    'Corinaldesi_2024/csv_paired/SRR25557617_1_Paired_All.csv.gz'), which is
    what Airrflow.buildSamplesheet records into the sample_name column.

    Arguments:
      text (str): text to search.

    Returns:
      str: the accession, or None if none was found.
    """
    match = ACCESSION_RE.search(text or '')
    return match.group(1) if match else None


def poolCodes(text):
    """
    Read off the donor codes named in a pooled/hashed sample's description.

    Arguments:
      text (str): a BioSample sample name or SRA title, e.g.
        'Hashed scBCR sample (FA007, FA048)'.

    Returns:
      tuple: the codes found, e.g. ('FA007', 'FA048'); empty if the text does
        not read as a pooled sample or names only one code.
    """
    if not POOLED_RE.search(text or ''):
        return ()

    match = re.search(r'\(([^)]+)\)', text or '')
    if not match:
        return ()

    codes = tuple(code.strip() for code in match.group(1).split(',') if code.strip())
    return codes if len(codes) > 1 else ()


def suggestSubject(text):
    """
    Best-effort, study-agnostic guess at the subject-identifying part of a
    BioSample sample name.

    Strips only patterns generic enough to be safe regardless of which
    study's naming convention produced the text: a trailing assay/locus token
    ('BL-110_VDJ' -> 'BL-110') or a trailing visit marker
    ('TT04_subj6_V3' -> 'TT04_subj6'). It does not attempt to also collapse a
    leading numeric or prefixed subject code down to a canonical form, because
    that requires knowing the study: 'P05_FNA_d0_1_Y1' and
    '32_CSF_uns_5pIGSEQ_2' both start with what should become the subject id,
    but where that field ends is a convention only the specific submission
    follows, not a pattern this function can read off the string alone.

    This is a hint for `sourcerer oas verify`'s evidence report, not a value
    ever written automatically into biosample_id; see the module docstring.

    Arguments:
      text (str): a raw BioSample sample name.

    Returns:
      str: the suggestion, or the input stripped, unchanged, if no generic
        pattern matched.
    """
    text = (text or '').strip()
    stripped = _LOCUS_SUFFIX_RE.sub('', text)
    stripped = _VISIT_SUFFIX_RE.sub('', stripped)
    return stripped or text


def _chunks(items, size):
    """Yield items in slices of at most size, preserving order."""
    for start in range(0, len(items), size):
        yield items[start:start + size]


def _eutilsUrl(endpoint, **params):
    """Build a full E-utilities URL, params baked into the query string."""
    return '%s/%s.fcgi?%s' % (EUTILS_BASE, endpoint, urlencode(params))


def resolveAccessions(client, accessions, api_key=None, chunk_size=CHUNK_SIZE):
    """
    Look up SRA run/sample accessions and return each one's BioSample.

    Batches esearch (accessions OR-joined, chunked) and esummary (the
    resulting UIDs, comma joined), then confirms every docsum against the
    accessions actually in its batch before trusting it, rather than relying
    on NCBI to return results in request order. In practice E-utilities does
    preserve order for an OR'd term search, but that is not documented
    behavior worth depending on for identity: a docsum is attributed to an
    accession only when that accession appears as the docsum's own Run acc
    (SRA-native submissions) or as a GSM token in its own title (GEO-mediated
    submissions), never by position.

    Pacing between requests is not this function's concern: it is whatever
    the client itself was built with (see HttpClient's own delay argument,
    and the module docstring's note on DEFAULT_DELAY/KEYED_DELAY) --
    HttpClient already enforces politeness once, at the one seam every
    request goes through, and a second delay layered on top here would just
    slow every run down for nothing.

    Arguments:
      client (HttpClient): the shared HTTP client.
      accessions (list): SRR/ERR/DRR/GSM accessions; deduplicated internally.
      api_key (str): an NCBI API key, if available (see module docstring).
      chunk_size (int): accessions per esearch/esummary batch.

    Returns:
      dict: accession to (biosample_accession, description), description
        being the free text between the accession and '; Homo sapiens...' in
        the SRA experiment title, e.g. 'BL-110_VDJ'. Accessions with no SRA
        hit are absent from the result; the caller decides how to report that.
    """
    key_param = {'api_key': api_key} if api_key else {}

    resolved = {}
    accessions = sorted(set(a for a in accessions if a))

    for batch in _chunks(accessions, chunk_size):
        term = ' OR '.join(batch)
        response = client.get(_eutilsUrl('esearch', db='sra', term=term,
                                         retmax=500, **key_param))

        ids = [node.text for node in ET.fromstring(response.text).findall('.//IdList/Id')]
        if not ids:
            continue

        for id_batch in _chunks(ids, 100):
            response = client.get(_eutilsUrl('esummary', db='sra',
                                             id=','.join(id_batch), **key_param))

            for docsum in ET.fromstring(response.text).findall('.//DocSum'):
                exp_xml, runs_xml = '', ''
                for item in docsum.findall('Item'):
                    if item.get('Name') == 'ExpXml':
                        exp_xml = unescape(item.text or '')
                    elif item.get('Name') == 'Runs':
                        runs_xml = unescape(item.text or '')

                title_match = re.search(r'<Title>(.*?)</Title>', exp_xml, re.DOTALL)
                title = title_match.group(1) if title_match else ''
                biosample_match = re.search(r'<Biosample>(SAM[A-Z]\w+)</Biosample>',
                                            exp_xml)
                biosample_accession = biosample_match.group(1) if biosample_match else ''
                run_accessions = set(re.findall(r'Run acc="([^"]+)"', runs_xml))

                matched = run_accessions & set(batch)
                gsm_match = re.search(r'\b(GSM\d+)\b', title)
                if gsm_match and gsm_match.group(1) in batch:
                    matched.add(gsm_match.group(1))

                description = re.sub(r'^GSM\d+:\s*', '', title)
                description = re.sub(r';\s*Homo sapiens.*$', '', description).strip()

                for accession in matched:
                    resolved[accession] = (biosample_accession, description)

    return resolved


def fetchBiosamples(client, biosample_accessions, api_key=None, chunk_size=100):
    """
    Fetch the sample name NCBI shows for a batch of BioSample accessions.

    Prefers the 'Sample name' identifier the BioSample page displays when the
    submitter set one; falls back to the record's own Description/Title
    otherwise. This is the field the OAS submission for a run's BioSample
    shows as its sample name when you open the page by hand.

    Arguments:
      client (HttpClient): the shared HTTP client.
      biosample_accessions (list): SAMN/SAME/SAMD accessions; deduplicated
        internally.
      api_key (str): an NCBI API key, if available.
      chunk_size (int): accessions per efetch batch.

    Returns:
      dict: biosample_accession to sample name (str; '' if the record carries
        neither a Sample name nor a Title, which efetch itself would have to
        be broken for).
    """
    key_param = {'api_key': api_key} if api_key else {}

    names = {}
    accessions = sorted(set(a for a in biosample_accessions if a))

    for batch in _chunks(accessions, chunk_size):
        response = client.get(_eutilsUrl('efetch', db='biosample',
                                         id=','.join(batch), rettype='full',
                                         retmode='xml', **key_param))

        for sample in ET.fromstring(response.text).findall('.//BioSample'):
            accession = sample.get('accession')
            sample_name = ''
            for id_node in sample.findall('./Ids/Id'):
                if id_node.get('db_label') == 'Sample name':
                    sample_name = (id_node.text or '').strip()

            if not sample_name:
                sample_name = (sample.findtext('./Description/Title') or '').strip()

            names[accession] = sample_name

    return names


def gatherEvidence(client, accessions, api_key=None):
    """
    Resolve a batch of run/sample accessions to NCBI evidence in one pass.

    Combines resolveAccessions and fetchBiosamples, then classifies each
    result: a sample name that reads as a multi-donor pool (see poolCodes)
    gets status 'pooled' rather than a guessed single subject.

    Arguments:
      client (HttpClient): the shared HTTP client.
      accessions (list): SRR/ERR/DRR/GSM accessions; deduplicated internally.
      api_key (str): an NCBI API key, if available.

    Returns:
      dict: accession to Evidence, one entry per input accession (including
        ones NCBI could not resolve, so callers never have to guard a missing
        key).
    """
    accessions = sorted(set(a for a in accessions if a))
    resolved = resolveAccessions(client, accessions, api_key=api_key)

    biosample_accessions = [biosample for biosample, _ in resolved.values()]
    names = fetchBiosamples(client, biosample_accessions, api_key=api_key)

    evidence = {}
    for accession in accessions:
        if accession not in resolved:
            evidence[accession] = Evidence(accession=accession, status='not_found')
            continue

        biosample_accession, description = resolved[accession]
        # fetchBiosamples already prefers the Sample name Id over a
        # BioSample's own Title (see its docstring), so this only falls back
        # further when the record contributed no usable text at all -- no
        # Sample name, no Title either, or efetch simply never returned that
        # accession. The SRA experiment's own description is what is left:
        # the same text a person would see on the run's trace page without
        # ever following the BioSample link at all.
        sample_name = names.get(biosample_accession) or description

        codes = poolCodes(sample_name)
        evidence[accession] = Evidence(
            accession=accession,
            status='pooled' if codes else 'ok',
            biosample_accession=biosample_accession,
            sample_name=sample_name,
            url=BIOSAMPLE_URL % biosample_accession if biosample_accession else '',
            pooled_codes=codes)

    return evidence
