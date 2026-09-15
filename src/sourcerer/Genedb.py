"""
IMGT/GENE-DB historical releases (JamieHeather/genedb-releases)

IMGT's GENElect serves only the current release, so a reference pinned to an
older release cannot be re-fetched from IMGT itself. Jamie Heather's
genedb-releases repository fills that gap: it archives GENE-DB release by
release, each a directory named ``<access-date>_<source>_<release>`` holding the
bulk ReferenceSequences FASTAs, which carry the same pipe-delimited IMGT headers
as GENElect.

This module is what ``download --from`` uses to reconstruct an IMGT release. It
resolves a release tag to the archived directory -- exactly when the archive
holds it, or the closest release with a warning when it does not -- fetches the
bulk FASTA, and filters it down to one species' reference chains. Because IMGT
serves only the current build, a pinned mouse reference in particular is a
best-effort reconstruction; ``download --compare`` against the original is how a
caller confirms it matches.
"""

# Info
__author__ = 'Ayelet Peres'

# Imports
import logging
from urllib.parse import quote

# Sourcerer imports
from sourcerer.Reference import (
    KIND_AA,
    KIND_CONSTANT,
    KIND_VDJ,
    KNOWN_CHAINS,
    parseFasta,
)

log = logging.getLogger(__name__)

#: The GitHub repository that archives the releases, and how to read it: the
#: contents API lists the release directories, the raw host serves their files.
CONTENTS_API = ('https://api.github.com/repos/JamieHeather/genedb-releases/'
                'contents/releases')
RAW_ROOT = ('https://raw.githubusercontent.com/JamieHeather/genedb-releases/'
            'main/releases/')

#: The bulk files taken from each release directory: the gapped nucleotide set,
#: whose V keeps its IMGT numbering, drives vdj and constant, and the gapped
#: amino acid set drives the translated V. Both are the F+ORF+inframeP
#: functionality set, the closest archived equivalent of a working reference.
BULK_NT = 'IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP'
BULK_AA = 'IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F+ORF+inframeP'

#: The bulk file for each unit group.
GROUPS = {'nt': BULK_NT, 'aa': BULK_AA}

#: Species as they appear in the IMGT headers. A header may carry a strain suffix
#: (``Mus musculus_C57BL/6``), so a species is matched on this binomial prefix.
SPECIES_LABEL = {'human': 'Homo sapiens', 'mouse': 'Mus musculus'}

#: A constant region under this many nucleotides is the delta diversity segment
#: sharing the IGHD name, told apart from the delta constant by length, exactly
#: as the OGRDB source does.
CONSTANT_MIN_LENGTH = 100


def releaseKey(tag):
    """
    Turn a release tag into a sortable key.

    A GENE-DB release is ``YYYYWW-R`` (year, week, revision), e.g. '202631-7',
    which orders by date once each part is read as a number.

    Arguments:
      tag (str): the release tag.

    Returns:
      tuple: (build, revision) integers, or (0, 0) for an unparseable tag.
    """
    try:
        build, revision = tag.split('-')
        return int(build), int(revision)
    except (ValueError, AttributeError):
        return (0, 0)


def listReleases(client):
    """
    List every archived release, newest last.

    Arguments:
      client (HttpClient): the shared HTTP client.

    Returns:
      list: (release_tag, dirname) tuples, sorted by release.
    """
    entries = client.get(CONTENTS_API).json()
    releases = []
    for entry in entries:
        if entry.get('type') != 'dir':
            continue
        name = entry['name']
        # <access-date>_<source>_<release>; the release is the last token.
        tag = name.rsplit('_', 1)[-1]
        releases.append((tag, name))

    return sorted(releases, key=lambda item: releaseKey(item[0]))


def resolveRelease(client, wanted):
    """
    Find the archived directory for a release, or the closest one.

    IMGT builds are not all archived, so an exact match is not guaranteed; when
    the wanted release is missing the nearest by release number is used and the
    substitution is logged, loud enough that a caller notices the reference is a
    neighbour rather than the one asked for.

    Arguments:
      client (HttpClient): the shared HTTP client.
      wanted (str): the release tag to reconstruct, e.g. '202631-7'.

    Returns:
      tuple: (dirname, resolved_tag, exact) where exact says the archive held the
      wanted release.
    """
    releases = listReleases(client)
    if not releases:
        raise RuntimeError('genedb-releases listed no archived releases')

    for tag, name in releases:
        if tag == wanted:
            return name, tag, True

    target = releaseKey(wanted)
    tag, name = min(releases, key=lambda item: abs(releaseKey(item[0])[0] - target[0]))
    log.warning('IMGT release %s is not archived; using the closest, %s. '
                'Compare against the original reference to confirm it matches.',
                wanted, tag)

    return name, tag, False


def bulkUrl(dirname, group):
    """
    Build the raw URL of a release directory's bulk file for a group.

    Arguments:
      dirname (str): the release directory name.
      group (str): 'nt' or 'aa'.

    Returns:
      str: the absolute URL.
    """
    filename = GROUPS[group]

    return '%s%s/%s' % (RAW_ROOT, quote(dirname, safe=''), quote(filename, safe=''))


def _classify(allele, sequence, is_aa):
    """
    Decide which reference chain an IMGT allele belongs to.

    Arguments:
      allele (str): the allele name, e.g. 'IGHV1-2*02'.
      sequence (str): its sequence, used only to tell the delta segments apart.
      is_aa (bool): whether this came from the amino acid set.

    Returns:
      tuple: (chain, kind), or None for a name too short or an amino acid non-V.
    """
    if len(allele) < 4:
        return None

    segment = allele[3]
    if is_aa:
        return (allele[:4], KIND_AA) if segment == 'V' else None
    if segment in ('V', 'J'):
        return allele[:4], KIND_VDJ
    if segment == 'D':
        if len(sequence.replace('.', '')) < CONSTANT_MIN_LENGTH:
            return allele[:4], KIND_VDJ
        return allele[:3] + 'C', KIND_CONSTANT

    return allele[:3] + 'C', KIND_CONSTANT


def selectChains(text, species, is_aa):
    """
    Filter a bulk IMGT FASTA to one species' reference chains.

    Only records for the species are kept, matched on the binomial so every
    archived strain of the mouse counts. An allele name seen more than once, as
    it is when strains share a gene, is kept from its first occurrence, so the
    written reference has one sequence per allele. The header is left as IMGT
    wrote it, with the species binomial underscored to match the live source.

    Arguments:
      text (str): the bulk FASTA.
      species (str): the species key, e.g. 'human'.
      is_aa (bool): whether the text is the amino acid set.

    Returns:
      dict: (chain, kind) to a list of (header, sequence) tuples.
    """
    label = SPECIES_LABEL[species]
    underscored = label.replace(' ', '_')

    chains, seen = {}, {}
    for header, sequence in parseFasta(text):
        fields = header.split('|')
        if len(fields) < 3:
            continue
        allele = fields[1].strip()
        source_species = fields[2].strip()
        if not (source_species == label or source_species.startswith(label + '_')):
            continue

        target = _classify(allele, sequence, is_aa)
        # The bulk file carries genes beyond IG and TR (MICA, VPREB, ...); the
        # reference is only the immunoglobulin and T-cell receptor chains, the
        # same set the live GENElect path fetches.
        if target is None or target[0] not in KNOWN_CHAINS:
            continue

        taken = seen.setdefault(target, set())
        if allele in taken:
            continue
        taken.add(allele)
        chains.setdefault(target, []).append(
            (header.replace(label, underscored), sequence.upper()))

    return chains
