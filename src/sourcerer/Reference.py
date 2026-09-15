"""
Germline reference output

The airrflow artifact for a germline source is not a samplesheet: it is a
reference tree that IgBLAST and Change-O read. This module is to the reference
sources what Airrflow.py is to the dataset sources -- the one place that knows
the shape of the output nf-core/airrflow expects, kept out of the sources
themselves so a second reference source inherits it unchanged.

This module holds the builders that turn downloaded germline FASTAs into the two
directory layouts airrflow consumes: a reference_base of per-chain FASTAs, and --
only when asked, because it needs the BLAST+ binary -- an igblast_base of BLAST
databases plus the internal_data and optional_file trees mirrored from NCBI. The
base class germline sources extend, ReferenceSource, lives in
sourcerer.Sources.Germline, kept there rather than here so this module never
imports from sourcerer.Sources and the two stay free of an import cycle.

The reference_base keeps the source's own FASTA verbatim, gaps and all, exactly
as airrflow's bin/fetch_references.sh leaves it. Cleaning (gap removal, dedup)
happens only when the BLAST database is built, so nothing that reads the
reference for its IMGT numbering, such as Change-O's germline reconstruction,
loses it.
"""

# Info
__author__ = 'Ayelet Peres'

# Imports
import hashlib
import logging
import os
import shutil
import subprocess
import tarfile
import textwrap
from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import urljoin, urlparse

import yaml
from bs4 import BeautifulSoup

# Sourcerer imports
from sourcerer.Exceptions import SourcererError

log = logging.getLogger(__name__)

#: Species airrflow builds references for, and the leading directory in the
#: reference tree. New species are added here and in each source's SETS/CHAINS.
SPECIES = ('human', 'mouse')

#: Receptor classes, matching airrflow's canonical database basenames.
LOCI = ('ig', 'tr')

#: Gene segments, in the order IgBLAST names its databases.
SEGMENTS = ('v', 'd', 'j', 'c')

#: The chains that make up each canonical (locus, segment) database. This is the
#: aggregation airrflow's ref2igblast.sh performs: one BLAST database per class
#: and segment, built from every locus in that class.
LOCUS_CHAINS = {
    ('ig', 'v'): ('IGHV', 'IGKV', 'IGLV'),
    ('ig', 'd'): ('IGHD',),
    ('ig', 'j'): ('IGHJ', 'IGKJ', 'IGLJ'),
    ('ig', 'c'): ('IGHC', 'IGKC', 'IGLC'),
    ('tr', 'v'): ('TRAV', 'TRBV', 'TRDV', 'TRGV'),
    ('tr', 'd'): ('TRBD', 'TRDD'),
    ('tr', 'j'): ('TRAJ', 'TRBJ', 'TRDJ', 'TRGJ'),
    ('tr', 'c'): ('TRAC', 'TRBC', 'TRDC', 'TRGC'),
}

#: Every chain a reference FASTA may be named for, flattened from LOCUS_CHAINS.
KNOWN_CHAINS = frozenset(chain for chains in LOCUS_CHAINS.values()
                         for chain in chains)

#: Which reference_base subdirectory a chain's FASTA lives in. Constant regions
#: are kept apart from V/D/J because airrflow's tree does, and amino acid V has
#: its own directory because it becomes a protein database rather than a
#: nucleotide one.
KIND_VDJ = 'vdj'
KIND_CONSTANT = 'constant'
KIND_AA = 'vdj_aa'

#: The longest sequence identifier ``makeblastdb -parse_seqids`` accepts. A
#: longer name is not a warning: makeblastdb refuses the whole database with
#: "the local id is too long". Dropping -parse_seqids would lift the limit, but
#: igblastn then reports gnl|BL_ORD_ID|0 instead of the gene name, which MakeDb
#: cannot parse, so the name has to give instead.
BLAST_NAME_LIMIT = 50

#: Where a build records the names it had to shorten.
RENAME_LOG = 'shortened_alleles.tsv'

#: Provenance sidecars written inside reference_base, so they travel with the
#: reference when airrflow moves the tree. IMGT records the release the FASTAs
#: came from; AIRRC records each OGRDB set's version. `download --from <ref>`
#: reads them back to re-download the same versions.
IMGT_METADATA = 'IMGT.yaml'
AIRRC_METADATA = 'AIRRC.yaml'

#: NCBI's IgBLAST release trees, mirrored into igblast_base so that igblastn has
#: the auxiliary data it cannot derive from the germline FASTAs alone. The
#: old_* directories are the layout airrflow's fetch_igblastdb.sh already tracks.
NCBI_IGBLAST_ROOT = ('https://ftp.ncbi.nlm.nih.gov/blast/executables/igblast/'
                     'release/')
NCBI_DATABASE_URL = urljoin(NCBI_IGBLAST_ROOT, 'database/')
NCBI_INTERNAL_URL = urljoin(NCBI_IGBLAST_ROOT, 'old_internal_data/')
NCBI_OPTIONAL_URL = urljoin(NCBI_IGBLAST_ROOT, 'old_optional_file/')

#: Archives NCBI ships inside database/ that have to be unpacked in place for the
#: mirrored tree to be usable.
NCBI_TAR_ARCHIVES = ('mouse_gl_VDJ.tar', 'rhesus_monkey_VJ.tar')


@dataclass
class ReferenceReport:
    """
    Summary of a reference build, in the same spirit as OAS's conversion report.

    Arguments:
      written (list): reference_base FASTAs written, as (chain, path) or basename.
      built (list): canonical BLAST database basenames created.
      renamed (dict): original allele name to the shortened one, for names too
        long for makeblastdb.
      aux_missing (dict): species to the J alleles the mirrored NCBI auxiliary
        file does not name.
      skipped_empty (list): canonical databases skipped because no chain in them
        had any sequence. This is the normal outcome for what a source does not
        cover, such as TR from OGRDB.
    """
    written: list = field(default_factory=list)
    built: list = field(default_factory=list)
    skipped_empty: list = field(default_factory=list)
    renamed: dict = field(default_factory=dict)
    aux_missing: dict = field(default_factory=dict)

    def logSummary(self):
        """Log a one-line summary of what the build produced."""
        log.info('reference: %d files written, %d databases built, %d skipped',
                 len(self.written), len(self.built), len(self.skipped_empty))


@dataclass
class ReferencePlan:
    """
    What building an IgBLAST base from a reference folder would produce.

    Computed without running makeblastdb, so it doubles as the format check: it
    says which databases would build, which come up empty, which files were not
    recognised, and where duplicate allele names were dropped.

    Arguments:
      found_species (list): species seen in the folder.
      databases (list): (basename, dbtype, records) that will build.
      empty (list): canonical basenames with no sequence to build from.
      unrecognized (list): FASTA paths whose names are not in the reference format.
      empty_files (list): recognised FASTA paths that held no sequence.
      duplicates (dict): basename to the number of duplicate names dropped.
      renamed (dict): original allele name to the shortened one, for every name
        too long for makeblastdb.
    """
    found_species: list = field(default_factory=list)
    databases: list = field(default_factory=list)
    empty: list = field(default_factory=list)
    unrecognized: list = field(default_factory=list)
    empty_files: list = field(default_factory=list)
    duplicates: dict = field(default_factory=dict)
    renamed: dict = field(default_factory=dict)

    @property
    def ok(self):
        """bool: True if at least one database can be built."""
        return bool(self.databases)

    def summary(self):
        """
        Render the plan as a human-readable report.

        Returns:
          str: the report.
        """
        lines = ['species found: %s' % (', '.join(self.found_species) or 'none')]
        if self.databases:
            lines.append('databases to build (%d):' % len(self.databases))
            for basename, dbtype, records in self.databases:
                dropped = self.duplicates.get(basename)
                note = '  (%d duplicate name(s) dropped)' % dropped if dropped else ''
                lines.append('  %-16s %5d seq  %s%s'
                             % (basename, len(records), dbtype, note))
        if self.renamed:
            lines.append('%d allele name(s) shortened to the %d-character BLAST '
                         'limit' % (len(self.renamed), BLAST_NAME_LIMIT))
        if self.empty:
            lines.append('empty, nothing to build: %s' % ', '.join(self.empty))
        for path in self.empty_files:
            lines.append('warning: %s held no sequence' % path)
        for path in self.unrecognized:
            lines.append('warning: %s is not in the reference naming format, '
                         'skipped' % path.name)

        return '\n'.join(lines)


# ---------------------------------------------------------------------------
# FASTA helpers
# ---------------------------------------------------------------------------

def parseFasta(text):
    """
    Read FASTA text into (header, sequence) pairs, in source order.

    Whitespace inside a sequence is collapsed; the header keeps everything after
    the '>' verbatim, because a source's own header, IMGT's pipe-delimited line
    for one and OGRDB's allele name for another, is what identifies the allele.

    Arguments:
      text (str): FASTA text.

    Returns:
      list: (header, sequence) tuples.
    """
    records = []
    header, seq = None, []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if header is not None:
                records.append((header, ''.join(seq)))
            header = line[1:]
            seq = []
        else:
            seq.append(''.join(line.split()))
    if header is not None:
        records.append((header, ''.join(seq)))

    return records


def alleleName(header):
    """
    Take the allele name from a FASTA header.

    IMGT headers are pipe-delimited and put the name in the second field
    (``>X02897|IGHV1-2*02|Homo sapiens|F|...``); OGRDB writes the bare name. The
    first token is used when there is no pipe, so both are handled by one rule.

    Arguments:
      header (str): the header line without its leading '>'.

    Returns:
      str: the allele name.
    """
    if '|' in header:
        return header.split('|')[1].strip()

    return header.split()[0]


def writeFastaText(path, records):
    """
    Write (header, sequence) pairs to a FASTA file verbatim, one line per part.

    Arguments:
      path (Path): output path. Parent directories are created.
      records (iterable): (header, sequence) tuples.

    Returns:
      int: the number of records written.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    written = 0
    with open(path, 'w') as handle:
        for header, sequence in records:
            handle.write('>%s\n%s\n' % (header, sequence))
            written += 1

    return written


def shortenName(name):
    """
    Fit an allele name inside the BLAST identifier limit.

    The head of the name is kept and a short digest of the whole name replaces
    the tail, so the gene stays readable and two names that share a head still
    differ. The digest is taken over the original name rather than the discarded
    tail, so the same allele always shortens to the same thing whichever database
    it is built into -- which matters because IgBLAST's auxiliary and delineation
    files are keyed by name.

    Only VDJbase-style novel allele names reach the limit in practice
    (``IGHV3-20*01_a123g_t456c_...``); nothing IMGT or OGRDB publishes is close.

    Six hex characters of digest keeps the odds of two different names
    colliding negligible while leaving the readable head almost intact.

    Arguments:
      name (str): the allele name.

    Returns:
      str: the name, shortened only if it had to be.
    """
    if len(name) <= BLAST_NAME_LIMIT:
        return name

    digest = hashlib.sha1(name.encode('utf-8')).hexdigest()[:6]

    return '%s_%s' % (name[:BLAST_NAME_LIMIT - 7], digest)


def shortenForBlast(records):
    """
    Shorten any allele name makeblastdb would refuse, and say which.

    Arguments:
      records (list): (name, sequence) tuples from cleanForBlast.

    Returns:
      tuple: (records, renamed) where renamed maps original name to the
      shortened one, empty when nothing had to change.
    """
    renamed = {}
    shortened = []
    for name, sequence in records:
        short = shortenName(name)
        if short != name:
            renamed[name] = short
        shortened.append((short, sequence))

    return shortened, renamed


def cleanForBlast(records):
    """
    Prepare germline records for makeblastdb.

    Gaps are removed, sequences upper-cased and duplicate names dropped, keeping
    the first. Deduplication is not cosmetic: makeblastdb -parse_seqids refuses a
    database with a repeated identifier, so a duplicate allele name is a hard
    failure rather than a warning. The name is taken with alleleName so an IMGT
    pipe header collapses to just the allele, which is what IgBLAST reports.

    Arguments:
      records (iterable): (header, sequence) tuples.

    Returns:
      list: (name, sequence) tuples, cleaned and de-duplicated.
    """
    seen = set()
    cleaned = []
    for header, sequence in records:
        name = alleleName(header)
        if name in seen:
            continue
        seen.add(name)
        cleaned.append((name, sequence.replace('.', '').upper()))

    return cleaned


# ---------------------------------------------------------------------------
# reference_base
# ---------------------------------------------------------------------------

def referenceFastaPath(reference_dir, prefix, species, kind, chain):
    """
    Locate one chain's FASTA in the reference tree.

    The layout matches airrflow's: ``<species>/<kind>/<prefix>_<species>_<chain>``,
    with amino acid V spelled ``<prefix>_aa_<species>_<chain>`` so a nucleotide
    and a protein file for the same chain do not collide.

    Arguments:
      reference_dir (Path): the reference_base root.
      prefix (str): the source tag, 'imgt' or 'airrc'.
      species (str): the species.
      kind (str): the subdirectory, one of KIND_VDJ, KIND_CONSTANT, KIND_AA.
      chain (str): the chain, e.g. 'IGHV'.

    Returns:
      Path: where the chain's FASTA belongs.
    """
    if kind == KIND_AA:
        name = '%s_aa_%s_%s.fasta' % (prefix, species, chain)
    else:
        name = '%s_%s_%s.fasta' % (prefix, species, chain)

    return Path(reference_dir) / species / kind / name


#: Filename extensions read as germline FASTA. A reference written by sourcerer
#: is always .fasta, but a folder a user assembled themselves may not be, and a
#: file the manifest names is no use if it is never looked at.
FASTA_SUFFIXES = ('.fasta', '.fa', '.fna')


def parseReferenceName(filename):
    """
    Read (species, chain, is_aa) from a reference FASTA's name.

    The accepted form is ``[<prefix>_][aa_]<species>_<CHAIN>.fasta``: an optional
    source prefix (``imgt_``, ``airrc_``, ...) that is ignored, an optional
    ``aa_`` marking translated V, the species, and a known chain. Only the name is
    read, never the directory, so a file nested in a reference_base and a file in a
    flat folder are recognised the same way -- which is what lets both layouts
    build.

    A name that does not follow this form is not a failure: a manifest passed to
    discoverReference can say what such a file holds. See loadReferenceMap.

    Arguments:
      filename (str): a FASTA file's basename.

    Returns:
      tuple: (species, chain, is_aa), or None if the name does not match.
    """
    stem, extension = os.path.splitext(filename)
    if extension not in FASTA_SUFFIXES:
        return None

    tokens = stem.split('_')
    for index, token in enumerate(tokens):
        if token in SPECIES and index + 1 < len(tokens):
            chain = tokens[index + 1]
            if chain in KNOWN_CHAINS:
                return token, chain, 'aa' in tokens[:index]

    return None


def loadReferenceMap(path):
    """
    Read a manifest declaring what each FASTA in a reference folder holds.

    A reference someone assembled themselves rarely follows sourcerer's naming --
    an OGRDB set downloaded as ``IGH_VDJ_V.fasta`` says nothing about species or
    chain that a filename rule could read. The manifest says it outright rather
    than guessing, because a chain inferred wrongly does not fail: it builds a
    database with the alleles filed under the wrong locus.

    The format is tab- or whitespace-separated, ``#`` comments and blank lines
    ignored::

        #file                species  chain  [aa]
        IGH_VDJ_V.fasta      human    IGHV
        C57BL-6_IGH_V.fasta  mouse    IGHV
        my_translated.fasta  human    IGHV   aa

    A file is matched on its path relative to the reference folder, or failing
    that on its basename, so a manifest can stay short when names are unique and
    still disambiguate two files that share one.

    Arguments:
      path (Path): the manifest file.

    Returns:
      dict: file key to (species, chain, is_aa).

    Raises:
      SourcererError: if the manifest is not there, or on an unreadable line, an
        unknown species or an unknown chain, none of which are worth silently
        skipping in a file whose whole purpose is to be explicit.
    """
    path = Path(path)
    if not path.is_file():
        raise SourcererError('no such manifest: %s' % path)

    mapping = {}
    for number, raw in enumerate(path.read_text().splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith('#'):
            continue

        fields = line.split('\t') if '\t' in line else line.split()
        fields = [f.strip() for f in fields if f.strip()]
        if len(fields) < 3:
            raise SourcererError('%s line %d: expected at least file, species '
                                 'and chain, got %r' % (path, number, line))

        name, species, chain = fields[0], fields[1], fields[2]
        if species not in SPECIES:
            raise SourcererError('%s line %d: unknown species %r; known: %s'
                                 % (path, number, species, ', '.join(SPECIES)))
        if chain not in KNOWN_CHAINS:
            raise SourcererError('%s line %d: unknown chain %r; known: %s'
                                 % (path, number, chain,
                                    ', '.join(sorted(KNOWN_CHAINS))))

        is_aa = len(fields) > 3 and fields[3].lower() == 'aa'
        mapping[name] = (species, chain, is_aa)

    return mapping


def discoverReference(reference_dir, mapping=None):
    """
    Find every reference FASTA under a folder, by filename, in any layout.

    The folder is searched recursively, so a nested reference_base and a flat
    folder of FASTAs are both handled; classification is by name alone, unless a
    manifest declares it. The manifest wins over the filename rule, so it can
    correct a file whose name would otherwise parse to the wrong thing as well as
    one that does not parse at all.

    Arguments:
      reference_dir (Path): a reference_base tree or a flat folder of FASTAs.
      mapping (dict): a manifest from loadReferenceMap, or None.

    Returns:
      tuple: (files, unrecognized) where files is a list of
      (species, chain, is_aa, Path), and unrecognized is the list of FASTA paths
      that neither the manifest nor the naming rule could place.
    """
    reference_dir = Path(reference_dir)
    mapping = mapping or {}

    paths = sorted({path for suffix in FASTA_SUFFIXES
                    for path in reference_dir.rglob('*' + suffix)})

    files, unrecognized = [], []
    for path in paths:
        try:
            relative = str(path.relative_to(reference_dir))
        except ValueError:
            relative = path.name
        parsed = (mapping.get(relative) or mapping.get(path.name)
                  or parseReferenceName(path.name))
        if parsed is None:
            unrecognized.append(path)
        else:
            species, chain, is_aa = parsed
            files.append((species, chain, is_aa, path))

    return files, unrecognized


# ---------------------------------------------------------------------------
# reference_base provenance
# ---------------------------------------------------------------------------

def writeImgtMetadata(reference_dir, species, release, date, generated_by,
                      requested=None, exact=True):
    """
    Record which IMGT release each species in a reference_base was built from.

    Written inside reference_base as IMGT.yaml so the release travels with the
    FASTAs. The release is the one datum a download date does not carry and the
    one `download --from` needs to fetch the same version again.

    A download fetches one species, but a reference_base holds as many as were
    downloaded into it, so this merges rather than replaces: writing mouse into
    a folder that already has human leaves human alone. The release is recorded
    per species because the two can be downloaded weeks apart, from different
    IMGT builds, and a single release field would then misname one of them.

    When a release was reconstructed from the archive, what was asked for is
    recorded beside what was found. The archive does not hold every IMGT build,
    so a pinned re-download can land on a neighbouring release; recording only
    the release used would present that substitute as the one requested, and a
    later `--from` would pin the substitute in turn.

    Arguments:
      reference_dir (Path): the reference_base root.
      species (iterable): the species this download wrote, e.g. ('human',).
      release (str): the IMGT GENE-DB release tag used, e.g. '202631-7', or ''.
      date (str): the download date, YYYY-MM-DD.
      generated_by (str): the tool and version that wrote it.
      requested (str): the release asked for, when this was a pinned
        re-download, or None for a plain download of the current release.
      exact (bool): whether the archive held the requested release.

    Returns:
      Path: the file written.
    """
    path = Path(reference_dir) / IMGT_METADATA
    record = _loadMetadata(path) or {'source': 'https://www.imgt.org/genedb'}
    record['generated_by'] = generated_by
    known = dict(record.get('species') or {})

    for name in species:
        entry = {'release': release or None, 'date': date}
        if requested is not None:
            entry['requested'] = requested
            entry['exact'] = bool(exact)
        known[name] = entry
    record['species'] = {name: known[name] for name in sorted(known)}

    return _writeMetadata(path, record)


def writeAirrcMetadata(reference_dir, sets, date, generated_by):
    """
    Record the OGRDB sets, and their versions, a reference_base was built from.

    Written inside reference_base as AIRRC.yaml. Each set keeps its version and
    release date, and its DOI when it was resolved, so the exact sets can be
    fetched again with `download --from` and cited.

    A download fetches one species, but a reference_base holds as many as were
    downloaded into it, so this merges on (species, set): a species downloaded
    later joins the ones already recorded, and re-downloading a set replaces its
    entry rather than adding a second one.

    Arguments:
      reference_dir (Path): the reference_base root.
      sets (list): per-set dicts (species, locus, set, version, release_date,
        and optionally doi, zenodo_record_id, zenodo_url).
      date (str): the download date, YYYY-MM-DD.
      generated_by (str): the tool and version that wrote it.

    Returns:
      Path: the file written.
    """
    path = Path(reference_dir) / AIRRC_METADATA
    record = _loadMetadata(path) or {'source': 'https://ogrdb.airr-community.org'}
    record['generated_by'] = generated_by
    record['date'] = date

    known = {(item.get('species'), item.get('set')): item
             for item in record.get('sets') or []}
    known.update({(item.get('species'), item.get('set')): item for item in sets})
    record['sets'] = [known[key] for key in sorted(known, key=lambda k: (k[0] or '',
                                                                        k[1] or ''))]

    return _writeMetadata(path, record)


def _loadMetadata(path):
    """
    Read a provenance sidecar if it is there, so a write can merge into it.

    Arguments:
      path (Path): the sidecar.

    Returns:
      dict: the record, or None if there is nothing readable to merge with.
    """
    path = Path(path)
    if not path.is_file():
        return None
    try:
        with open(path) as handle:
            return yaml.safe_load(handle) or None
    except (OSError, yaml.YAMLError):
        log.warning('%s is unreadable; rewriting it from this download alone',
                    path.name)
        return None


def _writeMetadata(path, record):
    """Write a provenance record as YAML, keys in the order given."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w') as handle:
        yaml.safe_dump(record, handle, sort_keys=False, default_flow_style=False)

    return path


#: Provenance sidecar written into an igblast_base by `reference build`, so a
#: user-supplied reference records what was built and where it came from, the
#: same way a download records its release. Named apart from the source sidecars
#: because it describes the build, not an upstream release.
BUILD_METADATA = 'sourcerer_build.yaml'


def writeBuildMetadata(out_dir, report, source_dir, date, generated_by):
    """
    Record what `reference build` produced, and carry the source's provenance.

    A custom reference is built into an igblast_base with no upstream release to
    name, so this records the build itself: what was built, from where, what came
    up empty, which allele names had to be shortened, and which J alleles the
    mirrored auxiliary file does not cover. When the source reference already carries an IMGT.yaml or
    AIRRC.yaml, from an earlier download, those are copied in beside it so the
    version provenance travels with the built databases.

    Arguments:
      out_dir (Path): the igblast_base that was built.
      report (ReferenceReport): the build outcome.
      source_dir (Path): the reference folder the build read from.
      date (str): the build date, YYYY-MM-DD.
      generated_by (str): the tool and version that built it.

    Returns:
      list: the files written into out_dir.
    """
    out_dir, source_dir = Path(out_dir), Path(source_dir)
    record = {'generated_by': generated_by,
              'date': date,
              'built_from': _relativeTo(source_dir, out_dir),
              'databases': list(report.built),
              'skipped_empty': list(report.skipped_empty)}
    # Both of these change how the databases behave, so they travel with them
    # rather than living only in the log of the run that built them.
    if report.renamed:
        record['shortened_alleles'] = len(report.renamed)
    if report.aux_missing:
        record['aux_not_covered'] = {species: list(names) for species, names
                                     in sorted(report.aux_missing.items())}
    written = [_writeMetadata(out_dir / BUILD_METADATA, record)]

    for name in (IMGT_METADATA, AIRRC_METADATA):
        found = _findSidecar(source_dir, name)
        if found is not None:
            destination = out_dir / name
            shutil.copyfile(found, destination)
            written.append(destination)

    return written


def _relativeTo(path, root):
    """
    Express a path relative to the built directory, so a move does not break it.

    The reference a build reads is usually a sibling of the igblast_base it
    writes, so this walks up as well as down, unlike Provenance.relativize,
    which records paths inside one download root.

    Arguments:
      path (Path): the path to express.
      root (Path): the directory the record is written into.

    Returns:
      str: a relative path where one exists, else the absolute path.
    """
    try:
        return os.path.relpath(Path(path).resolve(), Path(root).resolve())
    except (ValueError, OSError):
        return str(path)


def _findSidecar(source_dir, name):
    """
    Find a provenance sidecar at the reference root or one level up.

    `reference build` is often pointed at a reference_base whose sidecar sits
    beside it, or at that reference_base itself, so both are checked.

    Arguments:
      source_dir (Path): the reference folder given to the build.
      name (str): the sidecar filename.

    Returns:
      Path: the sidecar, or None.
    """
    for candidate in (source_dir / name, source_dir.parent / name):
        if candidate.is_file():
            return candidate

    return None


def loadReferencePins(path):
    """
    Read the IMGT.yaml / AIRRC.yaml a `download --from` should re-download.

    Accepts a reference_base directory or a single metadata file. A directory is
    searched for both sidecars, so a blended reference pins both sources at once.

    Arguments:
      path (Path): a reference_base directory or an IMGT.yaml/AIRRC.yaml file.

    Returns:
      dict: {'imgt': record or None, 'airrc': record or None}.

    Raises:
      SourcererError: if nothing readable is found.
    """
    path = Path(path)
    pins = {'imgt': None, 'airrc': None}

    candidates = []
    if path.is_dir():
        candidates = [path / IMGT_METADATA, path / AIRRC_METADATA]
    else:
        candidates = [path]

    for candidate in candidates:
        if not candidate.exists():
            continue
        with open(candidate) as handle:
            record = yaml.safe_load(handle) or {}
        if 'sets' in record or candidate.name == AIRRC_METADATA:
            pins['airrc'] = record
        else:
            pins['imgt'] = record

    if pins['imgt'] is None and pins['airrc'] is None:
        raise SourcererError('no IMGT.yaml or AIRRC.yaml found at %s; point '
                             '--from at a reference_base or one of those files'
                             % path)

    return pins


def describeReference(path, mapping=None):
    """
    Render what a reference folder is and where it came from.

    The provenance sidecars are written so that a reference can be re-downloaded
    and cited without sourcerer's help, but reading three YAML files by hand to
    answer "what is this folder" is friction enough that it does not get done.
    This is that answer: the release and sets it was built from, what was built,
    and what the folder actually holds.

    Both a reference_base and a directory containing one are accepted, as are an
    igblast_base, which carries copies of the same sidecars.

    Arguments:
      path (Path): a reference folder, or a directory holding one.
      mapping (dict): a manifest from loadReferenceMap, or None.

    Returns:
      str: the report.

    Raises:
      SourcererError: if the folder does not exist.
    """
    path = Path(path)
    if not path.is_dir():
        raise SourcererError('no such reference folder: %s' % path)

    root = path / 'reference_base' if (path / 'reference_base').is_dir() else path
    lines = ['reference: %s' % root]

    found = False
    for name, title in ((IMGT_METADATA, 'IMGT'), (AIRRC_METADATA, 'AIRR-C'),
                        (BUILD_METADATA, 'build')):
        sidecar = _findSidecar(root, name)
        if sidecar is None:
            continue
        found = True
        with open(sidecar) as handle:
            record = yaml.safe_load(handle) or {}
        lines.append('')
        lines.append('%s (%s)' % (title, sidecar.name))
        # Rendered as the YAML it is, so what is shown matches the file on disk.
        lines.append(textwrap.indent(
            yaml.safe_dump(record, sort_keys=False).rstrip(), '  '))

    if not found:
        lines.append('')
        lines.append('no provenance sidecars; this folder was not written by a '
                     'sourcerer download or build')

    files, unrecognized = discoverReference(root, mapping=mapping)
    lines.append('')
    if files:
        counts = {}
        for species, chain, is_aa, _path in files:
            counts.setdefault(species, []).append(chain + (' (aa)' if is_aa else ''))
        lines.append('contents: %d FASTA file(s)' % len(files))
        for species in sorted(counts):
            lines.append('  %-8s %s' % (species, ', '.join(sorted(counts[species]))))
    else:
        lines.append('contents: no reference FASTAs found')
    for item in unrecognized:
        lines.append('  warning: %s is not in the reference naming format and no '
                     'manifest names it' % item.name)

    return '\n'.join(lines)


# ---------------------------------------------------------------------------
# comparing two references
# ---------------------------------------------------------------------------

@dataclass
class ChainDiff:
    """
    How one chain differs between two references.

    Arguments:
      species (str): the species.
      chain (str): the chain, e.g. 'IGHV'.
      is_aa (bool): whether this is the translated (amino acid) V.
      identical (int): alleles present in both with the same sequence.
      added (list): allele names only in the second reference.
      removed (list): allele names only in the first reference.
      changed (list): allele names in both whose sequence differs.
    """
    species: str
    chain: str
    is_aa: bool
    identical: int = 0
    added: list = field(default_factory=list)
    removed: list = field(default_factory=list)
    changed: list = field(default_factory=list)

    @property
    def same(self):
        """bool: True if nothing was added, removed or changed."""
        return not (self.added or self.removed or self.changed)


@dataclass
class ReferenceDiff:
    """
    An allele-level comparison of two reference folders.

    Arguments:
      chains (list): one ChainDiff per (species, chain) seen in either folder.
    """
    chains: list = field(default_factory=list)

    @property
    def same(self):
        """bool: True if every chain is identical between the two folders."""
        return all(chain.same for chain in self.chains)

    def summary(self):
        """
        Render the comparison, one line per chain, differences called out.

        Returns:
          str: the report.
        """
        lines = []
        for diff in self.chains:
            label = '%s %s%s' % (diff.species, diff.chain,
                                 ' (aa)' if diff.is_aa else '')
            if diff.same:
                lines.append('%-22s %d identical' % (label, diff.identical))
                continue
            parts = ['%d identical' % diff.identical]
            if diff.changed:
                parts.append('%d changed' % len(diff.changed))
            if diff.added:
                parts.append('%d added' % len(diff.added))
            if diff.removed:
                parts.append('%d removed' % len(diff.removed))
            lines.append('%-22s %s' % (label, ', '.join(parts)))
            for name in diff.changed:
                lines.append('    changed: %s' % name)
            for name in diff.added:
                lines.append('    added:   %s' % name)
            for name in diff.removed:
                lines.append('    removed: %s' % name)

        if not lines:
            return 'no reference FASTAs found to compare'
        verdict = 'identical' if self.same else 'differences found'
        return '\n'.join(lines + ['', verdict])


def _alleleMap(files):
    """
    Collapse discovered FASTAs into {(species, chain, is_aa): {name: sequence}}.

    Sequences are compared without gaps and in upper case, so the same allele
    filed gapped in one reference and ungapped in another is not a false change;
    only a real base difference is. The name is taken with alleleName so an IMGT
    pipe header and a bare OGRDB name for the same allele line up.

    Arguments:
      files (list): (species, chain, is_aa, Path) tuples from discoverReference.

    Returns:
      dict: (species, chain, is_aa) to {allele name: cleaned sequence}.
    """
    grouped = {}
    for species, chain, is_aa, path in files:
        target = grouped.setdefault((species, chain, is_aa), {})
        for header, sequence in parseFasta(path.read_text()):
            target[alleleName(header)] = sequence.replace('.', '').upper()

    return grouped


def diffReference(dir_a, dir_b, species=None, mapping=None):
    """
    Compare two reference folders allele by allele.

    Both folders are read by filename in any layout, so a nested reference_base
    and a flat folder compare the same way. This is what turns a re-download of
    the closest archived IMGT release into a checked result: run it against the
    reference you started with and it says whether the sets match or exactly what
    drifted.

    Arguments:
      dir_a (Path): the first reference folder (the baseline).
      dir_b (Path): the second reference folder.
      species (iterable): limit to these species, or None for every species seen.
      mapping (dict): a manifest from loadReferenceMap, applied to both folders,
        or None.

    Returns:
      ReferenceDiff: the per-chain comparison.
    """
    files_a, _ = discoverReference(dir_a, mapping=mapping)
    files_b, _ = discoverReference(dir_b, mapping=mapping)
    if species is not None:
        wanted = set(species)
        files_a = [f for f in files_a if f[0] in wanted]
        files_b = [f for f in files_b if f[0] in wanted]

    map_a = _alleleMap(files_a)
    map_b = _alleleMap(files_b)

    diff = ReferenceDiff()
    for key in sorted(set(map_a) | set(map_b)):
        species_, chain, is_aa = key
        left = map_a.get(key, {})
        right = map_b.get(key, {})
        chain_diff = ChainDiff(species=species_, chain=chain, is_aa=is_aa,
                               added=sorted(set(right) - set(left)),
                               removed=sorted(set(left) - set(right)))
        for name in sorted(set(left) & set(right)):
            if left[name] == right[name]:
                chain_diff.identical += 1
            else:
                chain_diff.changed.append(name)
        diff.chains.append(chain_diff)

    return diff


# ---------------------------------------------------------------------------
# igblast_base
# ---------------------------------------------------------------------------

def runMakeblastdb(fasta, out_base, dbtype):
    """
    Build one BLAST database from a cleaned FASTA.

    Arguments:
      fasta (Path): the input FASTA.
      out_base (Path): the database basename, without an extension.
      dbtype (str): 'nucl' or 'prot'.

    Raises:
      SourcererError: if makeblastdb is not on PATH or exits non-zero.
    """
    if shutil.which('makeblastdb') is None:
        raise SourcererError(
            'makeblastdb not found on PATH; install NCBI BLAST+ (for example '
            'conda install -c bioconda blast) or drop --igblast to write only '
            'the reference FASTAs')

    result = subprocess.run(
        ['makeblastdb', '-parse_seqids', '-dbtype', dbtype,
         '-in', str(fasta), '-out', str(out_base)],
        capture_output=True, text=True)
    if result.returncode != 0:
        raise SourcererError('makeblastdb failed for %s: %s'
                             % (Path(fasta).name, result.stderr.strip()))


def planReference(reference_dir, species=None, mapping=None):
    """
    Work out which IgBLAST databases a reference folder would produce.

    Files are grouped by name into the canonical (species, locus, segment)
    databases airrflow expects, whatever layout they came in and whatever prefix
    wrote them, so a nested reference_base and a flat folder both plan. No
    makeblastdb is run, so this is also the format check: the returned plan says
    what would build, what is empty, what was not recognised, and where duplicate
    names were dropped.

    Arguments:
      reference_dir (Path): a reference_base tree or a flat folder of FASTAs.
      species (iterable): limit to these species, or None for every species found.
      mapping (dict): a manifest from loadReferenceMap, for files whose names do
        not say what they hold, or None.

    Returns:
      ReferencePlan: the databases that would build and the diagnostics.
    """
    files, unrecognized = discoverReference(reference_dir, mapping=mapping)
    found = sorted({item[0] for item in files})
    wanted = list(species) if species else found

    contents = {path: parseFasta(path.read_text())
                for _sp, _chain, _aa, path in files}
    empty_files = [path for path, records in contents.items() if not records]

    def collect(sp, chains, is_aa):
        records = []
        for f_sp, f_chain, f_aa, f_path in files:
            if f_sp == sp and f_aa == is_aa and f_chain in chains:
                records.extend(contents[f_path])
        return records

    plan = ReferencePlan(found_species=found, unrecognized=unrecognized,
                         empty_files=empty_files)
    for sp in wanted:
        for locus in LOCI:
            for segment in SEGMENTS:
                _addToPlan(plan, '%s_%s_%s' % (sp, locus, segment), 'nucl',
                           collect(sp, LOCUS_CHAINS[(locus, segment)], False))
            # Amino acid V is a protein database, built only when the folder
            # actually carries translated V (OGRDB, for one, does not), so an
            # absent one is not reported as an empty gap.
            _addToPlan(plan, 'aa_%s_%s_v' % (sp, locus), 'prot',
                       collect(sp, LOCUS_CHAINS[(locus, 'v')], True),
                       keep_empty=False)

    return plan


def _addToPlan(plan, basename, dbtype, records, keep_empty=True):
    """
    Clean one canonical database's records and record it on the plan.

    Arguments:
      plan (ReferencePlan): the plan to add to.
      basename (str): the canonical database basename.
      dbtype (str): 'nucl' or 'prot'.
      records (list): the raw (header, sequence) tuples gathered for it.
      keep_empty (bool): whether an empty database is worth reporting as a gap.
    """
    if not records:
        if keep_empty:
            plan.empty.append(basename)
        return

    cleaned = cleanForBlast(records)
    dropped = len(records) - len(cleaned)
    if dropped:
        plan.duplicates[basename] = dropped

    cleaned, renamed = shortenForBlast(cleaned)
    plan.renamed.update(renamed)
    plan.databases.append((basename, dbtype, cleaned))


def buildIgblastBase(reference_dir, out_dir, client, species=None, mapping=None):
    """
    Build the IgBLAST database tree airrflow expects from a reference folder.

    A thin wrapper over planReference and buildFromPlan, kept so the download
    path and the standalone `reference build` command share one code path.

    Arguments:
      reference_dir (Path): a reference_base tree or a flat folder of FASTAs.
      out_dir (Path): the igblast_base to write.
      client (HttpClient): used to mirror the NCBI support trees.
      species (iterable): limit to these species, or None for every species found.
      mapping (dict): a manifest from loadReferenceMap, or None.

    Returns:
      ReferenceReport: what was built and what was skipped.

    Raises:
      SourcererError: if makeblastdb is unavailable.
    """
    plan = planReference(reference_dir, species=species, mapping=mapping)
    if plan.unrecognized:
        log.warning('%d file(s) skipped: names not in the reference format',
                    len(plan.unrecognized))

    return buildFromPlan(plan, out_dir, client)


def writeRenameLog(out_dir, renamed):
    """
    Record every allele name a build had to shorten.

    A shortened name is what IgBLAST reports and what lands in a v_call, so the
    mapping back to the original has to be written down or the results cannot be
    traced to the reference they came from.

    Arguments:
      out_dir (Path): the igblast_base being built.
      renamed (dict): original name to shortened name.

    Returns:
      Path: the file written.
    """
    path = Path(out_dir) / RENAME_LOG
    with open(path, 'w') as handle:
        handle.write('#original\tshortened\n')
        for original in sorted(renamed):
            handle.write('%s\t%s\n' % (original, renamed[original]))

    return path


def checkAuxCoverage(out_dir, plan):
    """
    Find J alleles the mirrored auxiliary file does not name.

    IgBLAST looks a J germline up in the auxiliary file **by name**, unlike the
    V delineation file, which it transfers by alignment. A J allele with no row
    there gets no CDR3 and no productivity call, silently -- igblastn reports
    ``V-J frame: N/A`` and no CDR3 sub-region rather than failing.

    The auxiliary file mirrored here is NCBI's, so it covers the names IMGT
    publishes. A reference naming its J alleles otherwise -- OGRDB's mouse sets
    call them ``IGKJ0-4JXG*00`` -- is not covered, and nothing else would say so.

    The remedy is to build an auxiliary file from the reference itself and pass
    it to igblastn with ``-auxiliary_data``. Sourcerer does not build one: that
    belongs with the pipeline that runs IgBLAST. Reporting the gap is what it can
    usefully do, so the names needing rows are known rather than guessed at.

    Arguments:
      out_dir (Path): the igblast_base being built.
      plan (ReferencePlan): the databases that were built.

    Returns:
      dict: species to the sorted J allele names with no auxiliary row.
    """
    missing = {}
    for basename, _dbtype, records in plan.databases:
        species, *_rest, segment = basename.split('_')
        if segment != 'j' or species == 'aa':
            continue

        path = Path(out_dir) / 'optional_file' / ('%s_gl.aux' % species)
        if not path.is_file():
            continue
        named = {line.split()[0] for line in path.read_text().splitlines()
                 if line.strip() and not line.startswith('#')}

        absent = {name for name, _sequence in records if name not in named}
        if absent:
            missing.setdefault(species, set()).update(absent)

    return {species: sorted(names) for species, names in missing.items()}


def buildFromPlan(plan, out_dir, client):
    """
    Write the databases a plan describes, then mirror the NCBI support trees.

    Arguments:
      plan (ReferencePlan): the databases to build, from planReference.
      out_dir (Path): the igblast_base to write; fasta/ and database/ are created
        inside it, alongside the mirrored internal_data/ and optional_file/.
      client (HttpClient): used to mirror the NCBI support trees.

    Returns:
      ReferenceReport: what was built and what was skipped.

    Raises:
      SourcererError: if makeblastdb is unavailable.
    """
    out_dir = Path(out_dir)
    fasta_out = out_dir / 'fasta'
    db_out = out_dir / 'database'
    fasta_out.mkdir(parents=True, exist_ok=True)
    db_out.mkdir(parents=True, exist_ok=True)

    report = ReferenceReport(skipped_empty=list(plan.empty))
    for basename, dbtype, records in plan.databases:
        fasta = fasta_out / ('%s.fasta' % basename)
        with open(fasta, 'w') as handle:
            for name, sequence in records:
                handle.write('>%s\n%s\n' % (name, sequence))
        runMakeblastdb(fasta, db_out / basename, dbtype)
        report.built.append(basename)

    mirrorSupport(out_dir, client)

    report.renamed = dict(plan.renamed)
    if plan.renamed:
        log.warning('%d allele name(s) were longer than the %d-character BLAST '
                    'limit and were shortened; the mapping is in %s',
                    len(plan.renamed), BLAST_NAME_LIMIT,
                    writeRenameLog(out_dir, plan.renamed))

    report.aux_missing = checkAuxCoverage(out_dir, plan)
    for species, names in sorted(report.aux_missing.items()):
        log.warning('%s: %d J allele(s) are not named in the mirrored NCBI '
                    'auxiliary file (%s_gl.aux), so IgBLAST will call no CDR3 '
                    'and no productivity for reads assigned to them. Build an '
                    'auxiliary file from this reference and pass it to igblastn '
                    'with -auxiliary_data. Uncovered: %s',
                    species, len(names), species, ', '.join(names[:5])
                    + (', ...' if len(names) > 5 else ''))

    report.logSummary()

    return report


def mirrorSupport(out_dir, client):
    """
    Mirror the NCBI IgBLAST support trees into an igblast_base.

    database/, internal_data/ and optional_file/ are copied from NCBI's release
    directory, and the tar archives NCBI ships inside database/ are unpacked in
    place, so the result matches what airrflow's fetch_igblastdb.sh produces.

    Arguments:
      out_dir (Path): the igblast_base root.
      client (HttpClient): the shared HTTP client.
    """
    out_dir = Path(out_dir)
    database_dir = out_dir / 'database'
    mirrorTree(NCBI_DATABASE_URL, database_dir, client)
    for name in NCBI_TAR_ARCHIVES:
        archive = database_dir / name
        if archive.exists():
            extractTar(archive, database_dir)

    mirrorTree(NCBI_INTERNAL_URL, out_dir / 'internal_data', client)
    mirrorTree(NCBI_OPTIONAL_URL, out_dir / 'optional_file', client)


def mirrorTree(url, dest_dir, client, seen=None):
    """
    Recursively mirror an Apache/NCBI directory index into a local tree.

    Only links that stay under the starting URL are followed, so a parent link
    or an absolute link elsewhere on the host cannot walk the mirror out of the
    subtree it was pointed at.

    Arguments:
      url (str): the directory index URL, ending in '/'.
      dest_dir (Path): where to mirror it.
      client (HttpClient): the shared HTTP client.
      seen (set): URLs already visited, to guard against a self-referential index.
    """
    seen = seen if seen is not None else set()
    if url in seen:
        return
    seen.add(url)

    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)

    soup = BeautifulSoup(client.get(url).text, 'html.parser')
    for anchor in soup.find_all('a'):
        href = anchor.get('href')
        if not href or href in ('../', './', '/'):
            continue

        child = urljoin(url, href)
        if urlparse(child).scheme not in ('http', 'https'):
            continue
        if not child.startswith(url):
            continue

        name = Path(urlparse(child).path).name
        if not name:
            continue

        if href.endswith('/'):
            mirrorTree(child, dest_dir / name, client, seen)
        else:
            client.fetch(child, dest_dir / name, progress=False)


def extractTar(archive, dest_dir):
    """
    Unpack a tar archive, refusing any member that would escape the destination.

    Arguments:
      archive (Path): the tar file.
      dest_dir (Path): where to extract it.

    Raises:
      SourcererError: if a member path points outside dest_dir.
    """
    dest_dir = Path(dest_dir).resolve()
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            target = (dest_dir / member.name).resolve()
            if target != dest_dir and dest_dir not in target.parents:
                raise SourcererError('unsafe path %s in %s'
                                     % (member.name, Path(archive).name))
        tar.extractall(dest_dir)
