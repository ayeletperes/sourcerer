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
import logging
import shutil
import subprocess
import tarfile
from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import urljoin, urlparse

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
      skipped_empty (list): canonical databases skipped because no chain in them
        had any sequence. This is the normal outcome for what a source does not
        cover, such as TR from OGRDB.
    """
    written: list = field(default_factory=list)
    built: list = field(default_factory=list)
    skipped_empty: list = field(default_factory=list)

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
    """
    found_species: list = field(default_factory=list)
    databases: list = field(default_factory=list)
    empty: list = field(default_factory=list)
    unrecognized: list = field(default_factory=list)
    empty_files: list = field(default_factory=list)
    duplicates: dict = field(default_factory=dict)

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


def parseReferenceName(filename):
    """
    Read (species, chain, is_aa) from a reference FASTA's name.

    The accepted form is ``[<prefix>_][aa_]<species>_<CHAIN>.fasta``: an optional
    source prefix (``imgt_``, ``airrc_``, ...) that is ignored, an optional
    ``aa_`` marking translated V, the species, and a known chain. Only the name is
    read, never the directory, so a file nested in a reference_base and a file in a
    flat folder are recognised the same way -- which is what lets both layouts
    build.

    Arguments:
      filename (str): a FASTA file's basename.

    Returns:
      tuple: (species, chain, is_aa), or None if the name does not match.
    """
    if not filename.endswith('.fasta'):
        return None

    tokens = filename[:-len('.fasta')].split('_')
    for index, token in enumerate(tokens):
        if token in SPECIES and index + 1 < len(tokens):
            chain = tokens[index + 1]
            if chain in KNOWN_CHAINS:
                return token, chain, 'aa' in tokens[:index]

    return None


def discoverReference(reference_dir):
    """
    Find every reference FASTA under a folder, by filename, in any layout.

    The folder is searched recursively, so a nested reference_base and a flat
    folder of FASTAs are both handled; classification is by name alone.

    Arguments:
      reference_dir (Path): a reference_base tree or a flat folder of FASTAs.

    Returns:
      tuple: (files, unrecognized) where files is a list of
      (species, chain, is_aa, Path), and unrecognized is the list of .fasta paths
      whose names are not in the reference format.
    """
    files, unrecognized = [], []
    for path in sorted(Path(reference_dir).rglob('*.fasta')):
        parsed = parseReferenceName(path.name)
        if parsed is None:
            unrecognized.append(path)
        else:
            species, chain, is_aa = parsed
            files.append((species, chain, is_aa, path))

    return files, unrecognized


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


def planReference(reference_dir, species=None):
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

    Returns:
      ReferencePlan: the databases that would build and the diagnostics.
    """
    files, unrecognized = discoverReference(reference_dir)
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
    plan.databases.append((basename, dbtype, cleaned))


def buildIgblastBase(reference_dir, out_dir, client, species=None):
    """
    Build the IgBLAST database tree airrflow expects from a reference folder.

    A thin wrapper over planReference and buildFromPlan, kept so the download
    path and the standalone `reference build` command share one code path.

    Arguments:
      reference_dir (Path): a reference_base tree or a flat folder of FASTAs.
      out_dir (Path): the igblast_base to write.
      client (HttpClient): used to mirror the NCBI support trees.
      species (iterable): limit to these species, or None for every species found.

    Returns:
      ReferenceReport: what was built and what was skipped.

    Raises:
      SourcererError: if makeblastdb is unavailable.
    """
    plan = planReference(reference_dir, species=species)
    if plan.unrecognized:
        log.warning('%d file(s) skipped: names not in the reference format',
                    len(plan.unrecognized))

    return buildFromPlan(plan, out_dir, client)


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
