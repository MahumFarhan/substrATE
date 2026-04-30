"""
dbCAN annotation wrapper.

Handles two input cases:
  1. Protein FASTA (.faa) — runs dbCAN directly in protein mode
  2. Nucleotide FASTA (.fna/.fasta/.fa) — runs Prodigal first to
     predict proteins, then runs dbCAN on the predicted .faa

Uses run_dbcan easy_substrate which performs CAZyme annotation,
GFF processing, CGC identification, and substrate prediction in
one step.

Prodigal output is written to <output_dir>/prodigal/
dbCAN output is written to <output_dir>/cgc_output/output_<sample>/
"""
import os
import subprocess


# ── Constants ─────────────────────────────────────────────────────────────────

PROTEIN_EXTENSIONS    = {'.faa'}
NUCLEOTIDE_EXTENSIONS = {'.fna', '.fasta', '.fa'}


# ── Exceptions ────────────────────────────────────────────────────────────────

class ToolNotFoundError(Exception):
    """Raised when a required external tool is not found on PATH."""
    pass


# ── Tool availability ─────────────────────────────────────────────────────────

def check_dbcan():
    """
    Check that run_dbcan is available and return its version string.

    Raises:
        ToolNotFoundError if run_dbcan is not found on PATH
    """
    try:
        result = subprocess.run(
            ['run_dbcan', 'version'],
            capture_output=True, text=True
        )
        version = result.stdout.strip() or result.stderr.strip()
        return version
    except FileNotFoundError:
        raise ToolNotFoundError(
            "run_dbcan not found on PATH. "
            "Install with: conda install -c bioconda dbcan"
        )


def check_prodigal():
    """
    Check that Prodigal is available and return its version string.

    Raises:
        ToolNotFoundError if prodigal is not found on PATH
    """
    try:
        result = subprocess.run(
            ['prodigal', '-v'],
            capture_output=True, text=True
        )
        version = result.stderr.strip() or result.stdout.strip()
        return version
    except FileNotFoundError:
        raise ToolNotFoundError(
            "Prodigal not found on PATH. "
            "Install with: conda install -c bioconda prodigal"
        )


# ── Input file detection ──────────────────────────────────────────────────────

def detect_input_type(filepath):
    """
    Detect whether a file is a protein or nucleotide FASTA.

    Args:
        filepath: path to FASTA file

    Returns:
        'protein' or 'nucleotide'

    Raises:
        ValueError if extension is not recognised
    """
    for ext in PROTEIN_EXTENSIONS:
        if filepath.endswith(ext):
            return 'protein'
    for ext in NUCLEOTIDE_EXTENSIONS:
        if filepath.endswith(ext):
            return 'nucleotide'
    raise ValueError(
        f"Unrecognised file extension for "
        f"{os.path.basename(filepath)}. "
        f"Protein inputs: {sorted(PROTEIN_EXTENSIONS)}. "
        f"Nucleotide inputs: {sorted(NUCLEOTIDE_EXTENSIONS)}."
    )


def find_input_fastas(input_dir):
    """
    Scan input_dir for protein and nucleotide FASTA files.

    Args:
        input_dir: path to directory containing input FASTA files

    Returns:
        list of (sample_name, filepath, input_type) tuples,
        sorted by sample name
    """
    all_extensions = PROTEIN_EXTENSIONS | NUCLEOTIDE_EXTENSIONS
    samples = []

    for fname in sorted(os.listdir(input_dir)):
        filepath = os.path.join(input_dir, fname)
        if not os.path.isfile(filepath):
            continue
        for ext in all_extensions:
            if fname.endswith(ext):
                sample_name = fname[:-len(ext)]
                input_type  = detect_input_type(filepath)
                samples.append((sample_name, filepath, input_type))
                break

    return samples


# ── Prodigal gene prediction ──────────────────────────────────────────────────

def run_prodigal(sample, nucleotide_fasta, prodigal_dir,
                 log_path=None):
    """
    Run Prodigal on a nucleotide FASTA to predict protein sequences.

    Args:
        sample:           sample name string
        nucleotide_fasta: path to input nucleotide FASTA
        prodigal_dir:     directory to write Prodigal output
        log_path:         path to append Prodigal log output (optional)

    Returns:
        path to predicted protein FASTA (.faa)
    """
    check_prodigal()
    os.makedirs(prodigal_dir, exist_ok=True)

    faa_path = os.path.join(prodigal_dir, f"{sample}.faa")
    gbk_path = os.path.join(prodigal_dir, f"{sample}.gbk")

    cmd = [
        'prodigal',
        '-i', nucleotide_fasta,
        '-a', faa_path,
        '-o', gbk_path,
        '-p', 'single',
        '-q',
    ]

    print(f"  Running Prodigal on {sample}...")

    if log_path:
        with open(log_path, 'a') as log:
            subprocess.run(cmd, stderr=log, check=True)
    else:
        subprocess.run(cmd, stderr=subprocess.DEVNULL, check=True)

    print(f"  Prodigal done -> {faa_path}")
    return faa_path


# ── dbCAN annotation ──────────────────────────────────────────────────────────

def run_dbcan_sample(sample, faa_path, cgc_output_dir, db_dir,
                     threads=8, log_path=None):
    """
    Run dbCAN easy_substrate on a single protein FASTA.

    Uses easy_substrate which performs CAZyme annotation, GFF
    processing, CGC identification, and substrate prediction in
    one step. This matches the output format of dbCAN 5.x.

    Output is written to <cgc_output_dir>/output_<sample>/

    Args:
        sample:         sample name string
        faa_path:       path to protein FASTA (.faa)
        cgc_output_dir: directory to write dbCAN output
        db_dir:         path to dbCAN database directory
        threads:        number of threads (default: 8)
        log_path:       path to append dbCAN log output (optional)

    Returns:
        path to dbCAN output directory for this sample
    """
    check_dbcan()
    os.makedirs(cgc_output_dir, exist_ok=True)

    sample_output_dir = os.path.join(
        cgc_output_dir, f"output_{sample}")

    cmd = [
        'run_dbcan', 'easy_substrate',
        '--mode',           'protein',
        '--input_raw_data', faa_path,
        '--output_dir',     sample_output_dir,
        '--db_dir',         db_dir,
        '--threads',        str(threads),
    ]

    print(f"  Running dbCAN on {sample}...")

    if log_path:
        with open(log_path, 'a') as log:
            subprocess.run(cmd, stderr=log, check=True)
    else:
        subprocess.run(cmd, stderr=subprocess.DEVNULL, check=True)

    print(f"  dbCAN done -> {sample_output_dir}")
    return sample_output_dir


# ── Main entry point ──────────────────────────────────────────────────────────

def annotate_genomes(input_dir, output_dir, db_dir, threads=8,
                     force=False):
    """
    Run gene prediction (if needed) and dbCAN annotation on all
    input FASTA files in input_dir.

    Args:
        input_dir:  path to directory containing input FASTA files
        output_dir: base output directory
        db_dir:     path to dbCAN database directory
        threads:    number of threads for dbCAN and Prodigal
        force:      if True, rerun even if output already exists

    Returns:
        path to cgc_output directory containing all sample outputs
    """
    check_dbcan()

    prodigal_dir   = os.path.join(output_dir, 'prodigal')
    cgc_output_dir = os.path.join(output_dir, 'cgc_output')
    log_dir        = os.path.join(output_dir, 'logs')
    os.makedirs(log_dir, exist_ok=True)

    prodigal_log = os.path.join(log_dir, 'prodigal.log')
    dbcan_log    = os.path.join(log_dir, 'dbcan.log')

    samples = find_input_fastas(input_dir)

    if not samples:
        raise FileNotFoundError(
            f"No input FASTA files found in {input_dir}. "
            f"Expected files with extensions: "
            f"{sorted(PROTEIN_EXTENSIONS | NUCLEOTIDE_EXTENSIONS)}"
        )

    print(f"Found {len(samples)} input files in {input_dir}")

    n_protein    = sum(1 for _, _, t in samples if t == 'protein')
    n_nucleotide = sum(1 for _, _, t in samples if t == 'nucleotide')
    if n_protein > 0:
        print(f"  {n_protein} protein FASTA(s) — dbCAN direct")
    if n_nucleotide > 0:
        print(f"  {n_nucleotide} nucleotide FASTA(s) — "
              f"Prodigal + dbCAN")
        check_prodigal()

    for sample, filepath, input_type in samples:
        sample_output = os.path.join(
            cgc_output_dir, f"output_{sample}")

        if os.path.exists(sample_output) and not force:
            print(f"  Skipping {sample} — output already exists "
                  f"(use --force to rerun)")
            continue

        print(f"\nProcessing {sample} ({input_type})...")

        if input_type == 'nucleotide':
            faa_path = run_prodigal(
                sample, filepath, prodigal_dir,
                log_path=prodigal_log)
        else:
            faa_path = filepath

        run_dbcan_sample(
            sample, faa_path, cgc_output_dir, db_dir,
            threads=threads, log_path=dbcan_log
        )

    print(f"\nAnnotation complete. dbCAN output: {cgc_output_dir}")
    return cgc_output_dir
