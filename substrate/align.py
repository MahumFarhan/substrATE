"""
MAFFT alignment wrapper.

Checks for MAFFT availability at import time and provides a single
align() function for use by the pipeline.
"""
import os
import subprocess


# ── Exceptions ────────────────────────────────────────────────────────────────

class TooFewSequencesError(Exception):
    """Raised when a FASTA file contains too few sequences to align."""
    pass


class ToolNotFoundError(Exception):
    """Raised when a required external tool is not found on PATH."""
    pass


# ── Tool availability ─────────────────────────────────────────────────────────

def check_mafft():
    """
    Check that MAFFT is available and return its version string.

    Raises:
        ToolNotFoundError if mafft is not found on PATH
    """
    try:
        result = subprocess.run(
            ['mafft', '--version'],
            capture_output=True, text=True
        )
        # mafft --version writes to stderr
        version = result.stderr.strip() or result.stdout.strip()
        return version
    except FileNotFoundError:
        raise ToolNotFoundError(
            "MAFFT not found on PATH. "
            "Install with: conda install -c bioconda mafft"
        )


# ── Sequence counting ─────────────────────────────────────────────────────────

def count_sequences(fasta_path):
    """Count sequences in a FASTA file by counting '>' lines."""
    with open(fasta_path) as f:
        return sum(1 for line in f if line.startswith('>'))


# ── Alignment ─────────────────────────────────────────────────────────────────

MIN_SEQUENCES = 3


def add_fragments(query_path, reference_aln_path, output_path,
                  threads=8, log_path=None):
    """
    Add query sequences into an existing reference alignment using
    MAFFT --addfragments.

    Query sequences are inserted into the reference alignment without
    changing the alignment columns of the reference sequences. This is
    used for the place tree mode, where genomic sequences are added
    to the pre-built reference alignment before IQ-TREE2 --tree-fix.

    Args:
        query_path:         path to query sequences FASTA (genomic seqs)
        reference_aln_path: path to existing reference alignment FASTA
                            (the .ref.trim file from build-reference-trees)
        output_path:        path to write combined aligned FASTA
        threads:            number of threads (default: 8)
        log_path:           path to append MAFFT log output (optional)

    Returns:
        output_path on success

    Raises:
        TooFewSequencesError if query_path has fewer than MIN_SEQUENCES
        ToolNotFoundError if MAFFT is not on PATH
        subprocess.CalledProcessError if MAFFT exits with non-zero status
    """
    n = count_sequences(query_path)
    if n < MIN_SEQUENCES:
        raise TooFewSequencesError(
            f"{os.path.basename(query_path)}: {n} sequence(s) "
            f"(minimum {MIN_SEQUENCES} required for alignment)"
        )

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    cmd = [
        'mafft',
        '--add', query_path,
        '--thread', str(threads),
        '--quiet',
        reference_aln_path,
    ]

    with open(output_path, 'w') as out:
        if log_path:
            os.makedirs(os.path.dirname(log_path), exist_ok=True)
            with open(log_path, 'a') as log:
                subprocess.run(cmd, stdout=out, stderr=log, check=True)
        else:
            subprocess.run(cmd, stdout=out,
                           stderr=subprocess.DEVNULL, check=True)

    return output_path


def align(fasta_path, output_path, threads=8, log_path=None):
    """
    Align sequences in fasta_path using MAFFT --auto.

    Args:
        fasta_path:  path to input FASTA file
        output_path: path to write aligned FASTA
        threads:     number of threads (default: 8)
        log_path:    path to append MAFFT log output (optional)

    Returns:
        output_path on success

    Raises:
        TooFewSequencesError if fasta_path has fewer than MIN_SEQUENCES
        ToolNotFoundError if MAFFT is not on PATH
        subprocess.CalledProcessError if MAFFT exits with non-zero status
    """
    n = count_sequences(fasta_path)
    if n < MIN_SEQUENCES:
        raise TooFewSequencesError(
            f"{os.path.basename(fasta_path)}: {n} sequence(s) "
            f"(minimum {MIN_SEQUENCES} required for alignment)"
        )

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    cmd = ['mafft', '--auto', '--thread', str(threads), fasta_path]

    with open(output_path, 'w') as out:
        if log_path:
            os.makedirs(os.path.dirname(log_path), exist_ok=True)
            with open(log_path, 'a') as log:
                subprocess.run(cmd, stdout=out, stderr=log, check=True)
        else:
            subprocess.run(cmd, stdout=out,
                           stderr=subprocess.DEVNULL, check=True)

    return output_path
