"""
trimAl trimming wrapper.

Checks for trimAl availability at import time and provides a single
trim() function for use by the pipeline.
"""
import os
import shutil
import subprocess


# ── Exceptions ────────────────────────────────────────────────────────────────

class TooFewSequencesError(Exception):
    """Raised when a FASTA file contains too few sequences to trim."""
    pass


class ToolNotFoundError(Exception):
    """Raised when a required external tool is not found on PATH."""
    pass


# ── Tool availability ─────────────────────────────────────────────────────────

def check_trimal():
    """
    Check that trimAl is available and return its version string.

    Raises:
        ToolNotFoundError if trimal is not found on PATH
    """
    try:
        result = subprocess.run(
            ['trimal', '--version'],
            capture_output=True, text=True
        )
        version = result.stdout.strip() or result.stderr.strip()
        return version
    except FileNotFoundError:
        raise ToolNotFoundError(
            "trimAl not found on PATH. "
            "Install with: conda install -c bioconda trimal"
        )


# ── Sequence counting ─────────────────────────────────────────────────────────

def count_sequences(fasta_path):
    """Count sequences in a FASTA file by counting '>' lines."""
    with open(fasta_path) as f:
        return sum(1 for line in f if line.startswith('>'))


# ── Trimming ──────────────────────────────────────────────────────────────────

MIN_SEQUENCES = 3


def trim(alignment_path, output_path, log_path=None):
    """
    Trim an alignment using trimAl -automated1.

    If trimAl produces an empty output file (can happen with very
    sparse alignments), falls back to copying the untrimmed alignment
    and logs a warning.

    Args:
        alignment_path: path to input aligned FASTA
        output_path:    path to write trimmed FASTA
        log_path:       path to append trimAl log output (optional)

    Returns:
        output_path on success

    Raises:
        TooFewSequencesError if alignment_path has fewer than MIN_SEQUENCES
        ToolNotFoundError if trimAl is not on PATH
        subprocess.CalledProcessError if trimAl exits with non-zero status
    """
    n = count_sequences(alignment_path)
    if n < MIN_SEQUENCES:
        raise TooFewSequencesError(
            f"{os.path.basename(alignment_path)}: {n} sequence(s) "
            f"(minimum {MIN_SEQUENCES} required for trimming)"
        )

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    cmd = ['trimal', '-in', alignment_path, '-out', output_path, '-automated1']

    if log_path:
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        with open(log_path, 'a') as log:
            subprocess.run(cmd, stderr=log, check=True)
    else:
        subprocess.run(cmd, stderr=subprocess.DEVNULL, check=True)

    # Fall back to untrimmed if output is empty
    if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
        print(f"  WARNING: trimAl produced empty output for "
              f"{os.path.basename(alignment_path)}, using untrimmed")
        shutil.copy(alignment_path, output_path)

    return output_path
