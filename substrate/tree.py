"""
IQ-TREE2 phylogenetic tree wrapper.

Checks for IQ-TREE2 availability at import time and provides a single
build_tree() function for use by the pipeline.

Note: IQ-TREE2 may be installed as either 'iqtree' or 'iqtree2'
depending on the conda package version. Both are checked at startup.
"""
import os
import subprocess


# ── Exceptions ────────────────────────────────────────────────────────────────

class TooFewSequencesError(Exception):
    """Raised when a FASTA file contains too few sequences to build a tree."""
    pass


class ToolNotFoundError(Exception):
    """Raised when a required external tool is not found on PATH."""
    pass


# ── Tool availability ─────────────────────────────────────────────────────────

def _find_iqtree():
    """
    Find the IQ-TREE2 binary name on PATH.

    Tries 'iqtree2' first, then 'iqtree', since the binary name varies
    between conda package versions.

    Returns:
        Binary name string ('iqtree2' or 'iqtree')

    Raises:
        ToolNotFoundError if neither binary is found
    """
    for binary in ['iqtree2', 'iqtree']:
        try:
            subprocess.run(
                [binary, '--version'],
                capture_output=True, check=True
            )
            return binary
        except FileNotFoundError:
            continue
        except subprocess.CalledProcessError:
            # Binary exists but returned non-zero — still found
            return binary
    raise ToolNotFoundError(
        "IQ-TREE2 not found on PATH (tried 'iqtree2' and 'iqtree'). "
        "Install with: conda install -c bioconda iqtree"
    )


def check_iqtree():
    """
    Check that IQ-TREE2 is available and return its version string.

    Raises:
        ToolNotFoundError if neither iqtree2 nor iqtree is found on PATH
    """
    binary = _find_iqtree()
    result = subprocess.run(
        [binary, '--version'],
        capture_output=True, text=True
    )
    version = result.stdout.strip().split('\n')[0]
    return binary, version


# ── Sequence counting ─────────────────────────────────────────────────────────

def count_sequences(fasta_path):
    """Count sequences in a FASTA file by counting '>' lines."""
    with open(fasta_path) as f:
        return sum(1 for line in f if line.startswith('>'))


# ── Tree building ─────────────────────────────────────────────────────────────

MIN_SEQUENCES = 4


def place_sequences(combined_aln_path, ref_treefile, output_prefix,
                    threads=8, log_path=None):
    """
    Build a tree using a reference tree as a starting topology.

    Uses the reference tree (-t) as a starting point for a fast
    IQ-TREE2 search (LG+G4 + --fast), which is much faster than a
    full de novo build while producing a tree rooted in the reference
    phylogeny. New genomic sequences are added alongside reference
    sequences in the combined alignment.

    Args:
        combined_aln_path: path to combined alignment FASTA containing
                           both reference and query sequences (output of
                           align.add_fragments())
        ref_treefile:      path to reference tree (.ref.treefile from
                           build-reference-trees)
        output_prefix:     prefix for IQ-TREE2 output files
        threads:           number of threads (default: 8)
        log_path:          path to append IQ-TREE2 log output (optional)

    Returns:
        path to the output treefile on success

    Raises:
        TooFewSequencesError if combined_aln_path has fewer than
        MIN_SEQUENCES
        ToolNotFoundError if IQ-TREE2 is not on PATH
        subprocess.CalledProcessError if IQ-TREE2 exits non-zero
    """
    n = count_sequences(combined_aln_path)
    if n < MIN_SEQUENCES:
        raise TooFewSequencesError(
            f"{os.path.basename(combined_aln_path)}: {n} sequence(s) "
            f"(minimum {MIN_SEQUENCES} required)"
        )

    os.makedirs(os.path.dirname(output_prefix), exist_ok=True)

    binary, _ = check_iqtree()

    cmd = [
        binary,
        '-s',        combined_aln_path,
        '-t',        ref_treefile,
        '--prefix',  output_prefix,
        '-m',        'LG+G4',
        '--fast',
        '-T',        str(threads),
        '--quiet',
        '--redo',
    ]

    if log_path:
        with open(log_path, 'a') as log:
            subprocess.run(cmd, stderr=log, check=True)
    else:
        subprocess.run(cmd, stderr=subprocess.DEVNULL, check=True)

    treefile = f"{output_prefix}.treefile"
    if not os.path.exists(treefile):
        raise FileNotFoundError(
            f"IQ-TREE2 completed but treefile not found: {treefile}"
        )

    return treefile


def build_tree(trimmed_path, output_prefix, threads=8,
               bootstrap=1000, log_path=None, fast=False):
    """
    Build a maximum likelihood tree using IQ-TREE2.

    Uses -m TEST for automatic model selection and ultrafast bootstrap.

    Args:
        trimmed_path:  path to trimmed alignment FASTA
        output_prefix: prefix for IQ-TREE2 output files
                       (treefile will be output_prefix.treefile)
        threads:       number of threads (default: 8)
        bootstrap:     number of ultrafast bootstrap replicates (default: 1000)
        log_path:      path to append IQ-TREE2 log output (optional)

    Returns:
        path to the output treefile on success

    Raises:
        TooFewSequencesError if trimmed_path has fewer than MIN_SEQUENCES
        ToolNotFoundError if IQ-TREE2 is not on PATH
        subprocess.CalledProcessError if IQ-TREE2 exits with non-zero status
    """
    n = count_sequences(trimmed_path)
    if n < MIN_SEQUENCES:
        raise TooFewSequencesError(
            f"{os.path.basename(trimmed_path)}: {n} sequence(s) "
            f"(minimum {MIN_SEQUENCES} required for a meaningful tree)"
        )

    os.makedirs(os.path.dirname(output_prefix), exist_ok=True)

    binary, _ = check_iqtree()

    if fast:
        cmd = [
            binary,
            '-s',       trimmed_path,
            '--prefix', output_prefix,
            '-m',       'LG+G4',
            '--fast',
            '-T',       str(threads),
            '--quiet',
            '--redo',
        ]
    else:
        cmd = [
            binary,
            '-s',       trimmed_path,
            '--prefix', output_prefix,
            '-m',       'TEST',
            '-B',       str(bootstrap),
            '-T',       str(threads),
            '--quiet',
            '--redo',
        ]

    if log_path:
        with open(log_path, 'a') as log:
            subprocess.run(cmd, stderr=log, check=True)
    else:
        subprocess.run(cmd, stderr=subprocess.DEVNULL, check=True)

    treefile = f"{output_prefix}.treefile"
    if not os.path.exists(treefile):
        raise FileNotFoundError(
            f"IQ-TREE2 completed but treefile not found: {treefile}"
        )

    return treefile
