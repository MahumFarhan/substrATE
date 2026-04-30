"""
EPA-RAxML sequence placement module.

For each CAZyme family with genomic hits, places query sequences
onto a pre-built reference tree using EPA-ng (evolutionary placement
algorithm). Falls back to de novo IQ-TREE2 tree building if no
reference tree exists or if the family has too few reference sequences.

Workflow per family:
  1. Check if reference tree exists in reference_trees/{family}/
  2. If yes:
       a. Build HMM profile from reference alignment (hmmbuild)
       b. Align query sequences to profile (hmmalign)
       c. Run epa-ng placement onto reference tree
       d. Write {family}.jplace and convert to {family}.newick
       e. Write placement_scores.tsv
  3. If no:
       Fall back to de novo: MAFFT -> trimAl -> IQ-TREE2
       Print WARNING with reason

Output:
  {output_dir}/placements/{family}/
      {family}.jplace              <- full EPA placement data
      {family}.newick              <- best placement as Newick
      {family}.placement_scores.tsv
  {output_dir}/trees/{family}/     <- de novo fallback trees
"""
import os
import re
import json
import subprocess
import pandas as pd
from Bio import SeqIO


# ── Constants ─────────────────────────────────────────────────────────────────

_PKG_DIR      = os.path.dirname(os.path.abspath(__file__))
REF_TREES_DIR = os.path.join(_PKG_DIR, 'data', 'reference_trees')
MIN_REF_SEQS  = 4


# ── Exceptions ────────────────────────────────────────────────────────────────

class ToolNotFoundError(Exception):
    """Raised when a required external tool is not found on PATH."""
    pass


class TooFewSequencesError(Exception):
    """Raised when too few sequences are available for placement."""
    pass


# ── Tool availability ─────────────────────────────────────────────────────────

def check_epang():
    """
    Check that epa-ng is available and return its version string.

    Raises:
        ToolNotFoundError if epa-ng is not found on PATH
    """
    try:
        result = subprocess.run(
            ['epa-ng', '--version'],
            capture_output=True, text=True
        )
        version = result.stdout.strip() or result.stderr.strip()
        return version
    except FileNotFoundError:
        raise ToolNotFoundError(
            "epa-ng not found on PATH. "
            "Install with: conda install -c bioconda epa-ng"
        )


def check_hmmer():
    """
    Check that HMMER (hmmbuild + hmmalign) is available.

    Raises:
        ToolNotFoundError if hmmbuild or hmmalign not found on PATH
    """
    for tool in ['hmmbuild', 'hmmalign']:
        try:
            subprocess.run(
                [tool, '-h'],
                capture_output=True
            )
        except FileNotFoundError:
            raise ToolNotFoundError(
                f"{tool} not found on PATH. "
                f"Install with: conda install -c bioconda hmmer"
            )
    return "HMMER available"


# ── Reference tree lookup ─────────────────────────────────────────────────────

def get_reference_files(family, ref_trees_dir=None):
    """
    Find reference tree and alignment files for a family.

    Args:
        family:        family name string (e.g. 'GH16')
        ref_trees_dir: path to reference_trees directory
                       (uses default if None)

    Returns:
        dict with keys: treefile, trimmed_aln, ref_faa
        or None if reference tree does not exist
    """
    if ref_trees_dir is None:
        ref_trees_dir = REF_TREES_DIR

    family_dir = os.path.join(ref_trees_dir, family)
    treefile   = os.path.join(family_dir, f'{family}.ref.treefile')
    trimmed    = os.path.join(family_dir, f'{family}.ref.trim')
    ref_faa    = os.path.join(family_dir, f'{family}.ref.faa')

    if not os.path.exists(treefile):
        return None
    if not os.path.exists(trimmed):
        return None

    n_ref = sum(1 for _ in SeqIO.parse(ref_faa, 'fasta')) \
        if os.path.exists(ref_faa) else 0

    if n_ref < MIN_REF_SEQS:
        return None

    return {
        'treefile':    treefile,
        'trimmed_aln': trimmed,
        'ref_faa':     ref_faa,
        'n_ref_seqs':  n_ref,
    }


# ── HMM profile alignment ─────────────────────────────────────────────────────

def build_hmm_profile(trimmed_aln, hmm_path, log_path=None):
    """
    Build an HMM profile from the reference alignment using hmmbuild.

    Args:
        trimmed_aln: path to trimmed reference alignment
        hmm_path:    path to write HMM profile
        log_path:    path to append hmmbuild log (optional)

    Returns:
        hmm_path on success
    """
    cmd = ['hmmbuild', '--amino', hmm_path, trimmed_aln]

    if log_path:
        with open(log_path, 'a') as log:
            subprocess.run(cmd, stdout=log, stderr=log, check=True)
    else:
        subprocess.run(cmd, capture_output=True, check=True)

    return hmm_path


def align_queries_to_profile(hmm_path, query_faa, query_aln_path,
                              log_path=None):
    """
    Align query sequences to the HMM profile using hmmalign.

    Uses --trim to remove columns not in the reference alignment,
    producing an alignment compatible with the reference tree.

    Args:
        hmm_path:       path to HMM profile
        query_faa:      path to query sequences FASTA
        query_aln_path: path to write aligned query sequences
        log_path:       path to append hmmalign log (optional)

    Returns:
        query_aln_path on success
    """
    cmd = [
        'hmmalign',
        '--amino',
        '--trim',
        '--outformat', 'afa',
        '-o', query_aln_path,
        hmm_path,
        query_faa,
    ]

    if log_path:
        with open(log_path, 'a') as log:
            subprocess.run(cmd, stdout=log, stderr=log, check=True)
    else:
        subprocess.run(cmd, capture_output=True, check=True)

    return query_aln_path


# ── EPA-ng placement ──────────────────────────────────────────────────────────

def run_epang(query_aln, ref_aln, treefile, output_dir,
              threads=8, log_path=None):
    """
    Run epa-ng to place query sequences onto the reference tree.

    Args:
        query_aln:  path to hmmalign-aligned query sequences
        ref_aln:    path to reference alignment (trimmed)
        treefile:   path to reference treefile
        output_dir: directory to write epa-ng output
        threads:    number of threads
        log_path:   path to append epa-ng log (optional)

    Returns:
        path to output .jplace file
    """
    os.makedirs(output_dir, exist_ok=True)

    cmd = [
        'epa-ng',
        '--tree',     treefile,
        '--ref-msa',  ref_aln,
        '--query',    query_aln,
        '--outdir',   output_dir,
        '--threads',  str(threads),
        '--redo',
    ]

    if log_path:
        with open(log_path, 'a') as log:
            subprocess.run(cmd, stdout=log, stderr=log, check=True)
    else:
        subprocess.run(cmd, capture_output=True, check=True)

    # epa-ng names the output file epa_result.jplace
    jplace = os.path.join(output_dir, 'epa_result.jplace')
    if not os.path.exists(jplace):
        raise FileNotFoundError(
            f"epa-ng completed but jplace not found: {jplace}")

    return jplace


# ── jplace processing ─────────────────────────────────────────────────────────

def parse_jplace(jplace_path):
    """
    Parse a .jplace file and return placement data.

    Args:
        jplace_path: path to .jplace file

    Returns:
        dict with keys:
          tree:       reference tree with edge numbers
          placements: list of placement dicts, each with:
            name:     query sequence name
            edge:     best placement edge number
            like_weight: likelihood weight ratio (confidence)
            distal:   distal length
            pendant:  pendant length
    """
    with open(jplace_path) as f:
        data = json.load(f)

    # Field order: edge_num, likelihood, like_weight_ratio,
    #              distal_length, pendant_length
    fields = data.get('fields', [])
    edge_idx   = fields.index('edge_num')
    lwr_idx    = fields.index('like_weight_ratio')
    distal_idx = fields.index('distal_length')
    pendant_idx = fields.index('pendant_length')

    placements = []
    for p in data.get('placements', []):
        # Take best placement (highest like_weight_ratio)
        best = max(p['p'], key=lambda x: x[lwr_idx])
        name = p.get('n', p.get('name', ['unknown']))[0]
        placements.append({
            'name':         name,
            'edge':         best[edge_idx],
            'like_weight':  best[lwr_idx],
            'distal':       best[distal_idx],
            'pendant':      best[pendant_idx],
        })

    return {
        'tree':       data.get('tree', ''),
        'placements': placements,
    }


def jplace_to_newick(jplace_data, output_path):
    """
    Convert jplace best placements to a Newick tree by inserting
    query sequences at their best placement positions.

    Each query sequence is added as a leaf node attached to its
    best placement edge with the pendant length as branch length.

    Args:
        jplace_data: dict from parse_jplace()
        output_path: path to write Newick file

    Returns:
        output_path on success
    """
    tree_str = jplace_data['tree']

    # Replace edge number annotations {N} with placement info
    for placement in jplace_data['placements']:
        edge    = placement['edge']
        name    = placement['name'].replace('|', '_')
        pendant = placement['pendant']

        # Find the edge in the tree string and insert query as sister
        # Edge numbers appear as {N} in the epa-ng tree format
        pattern     = r'\{' + str(edge) + r'\}'
        replacement = f'{{{edge},[{name}:{pendant:.6f}]}}'
        tree_str    = re.sub(pattern, replacement, tree_str, count=1)

    # Clean up the tree string — remove edge number annotations
    # and format query insertions as proper Newick
    tree_str = re.sub(r'\{[0-9]+\}', '', tree_str)
    tree_str = re.sub(r'\[([^\]]+)\]', r'\1', tree_str)

    with open(output_path, 'w') as f:
        f.write(tree_str.strip() + '\n')

    return output_path


def write_placement_scores(jplace_data, output_path, family):
    """
    Write a TSV of placement confidence scores for all query sequences.

    Args:
        jplace_data: dict from parse_jplace()
        output_path: path to write TSV
        family:      family name for the family column
    """
    rows = []
    for p in jplace_data['placements']:
        rows.append({
            'sequence':         p['name'],
            'family':           family,
            'best_edge':        p['edge'],
            'like_weight_ratio': p['like_weight'],
            'distal_length':    p['distal'],
            'pendant_length':   p['pendant'],
            'placement_method': 'epa-ng',
        })

    df = pd.DataFrame(rows)
    df.to_csv(output_path, sep='\t', index=False)
    return output_path


# ── Main placement function ───────────────────────────────────────────────────

def place_or_build(family, query_faa, output_dir, threads=8,
                   ref_trees_dir=None, log_dir=None):
    """
    Place query sequences onto reference tree, or fall back to
    de novo tree building if no reference tree exists.

    Args:
        family:        family name string (e.g. 'GH16')
        query_faa:     path to query sequences FASTA
        output_dir:    base output directory for this substrate
        threads:       number of threads
        ref_trees_dir: path to reference_trees directory
                       (uses default if None)
        log_dir:       directory for tool log files (optional)

    Returns:
        dict with keys:
          method:    'placement' or 'denovo'
          treefile:  path to output treefile or newick
          jplace:    path to jplace file (None if de novo)
          n_query:   number of query sequences placed
          family:    family name
    """
    from substrate.align import align, TooFewSequencesError as AlignErr
    from substrate.trim  import trim
    from substrate.tree  import build_tree

    result = {
        'family':   family,
        'method':   None,
        'treefile': None,
        'jplace':   None,
        'n_query':  0,
    }

    # Count query sequences
    n_query = sum(1 for _ in SeqIO.parse(query_faa, 'fasta'))
    result['n_query'] = n_query

    if n_query == 0:
        print(f"  {family:<12} WARNING: no query sequences")
        return result

    # Set up log paths
    hmm_log   = os.path.join(log_dir, f'{family}_hmmer.log') \
        if log_dir else None
    epa_log   = os.path.join(log_dir, f'{family}_epang.log') \
        if log_dir else None
    align_log = os.path.join(log_dir, f'{family}_mafft.log') \
        if log_dir else None
    trim_log  = os.path.join(log_dir, f'{family}_trimal.log') \
        if log_dir else None
    tree_log  = os.path.join(log_dir, f'{family}_iqtree.log') \
        if log_dir else None

    # Check for reference tree
    ref_files = get_reference_files(family, ref_trees_dir)

    if ref_files:
        # ── EPA placement ─────────────────────────────────────────
        print(f"  {family:<12} placing {n_query} sequences onto "
              f"reference tree "
              f"({ref_files['n_ref_seqs']} ref seqs)...")

        placement_dir = os.path.join(output_dir, 'placements', family)
        os.makedirs(placement_dir, exist_ok=True)

        hmm_path       = os.path.join(placement_dir,
                                      f'{family}.ref.hmm')
        query_aln_path = os.path.join(placement_dir,
                                      f'{family}.query.aln')
        jplace_path    = os.path.join(placement_dir,
                                      f'{family}.jplace')
        newick_path    = os.path.join(placement_dir,
                                      f'{family}.newick')
        scores_path    = os.path.join(placement_dir,
                                      f'{family}.placement_scores.tsv')

        try:
            # Build HMM profile from reference alignment
            build_hmm_profile(
                ref_files['trimmed_aln'], hmm_path,
                log_path=hmm_log)

            # Align queries to profile
            align_queries_to_profile(
                hmm_path, query_faa, query_aln_path,
                log_path=hmm_log)

            # Run EPA placement
            epa_output_dir = os.path.join(placement_dir, 'epa_output')
            raw_jplace = run_epang(
                query_aln_path,
                ref_files['trimmed_aln'],
                ref_files['treefile'],
                epa_output_dir,
                threads=threads,
                log_path=epa_log,
            )

            # Parse and convert jplace
            jplace_data = parse_jplace(raw_jplace)

            # Copy jplace to placement dir
            import shutil
            shutil.copy(raw_jplace, jplace_path)

            # Convert to Newick
            jplace_to_newick(jplace_data, newick_path)

            # Write placement scores
            write_placement_scores(jplace_data, scores_path, family)

            n_placed = len(jplace_data['placements'])
            print(f"  {family:<12} placed {n_placed}/{n_query} "
                  f"sequences -> {newick_path}")

            result.update({
                'method':   'placement',
                'treefile': newick_path,
                'jplace':   jplace_path,
            })

        except Exception as e:
            print(f"  {family:<12} WARNING: EPA placement failed "
                  f"({e}), falling back to de novo")
            ref_files = None  # trigger de novo fallback

    if not ref_files:
        # ── De novo fallback ──────────────────────────────────────
        reason = "no reference tree" if \
            get_reference_files(family, ref_trees_dir) is None \
            else "placement failed"
        print(f"  {family:<12} WARNING: {reason}, "
              f"building de novo tree ({n_query} sequences)...")

        tree_dir    = os.path.join(output_dir, 'trees', family)
        align_dir   = os.path.join(output_dir, 'alignments')
        trimmed_dir = os.path.join(output_dir, 'trimmed')
        os.makedirs(tree_dir,    exist_ok=True)
        os.makedirs(align_dir,   exist_ok=True)
        os.makedirs(trimmed_dir, exist_ok=True)

        aligned_path = os.path.join(align_dir,
                                    f'{family}.aln')
        trimmed_path = os.path.join(trimmed_dir,
                                    f'{family}.trim')
        tree_prefix  = os.path.join(tree_dir, family)

        try:
            align(query_faa, aligned_path, threads=threads,
                  log_path=align_log)
            trim(aligned_path, trimmed_path, log_path=trim_log)
            treefile = build_tree(trimmed_path, tree_prefix,
                                  threads=threads, log_path=tree_log)
            print(f"  {family:<12} de novo tree -> {treefile}")
            result.update({
                'method':   'denovo',
                'treefile': treefile,
            })
        except AlignErr as e:
            print(f"  {family:<12} skipping: {e}")

    return result


def place_all_families(fasta_paths, output_dir, threads=8,
                       ref_trees_dir=None, log_dir=None):
    """
    Run placement or de novo tree building for all families
    in a substrate run.

    Args:
        fasta_paths:   dict mapping family -> query FASTA path
                       (from extract_seqs.extract_sequences())
        output_dir:    substrate output directory
        threads:       number of threads
        ref_trees_dir: path to reference_trees directory
        log_dir:       directory for tool log files

    Returns:
        dict mapping family -> result dict from place_or_build()
    """
    check_epang()
    check_hmmer()

    if log_dir:
        os.makedirs(log_dir, exist_ok=True)

    results          = {}
    n_placement      = 0
    n_denovo         = 0
    skipped_families = []

    for family, fasta_path in sorted(fasta_paths.items()):
        result = place_or_build(
            family=family,
            query_faa=fasta_path,
            output_dir=output_dir,
            threads=threads,
            ref_trees_dir=ref_trees_dir,
            log_dir=log_dir,
        )
        results[family] = result

        if result['method'] == 'placement':
            n_placement += 1
        elif result['method'] == 'denovo':
            n_denovo += 1
        elif result['method'] is None:
            skipped_families.append(family)

    print(f"\nTree/placement summary:")
    print(f"  EPA placement:  {n_placement} families")
    print(f"  De novo trees:  {n_denovo} families")
    if skipped_families:
        print(f"  Skipped:        {len(skipped_families)} families "
              f"({', '.join(skipped_families)})")

    return results
