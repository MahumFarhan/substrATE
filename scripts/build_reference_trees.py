"""
Build IQ-TREE2 reference trees from CAZy characterised sequences.

For each family in substrate/data/reference_seqs/by_family/, runs:
  MAFFT alignment -> trimAl trimming -> IQ-TREE2 tree building

Output is written to substrate/data/reference_trees/{family}/
and used by place_sequences.py for EPA-RAxML placement.

Usage:
    python scripts/build_reference_trees.py \\
        [--families GH16 GH5 GH3] \\
        [--max_seqs_per_family 200] \\
        [--threads 8] \\
        [--min_seqs 4] \\
        [--force]

If --families is omitted, all families in by_family/ are processed.
Use --force to rebuild trees that already exist.

A build_log.tsv is written to reference_trees/ recording what was
built, how many sequences were used, and when.
"""
import os
import sys
import argparse
import datetime
import pandas as pd

# Add parent directory to path so we can import substrate modules
sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))

from Bio import SeqIO
from substrate.align import align, check_mafft, TooFewSequencesError
from substrate.trim  import trim,  check_trimal
from substrate.tree  import build_tree, check_iqtree


# ── Constants ─────────────────────────────────────────────────────────────────

_SCRIPTS_DIR   = os.path.dirname(os.path.abspath(__file__))
_PKG_DIR       = os.path.dirname(_SCRIPTS_DIR)
REF_SEQS_DIR   = os.path.join(_PKG_DIR, 'substrate', 'data',
                               'reference_seqs', 'by_family')
REF_TREES_DIR  = os.path.join(_PKG_DIR, 'substrate', 'data',
                               'reference_trees')
MIN_SEQS_DEFAULT = 4


# ── Subsampling (Option C) ────────────────────────────────────────────────────

def subsample_by_subfamily_diversity(records, max_seqs):
    """
    Subsample sequences to max_seqs using proportional subfamily
    representation with organism diversity within each subfamily.

    Sequence headers are expected in the format:
        accession|family|subfamily|organism

    Args:
        records:  list of SeqRecord objects
        max_seqs: maximum number of sequences to keep

    Returns:
        list of subsampled SeqRecord objects
    """
    if len(records) <= max_seqs:
        return records

    # Group by subfamily
    subfamily_groups = {}
    for record in records:
        parts  = record.id.split('|')
        subfam = parts[2] if len(parts) > 2 else 'unknown'
        org    = parts[3] if len(parts) > 3 else 'unknown'
        if subfam not in subfamily_groups:
            subfamily_groups[subfam] = {}
        if org not in subfamily_groups[subfam]:
            subfamily_groups[subfam][org] = []
        subfamily_groups[subfam][org].append(record)

    n_subfamilies = len(subfamily_groups)
    total_seqs    = len(records)

    # Allocate slots proportionally across subfamilies
    subfam_sizes = {sf: sum(len(orgs) for orgs in org_dict.values())
                    for sf, org_dict in subfamily_groups.items()}
    slots = {}
    for sf, size in subfam_sizes.items():
        slots[sf] = max(1, int(max_seqs * size / total_seqs))

    # Adjust to hit exact max_seqs
    while sum(slots.values()) > max_seqs:
        largest = max(slots, key=lambda s: slots[s])
        slots[largest] -= 1
    while sum(slots.values()) < max_seqs:
        for sf in sorted(slots, key=lambda s: -subfam_sizes[s]):
            if slots[sf] < subfam_sizes[sf]:
                slots[sf] += 1
                if sum(slots.values()) >= max_seqs:
                    break

    # Within each subfamily, round-robin across organisms
    selected = []
    for subfam, n_slots in slots.items():
        org_lists = list(subfamily_groups[subfam].values())
        i = 0
        n_added = 0
        while n_added < n_slots and any(org_lists):
            org_idx = i % len(org_lists)
            if org_lists[org_idx]:
                selected.append(org_lists[org_idx].pop(0))
                n_added += 1
            i += 1
            if all(not ol for ol in org_lists):
                break

    print(f"    Subsampled {len(records)} -> {len(selected)} sequences "
          f"({n_subfamilies} subfamilies represented)")
    return selected


# ── Tree building ─────────────────────────────────────────────────────────────

def build_reference_tree(family, faa_path, output_dir, threads=8,
                         max_seqs=None, min_seqs=MIN_SEQS_DEFAULT,
                         force=False, log_dir=None):
    """
    Build a reference tree for one CAZyme family.

    Args:
        family:     family name string (e.g. 'GH16')
        faa_path:   path to input FASTA file
        output_dir: base reference_trees directory
        threads:    number of threads for MAFFT and IQ-TREE2
        max_seqs:   maximum sequences after subsampling (None = no limit)
        min_seqs:   minimum sequences required to build tree
        force:      if True, rebuild existing trees
        log_dir:    directory for tool log files

    Returns:
        dict with keys: family, n_seqs, status, treefile, message
    """
    family_dir = os.path.join(output_dir, family)
    treefile   = os.path.join(family_dir, f'{family}.ref.treefile')

    result = {
        'family':   family,
        'n_seqs':   0,
        'status':   'unknown',
        'treefile': None,
        'message':  '',
    }

    # Check if tree already exists
    if os.path.exists(treefile) and not force:
        n = sum(1 for _ in SeqIO.parse(faa_path, 'fasta'))
        print(f"  {family:<12} skipping — tree exists "
              f"({n} seqs, use --force to rebuild)")
        result.update({'n_seqs': n, 'status': 'skipped',
                       'treefile': treefile,
                       'message': 'tree already exists'})
        return result

    os.makedirs(family_dir, exist_ok=True)

    # Load sequences
    records = list(SeqIO.parse(faa_path, 'fasta'))
    n_input = len(records)

    if n_input < min_seqs:
        print(f"  {family:<12} skipping — only {n_input} sequences "
              f"(minimum {min_seqs}), will use de novo approach")
        result.update({'n_seqs': n_input, 'status': 'skipped_too_few',
                       'message': f'only {n_input} sequences, '
                                  f'minimum {min_seqs}'})
        return result

    # Subsample if needed
    if max_seqs and n_input > max_seqs:
        records = subsample_by_subfamily_diversity(records, max_seqs)

    n_seqs = len(records)
    result['n_seqs'] = n_seqs

    # Write working FASTA (post-subsampling)
    ref_faa  = os.path.join(family_dir, f'{family}.ref.faa')
    ref_aln  = os.path.join(family_dir, f'{family}.ref.aln')
    ref_trim = os.path.join(family_dir, f'{family}.ref.trim')
    tree_pfx = os.path.join(family_dir, f'{family}.ref')

    SeqIO.write(records, ref_faa, 'fasta')

    # Log paths
    mafft_log  = os.path.join(log_dir, f'{family}_mafft.log') \
        if log_dir else None
    trimal_log = os.path.join(log_dir, f'{family}_trimal.log') \
        if log_dir else None
    iqtree_log = os.path.join(log_dir, f'{family}_iqtree.log') \
        if log_dir else None

    try:
        # Step 1: Align
        print(f"  {family:<12} aligning {n_seqs} sequences "
              f"(MAFFT)...")
        align(ref_faa, ref_aln, threads=threads,
              log_path=mafft_log)

        # Step 2: Trim
        print(f"  {family:<12} trimming (trimAl)...")
        trim(ref_aln, ref_trim, log_path=trimal_log)

        # Step 3: Build tree
        print(f"  {family:<12} building tree (IQ-TREE2)...")
        treefile_out = build_tree(
            ref_trim, tree_pfx,
            threads=threads,
            log_path=iqtree_log
        )

        print(f"  {family:<12} done -> {treefile_out}")
        result.update({
            'status':   'success',
            'treefile': treefile_out,
            'message':  f'{n_seqs} sequences used',
        })

    except TooFewSequencesError as e:
        print(f"  {family:<12} WARNING: {e}")
        result.update({'status':  'skipped_too_few',
                       'message': str(e)})
    except Exception as e:
        print(f"  {family:<12} ERROR: {e}")
        result.update({'status': 'failed', 'message': str(e)})

    return result


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Build IQ-TREE2 reference trees from CAZy '
                    'characterised sequences'
    )
    parser.add_argument(
        '--families', nargs='+', default=None,
        help='Families to build trees for (default: all in by_family/)'
    )
    parser.add_argument(
        '--max_seqs_per_family', type=int, default=None,
        help='Maximum sequences per family tree after subsampling '
             'using proportional subfamily representation. '
             'Default: no limit.'
    )
    parser.add_argument(
        '--threads', type=int, default=8,
        help='Number of threads for MAFFT and IQ-TREE2 (default: 8)'
    )
    parser.add_argument(
        '--min_seqs', type=int, default=MIN_SEQS_DEFAULT,
        help=f'Minimum sequences required to build a reference tree '
             f'(default: {MIN_SEQS_DEFAULT}). Families below this '
             f'threshold fall back to de novo tree building.'
    )
    parser.add_argument(
        '--force', action='store_true',
        help='Rebuild trees that already exist'
    )
    parser.add_argument(
        '--ref_seqs_dir', default=REF_SEQS_DIR,
        help=f'Path to by_family/ directory '
             f'(default: {REF_SEQS_DIR})'
    )
    parser.add_argument(
        '--output', default=REF_TREES_DIR,
        help=f'Output directory for reference trees '
             f'(default: {REF_TREES_DIR})'
    )
    args = parser.parse_args()

    # Check tools
    print("Checking required tools...")
    try:
        v = check_mafft()
        print(f"  ✓ MAFFT: {v}")
    except Exception as e:
        print(f"  ✗ MAFFT: {e}")
        sys.exit(1)
    try:
        v = check_trimal()
        print(f"  ✓ trimAl: {v}")
    except Exception as e:
        print(f"  ✗ trimAl: {e}")
        sys.exit(1)
    try:
        binary, v = check_iqtree()
        print(f"  ✓ IQ-TREE2 ({binary}): {v}")
    except Exception as e:
        print(f"  ✗ IQ-TREE2: {e}")
        sys.exit(1)

    # Find families to process
    if not os.path.exists(args.ref_seqs_dir):
        print(f"\nERROR: Reference sequences directory not found: "
              f"{args.ref_seqs_dir}")
        print(f"Run scripts/build_reference_db.py first.")
        sys.exit(1)

    available = [
        f[:-4] for f in os.listdir(args.ref_seqs_dir)
        if f.endswith('.faa')
    ]

    if args.families:
        families = [f for f in args.families if f in available]
        missing  = [f for f in args.families if f not in available]
        if missing:
            print(f"\nWARNING: Families not found in {args.ref_seqs_dir}: "
                  f"{missing}")
    else:
        families = sorted(available)

    if not families:
        print(f"\nNo families to process.")
        sys.exit(0)

    print(f"\nBuilding reference trees for {len(families)} families")
    if args.max_seqs_per_family:
        print(f"Max sequences per family: {args.max_seqs_per_family} "
              f"(Option C subsampling)")
    print(f"Minimum sequences required: {args.min_seqs}")
    print(f"Threads: {args.threads}")
    print(f"Output: {args.output}\n")

    os.makedirs(args.output, exist_ok=True)
    log_dir = os.path.join(args.output, 'logs')
    os.makedirs(log_dir, exist_ok=True)

    # Build trees
    results = []
    for family in families:
        faa_path = os.path.join(args.ref_seqs_dir, f'{family}.faa')
        result   = build_reference_tree(
            family=family,
            faa_path=faa_path,
            output_dir=args.output,
            threads=args.threads,
            max_seqs=args.max_seqs_per_family,
            min_seqs=args.min_seqs,
            force=args.force,
            log_dir=log_dir,
        )
        results.append(result)

    # Write build log
    log_df   = pd.DataFrame(results)
    log_df['timestamp'] = datetime.datetime.now().isoformat()
    log_path = os.path.join(args.output, 'build_log.tsv')

    if os.path.exists(log_path) and not args.force:
        existing = pd.read_csv(log_path, sep='\t')
        log_df   = pd.concat([existing, log_df],
                              ignore_index=True)

    log_df.to_csv(log_path, sep='\t', index=False)

    # Print summary
    n_success     = sum(1 for r in results if r['status'] == 'success')
    n_skipped_few = sum(1 for r in results
                        if r['status'] == 'skipped_too_few')
    n_skipped_ex  = sum(1 for r in results
                        if r['status'] == 'skipped')
    n_failed      = sum(1 for r in results if r['status'] == 'failed')

    print(f"\n{'='*50}")
    print(f"Reference tree build complete")
    print(f"{'='*50}")
    print(f"  Built:              {n_success}")
    print(f"  Skipped (exists):   {n_skipped_ex}")
    print(f"  Skipped (too few):  {n_skipped_few} "
          f"(will use de novo)")
    print(f"  Failed:             {n_failed}")
    print(f"  Log:                {log_path}")

    if n_skipped_few > 0:
        few_families = [r['family'] for r in results
                        if r['status'] == 'skipped_too_few']
        print(f"\nFamilies falling back to de novo: "
              f"{', '.join(few_families)}")

    if n_failed > 0:
        failed = [(r['family'], r['message']) for r in results
                  if r['status'] == 'failed']
        print(f"\nFailed families:")
        for fam, msg in failed:
            print(f"  {fam}: {msg}")


if __name__ == '__main__':
    main()
