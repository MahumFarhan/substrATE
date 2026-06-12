"""
Extract sequences for target substrate families from each genome's
uniInput.faa, grouped by top-level family (e.g. all GH16 subfamilies
together). Adds sample name and localisation to sequence headers.

Also appends characterised reference sequences from ref_seq_dir,
using ref_metadata to assign family, substrate and localisation.
Only reference sequences matching the current substrate are included.
"""
import os
import re
import pandas as pd
from Bio import SeqIO


# ── Constants ─────────────────────────────────────────────────────────────────

# dbCAN always produces uniInput.faa — this is fixed, not user-controlled
DBCAN_FAA = 'uniInput.faa'

# Accepted extensions for user-provided protein FASTA files
# (reference sequences)
PROTEIN_EXTENSIONS = {'.faa', '.fasta', '.fa'}


# ── Helpers ───────────────────────────────────────────────────────────────────

def is_protein_fasta(filename):
    """Return True if filename has a recognised protein FASTA extension."""
    return any(filename.endswith(ext) for ext in PROTEIN_EXTENSIONS)


def top_family(subfamily):
    """
    Extract top-level family from subfamily annotation.

    Examples:
        'GH16_3'         -> 'GH16'
        'GH16_3|GH16_5'  -> 'GH16'
        '-'              -> 'unknown'

    Args:
        subfamily: subfamily annotation string from family_hits.tsv

    Returns:
        Top-level family string, or 'unknown' if not parseable
    """
    if pd.isna(subfamily) or str(subfamily).strip() in ['-', 'nan', '']:
        return 'unknown'
    first = str(subfamily).split('|')[0].strip()
    match = re.match(r'([A-Za-z]+[0-9]+)', first)
    return match.group(1) if match else first


def _clean_id(seq_id):
    """Replace characters illegal in FASTA headers with underscores."""
    return seq_id.replace('|', '_').replace('=', '_')


# ── Reference metadata loading ────────────────────────────────────────────────

def load_ref_metadata(ref_metadata, families):
    """
    Load reference sequence metadata filtered to the current substrate.

    Args:
        ref_metadata: path to reference_metadata.tsv
        families:     collection of family names to include (e.g. ['GH16', 'GH3'])

    Returns:
        dict mapping accession string to metadata row (pandas Series),
        or empty dict if file does not exist
    """
    if not os.path.exists(ref_metadata):
        print(f"WARNING: No reference metadata found at {ref_metadata}")
        return {}

    ref_meta = pd.read_csv(ref_metadata, sep='\t')
    ref_meta = ref_meta[ref_meta['family'].astype(str).isin(families)].copy()
    ref_lookup = {str(row['accession']): row
                  for _, row in ref_meta.iterrows()}
    print(f"Loaded {len(ref_meta)} reference sequences for families: {sorted(families)}")
    return ref_lookup


# ── Genome sequence extraction ────────────────────────────────────────────────

def extract_genome_sequences(cgc_output_dir, hits_df):
    """
    Extract sequences for all hits from each genome's uniInput.faa.

    Uses SeqIO.index() for memory-efficient access to large FAA files.
    Index is always closed cleanly via try/finally even if an error occurs.

    Args:
        cgc_output_dir: path to directory containing per-sample dbCAN output
        hits_df:        DataFrame of family hits for the current substrate

    Returns:
        dict mapping top-level family string to list of SeqRecord objects
    """
    family_seqs = {}

    for sample_dir in sorted(os.listdir(cgc_output_dir)):
        full_path = os.path.join(cgc_output_dir, sample_dir)
        faa_file  = os.path.join(full_path, DBCAN_FAA)

        if not os.path.isdir(full_path) or not os.path.exists(faa_file):
            continue

        sample      = sample_dir.replace('output_', '')
        sample_hits = hits_df[hits_df['sample'] == sample]

        if sample_hits.empty:
            print(f"  No hits for {sample}")
            continue

        seq_index = SeqIO.index(faa_file, 'fasta')
        try:
            for _, row in sample_hits.iterrows():
                gene_id      = str(row['Gene ID'])
                family       = str(row['top_family'])
                subfam       = str(row['subfamily_annotation']).replace('|', '_')
                localisation = str(row['localisation'])

                if gene_id not in seq_index:
                    print(f"  WARNING: {gene_id} not found in {faa_file}")
                    continue

                record             = seq_index[gene_id]
                record.id          = (f"{sample}__{_clean_id(gene_id)}"
                                      f"__{subfam}__{localisation}")
                record.description = ''

                if family not in family_seqs:
                    family_seqs[family] = []
                family_seqs[family].append(record)
        finally:
            seq_index.close()

    return family_seqs


# ── Reference sequence appending ──────────────────────────────────────────────


def subsample_by_relevance(genomic_records, ref_records, max_ref_seqs,
                           threads=8):
    """
    Subsample reference sequences by phylogenetic relevance to genomic
    sequences. Runs a fast MAFFT alignment, then keeps the top
    max_ref_seqs references ranked by maximum pairwise identity to
    any genomic sequence.
    Args:
        genomic_records: list of SeqRecord objects (genomic sequences)
        ref_records:     list of SeqRecord objects (reference sequences)
        max_ref_seqs:    maximum number of references to keep
        threads:         number of threads for MAFFT
    Returns:
        list of subsampled SeqRecord objects
    """
    import tempfile
    import subprocess
    if len(ref_records) <= max_ref_seqs:
        return ref_records
    # Write combined FASTA to temp file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.faa',
                                    delete=False) as tmp_in:
        tmp_in_path = tmp_in.name
        for rec in genomic_records + ref_records:
            tmp_in.write(f'>{rec.id}\n{str(rec.seq)}\n')
    tmp_aln_path = tmp_in_path + '.aln'
    try:
        # Fast MAFFT alignment
        cmd = ['mafft', '--retree', '1', '--maxiterate', '0',
               '--quiet', '--thread', str(threads), tmp_in_path]
        with open(tmp_aln_path, 'w') as aln_out:
            subprocess.run(cmd, stdout=aln_out, stderr=subprocess.DEVNULL,
                           check=True)
        # Parse alignment
        aln = {rec.id: str(rec.seq).upper()
               for rec in SeqIO.parse(tmp_aln_path, 'fasta')}
        genomic_ids = {rec.id for rec in genomic_records}
        ref_ids     = [rec.id for rec in ref_records]
        # Compute max pairwise identity of each ref to any genomic seq
        def _identity(seq_a, seq_b):
            matches = sum(a == b for a, b in zip(seq_a, seq_b)
                          if a != '-' and b != '-')
            comparable = sum(1 for a, b in zip(seq_a, seq_b)
                             if a != '-' and b != '-')
            return matches / comparable if comparable > 0 else 0.0
        scores = {}
        for ref_id in ref_ids:
            if ref_id not in aln:
                scores[ref_id] = 0.0
                continue
            ref_seq = aln[ref_id]
            scores[ref_id] = max(
                (_identity(ref_seq, aln[g]) for g in genomic_ids
                 if g in aln),
                default=0.0
            )
        # Rank and select top max_ref_seqs
        ranked = sorted(ref_ids, key=lambda r: scores[r], reverse=True)
        selected_ids = set(ranked[:max_ref_seqs])
        selected = [r for r in ref_records if r.id in selected_ids]
        print(f"    Relevance subsampling: kept {len(selected)}/{len(ref_records)} "
              f"reference sequences")
        return selected
    except Exception as e:
        print(f"    WARNING: relevance subsampling failed ({e}), "
              f"falling back to random subsampling")
        import random
        return random.sample(ref_records, max_ref_seqs)
    finally:
        for p in [tmp_in_path, tmp_aln_path]:
            if os.path.exists(p):
                os.remove(p)

def append_reference_sequences(family_seqs, ref_seq_dir, ref_lookup, max_ref_seqs=None, seed=None, ref_mode="diverse", threads=8):
    """
    Append characterised reference sequences to the family_seqs dict.

    Only sequences present in ref_lookup (i.e. matching the current
    substrate) are included. Others are skipped with a warning.

    Args:
        family_seqs: dict mapping family -> list of SeqRecords,
                     as returned by extract_genome_sequences().
                     Modified in place.
        ref_seq_dir: path to directory containing reference FASTA files
        ref_lookup:  dict from load_ref_metadata()

    Returns:
        family_seqs with reference sequences appended
    """
    if not ref_lookup or not os.path.exists(ref_seq_dir):
        print("\nNo reference sequences to append.")
        return family_seqs

    print("\nAppending reference sequences...")
    import random as _random
    ref_by_family = {}
    for ref_file in sorted(os.listdir(ref_seq_dir)):
        if not is_protein_fasta(ref_file):
            continue
        ref_path = os.path.join(ref_seq_dir, ref_file)
        for record in SeqIO.parse(ref_path, 'fasta'):
            accession = record.id.split()[0].split('|')[0]
            if accession not in ref_lookup:
                continue
            meta   = ref_lookup[accession]
            # Skip references with unknown activity (no EC number)
            _ec = str(meta.get('ec_numbers', '')).strip()
            if not _ec or _ec in ('nan', '-', ''):
                continue
            family = str(meta['family'])
            record.id = (f"Reference__{_clean_id(accession)}"
                         f"__{family}__characterised_reference")
            record.description = ''
            ref_by_family.setdefault(family, []).append(record)
    rng = _random.Random(seed)
    for family, records in ref_by_family.items():
        if max_ref_seqs is not None and len(records) > max_ref_seqs:
            if ref_mode == "relevant":
                genomic = family_seqs.get(family, [])
                records = subsample_by_relevance(
                    genomic, records, max_ref_seqs, threads=threads)
            else:
                records = rng.sample(records, max_ref_seqs)
            print(f"  {family}: subsampled to {len(records)} reference sequences [{ref_mode}]")
        if family not in family_seqs:
            family_seqs[family] = []
        family_seqs[family].extend(records)
        print(f"  {family}: added {len(records)} reference sequences")
    return family_seqs


# ── Deduplication and writing ─────────────────────────────────────────────────

def _deduplicate(records):
    """Remove duplicate sequences by ID, keeping first occurrence."""
    seen   = set()
    unique = []
    for r in records:
        if r.id not in seen:
            seen.add(r.id)
            unique.append(r)
    return unique


def write_family_fastas(family_seqs, seq_dir, substrate):
    """
    Write one deduplicated FASTA file per family to seq_dir.

    Args:
        family_seqs: dict mapping family -> list of SeqRecords
        seq_dir:     output directory for FASTA files
        substrate:   substrate name used as filename prefix

    Returns:
        dict mapping family string to output file path
    """
    os.makedirs(seq_dir, exist_ok=True)
    output_paths = {}

    print("\nSequences extracted per family:")
    for family, records in sorted(family_seqs.items()):
        unique   = _deduplicate(records)
        out_file = os.path.join(seq_dir, f"{substrate}_{family}.faa")
        SeqIO.write(unique, out_file, 'fasta')

        n_ref = sum(1 for r in unique if r.id.startswith('Reference__'))
        print(f"  {family}: {len(unique)} sequences "
              f"({len(unique) - n_ref} genomic + {n_ref} reference)"
              f" -> {out_file}")

        output_paths[family] = out_file

    return output_paths


# ── Main entry point ──────────────────────────────────────────────────────────

def extract_sequences(cgc_output_dir, hits_df, output_dir, substrate,
                      ref_metadata=None, ref_seq_dir=None, max_ref_seqs=None, seed=None, ref_mode="diverse", threads=8):
    """
    Full sequence extraction pipeline for one substrate.

    Extracts genomic sequences from dbCAN uniInput.faa files, optionally
    appends reference sequences, and writes one FASTA per family.

    Args:
        cgc_output_dir: path to dbCAN cgc_output directory
        hits_df:        family hits DataFrame (from classify_pul.process_samples)
        output_dir:     base output directory
        substrate:      substrate name
        ref_metadata:   path to reference_metadata.tsv (optional)
        ref_seq_dir:    path to reference sequences directory (optional)

    Returns:
        dict mapping family string to output FASTA path
    """
    hits_df = hits_df.copy()
    hits_df['top_family'] = hits_df['matched_family']

    print(f"Substrate: {substrate}")
    print(f"Total hits: {len(hits_df)}")

    # Extract genomic sequences first so we know which families are present
    family_seqs = extract_genome_sequences(cgc_output_dir, hits_df)
    # Load reference metadata filtered to families present in this run
    ref_lookup = {}
    if ref_metadata:
        ref_lookup = load_ref_metadata(ref_metadata, set(family_seqs.keys()))

    # Append reference sequences if provided
    if ref_seq_dir:
        family_seqs = append_reference_sequences(
            family_seqs, ref_seq_dir, ref_lookup, max_ref_seqs=max_ref_seqs, seed=seed, ref_mode=ref_mode, threads=threads)

    # Write output FASTAs
    seq_dir = os.path.join(output_dir, 'sequences')
    output_paths = write_family_fastas(family_seqs, seq_dir, substrate)

    print(f"\nDone. Sequences written to {seq_dir}")
    return output_paths
