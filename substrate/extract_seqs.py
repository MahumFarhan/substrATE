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

def load_ref_metadata(ref_metadata, substrate):
    """
    Load reference sequence metadata filtered to the current substrate.

    Args:
        ref_metadata: path to reference_metadata.tsv
        substrate:    substrate name to filter by

    Returns:
        dict mapping accession string to metadata row (pandas Series),
        or empty dict if file does not exist
    """
    if not os.path.exists(ref_metadata):
        print(f"WARNING: No reference metadata found at {ref_metadata}")
        return {}

    ref_meta = pd.read_csv(ref_metadata, sep='\t')
    ref_meta = ref_meta[ref_meta['substrate'] == substrate].copy()
    ref_lookup = {str(row['accession']): row
                  for _, row in ref_meta.iterrows()}
    print(f"Loaded {len(ref_meta)} reference sequences for {substrate}")
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

def append_reference_sequences(family_seqs, ref_seq_dir, ref_lookup):
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

    for ref_file in sorted(os.listdir(ref_seq_dir)):
        if not is_protein_fasta(ref_file):
            continue

        ref_path = os.path.join(ref_seq_dir, ref_file)
        for record in SeqIO.parse(ref_path, 'fasta'):
            accession = record.id.split()[0]

            if accession not in ref_lookup:
                print(f"  SKIPPING {accession} ({ref_file}) "
                      f"— not in metadata for current substrate")
                continue

            meta   = ref_lookup[accession]
            family = str(meta['family'])

            record.id = (f"Reference__{_clean_id(accession)}"
                         f"__{family}__characterised_reference")
            record.description = ''

            if family not in family_seqs:
                family_seqs[family] = []
            family_seqs[family].append(record)
            print(f"  Added {accession} -> {family} ({meta['label']})")

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
                      ref_metadata=None, ref_seq_dir=None):
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

    # Load reference metadata if provided
    ref_lookup = {}
    if ref_metadata:
        ref_lookup = load_ref_metadata(ref_metadata, substrate)

    # Extract genomic sequences
    family_seqs = extract_genome_sequences(cgc_output_dir, hits_df)

    # Append reference sequences if provided
    if ref_seq_dir:
        family_seqs = append_reference_sequences(
            family_seqs, ref_seq_dir, ref_lookup)

    # Write output FASTAs
    seq_dir = os.path.join(output_dir, 'sequences')
    output_paths = write_family_fastas(family_seqs, seq_dir, substrate)

    print(f"\nDone. Sequences written to {seq_dir}")
    return output_paths
