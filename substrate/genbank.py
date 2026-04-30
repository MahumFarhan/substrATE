"""
Extract CGC regions containing target substrate families and convert
to GenBank format for clinker comparison.

One GenBank file per CGC, organised into per-genome subdirectories
under <output_dir>/genbank/.

Only canonical_PUL CGCs with >= min_cazymes genes having
substrate-relevant activity annotations are included.

Activity relevance is determined by matching activity strings against
substrate-specific inclusion patterns defined in ACTIVITY_PATTERNS.
This keeps the approach data-driven and flexible across substrates.
To add a new substrate, add a new entry to ACTIVITY_PATTERNS only.
"""
import os
import shutil
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


# ── Constants ─────────────────────────────────────────────────────────────────

# Accepted extensions for user-provided nucleotide genome FASTA files
NUCLEOTIDE_EXTENSIONS = {'.fna', '.fasta', '.fa'}

# Path to activity patterns TSV relative to this file
_DATA_DIR              = os.path.join(os.path.dirname(__file__), 'data')
ACTIVITY_PATTERNS_FILE = os.path.join(_DATA_DIR, 'activity_patterns.tsv')


# ── Activity pattern loading ──────────────────────────────────────────────────

def load_activity_patterns(patterns_file=None):
    """
    Load activity relevance patterns from activity_patterns.tsv.

    Args:
        patterns_file: path to activity_patterns.tsv. If None, uses
                       the default file in substrate/data/

    Returns:
        tuple of (patterns_dict, reviewed_dict) where:
          patterns_dict maps substrate -> list of pattern strings
          reviewed_dict maps substrate -> bool (True if all reviewed)
    """
    if patterns_file is None:
        patterns_file = ACTIVITY_PATTERNS_FILE

    df = pd.read_csv(patterns_file, sep='\t')

    patterns_dict = {}
    reviewed_dict = {}

    for substrate, group in df.groupby('substrate'):
        patterns_dict[substrate] = group['pattern'].tolist()
        reviewed_dict[substrate] = bool(group['reviewed'].all())

    return patterns_dict, reviewed_dict


def check_patterns_reviewed(substrate, reviewed_dict, output_dir):
    """
    Check if activity patterns for a substrate have been reviewed.
    Prints a warning and notes the review report location if not.

    Args:
        substrate:     substrate name string
        reviewed_dict: dict from load_activity_patterns()
        output_dir:    substrate output directory for reminder
    """
    if reviewed_dict.get(substrate, False):
        return

    print(
        f"\nWARNING: Activity patterns for '{substrate}' were "
        f"auto-derived and have not been manually reviewed.\n"
        f"         CGC filtering results should be verified against "
        f"known biology.\n"
        f"         Inspect {substrate}_activity_annotated.tsv, then\n"
        f"         edit substrate/data/activity_patterns.tsv and "
        f"set reviewed=True.\n"
        f"         Pattern review report: "
        f"{output_dir}/{substrate}_pattern_review.tsv\n"
    )



# Activity relevance patterns per substrate.
# A gene is considered substrate-relevant if its activity string contains
# any of the listed patterns (case-insensitive substring match).
#
# To add a new substrate: add a new key with appropriate patterns.
# Patterns should be broad enough to catch variant naming in EXPASY/dbCAN
# but specific enough to exclude unrelated activities.
#
# NOTE: These will be replaced by dynamic pattern derivation in a future
# update (see refactoring notes — dynamic activity patterns).
ACTIVITY_PATTERNS = {
    'laminarin': [
        'glucan',
        'glucosidase',
        'laminarin',
        '1,3',
    ],
    'glycogen': [
        'glucan',
        'glucos',
        'glucosyltransferase',
        'dextrin',
        'glucosidase',
        'amylase',
        'glucanotransferase',
        'glycogen',
        'starch',
        '1,4',
        '1,6',
    ],
    'pullulan': [
        'pullulan',
        'glucan',
        'glucosidase',
        'amylase',
        'dextran',
        'glucanotransferase',
        '1,4',
        '1,6',
    ],
}


# ── TCDB label loading ────────────────────────────────────────────────────────

def load_tcdb_labels(tcdb_file):
    """
    Load TC family definitions from TCDB download.

    Args:
        tcdb_file: path to tc_family_definitions.tsv

    Returns:
        dict mapping TC ID string to description string
    """
    labels = {}
    with open(tcdb_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                tc_id   = parts[0].strip()
                tc_name = parts[1].strip()
                labels[tc_id] = tc_name
    return labels


def label_tc(tc_id, tcdb_labels):
    """
    Convert a TCDB ID to a readable label, matching at family level
    if an exact match is not found.

    Args:
        tc_id:       TC identifier string
        tcdb_labels: dict from load_tcdb_labels()

    Returns:
        Human-readable TC family name, or 'TC:<tc_id>' if not found
    """
    if tc_id in tcdb_labels:
        return tcdb_labels[tc_id]
    parts = tc_id.split('.')
    for length in [3, 2]:
        family = '.'.join(parts[:length])
        if family in tcdb_labels:
            return tcdb_labels[family]
    return f'TC:{tc_id}'


# ── Activity relevance ────────────────────────────────────────────────────────

def is_relevant_activity(activity, substrate, patterns_dict):
    """
    Check if an activity string is relevant for the given substrate.

    Args:
        activity:      activity label string from activity_annotated.tsv
        substrate:     substrate name string
        patterns_dict: dict from load_activity_patterns()

    Returns:
        True if any inclusion pattern matches the activity string
    """
    if not activity or str(activity).strip() in ['unknown', 'nan', '-', '']:
        return False
    patterns     = patterns_dict.get(substrate, [])
    activity_low = str(activity).lower()
    return any(p.lower() in activity_low for p in patterns)


# ── Genome FASTA auto-detection ───────────────────────────────────────────────

def find_genome_fasta(sample, genomes_dir):
    """
    Find the genome FASTA file for a sample by exact name match.

    Strips the file extension from each file in genomes_dir and checks
    for an exact match against the sample name. This handles filenames
    like 'Strain.final.assembly.fasta' where the sample name is
    'Strain.final.assembly'.

    Args:
        sample:      sample name string (from dbCAN output directory name)
        genomes_dir: path to directory containing genome FASTA files

    Returns:
        Full path to the matched FASTA file, or None if not found
    """
    for fname in os.listdir(genomes_dir):
        for ext in NUCLEOTIDE_EXTENSIONS:
            if fname.endswith(ext):
                stem = fname[:-len(ext)]
                if stem == sample:
                    return os.path.join(genomes_dir, fname)

    print(f"  WARNING: No genome FASTA found for sample '{sample}' "
          f"in {genomes_dir}.")
    print(f"           Expected a file named '{sample}' with extension "
          f"{', '.join(sorted(NUCLEOTIDE_EXTENSIONS))}. "
          f"Skipping this sample.")
    return None


# ── CGC filtering ─────────────────────────────────────────────────────────────

def filter_qualifying_cgcs(hits_df, activity_file, substrate,
                           output_dir, min_cazymes=2,
                           patterns_file=None):
    """
    Filter CGCs to those with >= min_cazymes substrate-relevant activity
    annotations, using canonical_PUL localisation only.

    Args:
        hits_df:       family hits DataFrame filtered to current substrate
        activity_file: path to <substrate>_activity_annotated.tsv
        substrate:     substrate name
        output_dir:    substrate output directory (for review report)
        min_cazymes:   minimum relevant activity genes per CGC
        patterns_file: optional path to activity_patterns.tsv override

    Returns:
        set of (sample, cgc_id) tuples for qualifying CGCs
    """
    patterns_dict, reviewed_dict = load_activity_patterns(patterns_file)
    check_patterns_reviewed(substrate, reviewed_dict, output_dir)

    hits = hits_df[hits_df['localisation'] == 'canonical_PUL'].copy()
    # If activity already present in hits_df (added by annotate_hits())
    # use it directly to avoid merge conflicts
    if 'activity' in hits.columns:
        hits['activity'] = hits['activity'].fillna('unknown')
        print(f"Using existing activity annotations "
              f"({len(hits)} canonical PUL genes)")
    elif os.path.exists(activity_file):
        act_df = pd.read_csv(activity_file, sep='\t')
        act_df = act_df[act_df['sample'] != 'Reference'].copy()
        act_df = act_df[['Gene ID', 'activity']].drop_duplicates()
        hits   = hits.merge(act_df, on='Gene ID', how='left')
        hits['activity'] = hits['activity'].fillna('unknown')
        print(f"Merged activity annotations for {len(hits)} "
              f"canonical PUL genes")
    else:
        hits['activity'] = 'unknown'
        print(f"WARNING: No activity file found at {activity_file}")
        print("Run activity annotation first — all CGCs will be excluded")
    hits['activity_relevant'] = hits['activity'].apply(
        lambda a: is_relevant_activity(a, substrate, patterns_dict)
    )

    cgc_counts = (hits.groupby(['sample', 'CGC_id'])['activity_relevant']
                  .sum()
                  .reset_index())
    cgc_counts.columns = ['sample', 'CGC_id', 'relevant_gene_count']

    qualifying = cgc_counts[cgc_counts['relevant_gene_count'] >= min_cazymes]

    print(f"\nSubstrate: {substrate}")
    print(f"Total canonical PUL CGCs before activity filter: "
          f"{hits.groupby(['sample', 'CGC_id']).ngroups}")
    print(f"CGCs with >= {min_cazymes} {substrate}-relevant activity genes: "
          f"{len(qualifying)}")

    # Report excluded CGCs for transparency
    all_cgcs  = hits.groupby(['sample', 'CGC_id'])['activity_relevant'].sum()
    excluded  = all_cgcs[all_cgcs < min_cazymes]
    if len(excluded) > 0:
        print(f"\nExcluded CGCs (< {min_cazymes} relevant activity genes):")
        for (sample, cgc_id), count in excluded.items():
            activities = hits[
                (hits['sample'] == sample) &
                (hits['CGC_id'] == cgc_id)
            ]['activity'].value_counts()
            print(f"  {sample} {cgc_id} "
                  f"(relevant={int(count)}, "
                  f"activities: {dict(activities)})")

    return set(zip(qualifying['sample'], qualifying['CGC_id']))


# ── GenBank record construction ───────────────────────────────────────────────

def _make_gene_label(gene, tcdb_labels):
    """
    Derive a clean display label for a CGC gene feature.

    Args:
        gene:        row from cgc_standard_out.tsv
        tcdb_labels: dict from load_tcdb_labels()

    Returns:
        Label string
    """
    annot     = str(gene['Gene Annotation'])
    gene_type = str(gene['Gene Type'])

    if 'CAZyme' in gene_type:
        raw   = annot.split('|')[-1] if '|' in annot else annot
        return raw.split('+')[0].split('_e')[0]
    elif 'TC' in gene_type:
        tc_id = annot.split('|')[-1] if '|' in annot else annot
        return label_tc(tc_id, tcdb_labels)
    elif 'TF' in gene_type:
        return 'TF'
    elif 'STP' in gene_type:
        return 'STP'
    elif annot in ['nan', ''] or pd.isna(gene['Gene Annotation']):
        return 'hypothetical'
    else:
        return annot.split('|')[-1] if '|' in annot else annot


def build_genbank_record(sample, cgc_id, cgc_genes, genome,
                         sample_label, substrate, tcdb_labels):
    """
    Build a BioPython SeqRecord for a single CGC genomic region.

    Args:
        sample:       sample name string
        cgc_id:       CGC identifier string
        cgc_genes:    DataFrame rows for this CGC from cgc_standard_out.tsv
        genome:       dict mapping contig ID to SeqRecord (from SeqIO.to_dict)
        sample_label: display name for the organism
        substrate:    substrate name
        tcdb_labels:  dict from load_tcdb_labels()

    Returns:
        SeqRecord with CDS features, or None if contig not found
    """
    contig_id = cgc_genes['Contig ID'].iloc[0]
    cgc_start = int(cgc_genes['Gene Start'].min()) - 1  # 0-based
    cgc_end   = int(cgc_genes['Gene Stop'].max())

    if contig_id not in genome:
        print(f"  WARNING: Contig {contig_id} not found in genome")
        return None

    region_seq = genome[contig_id].seq[cgc_start:cgc_end]
    clean_cgc  = cgc_id.replace('/', '_').replace('|', '_')
    safe_cgc   = clean_cgc.replace(' ', '_')
    org_short  = sample_label.replace(' ', '_')[:10]

    record = SeqRecord(
        region_seq,
        id=f"{sample}__{clean_cgc}",
        name=f"{org_short}_{safe_cgc}"[:16],
        description=f"{sample_label} {cgc_id} {substrate} CGC",
        annotations={
            'molecule_type': 'DNA',
            'organism':      sample_label,
        }
    )

    for _, gene in cgc_genes.iterrows():
        start  = int(int(gene['Gene Start']) - 1 - cgc_start)
        end    = int(int(gene['Gene Stop']) - cgc_start)
        strand = 1 if gene['Gene Strand'] == '+' else -1
        label  = _make_gene_label(gene, tcdb_labels)

        feature = SeqFeature(
            FeatureLocation(start, end, strand=strand),
            type='CDS',
            qualifiers={
                'gene':      [label],
                'product':   [str(gene['Gene Annotation'])],
                'locus_tag': [str(gene['Protein ID'])],
                'note':      [str(gene['Gene Type'])],
            }
        )
        record.features.append(feature)

    return record


# ── Main entry point ──────────────────────────────────────────────────────────

def make_genbank_files(cgc_output_dir, hits_df, genomes_dir, output_dir,
                       substrate, tcdb_file, activity_file,
                       sample_labels=None, min_cazymes=2,
                       patterns_file=None):
    """
    Extract qualifying CGC regions and write GenBank files.

    Clears the genbank output directory before writing to prevent old
    substrate files persisting between runs.

    Args:
        cgc_output_dir: path to dbCAN cgc_output directory
        hits_df:        family hits DataFrame for current substrate
        genomes_dir:    path to directory containing genome FASTA files
        output_dir:     substrate output directory (<base>/<substrate>/)
        substrate:      substrate name
        tcdb_file:      path to tc_family_definitions.tsv
        activity_file:  path to <substrate>_activity_annotated.tsv
        sample_labels:  dict mapping sample name to display name (optional).
                        If None, sample names are used as-is.
        min_cazymes:    minimum relevant activity genes per CGC
        patterns_file:  optional path to activity_patterns.tsv override

    Returns:
        Total number of GenBank files written
    """
    gbk_dir = os.path.join(output_dir, 'genbank')

    # Always clear genbank directory before writing to prevent old files
    # from previous substrate runs being picked up by clinker
    if os.path.exists(gbk_dir):
        shutil.rmtree(gbk_dir)
        print(f"Cleared old genbank directory: {gbk_dir}")
    os.makedirs(gbk_dir)

    tcdb_labels = load_tcdb_labels(tcdb_file)

    qualifying_cgcs = filter_qualifying_cgcs(
        hits_df, activity_file, substrate,
        output_dir=output_dir,
        min_cazymes=min_cazymes,
        patterns_file=patterns_file
    )

    if not qualifying_cgcs:
        print("No qualifying CGCs found — no GenBank files written")
        return 0

    total_written = 0

    for sample_dir in sorted(os.listdir(cgc_output_dir)):
        full_path = os.path.join(cgc_output_dir, sample_dir)
        if not os.path.isdir(full_path):
            continue

        sample = sample_dir.replace('output_', '')

        sample_qualifying = [
            cgc_id for (s, cgc_id) in qualifying_cgcs
            if s == sample
        ]
        if not sample_qualifying:
            continue

        print(f"\n{sample}: {len(sample_qualifying)} qualifying CGCs — "
              f"{sorted(sample_qualifying)}")

        # Load CGC gene coordinates
        cgc_file = os.path.join(full_path, 'cgc_standard_out.tsv')
        if not os.path.exists(cgc_file):
            print(f"  WARNING: No cgc_standard_out.tsv for {sample}")
            continue

        cgc_df = pd.read_csv(cgc_file, sep='\t')
        cgc_df.columns = cgc_df.columns.str.strip()
        cgc_df['Gene Annotation'] = cgc_df['Gene Annotation'].fillna('')

        # Auto-detect genome FASTA
        fasta_path = find_genome_fasta(sample, genomes_dir)
        if not fasta_path:
            continue

        genome = SeqIO.to_dict(SeqIO.parse(fasta_path, 'fasta'))

        sample_label  = (sample_labels.get(sample, sample)
                         if sample_labels else sample)
        sample_gbk_dir = os.path.join(gbk_dir, sample)
        os.makedirs(sample_gbk_dir, exist_ok=True)

        for cgc_id in sorted(sample_qualifying):
            cgc_genes = cgc_df[cgc_df['CGC#'] == cgc_id]
            if cgc_genes.empty:
                print(f"  WARNING: No genes found for {cgc_id}")
                continue

            record = build_genbank_record(
                sample, cgc_id, cgc_genes, genome,
                sample_label, substrate, tcdb_labels
            )
            if record is None:
                continue

            safe_cgc = cgc_id.replace('/', '_').replace(
                '|', '_').replace(' ', '_')
            out_file = os.path.join(
                sample_gbk_dir, f"{sample}__{safe_cgc}.gbk")

            with open(out_file, 'w') as f:
                SeqIO.write(record, f, 'genbank')
            total_written += 1

    print(f"\nTotal GenBank files written: {total_written}")
    print(f"Output directory: {gbk_dir}")
    return total_written
