"""
Clinker synteny plot wrapper.

Provides two functions:
  1. generate_clinker_inputs() — scans GenBank files and writes
     gene_functions.csv and colour_map.csv for clinker
  2. run_clinker() — wraps the clinker subprocess call

Colour assignment rules:
  - CAZymes: coloured by top-level family (GH16_3 -> GH16)
  - Known accessory enzyme classes (sulfatases, deacetylases etc):
    single colour per class regardless of subfamily
  - Transporters: fixed colours per category (SusC, SusD etc.)
  - Hypothetical: grey
  - TF, STP: fixed colours

All fixed category colours are loaded from default_colours.tsv
so they can be edited without touching code.
"""
import os
import csv
import glob
import re
import subprocess
import pandas as pd
from Bio import SeqIO


# ── Exceptions ────────────────────────────────────────────────────────────────

class ToolNotFoundError(Exception):
    """Raised when clinker is not found on PATH."""
    pass


# ── Accessory enzyme label groups ─────────────────────────────────────────────

# Maps display group name to list of patterns (case-insensitive
# substring match against the gene label).
# Extend this dict as new accessory enzyme classes are encountered.
LABEL_GROUPS = {
    'sulfatase':   ['sulfatase'],
    'deacetylase': ['deacetylase'],
    'epimerase':   ['epimerase'],
    'isomerase':   ['isomerase'],
    'oxidase':     ['oxidase'],
    'peroxidase':  ['peroxidase'],
}


# ── Tool availability ─────────────────────────────────────────────────────────

def check_clinker():
    """
    Check that clinker is available and return its version string.

    Raises:
        ToolNotFoundError if clinker is not found on PATH
    """
    try:
        result = subprocess.run(
            ['clinker', '--version'],
            capture_output=True, text=True
        )
        version = result.stdout.strip() or result.stderr.strip()
        return version
    except FileNotFoundError:
        raise ToolNotFoundError(
            "clinker not found on PATH. "
            "Install with: conda install -c bioconda clinker-py"
        )


# ── Colour loading ────────────────────────────────────────────────────────────

def load_clinker_colours(colours_file):
    """
    Load clinker category colours from default_colours.tsv.

    Loads rows with palette values starting with 'clinker_' for
    fixed categories (transporters, TF, STP, hypothetical) and
    rows starting with 'accessory_' for accessory enzyme groups.

    Args:
        colours_file: path to default_colours.tsv

    Returns:
        tuple of (category_colours, accessory_colours) where each
        is a dict mapping label to hex colour string
    """
    df = pd.read_csv(colours_file, sep='\t')

    category_colours = {
        row['palette'][len('clinker_'):]: row['hex']
        for _, row in df.iterrows()
        if row['palette'].startswith('clinker_')
    }

    accessory_colours = {
        row['palette'][len('accessory_'):]: row['hex']
        for _, row in df.iterrows()
        if row['palette'].startswith('accessory_')
    }

    return category_colours, accessory_colours


# ── Gene label extraction ─────────────────────────────────────────────────────

def _top_level_family(family_str):
    """
    Extract top-level CAZyme family from a subfamily annotation.

    Examples:
        'GH16_3'     -> 'GH16'
        'GH16_3_e12' -> 'GH16'
        'CBM48'      -> 'CBM48'
        'PL7_5'      -> 'PL7'

    Args:
        family_str: raw family/subfamily string

    Returns:
        Top-level family string
    """
    match = re.match(r'([A-Za-z]+[0-9]+)', str(family_str))
    return match.group(1) if match else str(family_str)


def _accessory_group(label):
    """
    Check if a label matches any known accessory enzyme group.

    Args:
        label: gene label string

    Returns:
        Group name string if matched, or None if no match
    """
    label_lower = label.lower()
    for group, patterns in LABEL_GROUPS.items():
        if any(p in label_lower for p in patterns):
            return group
    return None


def get_gene_label(feature):
    """
    Derive a clean display label for a CGC gene feature.

    Priority:
    1. CAZymes — collapse to top-level family (GH16_3 -> GH16)
    2. Known accessory enzyme classes — collapse to group name
    3. Transporters — use TC label as-is
    4. TF, STP — use gene type as label
    5. Hypothetical/unknown — 'hypothetical'

    Args:
        feature: BioPython SeqFeature object

    Returns:
        Label string
    """
    gene      = feature.qualifiers.get('gene',    [''])[0]
    note      = feature.qualifiers.get('note',    [''])[0]
    product   = feature.qualifiers.get('product', [''])[0]

    if 'CAZyme' in note:
        # Collapse subfamily to top-level family
        raw = gene.split('_e')[0] if gene else product
        return _top_level_family(raw)

    if 'SULFATLAS' in note or 'sulfatase' in product.lower():
        return 'sulfatase'

    # Check accessory enzyme groups from product annotation
    if product and product not in ['', 'nan', 'hypothetical']:
        group = _accessory_group(product)
        if group:
            return group

    # Check gene label for accessory groups
    if gene and gene not in ['', 'nan', 'hypothetical']:
        group = _accessory_group(gene)
        if group:
            return group

    if 'TC' in note:
        tc_id = product.split('|')[-1] if '|' in product else product
        return tc_id if tc_id else 'transporter'

    if 'TF' in note:
        return 'TF'

    if 'STP' in note:
        return 'STP'

    return 'hypothetical'


# ── CSV generation ────────────────────────────────────────────────────────────

def generate_clinker_inputs(gbk_dir, output_dir, substrate,
                            colours_file):
    """
    Scan GenBank files and write gene_functions.csv and
    colour_map.csv for use with clinker.

    CAZyme families get colours from the activity palette.
    Known accessory enzyme classes get colours from the accessory
    palette. Fixed categories use clinker palette colours.
    All colours are loaded from default_colours.tsv.

    Args:
        gbk_dir:      path to genbank output directory
        output_dir:   substrate output directory
        substrate:    substrate name
        colours_file: path to default_colours.tsv

    Returns:
        tuple of (gene_functions_path, colour_map_path),
        or (None, None) if no CDS features found
    """
    category_colours, accessory_colours = load_clinker_colours(
        colours_file)

    # Load activity palette for CAZyme family colours
    df_colours  = pd.read_csv(colours_file, sep='\t')
    cazyme_palette = df_colours[
        df_colours['palette'] == 'activity']['hex'].tolist()

    clin_dir = os.path.join(output_dir, 'clinker')
    os.makedirs(clin_dir, exist_ok=True)

    # Collect all unique labels and locus->label mappings
    all_labels     = set()
    locus_to_label = {}

    for root, dirs, files in os.walk(gbk_dir):
        for fname in sorted(files):
            if not fname.endswith('.gbk'):
                continue
            gbk_path = os.path.join(root, fname)
            for record in SeqIO.parse(gbk_path, 'genbank'):
                for feature in record.features:
                    if feature.type != 'CDS':
                        continue
                    locus = feature.qualifiers.get(
                        'locus_tag', [''])[0]
                    label = get_gene_label(feature)
                    locus_to_label[locus] = label
                    all_labels.add(label)

    if not locus_to_label:
        print(f"WARNING: No CDS features found in {gbk_dir}")
        return None, None

    # Assign colours
    # Fixed categories first
    label_colours = {}
    for label, colour in category_colours.items():
        label_colours[label] = colour
    for group, colour in accessory_colours.items():
        label_colours[group] = colour

    # Dynamic colours for CAZyme top-level families
    cazyme_labels = sorted([
        l for l in all_labels
        if l not in label_colours
        and re.match(r'^(GH|PL|CE|AA|CBM)[0-9]+$', l)
    ])
    for i, label in enumerate(cazyme_labels):
        label_colours[label] = cazyme_palette[
            i % len(cazyme_palette)]

    # Any remaining labels get grey
    for label in all_labels:
        if label not in label_colours:
            label_colours[label] = '#cccccc'

    # Write gene_functions.csv
    gf_file = os.path.join(
        clin_dir, f'{substrate}_gene_functions.csv')
    with open(gf_file, 'w', newline='') as f:
        writer = csv.writer(f)
        for locus, label in sorted(locus_to_label.items()):
            writer.writerow([locus, label])

    # Write colour_map.csv
    cm_file = os.path.join(
        clin_dir, f'{substrate}_colour_map.csv')
    with open(cm_file, 'w', newline='') as f:
        writer = csv.writer(f)
        for label, colour in sorted(label_colours.items()):
            writer.writerow([label, colour])

    print(f"Gene functions: {gf_file} ({len(locus_to_label)} genes)")
    print(f"Colour map:     {cm_file} ({len(label_colours)} categories)")
    print("\nColour assignments:")
    for label, colour in sorted(label_colours.items()):
        count = sum(1 for l in locus_to_label.values()
                    if l == label)
        if count > 0:
            print(f"  {colour}  {label} ({count} genes)")

    return gf_file, cm_file


# ── Clinker subprocess wrapper ────────────────────────────────────────────────

def run_clinker(gbk_dir, output_dir, substrate, gene_functions_file,
                colour_map_file, identity=0.3, jobs=8, log_path=None):
    """
    Run clinker on all GenBank files for the current substrate.

    Args:
        gbk_dir:              path to genbank output directory
        output_dir:           substrate output directory
        substrate:            substrate name
        gene_functions_file:  path to gene_functions.csv
        colour_map_file:      path to colour_map.csv
        identity:             minimum identity threshold (default: 0.3)
        jobs:                 number of parallel jobs (default: 8)
        log_path:             path to append clinker log (optional)

    Returns:
        path to output HTML file

    Raises:
        ToolNotFoundError if clinker is not on PATH
        FileNotFoundError if no GenBank files are found
        subprocess.CalledProcessError if clinker exits non-zero
    """
    check_clinker()

    gbk_files = sorted(glob.glob(
        os.path.join(gbk_dir, '**', '*.gbk'), recursive=True))
    if not gbk_files:
        raise FileNotFoundError(
            f"No GenBank files found in {gbk_dir}. "
            "Run genbank step first."
        )

    clin_dir = os.path.join(output_dir, 'clinker')
    os.makedirs(clin_dir, exist_ok=True)

    html_file = os.path.join(
        clin_dir, f'{substrate}_all_cgcs.html')
    tsv_file  = os.path.join(
        clin_dir, f'{substrate}_clinker.tsv')

    cmd = [
        'clinker',
        *gbk_files,
        '--plot',           html_file,
        '--output',         tsv_file,
        '--identity',       str(identity),
        '--jobs',           str(jobs),
        '--gene_functions', gene_functions_file,
        '--colour_map',     colour_map_file,
        '--force',
    ]

    print(f"Running clinker on {len(gbk_files)} GenBank files...")

    if log_path:
        with open(log_path, 'a') as log:
            subprocess.run(cmd, stderr=log, check=True)
    else:
        subprocess.run(cmd, stderr=subprocess.DEVNULL, check=True)

    print(f"Done. Output: {html_file}")
    return html_file
