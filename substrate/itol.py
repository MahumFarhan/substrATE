"""
Generate iTOL annotation files for all family trees.

Produces colour strip files for sample, localisation and enzymatic
activity. Supports a two-phase workflow:

  Phase 1 — assign colours, write <substrate>_colour_config.tsv
  Phase 2 — read colour_config.tsv, write iTOL annotation files

This separation allows users to edit colours between phases and
regenerate annotation files without rebuilding the tree.

Sample labels are auto-generated from sample names by default.
A metadata TSV can be provided to override with custom labels.
"""
import os
import re
import pandas as pd
from Bio import SeqIO


# ── Constants ─────────────────────────────────────────────────────────────────

# Fixed localisation colours — tied to biological categories used in code
# Steel blue / amber / dark red / dark grey
LOCALISATION_COLOURS = {
    'canonical_PUL':           '#1a6faf',
    'non_canonical_CGC':       '#e08214',
    'outside_CGC':             '#b2182b',
    'characterised_reference': '#555555',
}

# Suffixes to strip when auto-generating sample labels from directory names
STRIP_SUFFIXES = [
    '.final.assembly',
    '_wholegenome',
    '_assembly',
    '_genome',
]


# ── Sample label generation ───────────────────────────────────────────────────

def _strip_sample_suffix(name):
    """Strip common assembly suffixes from a sample name."""
    for suffix in STRIP_SUFFIXES:
        if name.endswith(suffix):
            name = name[:-len(suffix)]
    return name


def generate_sample_labels(samples, metadata_file=None, substrate=None):
    """
    Generate display labels for samples.

    If metadata_file is provided and contains 'sample' and 'label'
    columns, those labels are used. Any samples not found in the
    metadata fall back to auto-generation.

    Auto-generation replaces underscores with spaces and strips
    common assembly suffixes.

    Args:
        samples:       iterable of sample name strings
        metadata_file: optional path to metadata TSV with sample/label columns
        substrate:     optional substrate name for filtering metadata

    Returns:
        dict mapping sample name to display label
    """
    labels = {}

    # Load metadata labels if provided
    meta_labels = {}
    if metadata_file and os.path.exists(metadata_file):
        meta_df = pd.read_csv(metadata_file, sep='\t')
        if 'sample' in meta_df.columns and 'label' in meta_df.columns:
            if substrate and 'substrate' in meta_df.columns:
                meta_df = meta_df[meta_df['substrate'] == substrate]
            meta_labels = dict(zip(
                meta_df['sample'].astype(str),
                meta_df['label'].astype(str)
            ))

    for sample in samples:
        if sample in meta_labels:
            labels[sample] = meta_labels[sample]
        else:
            # Auto-generate: strip suffixes, replace underscores with spaces
            clean = _strip_sample_suffix(sample)
            labels[sample] = clean.replace('_', ' ')

    return labels


# ── Colour palette loading ────────────────────────────────────────────────────

def load_colour_palettes(colours_file):
    """
    Load sample and activity colour palettes from default_colours.tsv.

    Args:
        colours_file: path to default_colours.tsv

    Returns:
        tuple of (sample_palette, activity_palette) where each is a
        list of hex colour strings
    """
    df = pd.read_csv(colours_file, sep='\t')
    sample_palette   = df[df['palette'] == 'sample']['hex'].tolist()
    activity_palette = df[df['palette'] == 'activity']['hex'].tolist()
    return sample_palette, activity_palette


# ── Colour assignment ─────────────────────────────────────────────────────────

def _generate_hsl_palette(n):
    """
    Generate N visually distinct colours using evenly spaced HSL hues.

    Uses fixed saturation (0.65) and lightness (0.50) for colours that
    are both distinct and visible against white backgrounds in iTOL.

    Args:
        n: number of colours to generate

    Returns:
        list of hex colour strings
    """
    import colorsys
    colours = []
    for i in range(n):
        hue        = i / n
        r, g, b    = colorsys.hls_to_rgb(hue, 0.50, 0.65)
        colours.append('#{:02x}{:02x}{:02x}'.format(
            int(r * 255), int(g * 255), int(b * 255)))
    return colours


def assign_sample_colours(samples, sample_palette, max_colours=None):
    """
    Assign colours to samples from the sample palette.

    Samples are sorted alphabetically for reproducibility.
    'Reference' is always assigned the last colour in the palette
    (dark grey by convention).

    If the number of samples exceeds the palette size:
      - If max_colours is set: switch to programmatic HSL colour
        generation for all samples, producing max_colours evenly
        spaced hues. Use this for datasets with >19 samples.
      - If max_colours is not set: colours cycle and a warning is
        printed suggesting the user set --max_colours.

    Args:
        samples:        iterable of sample name strings
        sample_palette: list of hex colour strings
        max_colours:    optional int — if set and sample count exceeds
                        palette size, generate this many HSL colours
                        instead of cycling the palette

    Returns:
        dict mapping sample name to hex colour string
    """
    sorted_samples = sorted(s for s in samples if s != 'Reference')
    n_samples      = len(sorted_samples)

    # Reserve last palette colour for Reference
    available_palette = sample_palette[:-1]
    palette_size      = len(available_palette)

    if n_samples > palette_size:
        if max_colours is not None:
            n_colours   = max(max_colours, n_samples)
            hsl_palette = _generate_hsl_palette(n_colours)
            print(f"INFO: {n_samples} samples exceed palette size "
                  f"({palette_size}) — using HSL colour generation "
                  f"({n_colours} colours)")
            colours = {}
            for i, sample in enumerate(sorted_samples):
                colours[sample] = hsl_palette[i]
        else:
            print(f"WARNING: {n_samples} samples exceed the colour "
                  f"palette size ({palette_size}).")
            print(f"         Colours will repeat. Use --max_colours "
                  f"{n_samples} to generate distinct colours for all "
                  f"samples automatically.")
            colours = {}
            for i, sample in enumerate(sorted_samples):
                colours[sample] = available_palette[i % palette_size]
    else:
        colours = {}
        for i, sample in enumerate(sorted_samples):
            colours[sample] = available_palette[i]

    if 'Reference' in samples:
        colours['Reference'] = sample_palette[-1]

    return colours


def assign_activity_colours(activities, activity_palette):
    """
    Assign colours to activity labels from the activity palette.

    Activities are sorted alphabetically for reproducibility across
    substrates.

    Args:
        activities:       iterable of activity label strings
        activity_palette: list of hex colour strings

    Returns:
        dict mapping activity label to hex colour string
    """
    # Fixed muted palette for reference substrates — kept separate from
    # genomic activity colours to avoid palette collisions
    ref_palette = [
        '#767676', '#b8860b', '#6b8e23', '#8b6914', '#4682b4',
        '#8b4513', '#708090', '#9400d3', '#2e8b57', '#cd853f',
        '#4a708b', '#8b2252', '#6b8e6b', '#b8733a', '#5f9ea0',
        '#8b7355', '#7b6b8b', '#4f7942', '#8b3a62', '#6b7b5b',
    ]
    genomic = sorted(a for a in set(activities) if not a.startswith('reference: '))
    refs    = sorted(a for a in set(activities) if a.startswith('reference: '))
    result  = {act: activity_palette[i % len(activity_palette)]
               for i, act in enumerate(genomic)}
    result.update({act: ref_palette[i % len(ref_palette)]
                   for i, act in enumerate(refs)})
    return result


# ── Colour config TSV ─────────────────────────────────────────────────────────

def write_colour_config(colour_config_path, sample_colours,
                        activity_colours):
    """
    Write a colour configuration TSV that users can edit between runs.

    Columns: category, label, hex

    Args:
        colour_config_path: path to write the TSV
        sample_colours:     dict mapping sample name to hex colour
        activity_colours:   dict mapping activity label to hex colour
    """
    rows = []
    for label, hex_col in sorted(sample_colours.items()):
        rows.append({'category': 'sample', 'label': label, 'hex': hex_col})
    for label, hex_col in sorted(activity_colours.items()):
        rows.append({'category': 'activity', 'label': label, 'hex': hex_col})
    for label, hex_col in sorted(LOCALISATION_COLOURS.items()):
        rows.append({'category': 'localisation', 'label': label,
                     'hex': hex_col})

    df = pd.DataFrame(rows)
    df.to_csv(colour_config_path, sep='\t', index=False)
    print(f"Colour config written to {colour_config_path}")
    print("Edit this file to customise colours, then rerun visualise step.")


def load_colour_config(colour_config_path):
    """
    Load a colour configuration TSV, returning dicts for each category.

    Args:
        colour_config_path: path to colour config TSV

    Returns:
        tuple of (sample_colours, activity_colours, localisation_colours)
        where each is a dict mapping label to hex colour string
    """
    df = pd.read_csv(colour_config_path, sep='\t')

    sample_colours = dict(zip(
        df[df['category'] == 'sample']['label'],
        df[df['category'] == 'sample']['hex']
    ))
    activity_colours = dict(zip(
        df[df['category'] == 'activity']['label'],
        df[df['category'] == 'activity']['hex']
    ))
    localisation_colours = dict(zip(
        df[df['category'] == 'localisation']['label'],
        df[df['category'] == 'localisation']['hex']
    ))
    return sample_colours, activity_colours, localisation_colours


# ── iTOL file writers ─────────────────────────────────────────────────────────

def write_colour_strip(out_file, leaf_data, colour_map, title,
                       legend_title):
    """Write an iTOL DATASET_COLORSTRIP annotation file."""
    # Filter legend to only values present in this family's data
    present_values = set(leaf_data.values())
    filtered_map = {k: v for k, v in colour_map.items()
                    if k in present_values}
    if not filtered_map:
        filtered_map = colour_map
    with open(out_file, 'w') as f:
        f.write('DATASET_COLORSTRIP\n')
        f.write('SEPARATOR TAB\n')
        f.write(f'DATASET_LABEL\t{title}\n')
        f.write('COLOR\t#000000\n')
        f.write(f'LEGEND_TITLE\t{legend_title}\n')
        f.write('LEGEND_SHAPES\t' +
                '\t'.join(['1'] * len(filtered_map)) + '\n')
        f.write('LEGEND_COLORS\t' +
                '\t'.join(filtered_map.values()) + '\n')
        f.write('LEGEND_LABELS\t' +
                '\t'.join(filtered_map.keys()) + '\n')
        f.write('DATA\n')
        for leaf_id, value in leaf_data.items():
            colour = colour_map.get(value, '#cccccc')
            f.write(f'{leaf_id}\t{colour}\t{value}\n')


def write_branch_colours(out_file, leaf_data, colour_map, title):
    """Write an iTOL TREE_COLORS file to colour branches by localisation."""
    with open(out_file, 'w') as f:
        f.write('TREE_COLORS\n')
        f.write('SEPARATOR TAB\n')
        f.write(f'DATASET_LABEL\t{title}\n')
        f.write('DATA\n')
        for leaf_id, value in leaf_data.items():
            colour = colour_map.get(value, '#cccccc')
            f.write(f'{leaf_id}\tbranch\t{colour}\tnormal\t2\n')


def write_leaf_symbols(out_file, leaf_data, colour_map, title):
    """Write an iTOL DATASET_SYMBOL file for leaf tip shapes."""
    with open(out_file, 'w') as f:
        f.write('DATASET_SYMBOL\n')
        f.write('SEPARATOR TAB\n')
        f.write(f'DATASET_LABEL\t{title}\n')
        f.write('COLOR\t#000000\n')
        f.write('DATA\n')
        for leaf_id, value in leaf_data.items():
            colour = colour_map.get(value, '#cccccc')
            f.write(f'{leaf_id}\t1\t10\t{colour}\t1\t1\n')


def write_labels(out_file, leaf_data):
    """Write an iTOL LABELS annotation file with clean display names."""
    with open(out_file, 'w') as f:
        f.write('LABELS\n')
        f.write('SEPARATOR TAB\n')
        f.write('DATA\n')
        for leaf_id, label in leaf_data.items():
            f.write(f'{leaf_id}\t{label}\n')


def write_label_styles(out_file, leaf_data):
    """Write an iTOL DATASET_STYLE file to bold genomic sequences."""
    with open(out_file, 'w') as f:
        f.write('DATASET_STYLE\n')
        f.write('SEPARATOR TAB\n')
        f.write('DATASET_LABEL\tLabel styles\n')
        f.write('COLOR\t#000000\n')
        f.write('DATA\n')
        for leaf_id in leaf_data:
            is_ref = leaf_id.startswith('Reference__')
            if not is_ref:
                f.write(f'{leaf_id}\tlabel\tlabel\t#000000\t1\tbold\n')


# ── Reference label loading ───────────────────────────────────────────────────

def load_ref_labels(ref_metadata, families):
    """
    Build accession -> label lookup from reference_metadata.tsv.

    Args:
        ref_metadata: path to reference_metadata.tsv
        families:     collection of family names to include

    Returns:
        dict mapping cleaned accession string to label string
    """
    if not os.path.exists(ref_metadata):
        return {}
    ref_meta = pd.read_csv(ref_metadata, sep='\t')
    ref_meta = ref_meta[ref_meta['family'].astype(str).isin(families)].copy()
    return {
        str(row['accession']).replace('|', '_').replace('=', '_'):
        str(row['label'])
        for _, row in ref_meta.iterrows()
    }


def load_ref_substrate_map(ref_metadata, families):
    """Build accession -> cleaned protein_name lookup from reference_metadata.tsv.
    Uses protein_name rather than substrate because the substrate field is
    assigned at family level and is often inaccurate for multi-functional
    families (e.g. GH16 contains agarases, carrageenases, lichenases etc.
    all labelled 'agar'). Protein name is stripped of gene code parentheticals
    to give a clean activity label, e.g. 'κ-carrageenase (Ce384; CeCgkA)'
    becomes 'κ-carrageenase'.
    """
    if not os.path.exists(ref_metadata):
        return {}
    ref_meta = pd.read_csv(ref_metadata, sep='\t')
    ref_meta = ref_meta[ref_meta['family'].astype(str).isin(families)].copy()
    # EC number -> canonical activity label for GH16 and other common cases
    # where protein_name variants are too numerous to normalise by text alone
    EC_LABELS = {
        '1.1.3.10': 'pyranose oxidase',
        '1.1.3.13': 'alcohol oxidase',
        '1.1.3.16': 'ecdysone oxidase',
        '1.1.3.4': 'glucose oxidase',
        '1.1.3.7': 'aryl-alcohol oxidase',
        '1.1.5.9': 'glucose 1-dehydrogenase (FAD, quinone)',
        '1.1.99.18': 'cellobiose dehydrogenase (acceptor)',
        '1.1.99.29': 'pyranose dehydrogenase (acceptor)',
        '1.14.99.53': 'lytic chitin monooxygenase',
        '1.14.99.54': 'lytic cellulose monooxygenase (C1-hydroxylating)',
        '1.14.99.55': 'lytic starch monooxygenase',
        '1.14.99.56': 'lytic cellulose monooxygenase (C4-dehydrogenating)',
        '2.3.1.122': 'trehalose O-mycolyltransferase',
        '2.3.1.20': 'diacylglycerol O-acyltransferase',
        '2.4.1.10': 'levansucrase',
        '2.4.1.100': '2,1-fructan:2,1-fructan 1-fructosyltransferase',
        '2.4.1.140': 'alternansucrase',
        '2.4.1.161': 'oligosaccharide 4-alpha-D-glucosyltransferase',
        '2.4.1.18': '1,4-alpha-glucan branching enzyme',
        '2.4.1.19': 'cyclomaltodextrin glucanotransferase',
        '2.4.1.2': 'dextrin dextranase',
        '2.4.1.20': 'cellobiose phosphorylase',
        '2.4.1.207': 'xyloglucan:xyloglucosyl transferase',
        '2.4.1.24': '1,4-alpha-glucan 6-alpha-glucosyltransferase',
        '2.4.1.243': '6(G)-fructosyltransferase',
        '2.4.1.25': '4-alpha-glucanotransferase',
        '2.4.1.280': "N,N'-diacetylchitobiose phosphorylase",
        '2.4.1.281': '4-O-beta-D-mannosyl-D-glucose phosphorylase',
        '2.4.1.31': 'laminaribiose phosphorylase',
        '2.4.1.319': 'beta-1,4-mannooligosaccharide phosphorylase',
        '2.4.1.320': '1,4-beta-mannosyl-N-acetylglucosamine phosphorylase',
        '2.4.1.321': 'cellobionic acid phosphorylase',
        '2.4.1.329': 'sucrose 6(F)-phosphate phosphorylase',
        '2.4.1.333': '1,2-beta-oligoglucan phosphorylase',
        '2.4.1.339': 'beta-1,2-mannobiose phosphorylase',
        '2.4.1.340': '1,2-beta-oligomannan phosphorylase',
        '2.4.1.352': 'glucosylglycerate phosphorylase',
        '2.4.1.359': 'glucosylglycerol phosphorylase (configuration-retaining)',
        '2.4.1.387': 'isomaltosyltransferase',
        '2.4.1.389': 'solabiose phosphorylase',
        '2.4.1.391': 'beta-1,2-glucosyltransferase',
        '2.4.1.392': '3-O-beta-D-glucopyranosyl-beta-D-glucuronide phosphorylase',
        '2.4.1.4': 'amylosucrase',
        '2.4.1.49': 'cellodextrin phosphorylase',
        '2.4.1.5': 'dextransucrase',
        '2.4.1.67': 'galactinol--raffinose galactosyltransferase',
        '2.4.1.7': 'sucrose phosphorylase',
        '2.4.1.82': 'galactinol--sucrose galactosyltransferase',
        '2.4.1.9': 'inulosucrase',
        '2.4.1.99': 'sucrose:sucrose fructosyltransferase',
        '2.4.99.16': 'starch synthase (maltosyl-transferring)',
        '3.1.1.11': 'pectinesterase',
        '3.1.1.3': 'triacylglycerol lipase',
        '3.1.1.41': 'cephalosporin-C deacetylase',
        '3.1.1.6': 'acetylesterase',
        '3.1.1.72': 'acetylxylan esterase',
        '3.1.1.73': 'feruloyl esterase',
        '3.1.1.74': 'cutinase',
        '3.1.1.86': 'rhamnogalacturonan acetylesterase',
        '3.2.1.1': 'alpha-amylase',
        '3.2.1.10': 'oligo-1,6-glucosidase',
        '3.2.1.100': 'mannan 1,4-mannobiosidase',
        '3.2.1.101': 'mannan endo-1,6-alpha-mannosidase',
        '3.2.1.102': 'blood-group-substance endo-1,4-beta-galactosidase',
        '3.2.1.103': 'keratan-sulfate endo-1,4-beta-galactosidase',
        '3.2.1.104': 'steryl-beta-glucosidase',
        '3.2.1.108': 'lactase',
        '3.2.1.11': 'dextranase',
        '3.2.1.116': 'glucan 1,4-alpha-maltotriohydrolase',
        '3.2.1.120': 'oligoxyloglucan beta-glycosidase',
        '3.2.1.122': "maltose-6'-phosphate glucosidase",
        '3.2.1.123': 'endoglycosylceramidase',
        '3.2.1.124': '3-deoxy-2-octulosonidase',
        '3.2.1.126': 'coniferin beta-glucosidase',
        '3.2.1.128': 'glycyrrhizinate beta-glucuronidase',
        '3.2.1.131': 'xylan alpha-1,2-glucuronosidase',
        '3.2.1.132': 'chitosanase',
        '3.2.1.133': 'glucan 1,4-alpha-maltohydrolase',
        '3.2.1.135': 'neopullulanase',
        '3.2.1.136': 'glucuronoarabinoxylan endo-1,4-beta-xylanase',
        '3.2.1.139': 'alpha-glucuronidase',
        '3.2.1.14': 'chitinase',
        '3.2.1.140': 'lacto-N-biosidase',
        '3.2.1.141': '4-alpha-D-{(1->4)-alpha-D-glucano}trehalose trehalohydrolase',
        '3.2.1.145': 'galactan 1,3-beta-galactosidase',
        '3.2.1.146': 'beta-galactofuranosidase',
        '3.2.1.149': 'beta-primeverosidase',
        '3.2.1.15': 'endo-polygalacturonase',
        '3.2.1.150': 'oligoxyloglucan reducing-end-specific cellobiohydrolase',
        '3.2.1.151': 'xyloglucan-specific endo-beta-1,4-glucanase',
        '3.2.1.152': 'mannosylglycoprotein endo-beta-mannosidase',
        '3.2.1.153': 'fructan beta-(2,1)-fructosidase',
        '3.2.1.154': 'fructan beta-(2,6)-fructosidase',
        '3.2.1.156': 'oligosaccharide reducing-end xylanase',
        '3.2.1.157': 'iota-carrageenase',
        '3.2.1.158': 'alpha-agarase',
        '3.2.1.159': 'alpha-neoagaro-oligosaccharide hydrolase',
        '3.2.1.162': 'lambda-carrageenase',
        '3.2.1.163': '1,6-alpha-D-mannosidase',
        '3.2.1.164': 'galactan endo-1,6-beta-galactosidase',
        '3.2.1.165': 'exo-1,4-beta-D-glucosaminidase',
        '3.2.1.166': 'heparanase',
        '3.2.1.167': 'baicalin-beta-D-glucuronidase',
        '3.2.1.168': 'hesperidin 6-O-alpha-L-rhamnosyl-beta-D-glucosidase',
        '3.2.1.169': 'protein O-GlcNAcase',
        '3.2.1.17': 'lysozyme',
        '3.2.1.171': 'rhamnogalacturonan hydrolase',
        '3.2.1.172': 'unsaturated rhamnogalacturonyl hydrolase',
        '3.2.1.173': 'rhamnogalacturonan galacturonohydrolase',
        '3.2.1.174': 'rhamnogalacturonan rhamnohydrolase',
        '3.2.1.175': 'beta-D-glucopyranosyl abscisate beta-glucosidase',
        '3.2.1.176': 'cellulose 1,4-beta-cellobiosidase (reducing end)',
        '3.2.1.177': 'alpha-D-xyloside xylohydrolase',
        '3.2.1.178': 'beta-porphyranase',
        '3.2.1.179': 'gellan tetrasaccharide unsaturated glucuronosyl hydrolase',
        '3.2.1.18': 'exo-alpha-sialidase',
        '3.2.1.180': 'unsaturated chondroitin disaccharide hydrolase',
        '3.2.1.181': 'galactan endo-beta-1,3-galactanase',
        '3.2.1.185': 'non-reducing end beta-L-arabinofuranosidase',
        '3.2.1.186': 'protodioscin 26-O-beta-D-glucosidase',
        '3.2.1.197': 'beta-1,2-mannosidase',
        '3.2.1.199': 'sulfoquinovosidase',
        '3.2.1.2': 'beta-amylase',
        '3.2.1.20': 'alpha-glucosidase',
        '3.2.1.200': 'exo-chitinase (non-reducing end)',
        '3.2.1.201': 'exo-chitinase (reducing end)',
        '3.2.1.204': '1,3-alpha-isomaltosidase',
        '3.2.1.205': 'isomaltose glucohydrolase',
        '3.2.1.207': 'mannosyl-oligosaccharide alpha-1,3-glucosidase',
        '3.2.1.21': 'beta-glucosidase',
        '3.2.1.211': 'endo-(1->3)-fucoidanase',
        '3.2.1.212': 'endo-(1->4)-fucoidanase',
        '3.2.1.213': 'galactan exo-1,6-beta-galactobiohydrolase (non-reducing end)',
        '3.2.1.215': 'arabinofuranosidase (non-reducing end)',
        '3.2.1.217': 'exo-acting protein-alpha-N-acetylgalactosaminidase',
        '3.2.1.22': 'alpha-galactosidase',
        '3.2.1.222': 'funoran endo-beta-hydrolase',
        '3.2.1.223': 'arabinofuranosidase (non-reducing end)',
        '3.2.1.228': 'funoran endo-alpha-hydrolase',
        '3.2.1.23': 'beta-galactosidase',
        '3.2.1.24': 'alpha-mannosidase',
        '3.2.1.25': 'beta-mannosidase',
        '3.2.1.26': 'beta-fructofuranosidase',
        '3.2.1.28': 'alpha,alpha-trehalase',
        '3.2.1.3': 'glucan 1,4-alpha-glucosidase',
        '3.2.1.31': 'beta-glucuronidase',
        '3.2.1.32': 'endo-1,3-beta-xylanase',
        '3.2.1.35': 'hyaluronoglucosaminidase',
        '3.2.1.36': 'hyaluronoglucuronidase',
        '3.2.1.37': 'xylan 1,4-beta-xylosidase',
        '3.2.1.38': 'beta-D-fucosidase',
        '3.2.1.39': 'glucan endo-1,3-beta-D-glucosidase',
        '3.2.1.4': 'cellulase',
        '3.2.1.40': 'alpha-L-rhamnosidase',
        '3.2.1.41': 'pullulanase',
        '3.2.1.45': 'glucosylceramidase',
        '3.2.1.46': 'galactosylceramidase',
        '3.2.1.48': 'sucrose alpha-glucosidase',
        '3.2.1.49': 'alpha-N-acetylgalactosaminidase',
        '3.2.1.51': 'alpha-L-fucosidase',
        '3.2.1.52': 'beta-N-acetylhexosaminidase',
        '3.2.1.53': 'beta-N-acetylgalactosaminidase',
        '3.2.1.54': 'cyclomaltodextrinase',
        '3.2.1.55': 'non-reducing end alpha-L-arabinofuranosidase',
        '3.2.1.57': 'isopullulanase',
        '3.2.1.58': 'glucan 1,3-beta-glucosidase',
        '3.2.1.6': 'endo-1,3(4)-beta-glucanase',
        '3.2.1.60': 'glucan 1,4-alpha-maltotetraohydrolase',
        '3.2.1.63': '1,2-alpha-L-fucosidase',
        '3.2.1.64': '2,6-beta-fructan 6-levanbiohydrolase',
        '3.2.1.65': 'levanase',
        '3.2.1.67': 'galacturonan 1,4-alpha-galacturonidase',
        '3.2.1.68': 'isoamylase',
        '3.2.1.7': 'inulinase',
        '3.2.1.70': 'glucan 1,6-alpha-glucosidase',
        '3.2.1.72': 'xylan 1,3-beta-xylosidase',
        '3.2.1.73': 'licheninase',
        '3.2.1.74': 'glucan 1,4-beta-glucosidase',
        '3.2.1.75': 'glucan endo-1,6-beta-glucosidase',
        '3.2.1.76': 'L-iduronidase',
        '3.2.1.78': 'mannan endo-1,4-beta-mannosidase',
        '3.2.1.8': 'endo-1,4-beta-xylanase',
        '3.2.1.81': 'beta-agarase',
        '3.2.1.82': 'exo-poly-alpha-digalacturonosidase',
        '3.2.1.83': 'kappa-carrageenase',
        '3.2.1.84': 'glucan 1,3-alpha-glucosidase',
        '3.2.1.86': '6-phospho-beta-glucosidase',
        '3.2.1.88': 'non-reducing end beta-L-arabinopyranosidase',
        '3.2.1.89': 'arabinogalactan endo-beta-1,4-galactanase',
        '3.2.1.91': 'cellulose 1,4-beta-cellobiosidase (non-reducing end)',
        '3.2.1.93': 'alpha,alpha-phosphotrehalase',
        '3.2.1.94': 'glucan 1,6-alpha-isomaltosidase',
        '3.2.1.95': 'dextran 1,6-alpha-isomaltotriosidase',
        '3.2.1.96': 'mannosyl-glycoprotein endo-beta-N-acetylglucosaminidase',
        '3.2.1.97': 'endo-alpha-N-acetylgalactosaminidase',
        '3.2.1.98': 'glucan 1,4-alpha-maltohexaosidase',
        '3.2.1.99': 'arabinan endo-1,5-alpha-L-arabinosidase',
        '3.5.1.104': 'peptidoglycan-N-acetylglucosamine deacetylase',
        '3.5.1.105': 'chitin disaccharide deacetylase',
        '3.5.1.115': 'mycothiol S-conjugate amidase',
        '3.5.1.136': "N,N'-diacetylchitobiose non-reducing end deacetylase",
        '3.5.1.33': 'N-acetylglucosamine deacetylase',
        '3.5.1.41': 'chitin deacetylase',
        '3.5.1.89': 'N-acetylglucosaminylphosphatidylinositol deacetylase',
        '4.2.2.1': 'hyaluronate lyase',
        '4.2.2.10': 'pectin lyase',
        '4.2.2.11': 'guluronate-specific alginate lyase',
        '4.2.2.12': 'xanthan lyase',
        '4.2.2.13': 'exo-(1->4)-alpha-D-glucan lyase',
        '4.2.2.14': 'glucuronan lyase',
        '4.2.2.15': 'anhydrosialidase',
        '4.2.2.16': 'levan fructotransferase (DFA-IV-forming)',
        '4.2.2.17': 'inulin fructotransferase (DFA-I-forming)',
        '4.2.2.18': 'inulin fructotransferase (DFA-III-forming)',
        '4.2.2.19': 'chondroitin B lyase',
        '4.2.2.2': 'pectate lyase',
        '4.2.2.23': 'rhamnogalacturonan endolyase',
        '4.2.2.24': 'rhamnogalacturonan exolyase',
        '4.2.2.25': 'gellan lyase',
        '4.2.2.26': 'oligo-alginate lyase',
        '4.2.2.28': 'alpha-L-rhamnosyl-(1->4)-beta-D-glucuronate lyase',
        '4.2.2.29': 'peptidoglycan lytic transglycosylase',
        '4.2.2.3': 'mannuronate-specific alginate lyase',
        '4.2.2.5': 'chondroitin AC lyase',
        '4.2.2.6': 'oligogalacturonide lyase',
        '4.2.2.7': 'heparin lyase',
        '4.2.2.8': 'heparin-sulfate lyase',
        '4.2.2.9': 'pectate disaccharide-lyase',
        '5.4.99.11': 'isomaltulose synthase',
        '5.4.99.15': '(1->4)-alpha-D-glucan 1-alpha-D-glucosylmutase',
        '5.4.99.16': 'maltose alpha-D-glucosyltransferase',
    }

    # Keyword normalisation for protein names when EC is ambiguous (3.2.1.-)
    # Order matters — more specific patterns first
    KEYWORD_NORMS = [
        (r'κ-carrageen',                'κ-carrageenase'),
        (r'kappa-carrageen',            'κ-carrageenase'),
        (r'ι-carrageen|iota-carrageen', 'ι-carrageenase'),
        (r'β-carrageen',                'β-carrageenase'),
        (r'b-/k-carrageen|b-/κ',        'β/κ-carrageenase'),
        (r'porphyran',                  'β-porphyranase'),
        (r'agarase|agarohydrolase',     'β-agarase'),
        (r'laminarinase|laminarin',     'laminarinase'),
        (r'lichenase|lichenin',         'lichenase'),
        (r'xyloglucan',                 'xyloglucan endotransglycosylase'),
        (r'1,3-1,4-glucanase|1,3.1,4.glucanase|1,3\(4\)-glucanase|1,3\(4\)glucanase', 'endo-β-1,3(4)-glucanase'),
        (r'1,3-glucanase|1,3.glucanase', 'endo-β-1,3-glucanase'),
        (r'hyaluronidase',              'hyaluronidase'),
        (r'galactanase',                'endo-β-galactanase'),
        (r'galactosidase',              'β-galactosidase'),
        (r'glucosyltransferase',        'glucosyltransferase'),
        (r'glucanase',                  'β-glucanase'),
    ]

    import re as _re

    def _greek(s):
        """Standardise alpha/beta/kappa/iota to Greek symbols."""
        s = re.sub(r'\balpha\b', 'α', s)
        s = re.sub(r'\bbeta\b',  'β', s)
        s = re.sub(r'\bkappa\b', 'κ', s)
        s = re.sub(r'\biota\b',  'ι', s)
        s = re.sub(r'\blambda\b','λ', s)
        return s

    def _normalise_name(ec_str, protein_name):
        # Try primary EC keys first (use first EC if comma-separated)
        for ec in [e.strip() for e in ec_str.split(',')]:
            if ec in EC_LABELS:
                return _greek(EC_LABELS[ec])
        # Fall back to keyword matching on lowercased protein name
        pn_lower = protein_name.lower()
        for pattern, label in KEYWORD_NORMS:
            if _re.search(pattern, pn_lower):
                return label
        # Last resort: strip parentheticals and return cleaned name
        name = protein_name.split('(')[0].strip()
        return ' '.join(name.split()) or 'unknown'

    result = {}
    for _, row in ref_meta.iterrows():
        acc      = str(row['accession']).replace('|', '_').replace('=', '_')
        ec_str   = str(row.get('ec_numbers', ''))
        pname    = str(row.get('protein_name', row.get('substrate', 'unknown')))
        result[acc] = _normalise_name(ec_str, pname)
    return result

# ── Sequence parsing ──────────────────────────────────────────────────────────

def parse_faa_annotations(faa_path, activity_map, sample_labels,
                          ref_label_map, ref_substrate_map=None):
    """
    Parse a family FASTA file and extract per-leaf annotation data.

    Args:
        faa_path:      path to family .faa file
        activity_map:  dict mapping gene_id to activity string
        sample_labels: dict mapping sample name to display label
        ref_label_map: dict mapping accession to reference label

    Returns:
        tuple of (samples, localisations, activities, labels) dicts
        each mapping leaf_id to its annotation value
    """
    samples       = {}
    localisations = {}
    activities    = {}
    labels        = {}

    for record in SeqIO.parse(faa_path, 'fasta'):
        parts = record.id.split('__')
        if len(parts) != 4:
            continue
        sample, gene_id, subfam, localisation = parts
        leaf_id = record.id

        samples[leaf_id]       = sample
        localisations[leaf_id] = localisation
        if sample == 'Reference':
            _ref_sub = ref_substrate_map.get(gene_id) if ref_substrate_map else None
            activities[leaf_id] = (
                f'reference: {_ref_sub}' if _ref_sub
                else activity_map.get(gene_id, 'unknown'))
        else:
            activities[leaf_id] = activity_map.get(gene_id, 'unknown')

        if sample == 'Reference':
            labels[leaf_id] = ref_label_map.get(
                gene_id, f'Reference {gene_id}')
        else:
            labels[leaf_id] = sample_labels.get(sample, sample)

    return samples, localisations, activities, labels


def _get_tree_leaf_ids(treefile_path):
    """Extract all leaf IDs from a Newick treefile."""
    if not os.path.exists(treefile_path):
        return None
    with open(treefile_path) as f:
        tree_str = f.read()
    # Leaf IDs are strings before : that follow ( or ,
    import re
    return set(re.findall(r'[,(]([^,()\[\]:]+):', tree_str))


# ── Phase 1: assign colours and write config ──────────────────────────────────

def assign_colours(seq_dir, output_dir, substrate, colours_file,
                   activity_file, ref_metadata=None,
                   sample_metadata=None, max_colours=None):
    """
    Phase 1: scan sequence files, assign colours, write colour_config.tsv.

    Args:
        seq_dir:         path to sequences directory
        output_dir:      substrate output directory
        substrate:       substrate name
        colours_file:    path to default_colours.tsv
        activity_file:   path to <substrate>_activity_annotated.tsv
        ref_metadata:    optional path to reference_metadata.tsv
        sample_metadata: optional path to metadata TSV with sample labels
        max_colours:     optional int — if set and sample count exceeds
                         palette size, generate this many HSL colours

    Returns:
        path to written colour_config.tsv
    """
    sample_palette, activity_palette = load_colour_palettes(colours_file)

    # Load activity annotations
    activity_map = {}
    if os.path.exists(activity_file):
        act_df = pd.read_csv(activity_file, sep='\t')
        for _, row in act_df.iterrows():
            gene_id  = str(row['Gene ID'])
            activity = str(row['activity'])
            activity_map[gene_id] = activity
            activity_map[gene_id.replace('|', '_').replace('=', '_')] = activity

    # Collect all samples and activities across all family FASTAs
    all_samples    = set()
    all_activities = set()

    for faa_file in sorted(os.listdir(seq_dir)):
        if not faa_file.endswith('.faa'):
            continue
        faa_path = os.path.join(seq_dir, faa_file)
        for record in SeqIO.parse(faa_path, 'fasta'):
            parts = record.id.split('__')
            if len(parts) != 4:
                continue
            sample, gene_id, _, _ = parts
            all_samples.add(sample)
            all_activities.add(activity_map.get(gene_id, 'unknown'))

    # Add reference substrates to activity set so they get colours
    if ref_metadata and os.path.exists(ref_metadata):
        _families = {f.replace(".faa", "").split("_", 1)[-1]
                     for f in os.listdir(seq_dir) if f.endswith(".faa")}
        ref_sub_map = load_ref_substrate_map(ref_metadata, _families)
        for sub_val in ref_sub_map.values():
            all_activities.add(f'reference: {sub_val}')
    sample_colours   = assign_sample_colours(all_samples, sample_palette,
                                              max_colours=max_colours)
    activity_colours = assign_activity_colours(all_activities,
                                               activity_palette)

    colour_config_path = os.path.join(
        output_dir, f'{substrate}_colour_config.tsv')
    write_colour_config(colour_config_path, sample_colours, activity_colours)

    return colour_config_path


# ── Phase 2: write iTOL annotation files ─────────────────────────────────────

def write_itol_annotations(seq_dir, output_dir, substrate,
                           colour_config_path, activity_file,
                           ref_metadata=None, sample_metadata=None):
    """
    Phase 2: read colour_config.tsv and write iTOL annotation files
    for all family trees.

    Args:
        seq_dir:            path to sequences directory
        output_dir:         substrate output directory
        substrate:          substrate name
        colour_config_path: path to colour_config.tsv (from Phase 1
                            or edited by user)
        activity_file:      path to <substrate>_activity_annotated.tsv
        ref_metadata:       optional path to reference_metadata.tsv
        sample_metadata:    optional path to metadata TSV with sample labels
    """
    sample_colours, activity_colours, localisation_colours = (
        load_colour_config(colour_config_path))

    # Load activity annotations
    activity_map = {}
    if os.path.exists(activity_file):
        act_df = pd.read_csv(activity_file, sep='\t')
        for _, row in act_df.iterrows():
            gene_id  = str(row['Gene ID'])
            activity = str(row['activity'])
            activity_map[gene_id] = activity
            activity_map[gene_id.replace('|', '_').replace('=', '_')] = activity

    # Load reference labels
    ref_label_map = {}
    ref_substrate_map = {}
    if ref_metadata:
        _families = {f.replace(".faa", "").split("_", 1)[-1]
                     for f in os.listdir(seq_dir) if f.endswith(".faa")}
        ref_label_map = load_ref_labels(ref_metadata, _families)
        ref_substrate_map = load_ref_substrate_map(ref_metadata, _families)

    itol_dir = os.path.join(output_dir, 'itol_annotations')
    os.makedirs(itol_dir, exist_ok=True)

    for faa_file in sorted(os.listdir(seq_dir)):
        if not faa_file.endswith('.faa'):
            continue

        family   = faa_file.replace('.faa', '')
        faa_path = os.path.join(seq_dir, faa_file)

        # Generate sample labels for samples present in this family
        samples_in_file = set()
        for record in SeqIO.parse(faa_path, 'fasta'):
            parts = record.id.split('__')
            if len(parts) == 4:
                samples_in_file.add(parts[0])

        sample_labels = generate_sample_labels(
            samples_in_file,
            metadata_file=sample_metadata,
            substrate=substrate
        )

        sample_data, localisation_data, activity_data, label_data = (
            parse_faa_annotations(
                faa_path, activity_map, sample_labels, ref_label_map))

        if not sample_data:
            continue

        # Filter to sequences present in the treefile
        family_key = family[len(substrate)+1:] if family.startswith(f'{substrate}_') else family
        treefile = os.path.join(output_dir, 'trees', f'{family_key}.treefile')
        tree_ids = _get_tree_leaf_ids(treefile)
        if tree_ids:
            sample_data       = {k: v for k, v in sample_data.items()       if k in tree_ids}
            localisation_data = {k: v for k, v in localisation_data.items() if k in tree_ids}
            activity_data     = {k: v for k, v in activity_data.items()     if k in tree_ids}
            label_data        = {k: v for k, v in label_data.items()        if k in tree_ids}
            # Add Reference__ sequences in tree but not in FAA (place mode)
            for _leaf in tree_ids:
                if _leaf.startswith('Reference__') and _leaf not in sample_data:
                    _parts = _leaf.split('__')
                    _acc   = _parts[1] if len(_parts) > 1 else _leaf
                    sample_data[_leaf]       = 'Reference'
                    localisation_data[_leaf] = 'characterised_reference'
                    _ref_sub = ref_substrate_map.get(_acc)
                    activity_data[_leaf]     = (
                        f'reference: {_ref_sub}' if _ref_sub else 'unknown')
                    label_data[_leaf]        = ref_label_map.get(_acc, f'Reference {_acc}')

        if not sample_data:
            continue

        print(f'Writing iTOL annotations for {family}...')

        write_colour_strip(
            os.path.join(itol_dir, f'{family}_sample.txt'),
            sample_data, sample_colours, 'Genome', 'Genome'
        )
        write_colour_strip(
            os.path.join(itol_dir, f'{family}_localisation.txt'),
            localisation_data, localisation_colours,
            'Localisation', 'Localisation'
        )
        write_colour_strip(
            os.path.join(itol_dir, f'{family}_activity.txt'),
            activity_data, activity_colours,
            'Activity', 'Enzymatic Activity'
        )
        write_branch_colours(
            os.path.join(itol_dir, f'{family}_branch_localisation.txt'),
            localisation_data, localisation_colours,
            'Localisation (branch)'
        )
        write_leaf_symbols(
            os.path.join(itol_dir, f'{family}_localisation_symbols.txt'),
            localisation_data, localisation_colours,
            'Localisation (symbol)'
        )
        write_labels(
            os.path.join(itol_dir, f'{family}_labels.txt'),
            label_data
        )
        write_label_styles(
            os.path.join(itol_dir, f'{family}_label_styles.txt'),
            label_data
        )

    print(f'\nDone. iTOL annotation files written to {itol_dir}')


# ── Main entry point ──────────────────────────────────────────────────────────

def generate_itol_annotations(seq_dir, output_dir, substrate,
                              colours_file, activity_file,
                              ref_metadata=None, sample_metadata=None,
                              colour_config_path=None,
                              max_colours=None):
    """
    Full iTOL annotation generation for one substrate.

    If colour_config_path is provided and exists, skips Phase 1 and
    uses the existing config directly. This allows users to edit colours
    and regenerate annotations without rerunning colour assignment.

    Args:
        seq_dir:            path to sequences directory
        output_dir:         substrate output directory
        substrate:          substrate name
        colours_file:       path to default_colours.tsv
        activity_file:      path to <substrate>_activity_annotated.tsv
        ref_metadata:       optional path to reference_metadata.tsv
        sample_metadata:    optional path to metadata TSV with sample labels
        colour_config_path: optional path to existing colour_config.tsv.
                            If None, Phase 1 runs automatically.
        max_colours:        optional int — if set and sample count exceeds
                            the default palette size, generate this many
                            HSL colours instead of cycling. Recommended
                            for datasets with more than 19 samples.
    """
    default_config = os.path.join(
        output_dir, f'{substrate}_colour_config.tsv')

    if colour_config_path and os.path.exists(colour_config_path):
        print(f"Using existing colour config: {colour_config_path}")
        print("Skipping colour assignment (Phase 1)")
    elif os.path.exists(default_config):
        print(f"Using existing colour config: {default_config}")
        print("Skipping colour assignment (Phase 1)")
        colour_config_path = default_config
    else:
        print("Phase 1: assigning colours...")
        colour_config_path = assign_colours(
            seq_dir, output_dir, substrate, colours_file,
            activity_file, ref_metadata=ref_metadata,
            sample_metadata=sample_metadata,
            max_colours=max_colours,
        )

    print("\nPhase 2: writing iTOL annotation files...")
    write_itol_annotations(
        seq_dir, output_dir, substrate, colour_config_path,
        activity_file, ref_metadata=ref_metadata,
        sample_metadata=sample_metadata
    )
