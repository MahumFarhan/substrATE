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
    unique = sorted(set(activities))
    return {act: activity_palette[i % len(activity_palette)]
            for i, act in enumerate(unique)}


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
    with open(out_file, 'w') as f:
        f.write('DATASET_COLORSTRIP\n')
        f.write('SEPARATOR TAB\n')
        f.write(f'DATASET_LABEL\t{title}\n')
        f.write('COLOR\t#000000\n')
        f.write(f'LEGEND_TITLE\t{legend_title}\n')
        f.write('LEGEND_SHAPES\t' +
                '\t'.join(['1'] * len(colour_map)) + '\n')
        f.write('LEGEND_COLORS\t' +
                '\t'.join(colour_map.values()) + '\n')
        f.write('LEGEND_LABELS\t' +
                '\t'.join(colour_map.keys()) + '\n')
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


# ── Reference label loading ───────────────────────────────────────────────────

def load_ref_labels(ref_metadata, substrate):
    """
    Build accession -> label lookup from reference_metadata.tsv.

    Args:
        ref_metadata: path to reference_metadata.tsv
        substrate:    substrate name to filter by

    Returns:
        dict mapping cleaned accession string to label string
    """
    if not os.path.exists(ref_metadata):
        return {}
    ref_meta = pd.read_csv(ref_metadata, sep='\t')
    ref_meta = ref_meta[ref_meta['substrate'] == substrate].copy()
    return {
        str(row['accession']).replace('|', '_').replace('=', '_'):
        str(row['label'])
        for _, row in ref_meta.iterrows()
    }


# ── Sequence parsing ──────────────────────────────────────────────────────────

def parse_faa_annotations(faa_path, activity_map, sample_labels,
                          ref_label_map):
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
        activities[leaf_id]    = activity_map.get(gene_id, 'unknown')

        if sample == 'Reference':
            labels[leaf_id] = ref_label_map.get(
                gene_id, f'Reference {gene_id}')
        else:
            labels[leaf_id] = sample_labels.get(sample, sample)

    return samples, localisations, activities, labels


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
    if ref_metadata:
        ref_label_map = load_ref_labels(ref_metadata, substrate)

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
