"""
Activity annotation for CAZyme hits.

Assigns primary EC numbers and human-readable activity labels to each
gene in a hits DataFrame, using:
  - EXPASY enzyme.dat for EC number -> activity name lookup
  - dbCAN fam-substrate-mapping.tsv as fallback for genes without EC

Reference sequence annotations are appended from REF_METADATA using
the activity column directly (no EC lookup needed since these are
experimentally characterised enzymes).
"""
import os
import re
import pandas as pd


# ── Activity name normalisation ───────────────────────────────────────────────

# Lookup table for known non-standard prefix variants in EXPASY/dbCAN output.
# Keys are the non-standard forms, values are the normalised forms.
# Extend this during testing as new variants are encountered.
PREFIX_NORMALISATION = {
    'a-gluco':   'alpha-gluco',
    'a-galacto': 'alpha-galacto',
    'a-manno':   'alpha-manno',
    'a-fuco':    'alpha-fuco',
    'a-xylo':    'alpha-xylo',
    'b-gluco':   'beta-gluco',
    'b-galacto': 'beta-galacto',
    'b-manno':   'beta-manno',
    'b-xylo':    'beta-xylo',
    'b-fuco':    'beta-fuco',
}

# Family class prefixes and their human-readable names.
# Used to generate fallback labels for verbose descriptions.
FAMILY_CLASS_NAMES = {
    'GH':  'glycoside hydrolase',
    'PL':  'polysaccharide lyase',
    'CE':  'carbohydrate esterase',
    'AA':  'auxiliary activity enzyme',
    'CBM': 'carbohydrate-binding module',
}


def normalise_activity(activity, family):
    """
    Normalise an activity string from EXPASY or dbCAN output.

    Handles two classes of problem:
    1. Verbose module descriptions (> 60 chars) — replaced with a clean
       label derived from the family class, e.g.:
         'Modules of approx. 100 residues with glycogen-binding...'
         -> 'carbohydrate-binding module (CBM48)'
    2. Non-standard alpha/beta prefix variants, e.g.:
         'a-glucosidase' -> 'alpha-glucosidase'
         'b-galactosidase' -> 'beta-galactosidase'

    Args:
        activity: raw activity string from EXPASY or dbCAN
        family:   matched CAZyme family string (e.g. 'CBM48', 'GH16')

    Returns:
        Normalised activity string.
    """
    if not activity or str(activity).strip() in ['unknown', 'nan', '-', '']:
        return 'unknown'

    activity = str(activity).strip()

    # Replace verbose descriptions with clean family-derived label
    if len(activity) > 60:
        family_str = str(family).strip()
        for prefix, class_name in FAMILY_CLASS_NAMES.items():
            if family_str.startswith(prefix):
                return f'{class_name} ({family_str})'
        # Family prefix not recognised — truncate with ellipsis as last resort
        return activity[:57] + '...'

    # Normalise non-standard alpha/beta prefixes
    for variant, canonical in PREFIX_NORMALISATION.items():
        if variant in activity:
            activity = activity.replace(variant, canonical)

    return activity


# ── Data loading ──────────────────────────────────────────────────────────────

def load_ec_names(expasy_file):
    """
    Parse EXPASY enzyme.dat and return a dict of EC -> accepted name.
    Uses the DE (description) field as the primary name.

    Args:
        expasy_file: path to EXPASY enzyme.dat

    Returns:
        dict mapping EC string to activity name string
    """
    ec_names = {}
    current_ec = None
    with open(expasy_file) as f:
        for line in f:
            if line.startswith('ID'):
                current_ec = line.split()[1].strip()
            elif line.startswith('DE') and current_ec:
                name = line[5:].strip().rstrip('.')
                if not any(skip in name for skip in
                           ['Deleted', 'Transferred', 'Would be']):
                    ec_names[current_ec] = name
    return ec_names


def load_family_activities(fam_sub_map):
    """
    Parse dbCAN fam-substrate-mapping.tsv and return a dict of
    family -> list of (ec, activity) tuples.
    Used as fallback when no EC number is available for a gene.

    Args:
        fam_sub_map: path to dbCAN fam-substrate-mapping.tsv

    Returns:
        dict mapping family string to list of (ec, activity) tuples
    """
    fam_map = pd.read_csv(fam_sub_map, sep='\t', header=None,
                          names=['substrate_cat', 'specific_substrate',
                                 'family', 'activity', 'ec_number'])
    fam_map['family']    = fam_map['family'].astype(str).str.strip()
    fam_map['ec_number'] = fam_map['ec_number'].astype(str).str.strip()
    fam_map['activity']  = fam_map['activity'].astype(str).str.strip()

    family_ec_map = {}
    for _, row in fam_map.iterrows():
        fam = row['family']
        ec  = row['ec_number']
        act = row['activity']
        if fam not in family_ec_map:
            family_ec_map[fam] = []
        if ec not in ['nan', '', '-']:
            family_ec_map[fam].append((ec, act))
    return family_ec_map


# ── EC extraction ─────────────────────────────────────────────────────────────

def extract_primary_ec(ec_string):
    """
    Extract the most supported EC number from a dbCAN EC string.

    dbCAN EC strings can be complex, e.g.:
      '3.2.1.39:1'              -> single EC with count
      '3.2.1.39:1;3.2.1.6:2'   -> multiple ECs, take highest count
      '3.2.1.39:1|-|-'          -> multi-tool, take first non-empty
      '-'                       -> no EC assigned

    Returns the EC number with the highest count, or '-' if none found.

    Args:
        ec_string: raw EC string from dbCAN overview.tsv

    Returns:
        EC number string, or '-' if none found
    """
    if pd.isna(ec_string) or str(ec_string).strip() in ['-', '', 'nan']:
        return '-'

    segments = [s.strip() for s in str(ec_string).split('|')
                if s.strip() not in ['-', '', 'nan']]
    if not segments:
        return '-'
    ec_string = segments[0]

    best_ec    = '-'
    best_count = 0
    for part in ec_string.split(';'):
        part  = part.strip()
        match = re.match(r'([\d]+\.[\d]+\.[\d]+\.[\d\-]+)(?::(\d+))?', part)
        if match:
            ec    = match.group(1)
            count = int(match.group(2)) if match.group(2) else 1
            if ec not in ['-', ''] and count > best_count:
                best_ec    = ec
                best_count = count
    return best_ec


# ── Activity label resolution ─────────────────────────────────────────────────

def get_activity_label(primary_ec, matched_family, ec_activities,
                       family_activities):
    """
    Get a normalised human-readable activity label for a gene.

    Priority:
    1. Look up primary EC in EXPASY ec_activities dict
    2. For partial EC numbers ending in '.-', find the best matching
       name from ec_activities and return it normalised
    3. Fall back to family-level activity from dbCAN fam-substrate-mapping
    4. Return 'unknown' if no source has an entry

    Args:
        primary_ec:        EC string from extract_primary_ec()
        matched_family:    matched CAZyme family string (e.g. 'GH16')
        ec_activities:     dict from load_ec_names()
        family_activities: dict from load_family_activities()

    Returns:
        Normalised activity label string
    """
    # 1. Direct EC lookup
    if primary_ec != '-' and primary_ec in ec_activities:
        return normalise_activity(ec_activities[primary_ec], matched_family)

    # 2. Partial EC lookup — use matched names rather than hardcoding
    if primary_ec != '-' and primary_ec.endswith('.-'):
        base = primary_ec[:-2]
        matches = [name for ec, name in ec_activities.items()
                   if ec.startswith(base + '.')]
        if matches:
            return normalise_activity(matches[0], matched_family)

    # 3. Family-level fallback from dbCAN mapping
    top_family_match = re.match(r'([A-Za-z]+[0-9]+)', str(matched_family))
    if top_family_match:
        fam = top_family_match.group(1)
        if fam in family_activities and family_activities[fam]:
            raw_activity = family_activities[fam][0][1]
            return normalise_activity(raw_activity, fam)

    return 'unknown'


# ── Main annotation functions ─────────────────────────────────────────────────

def annotate_hits(hits_df, expasy_file, fam_sub_map):
    """
    Assign primary EC numbers and activity labels to a genomic hits DataFrame.

    Args:
        hits_df:     DataFrame of family hits filtered to a single substrate
        expasy_file: path to EXPASY enzyme.dat
        fam_sub_map: path to dbCAN fam-substrate-mapping.tsv

    Returns:
        hits_df with 'primary_ec' and 'activity' columns added
    """
    ec_activities     = load_ec_names(expasy_file)
    family_activities = load_family_activities(fam_sub_map)

    print(f"Loaded {len(ec_activities)} EC names from EXPASY")
    print(f"Loaded activity fallbacks for {len(family_activities)} "
          f"CAZyme families")

    hits_df = hits_df.copy()
    hits_df['primary_ec'] = hits_df['EC#'].apply(extract_primary_ec)
    hits_df['activity']   = hits_df.apply(
        lambda row: get_activity_label(
            row['primary_ec'], row['matched_family'],
            ec_activities, family_activities
        ),
        axis=1
    )
    return hits_df


def annotate_references(ref_metadata, substrate):
    """
    Build activity annotation rows for reference sequences.

    Reference sequences use the activity column from REF_METADATA directly
    since these are experimentally characterised enzymes — no EC lookup needed.

    Args:
        ref_metadata: path to reference_metadata.tsv
        substrate:    substrate name to filter by

    Returns:
        DataFrame of reference annotation rows, or empty DataFrame if
        ref_metadata does not exist or has no rows for this substrate
    """
    if not os.path.exists(ref_metadata):
        print(f"WARNING: No reference metadata found at {ref_metadata}")
        return pd.DataFrame()

    ref_meta = pd.read_csv(ref_metadata, sep='\t')
    ref_meta = ref_meta[ref_meta['substrate'] == substrate].copy()

    if ref_meta.empty:
        return pd.DataFrame()

    # Determine activity — support legacy schema (activity column)
    # and current schema (protein_name + ec_numbers)
    if 'activity' in ref_meta.columns:
        ref_meta['_activity'] = ref_meta['activity'].fillna('unknown')
    elif 'protein_name' in ref_meta.columns:
        def _make_activity(row):
            name = str(row.get('protein_name', '')).strip()
            ec   = str(row.get('ec_numbers', '')).strip()
            if ec and ec not in ('', 'nan', '-'):
                return f"{name} [{ec}]" if name else ec
            return name or 'unknown'
        ref_meta['_activity'] = ref_meta.apply(_make_activity, axis=1)
    else:
        ref_meta['_activity'] = 'unknown'

    if 'ec_numbers' in ref_meta.columns:
        ref_meta['_ec'] = ref_meta['ec_numbers'].fillna('-')
    else:
        ref_meta['_ec'] = '-'

    rows = []
    for _, row in ref_meta.iterrows():
        rows.append({
            'Gene ID':              str(row['accession']),
            'sample':               'Reference',
            'substrate_category':   substrate,
            'matched_family':       str(row['family']),
            'localisation':         'characterised_reference',
            'subfamily_annotation': str(row.get('subfamily',
                                                 row['family'])),
            'EC#':                  str(row['_ec']),
            'primary_ec':           str(row['_ec']),
            'activity':             str(row['_activity']),
        })

    return pd.DataFrame(rows)
