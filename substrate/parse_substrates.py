"""
Substrate and CAZyme family parsing from dbCAN fam-substrate-mapping.tsv.

Provides two functions:
  1. get_families_for_substrate() — derives the CAZyme family list for
     a substrate by searching fam-substrate-mapping.tsv. This replaces
     the one-time manual get_families.py step and runs automatically
     at pipeline startup for any substrate not in the built-in FAMILY_MAP.

  2. validate_substrate() — checks whether a substrate is built-in or
     can be derived from the mapping file, and raises an informative
     error if neither is true.

This module is the foundation for supporting user-defined substrates
via --substrate_terms in a future update.
"""
import os
import re
import pandas as pd
from substrate.classify_pul import FAMILY_MAP, SUBSTRATE_TERMS

# ── Path to data files ────────────────────────────────────────────────────────

_DATA_DIR              = os.path.join(os.path.dirname(__file__), 'data')
ACTIVITY_PATTERNS_FILE = os.path.join(_DATA_DIR, 'activity_patterns.tsv')

PATTERN_MODES = ('permissive', 'strict')


# ── Pattern loading ───────────────────────────────────────────────────────────

def load_patterns(substrate=None, patterns_file=None,
                  pattern_mode='permissive'):
    """
    Load activity patterns from activity_patterns.tsv.

    In permissive mode (default), all patterns are loaded.
    In strict mode, only patterns with mode='strict' are loaded,
    giving conservative, substrate-specific matching that avoids
    false positives from shared enzyme activities.

    The 'mode' column is optional for backwards compatibility — files
    without it are treated as all-permissive.

    Args:
        substrate:     optional substrate name to filter by
        patterns_file: path to activity_patterns.tsv (uses default
                       if None)
        pattern_mode:  'permissive' (default) or 'strict'

    Returns:
        DataFrame with columns: substrate, pattern, source, reviewed,
        mode.  Empty DataFrame if file does not exist.

    Raises:
        ValueError if pattern_mode is not 'permissive' or 'strict'
    """
    if pattern_mode not in PATTERN_MODES:
        raise ValueError(
            f"pattern_mode must be one of {PATTERN_MODES}, "
            f"got '{pattern_mode}'"
        )

    if patterns_file is None:
        patterns_file = ACTIVITY_PATTERNS_FILE

    if not os.path.exists(patterns_file):
        return pd.DataFrame(
            columns=['substrate', 'pattern', 'source',
                     'reviewed', 'mode']
        )

    df = pd.read_csv(patterns_file, sep='\t', dtype=str)

    # Back-compat: files without a mode column are treated as permissive
    if 'mode' not in df.columns:
        df['mode'] = 'permissive'

    if substrate is not None:
        df = df[df['substrate'] == substrate]

    if pattern_mode == 'strict':
        df = df[df['mode'] == 'strict']

    return df.reset_index(drop=True)


# ── Family derivation ─────────────────────────────────────────────────────────

def get_families_for_substrate(substrate, fam_sub_map,
                               substrate_terms=None):
    """
    Derive the CAZyme family list for a substrate from the dbCAN
    fam-substrate-mapping.tsv file.

    If the substrate is in the built-in FAMILY_MAP, that list is
    returned directly without reading the mapping file.

    For substrates not in FAMILY_MAP, the mapping file is searched
    using substrate_terms (if provided) or the substrate name itself
    as the search term.

    Args:
        substrate:       substrate name string
        fam_sub_map:     path to dbCAN fam-substrate-mapping.tsv
        substrate_terms: optional list of search terms to use instead
                         of the substrate name. Used for substrates
                         with multiple synonyms (e.g. ['fucoidan',
                         'sulfated fucan']).

    Returns:
        sorted list of CAZyme family strings

    Raises:
        ValueError if no families are found for the substrate
    """
    # Return built-in list directly if available
    if substrate in FAMILY_MAP:
        return FAMILY_MAP[substrate]

    # Load mapping file
    mapping = pd.read_csv(
        fam_sub_map, sep='\t', header=None,
        names=['substrate_cat', 'specific_substrate',
               'family', 'activity', 'ec_number']
    )

    # Use provided terms or fall back to substrate name
    terms   = substrate_terms or [substrate]
    pattern = '|'.join(re.escape(t) for t in terms)

    hits = mapping[
        mapping['specific_substrate'].str.contains(
            pattern, case=False, na=False) |
        mapping['activity'].str.contains(
            pattern, case=False, na=False)
    ]

    families = sorted(hits['family'].dropna().unique().tolist())

    if not families:
        raise ValueError(
            f"No CAZyme families found for substrate '{substrate}' "
            f"in {fam_sub_map}.\n"
            f"Search terms used: {terms}\n"
            f"For a custom substrate, provide search terms with "
            f"--substrate_terms."
        )

    print(f"Derived {len(families)} families for '{substrate}' "
          f"from dbCAN mapping: {families}")
    return families


def get_terms_for_substrate(substrate, substrate_terms=None):
    """
    Get the substrate search terms for a substrate.

    Returns built-in terms from SUBSTRATE_TERMS if available,
    otherwise returns the provided substrate_terms or the substrate
    name itself as a single-item list.

    Args:
        substrate:       substrate name string
        substrate_terms: optional list of custom search terms

    Returns:
        list of search term strings
    """
    if substrate in SUBSTRATE_TERMS:
        return SUBSTRATE_TERMS[substrate]
    if substrate_terms:
        return substrate_terms
    return [substrate]


# ── Substrate validation ──────────────────────────────────────────────────────

def validate_substrate(substrate, fam_sub_map, substrate_terms=None):
    """
    Validate that a substrate is usable by the pipeline.

    Checks in order:
      1. Built-in substrate (in FAMILY_MAP) — always valid
      2. Derivable from fam-substrate-mapping.tsv — valid if families found
      3. Neither — raises ValueError with helpful message

    Args:
        substrate:       substrate name string
        fam_sub_map:     path to dbCAN fam-substrate-mapping.tsv
        substrate_terms: optional list of search terms for custom substrates

    Returns:
        list of CAZyme family strings

    Raises:
        ValueError if substrate cannot be validated
    """
    if substrate in FAMILY_MAP:
        print(f"'{substrate}' is a built-in substrate "
              f"({len(FAMILY_MAP[substrate])} families)")
        return FAMILY_MAP[substrate]

    print(f"'{substrate}' is not a built-in substrate — "
          f"deriving families from dbCAN mapping file...")
    return get_families_for_substrate(
        substrate, fam_sub_map, substrate_terms=substrate_terms)


# ── Summary ───────────────────────────────────────────────────────────────────

def list_builtin_substrates():
    """
    Print a summary of all built-in substrates and their family counts.
    """
    print("Built-in substrates:")
    for substrate, families in sorted(FAMILY_MAP.items()):
        terms = SUBSTRATE_TERMS.get(substrate, [substrate])
        print(f"  {substrate:<12} "
              f"{len(families)} families, "
              f"search terms: {terms}")


# ── Pattern overlap checking ──────────────────────────────────────────────────

def check_pattern_overlap(substrates, patterns_file=None,
                          overlap_threshold=5,
                          pattern_mode='permissive'):
    """
    Check for significant pattern overlap between each substrate being
    analysed and ALL other substrates in activity_patterns.tsv.

    This runs even for single-substrate analyses, since a substrate's
    patterns may overlap with others not being analysed — which still
    affects how results should be interpreted.

    For each substrate being analysed, two classes of overlap are
    reported separately:
      - Overlap with another substrate also being analysed in this run
        (CGCs may appear in both outputs)
      - Overlap with a substrate NOT being analysed in this run
        (results may contain enzymes active on that substrate)

    Args:
        substrates:        list of substrate name strings being analysed
        patterns_file:     path to activity_patterns.tsv (uses default
                           if None)
        overlap_threshold: minimum shared patterns to trigger warning
                           (default: 5)
        pattern_mode:      'permissive' (default) or 'strict'

    Returns:
        list of (substrate1, substrate2, shared_patterns, in_session)
        tuples for pairs exceeding the threshold, where in_session is
        True if both substrates are being analysed in this run.
        Returns empty list if no overlaps found.
    """
    df = load_patterns(patterns_file=patterns_file,
                       pattern_mode=pattern_mode)
    if df.empty:
        return []

    session_set = set(substrates)

    # Build pattern sets for ALL substrates in the file
    all_patterns = {}
    for substrate in df['substrate'].unique():
        all_patterns[substrate] = set(
            df[df['substrate'] == substrate]['pattern'].tolist()
        )

    overlapping_pairs = []

    # For each substrate being analysed, check against all others
    for s1 in sorted(substrates):
        if s1 not in all_patterns:
            continue
        p1 = all_patterns[s1]

        for s2 in sorted(all_patterns):
            if s2 == s1:
                continue
            # Avoid reporting the same pair twice
            if s2 in session_set and s2 < s1:
                continue

            shared = p1 & all_patterns[s2]
            if len(shared) >= overlap_threshold:
                in_session = s2 in session_set
                overlapping_pairs.append(
                    (s1, s2, sorted(shared), in_session))

    if overlapping_pairs:
        # Separate into in-session and cross-session overlaps
        in_session_pairs = [
            (s1, s2, shared) for s1, s2, shared, ins
            in overlapping_pairs if ins
        ]
        cross_session_pairs = [
            (s1, s2, shared) for s1, s2, shared, ins
            in overlapping_pairs if not ins
        ]

        print(f"\n{'='*60}")
        print(f"WARNING: Activity pattern overlap detected "
              f"[mode: {pattern_mode}]")
        print(f"{'='*60}")

        if in_session_pairs:
            print(f"\nSubstrates being analysed in this run that share")
            print(f">= {overlap_threshold} activity patterns")
            print(f"(CGCs may appear in both outputs):\n")
            for s1, s2, shared in in_session_pairs:
                print(f"  {s1} / {s2}:")
                print(f"    {len(shared)} shared: "
                      f"{', '.join(shared[:8])}"
                      f"{'...' if len(shared) > 8 else ''}")

        if cross_session_pairs:
            print(f"\nSubstrates NOT in this run that share patterns")
            print(f"with your selected substrate(s)")
            print(f"(results may contain enzymes active on these")
            print(f"substrates):\n")
            for s1, s2, shared in cross_session_pairs:
                print(f"  {s1} overlaps with {s2}:")
                print(f"    {len(shared)} shared: "
                      f"{', '.join(shared[:8])}"
                      f"{'...' if len(shared) > 8 else ''}")

        print(f"\n--overlap_threshold {overlap_threshold} "
              f"(use 0 to suppress this warning)")
        print(f"--pattern_mode {pattern_mode} "
              f"(use 'strict' to reduce overlap)")
        print(f"{'='*60}\n")

    return overlapping_pairs


# ── Auto pattern derivation ───────────────────────────────────────────────────

GENERIC_TOKENS = {
    'enzyme', 'protein', 'activity', 'family', 'domain', 'module',
    'binding', 'function', 'related', 'like', 'type', 'class',
    'subunit', 'chain', 'component', 'factor', 'involved', 'putative',
    'hypothetical', 'unknown', 'uncharacterised', 'uncharacterized',
    'predicted', 'probable', 'possible', 'similar', 'homolog',
    'homologue', 'superfamily', 'fold', 'repeat', 'motif', 'site',
    'shown', 'demonstrated', 'case', 'cases', 'found', 'appended',
    'modules', 'approx', 'residues', 'terminus', 'structure',
    'crystal', 'pmid', 'doi', 'reference', 'strain', 'species',
    'bacteria', 'from', 'with', 'into', 'also', 'have', 'been',
    'that', 'this', 'which', 'these', 'their', 'other', 'some',
    'various', 'several', 'many', 'most', 'both', 'each', 'only',
    'mainly', 'usually', 'often', 'always', 'never', 'including',
    'derived', 'based', 'specific', 'general', 'broad', 'wide',
    'high', 'low', 'large', 'small', 'major', 'minor', 'novel',
    'new', 'first', 'second', 'wild', 'native', 'recombinant',
    'purified', 'mature', 'truncated', 'full', 'length', 'terminal',
    'internal', 'central', 'more', 'than', 'such', 'well', 'three',
}

MIN_TOKEN_LEN = 5
MAX_NAME_LENGTH = 60


def _extract_tokens(name_string):
    """Extract meaningful tokens from a Name string."""
    name_string = str(name_string).strip()
    if len(name_string) > MAX_NAME_LENGTH:
        return set()
    tokens = re.split(r'[\s\-,/();:]+', name_string.lower())
    result = set()
    for token in tokens:
        token = token.strip('.,\'\"[]{}')
        if (len(token) >= MIN_TOKEN_LEN
                and token not in GENERIC_TOKENS
                and not token.isnumeric()
                and not re.match(r'^\d+\.\d+', token)
                and not re.match(r'^[a-z]{1,2}\d+$', token)
                and not token.startswith('pmid')):
            result.add(token)
    return result


def get_all_substrates(fam_sub_map):
    """
    Extract all unique substrate names from fam-substrate-mapping.tsv.

    Args:
        fam_sub_map: path to dbCAN fam-substrate-mapping.tsv

    Returns:
        sorted list of substrate name strings
    """
    df = pd.read_csv(fam_sub_map, sep='\t', dtype=str)
    df.columns = [c.strip() for c in df.columns]
    substrates = df['Substrate_curated'].dropna().unique().tolist()
    return sorted(set(s.strip().lower() for s in substrates if s.strip()))


def auto_derive_patterns(substrate, fam_sub_map, substrate_terms=None):
    """
    Automatically derive activity patterns for a substrate from the
    dbCAN fam-substrate-mapping.tsv Name column.

    Args:
        substrate:       substrate name string
        fam_sub_map:     path to dbCAN fam-substrate-mapping.tsv
        substrate_terms: optional list of search terms

    Returns:
        list of pattern strings (lowercase, deduplicated, sorted)
    """
    df = pd.read_csv(fam_sub_map, sep='\t', dtype=str)
    df.columns = [c.strip() for c in df.columns]
    df = df.fillna('')

    terms   = substrate_terms or [substrate]
    pattern = '|'.join(re.escape(t) for t in terms)

    hits = df[
        df['Substrate_curated'].str.contains(
            pattern, case=False, na=False)
    ]

    if hits.empty:
        return []

    all_tokens = set()
    for term in terms:
        for part in re.split(r'[\s\-]+', term.lower()):
            if len(part) >= MIN_TOKEN_LEN:
                all_tokens.add(part)

    for name in hits['Name'].dropna():
        all_tokens.update(_extract_tokens(name))

    for sub in hits['Substrate_curated'].dropna():
        if sub.strip():
            all_tokens.update(_extract_tokens(sub))

    return sorted(all_tokens)


def write_pattern_review_report(substrate, hits_df, activity_file,
                                patterns, output_dir,
                                pattern_mode='permissive'):
    """
    Write a pattern review report TSV to help users curate patterns.

    Args:
        substrate:     substrate name string
        hits_df:       family hits DataFrame for this substrate
        activity_file: path to <substrate>_activity_annotated.tsv
        patterns:      list of pattern strings used
        output_dir:    substrate output directory
        pattern_mode:  'permissive' or 'strict' — recorded in report
    """
    if not os.path.exists(activity_file):
        return

    act_df = pd.read_csv(activity_file, sep='\t')
    act_df = act_df[act_df['sample'] != 'Reference'].copy()

    report_rows = []
    for pattern in sorted(patterns):
        matched = act_df[
            act_df['activity'].str.contains(
                pattern, case=False, na=False)
        ]
        report_rows.append({
            'pattern':       pattern,
            'genes_matched': len(matched),
            'activities':    '|'.join(matched['activity'].unique()),
            'pattern_mode':  pattern_mode,
            'reviewed':      False,
        })

    report_df   = pd.DataFrame(report_rows)
    report_path = os.path.join(
        output_dir, f'{substrate}_pattern_review.tsv')
    report_df.to_csv(report_path, sep='\t', index=False)
    print(f"Pattern review report written to: {report_path}")


def update_patterns_file(substrate, patterns, patterns_file=None,
                         source='auto_derived', reviewed=False,
                         mode='permissive'):
    """
    Add new patterns for a substrate to activity_patterns.tsv.

    Only adds patterns not already present. Existing patterns are
    not modified.

    Args:
        substrate:     substrate name string
        patterns:      list of pattern strings to add
        patterns_file: path to activity_patterns.tsv (uses default
                       if None)
        source:        'auto_derived' or 'curated'
        reviewed:      whether patterns have been manually reviewed
        mode:          'permissive' (default) or 'strict'
    """
    if patterns_file is None:
        patterns_file = ACTIVITY_PATTERNS_FILE

    if os.path.exists(patterns_file):
        existing = pd.read_csv(patterns_file, sep='\t')
        if 'mode' not in existing.columns:
            existing['mode'] = 'permissive'
        existing_patterns = set(
            existing[existing['substrate'] == substrate]['pattern']
        )
    else:
        existing = pd.DataFrame(
            columns=['substrate', 'pattern', 'source',
                     'reviewed', 'mode'])
        existing_patterns = set()

    new_rows = []
    for pattern in patterns:
        if pattern not in existing_patterns:
            new_rows.append({
                'substrate': substrate,
                'pattern':   pattern,
                'source':    source,
                'reviewed':  reviewed,
                'mode':      mode,
            })

    if new_rows:
        new_df  = pd.DataFrame(new_rows)
        updated = pd.concat([existing, new_df], ignore_index=True)
        updated.to_csv(patterns_file, sep='\t', index=False)
        print(f"Added {len(new_rows)} new patterns for '{substrate}'")
    else:
        print(f"No new patterns to add for '{substrate}'")


# ── Survey ────────────────────────────────────────────────────────────────────

def survey_substrates(cgc_output_dir, fam_sub_map,
                      pul_mode='bacteroidetes', min_cazymes=2):
    """
    Scan dbCAN output for all substrates in fam-substrate-mapping.tsv
    and report hit counts per substrate.

    Args:
        cgc_output_dir: path to dbCAN cgc_output directory
        fam_sub_map:    path to dbCAN fam-substrate-mapping.tsv
        pul_mode:       PUL classification mode (see classify_cgc)
        min_cazymes:    minimum substrate CAZymes required per CGC

    Returns:
        DataFrame with columns: substrate, n_genomes, n_canonical_pul,
        n_non_canonical, families. Sorted by n_canonical_pul descending.
        Empty substrates are excluded.
    """
    df = pd.read_csv(fam_sub_map, sep='\t', dtype=str)
    df.columns = [c.strip() for c in df.columns]
    df = df.fillna('')

    # Build substrate -> families mapping from Substrate_curated
    substrate_families = {}
    for substrate in df['Substrate_curated'].unique():
        substrate_clean = substrate.strip().lower()
        if not substrate_clean:
            continue
        families = df[
            df['Substrate_curated'].str.lower().str.strip()
            == substrate_clean
        ]['Family'].dropna().unique().tolist()
        if families:
            substrate_families[substrate_clean] = sorted(families)

    print(f"Surveying {len(substrate_families)} substrates across "
          f"dbCAN output in {cgc_output_dir}...")

    results = []

    for substrate, families in substrate_families.items():
        n_genomes       = 0
        n_canonical     = 0
        n_non_canonical = 0

        for sample_dir in sorted(os.listdir(cgc_output_dir)):
            full_path = os.path.join(cgc_output_dir, sample_dir)
            if not os.path.isdir(full_path):
                continue

            over_file = os.path.join(full_path, 'overview.tsv')
            cgc_file  = os.path.join(full_path, 'cgc_standard_out.tsv')

            if not os.path.exists(over_file):
                continue

            over_df = pd.read_csv(over_file, sep='\t', dtype=str)
            over_df.columns = over_df.columns.str.strip()
            over_df = over_df.fillna('')

            over_df['all_annotations'] = (
                over_df.get('dbCAN_hmm', pd.Series([''] * len(over_df))).astype(str) + '|' +
                over_df.get('dbCAN_sub', pd.Series([''] * len(over_df))).astype(str) + '|' +
                over_df.get('DIAMOND',   pd.Series([''] * len(over_df))).astype(str) + '|' +
                over_df.get('Recommend Results', pd.Series([''] * len(over_df))).astype(str)
            )

            has_hit = False
            for family in families:
                mask = over_df['all_annotations'].str.contains(
                    rf'(?<![0-9]){re.escape(family)}(?![0-9])',
                    case=False, na=False, regex=True
                )
                if mask.any():
                    has_hit = True
                    break

            if not has_hit:
                continue

            n_genomes += 1

            if os.path.exists(cgc_file):
                cgc_df = pd.read_csv(cgc_file, sep='\t', dtype=str)
                cgc_df.columns = cgc_df.columns.str.strip()
                cgc_df['Gene Annotation'] = (
                    cgc_df['Gene Annotation'].fillna('')
                )
                for cgc_id in cgc_df['CGC#'].unique():
                    classification = classify_cgc(
                        cgc_id, cgc_df, families,
                        pul_mode=pul_mode,
                        min_cazymes=min_cazymes,
                    )
                    if classification == 'canonical_PUL':
                        n_canonical += 1
                    elif classification == 'non_canonical_CGC':
                        n_non_canonical += 1

        if n_genomes > 0:
            results.append({
                'substrate':       substrate,
                'n_genomes':       n_genomes,
                'n_canonical_pul': n_canonical,
                'n_non_canonical': n_non_canonical,
                'families':        '|'.join(families),
            })

    if not results:
        return pd.DataFrame(columns=[
            'substrate', 'n_genomes', 'n_canonical_pul',
            'n_non_canonical', 'families'
        ])

    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values(
        ['n_canonical_pul', 'n_genomes'],
        ascending=False
    ).reset_index(drop=True)

    return results_df


def print_survey_results(results_df, total_substrates):
    """
    Print a formatted survey results table to the terminal.

    Args:
        results_df:       DataFrame from survey_substrates()
        total_substrates: total number of substrates checked
    """
    if results_df.empty:
        print("No substrate hits found in dbCAN output.")
        return

    n_no_hits = total_substrates - len(results_df)

    col_w = {
        'substrate':    max(20, results_df['substrate'].str.len().max() + 2),
        'n_genomes':    10,
        'canonical':    18,
        'noncanonical': 16,
    }

    header = (
        f"  {'Substrate':<{col_w['substrate']}}"
        f"{'Genomes':<{col_w['n_genomes']}}"
        f"{'Canonical PULs':<{col_w['canonical']}}"
        f"{'Non-canonical':<{col_w['noncanonical']}}"
    )
    divider = '  ' + '─' * (sum(col_w.values()))

    print(f"\n{header}")
    print(divider)

    for _, row in results_df.iterrows():
        print(
            f"  {row['substrate']:<{col_w['substrate']}}"
            f"{row['n_genomes']:<{col_w['n_genomes']}}"
            f"{row['n_canonical_pul']:<{col_w['canonical']}}"
            f"{row['n_non_canonical']:<{col_w['noncanonical']}}"
        )

    print(divider)
    print(f"  {len(results_df)} substrates with hits. "
          f"{n_no_hits} had no hits.")
