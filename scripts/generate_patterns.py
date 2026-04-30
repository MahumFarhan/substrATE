"""
Generate activity_patterns.tsv from dbCAN fam-substrate-mapping.tsv.

This script derives activity patterns for all built-in substrates
directly from the dbCAN database, replacing any manually curated
or hardcoded patterns.

Run this script when:
  - Setting up substrATE for the first time
  - Updating to a new version of the dbCAN database
  - Adding new built-in substrates

Usage:
    python scripts/generate_patterns.py \
        --fam_sub_map /path/to/fam-substrate-mapping.tsv \
        --output substrate/data/activity_patterns.tsv
"""
import os
import re
import argparse
import pandas as pd

SUBSTRATE_SEARCH_CONFIG = {
    'laminarin': {
        'terms': ['laminarin'],
        'curated_only': True,
    },
    'agar': {
        'terms': ['agarose', 'agarobiose'],
        'curated_only': False,
    },
    'carrageenan': {
        'terms': ['carrageenan'],
        'curated_only': True,
    },
    'alginate': {
        'terms': ['alginate'],
        'curated_only': True,
    },
    'fucoidan': {
        'terms': ['fucoidan', 'sulfated fucan'],
        'curated_only': False,
    },
    'ulvan': {
        'terms': ['ulvan'],
        'curated_only': True,
    },
    'porphyran': {
        'terms': ['porphyran'],
        'curated_only': True,
    },
    'xylan': {
        'terms': ['xylan', 'xylan and xylooligosaccharide',
                  'glucuronoxylan', 'beta-1,3-xylan'],
        'curated_only': True,
    },
    'arabinoxylan': {
        'terms': ['arabinoxylan'],
        'curated_only': True,
    },
    'pectin': {
        'terms': ['pectin', 'rhamnogalacturonan',
                  'polygalacturonic'],
        'curated_only': True,
    },
    'chitin': {
        'terms': ['chitin'],
        'curated_only': True,
    },
    'cellulose': {
        'terms': ['cellulose'],
        'curated_only': True,
    },
    'starch': {
        'terms': ['starch', 'amylopectin', 'amylose',
                  'granular starch'],
        'curated_only': True,
    },
    'beta_mannan': {
        'terms': ['mannan', 'glucomannan',
                  'mannooligosaccharide'],
        'curated_only': True,
    },
    'lichenan': {
        'terms': ['lichenin', 'lichenan',
                  'beta-1,3-1,4-glucan',
                  'mixed-linkage glucan'],
        'curated_only': True,
    },
    'mixed_linkage_glucan': {
        'terms': ['beta-1,3-1,4', 'mixed-linkage glucan',
                  'mixed linked glucan'],
        'curated_only': True,
    },
    'inulin': {
        'terms': ['inulin'],
        'curated_only': True,
    },
    'levan': {
        'terms': ['levan'],
        'curated_only': True,
    },
    'glycogen': {
        'terms': ['glycogen'],
        'curated_only': True,
    },
    'pullulan': {
        'terms': ['pullulan'],
        'curated_only': True,
    },
    'chondroitin_sulfate': {
        'terms': ['chondroitin sulfate', 'dermatan sulfate',
                  'glycosaminoglycan, hyaluronan, chondroitin',
                  'sulfated glycosaminoglycans'],
        'curated_only': True,
    },
    'heparan_sulfate': {
        'terms': ['heparan sulfate', 'heparin and heparan',
                  'heparin sulfate'],
        'curated_only': True,
    },
    'hyaluronic_acid': {
        'terms': ['hyaluronate', 'hyaluronan'],
        'curated_only': True,
    },
    'sucrose': {
        'terms': ['sucrose'],
        'curated_only': True,
    },
    'xyloglucan': {
        'terms': ['xyloglucan'],
        'curated_only': True,
    },
    'arabinogalactan': {
        'terms': ['arabinogalactan'],
        'curated_only': True,
    },
}

# Maximum length for Name strings to tokenise.
# Verbose CBM descriptions are long sentences containing organism names,
# journal references and other noise — skip tokenisation for these.
# Short Name strings (enzyme names) are always safe to tokenise.
MAX_NAME_LENGTH = 60

GENERIC_TOKENS = {
    # General English
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
    # Chemistry/biochemistry generics
    'acid', 'alpha', 'beta', 'bond', 'bonds', 'linked', 'linkage',
    'bind', 'binds', 'affinity', 'catalytic', 'reducing', 'oxidase',
    'monooxygenase', 'hydroxylating', 'dehydrogenating', 'lytic',
    'endo', 'exo', 'glycoside', 'glycosidase', 'hydrolase',
    'transferase', 'lyase', 'isomerase', 'ligase', 'phosphorylase',
    'deacetylase', 'esterase', 'oxidoreductase', 'mutase',
    'synthase', 'synthetase',
    # Structural descriptors
    'soluble', 'insoluble', 'crystalline', 'microcrystalline',
    'amorphous', 'granular', 'mixed', 'spectrum', 'range',
    'polysaccharide', 'polysaccharides', 'oligosaccharide',
    'oligosaccharides', 'disaccharide', 'monosaccharide',
    # Organism/source generics
    'plant', 'algae', 'fungi', 'fungal', 'bacterial', 'bacterium',
    'microbial', 'marine', 'soil', 'rumen', 'human', 'animal',
    'yeast', 'viral',
    # Literature noise
    'journal', 'world', 'microbiology', 'biotechnology', 'nature',
    'created', 'after', 'report', 'accord', 'pmid=26765840',
    'pmid=30635385',
    # Common words that pass length filter
    'number', 'members', 'proteins', 'enzymes', 'families',
    'tandem', 'repeats', 'conserved', 'appear', 'exhibit',
    'rather', 'between', 'sometimes', 'numerous',
}

MIN_TOKEN_LEN = 5


def load_mapping(fam_sub_map):
    """Load fam-substrate-mapping.tsv, stripping whitespace from headers."""
    df = pd.read_csv(fam_sub_map, sep='\t', dtype=str)
    df.columns = [c.strip() for c in df.columns]
    df = df.fillna('')
    return df


def find_matching_rows(df, search_terms, curated_only=False):
    """Find rows matching any search term."""
    pattern = '|'.join(re.escape(t) for t in search_terms)
    if curated_only:
        mask = df['Substrate_curated'].str.contains(
            pattern, case=False, na=False)
    else:
        mask = (
            df['Substrate_high_level'].str.contains(
                pattern, case=False, na=False) |
            df['Substrate_curated'].str.contains(
                pattern, case=False, na=False)
        )
    return df[mask]


def extract_tokens(name_string):
    """
    Extract meaningful tokens from a Name string.

    Skips strings longer than MAX_NAME_LENGTH to avoid tokenising
    verbose CBM descriptions that contain organism names and
    literature references.
    """
    name_string = str(name_string).strip()

    # Skip verbose descriptions — they contain too much noise
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


def derive_patterns_for_substrate(substrate, config, df):
    """Derive activity patterns for a substrate from matching rows."""
    search_terms = config['terms']
    curated_only = config.get('curated_only', False)

    hits = find_matching_rows(df, search_terms, curated_only)

    if hits.empty:
        print(f"  WARNING: No rows found for '{substrate}' "
              f"(terms: {search_terms})")
        return []

    print(f"  {substrate:<25} {len(hits):>3} rows matched "
          f"{'[curated only]' if curated_only else '[all columns]'}")

    all_tokens = set()

    # Always include search terms as patterns
    for term in search_terms:
        for part in re.split(r'[\s\-]+', term.lower()):
            if len(part) >= MIN_TOKEN_LEN:
                all_tokens.add(part)

    # Extract tokens from Name column (short entries only)
    for name in hits['Name'].dropna():
        all_tokens.update(extract_tokens(name))

    # Extract tokens from Substrate_curated (short entries only)
    for sub in hits['Substrate_curated'].dropna():
        if sub.strip():
            all_tokens.update(extract_tokens(sub))

    # Remove any tokens that are themselves substrate names of other
    # substrates — prevents cross-contamination
    # (applied after collection so search terms are protected)
    other_substrate_names = {
        s.replace('_', '') for s in SUBSTRATE_SEARCH_CONFIG
        if s != substrate
    }
    all_tokens = {
        t for t in all_tokens
        if t not in other_substrate_names
    }

    return sorted(all_tokens)


def generate_patterns(fam_sub_map, output_file):
    """Generate activity_patterns.tsv for all configured substrates."""
    print(f"Loading {fam_sub_map}...")
    df = load_mapping(fam_sub_map)
    print(f"Loaded {len(df)} rows, "
          f"{df['Substrate_high_level'].nunique()} high-level categories, "
          f"{df['Substrate_curated'].nunique()} curated substrates\n")

    rows = []
    print(f"Deriving patterns for "
          f"{len(SUBSTRATE_SEARCH_CONFIG)} substrates...\n")

    for substrate, config in sorted(SUBSTRATE_SEARCH_CONFIG.items()):
        patterns = derive_patterns_for_substrate(substrate, config, df)

        if not patterns:
            print(f"  SKIPPING {substrate} — no patterns derived")
            continue

        for pattern in patterns:
            rows.append({
                'substrate': substrate,
                'pattern':   pattern,
                'source':    'auto_derived',
                'reviewed':  False,
            })

    result_df = pd.DataFrame(rows)

    print(f"\nSummary:")
    print(f"  Total substrates: {result_df['substrate'].nunique()}")
    print(f"  Total patterns:   {len(result_df)}")
    print(f"\nPatterns per substrate:")
    for sub, count in result_df.groupby(
            'substrate').size().sort_values(
            ascending=False).items():
        print(f"  {sub:<25} {count:>3} patterns")

    # Report cross-substrate overlap
    print(f"\nCross-substrate pattern overlap:")
    substrates = result_df['substrate'].unique()
    for i, s1 in enumerate(sorted(substrates)):
        p1 = set(result_df[result_df['substrate'] == s1]['pattern'])
        for s2 in sorted(substrates)[i+1:]:
            p2 = set(result_df[result_df['substrate'] == s2]['pattern'])
            overlap = p1 & p2
            if len(overlap) > 5:
                print(f"  {s1} / {s2}: "
                      f"{len(overlap)} shared patterns")

    os.makedirs(os.path.dirname(os.path.abspath(output_file)),
                exist_ok=True)
    result_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nWritten to {output_file}")
    print(f"\nAll patterns written with reviewed=False.")
    print(f"Run the pipeline on your test dataset, then inspect")
    print(f"the output and set reviewed=True for correct patterns.")


def main():
    parser = argparse.ArgumentParser(
        description='Generate activity_patterns.tsv from dbCAN '
                    'fam-substrate-mapping.tsv'
    )
    parser.add_argument(
        '--fam_sub_map', required=True,
        help='Path to dbCAN fam-substrate-mapping.tsv'
    )
    parser.add_argument(
        '--output',
        default='substrate/data/activity_patterns.tsv',
        help='Output path for activity_patterns.tsv'
    )
    args = parser.parse_args()
    generate_patterns(args.fam_sub_map, args.output)


if __name__ == '__main__':
    main()
