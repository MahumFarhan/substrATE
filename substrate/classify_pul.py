"""
PUL classification and substrate family lookup.

Classifies CGCs as canonical_PUL, non_canonical_CGC, or outside_CGC
based on the presence of transporter genes and substrate-specific CAZymes.

Also provides substrate term and family lookups used to identify relevant
genes from dbCAN output.

NOTE: The SUBSTRATE_TERMS and FAMILY_MAP dicts are currently hardcoded
for the three substrates used in the original analysis. These will be
replaced by dynamic family derivation from dbCAN fam-substrate-mapping.tsv
at runtime (get_families.py logic) in a future update. To add a new
substrate in the meantime, add entries to both dicts here.
"""
import os
import pandas as pd


# ── Constants ─────────────────────────────────────────────────────────────────

# TCDB families considered SusC/SusD-type transporters in Bacteroidetes.
# 1.B.14 = TonB-dependent transporters (SusC-type)
# 8.A.46 = SusD-like outer membrane proteins
SUSC_FAMILIES = ['1.B.14', '8.A.46']

# Search terms for matching substrate_prediction.tsv columns.
# Each substrate maps to a list of terms used as case-insensitive
# substring matches against dbcan_pul_substrate and dbcansub_substrate.
SUBSTRATE_TERMS = {
    'laminarin': ['laminarin', 'beta-glucan'],
    'glycogen':  ['glycogen', 'alpha-glucan'],
    'pullulan':  ['pullulan', 'alpha-glucan'],
}

# CAZyme families per substrate, derived from dbCAN fam-substrate-mapping.tsv.
# Used to identify substrate-relevant genes in overview.tsv and to count
# substrate CAZymes per CGC for PUL classification.
FAMILY_MAP = {
    'agar':                ['CBM101', 'GH117', 'GH118', 'GH16', 'GH50',
                            'GH86', 'GH96'],
    'alginate':            ['CBM106', 'CBM96', 'PL14', 'PL15', 'PL17',
                            'PL18', 'PL20', 'PL31', 'PL32', 'PL34',
                            'PL36', 'PL38', 'PL39', 'PL41', 'PL5',
                            'PL6', 'PL7', 'PL8'],
    'arabinogalactan':     ['CBM62', 'GH115', 'GH127', 'GH154', 'GH16',
                            'GH27', 'GH30', 'GH35', 'GH42', 'GH43',
                            'GH5', 'GH53', 'GH62', 'GH78', 'GH79',
                            'PL42'],
    'arabinoxylan':        ['CBM42', 'GH10', 'GH3', 'GH42', 'GH43',
                            'GH5', 'GH51', 'GH54', 'GH62'],
    'beta_mannan':         ['CBM16', 'CBM23', 'CBM27', 'CBM29', 'CBM35',
                            'CBM59', 'CBM62', 'CBM72', 'CBM76', 'CBM80',
                            'CBM85', 'CBM88', 'CE17', 'CE2', 'GH113',
                            'GH125', 'GH130', 'GH134', 'GH148', 'GH2',
                            'GH26', 'GH27', 'GH36', 'GH45', 'GH5',
                            'GH76'],
    'carrageenan':         ['CBM92', 'GH127', 'GH129', 'GH150', 'GH16',
                            'GH167', 'GH82'],
    'cellulose':           ['AA10', 'AA15', 'AA16', 'AA18', 'AA3', 'AA9',
                            'CBM1', 'CBM10', 'CBM11', 'CBM16', 'CBM17',
                            'CBM2', 'CBM28', 'CBM3', 'CBM30', 'CBM37',
                            'CBM4', 'CBM44', 'CBM46', 'CBM49', 'CBM59',
                            'CBM6', 'CBM63', 'CBM64', 'CBM72', 'CBM8',
                            'CBM85', 'CBM9', 'GH10', 'GH12', 'GH124',
                            'GH44', 'GH45', 'GH48', 'GH5', 'GH51',
                            'GH6', 'GH7', 'GH74', 'GH8', 'GH9'],
    'chitin':              ['AA10', 'AA11', 'AA15', 'CBM1', 'CBM12',
                            'CBM14', 'CBM18', 'CBM19', 'CBM2', 'CBM3',
                            'CBM37', 'CBM5', 'CBM50', 'CBM54', 'CBM55',
                            'CBM73', 'CE14', 'CE4', 'GH16', 'GH179',
                            'GH18', 'GH19', 'GH2', 'GH20', 'GH23',
                            'GH3', 'GH35', 'GH48', 'GH5', 'GH94'],
    'chondroitin_sulfate': ['CBM100', 'GH39', 'GH56', 'GH88', 'PL6'],
    'fucoidan':            ['GH107', 'GH168', 'GH174', 'GH187', 'PL43'],
    'glycogen':            ['CBM48', 'GH126', 'GH13', 'GH133', 'GH14',
                            'GH31', 'GH57'],
    'heparan_sulfate':     ['GH39', 'GH79', 'PL12', 'PL13', 'PL21',
                            'PL37', 'PL8'],
    'hyaluronic_acid':     ['CBM70', 'GH16', 'GH56', 'GH79', 'GH84',
                            'PL16', 'PL30', 'PL33', 'PL8'],
    'inulin':              ['CBM38', 'GH32', 'GH91'],
    'laminarin':           ['CBM102', 'CBM103', 'GH131', 'GH16', 'GH17',
                            'GH3', 'GH5', 'GH55', 'GH70', 'GH8', 'GH9',
                            'GH94'],
    'levan':               ['GH32', 'GH68'],
    'lichenan':            ['CBM4', 'CBM6', 'GH12', 'GH131', 'GH148',
                            'GH16', 'GH17', 'GH3', 'GH44', 'GH45',
                            'GH5', 'GH7', 'GH70', 'GH8', 'GH9'],
    'pectin':              ['AA17', 'CBM41', 'CBM95', 'CBM97', 'CBM98',
                            'CE12', 'CE13', 'CE19', 'CE8', 'GH105',
                            'GH106', 'GH127', 'GH13', 'GH137', 'GH138',
                            'GH139', 'GH140', 'GH141', 'GH142', 'GH143',
                            'GH173', 'GH2', 'GH28', 'GH33', 'GH4',
                            'GH78', 'GH93', 'PL1', 'PL10', 'PL11',
                            'PL2', 'PL22', 'PL26', 'PL4', 'PL9'],
    'porphyran':           ['CBM99', 'GH16', 'GH50', 'GH86'],
    'pullulan':            ['CBM41', 'GH13', 'GH49', 'GH57'],
    'starch':              ['AA13', 'CBM20', 'CBM21', 'CBM25', 'CBM26',
                            'CBM34', 'CBM41', 'CBM45', 'CBM53', 'CBM69',
                            'CBM74', 'CBM82', 'CBM83', 'CBM98', 'GH119',
                            'GH122', 'GH126', 'GH13', 'GH14', 'GH15',
                            'GH31', 'GH49', 'GH57', 'GH70', 'GH97'],
    'sucrose':             ['GH100', 'GH13', 'GH31', 'GH32', 'GH4',
                            'GH68', 'GH70'],
    'ulvan':               ['CBM90', 'GH105', 'GH78'],
    'xylan':               ['AA10', 'AA14', 'CBM13', 'CBM15', 'CBM2',
                            'CBM22', 'CBM31', 'CBM35', 'CBM36', 'CBM37',
                            'CBM4', 'CBM42', 'CBM54', 'CBM59', 'CBM6',
                            'CBM60', 'CBM72', 'CBM85', 'CBM86', 'CBM89',
                            'CBM91', 'CE1', 'CE12', 'CE2', 'CE3', 'CE4',
                            'CE5', 'CE6', 'CE7', 'GH10', 'GH11',
                            'GH115', 'GH120', 'GH141', 'GH26', 'GH3',
                            'GH30', 'GH39', 'GH4', 'GH42', 'GH43',
                            'GH5', 'GH51', 'GH52', 'GH54', 'GH62',
                            'GH67', 'GH70', 'GH8', 'GH95', 'GH98'],
    'xyloglucan':          ['CBM44', 'CBM62', 'CBM65', 'CBM75', 'CBM76',
                            'CBM80', 'CBM81', 'CBM88', 'CE3', 'GH12',
                            'GH16', 'GH3', 'GH31', 'GH44', 'GH45',
                            'GH5', 'GH74', 'GH9'],
}


# ── CGC classification ────────────────────────────────────────────────────────

def classify_cgc(cgc_id, cgc_df, substrate_families,
                 pul_mode='bacteroidetes', min_cazymes=2):
    """
    Classify a CGC as canonical_PUL, non_canonical_CGC, or outside_CGC.

    Classification rules by pul_mode:

      bacteroidetes: SusC/SusD (TCDB 1.B.14 or 8.A.46) present
                     AND >= min_cazymes substrate CAZymes
                     -> canonical_PUL
                     No SusC/SusD BUT >= min_cazymes substrate CAZymes
                     -> non_canonical_CGC
                     Otherwise -> outside_CGC

      generic:       Any TC gene present AND >= min_cazymes substrate CAZymes
                     -> canonical_PUL
                     No TC gene BUT >= min_cazymes substrate CAZymes
                     -> non_canonical_CGC
                     Otherwise -> outside_CGC

      cazyme_only:   >= min_cazymes substrate CAZymes -> canonical_PUL
                     Otherwise -> outside_CGC
                     (non_canonical_CGC is not used in this mode)

    Args:
        cgc_id:            CGC identifier string (e.g. 'CGC44')
        cgc_df:            DataFrame from cgc_standard_out.tsv for this sample
        substrate_families: list of CAZyme family strings for the substrate
        pul_mode:          'bacteroidetes', 'generic', or 'cazyme_only'
        min_cazymes:       minimum substrate CAZymes required per CGC

    Returns:
        'canonical_PUL', 'non_canonical_CGC', or 'outside_CGC'
    """
    genes = cgc_df[cgc_df['CGC#'] == cgc_id].copy()
    genes['Gene Annotation'] = genes['Gene Annotation'].fillna('')

    all_annots = ' '.join(genes['Gene Annotation'].astype(str))

    # Count substrate-specific CAZymes
    sub_cazyme_count = sum(
        any(f in str(row['Gene Annotation']) for f in substrate_families)
        for _, row in genes.iterrows()
        if str(row['Gene Type']) == 'CAZyme'
    )

    if pul_mode == 'cazyme_only':
        if sub_cazyme_count >= min_cazymes:
            return 'canonical_PUL'
        return 'outside_CGC'

    elif pul_mode == 'generic':
        has_transporter = any(
            str(row['Gene Type']) == 'TC'
            for _, row in genes.iterrows()
        )

    else:  # bacteroidetes (default)
        has_transporter = any(tc in all_annots for tc in SUSC_FAMILIES)

    if has_transporter and sub_cazyme_count >= min_cazymes:
        return 'canonical_PUL'
    elif not has_transporter and sub_cazyme_count >= min_cazymes:
        return 'non_canonical_CGC'
    else:
        return 'outside_CGC'


# ── Sample processing ─────────────────────────────────────────────────────────

def process_samples(cgc_output_dir, substrate, pul_mode='bacteroidetes',
                    min_cazymes=2):
    """
    Iterate over all sample directories in cgc_output_dir and build
    substrate hit and family hit DataFrames.

    Args:
        cgc_output_dir: path to directory containing per-sample dbCAN output
        substrate:      substrate name (must be a key in SUBSTRATE_TERMS
                        and FAMILY_MAP)
        pul_mode:       PUL classification mode (see classify_cgc)
        min_cazymes:    minimum substrate CAZymes required per CGC

    Returns:
        tuple of (substrate_hits_df, family_hits_df, overview_df)
        Any of these may be an empty DataFrame if no hits were found.
    """
    if substrate not in FAMILY_MAP:
        raise ValueError(
            f"Unknown substrate '{substrate}'. "
            f"Known substrates: {sorted(FAMILY_MAP.keys())}"
        )

    families     = FAMILY_MAP[substrate]
    terms        = SUBSTRATE_TERMS[substrate]
    term_pattern = '|'.join(terms)

    all_substrate_hits = []
    all_family_hits    = []
    all_overview       = []

    for sample_dir in sorted(os.listdir(cgc_output_dir)):
        full_path = os.path.join(cgc_output_dir, sample_dir)
        if not os.path.isdir(full_path):
            continue

        sample    = sample_dir.replace('output_', '')
        sub_file  = os.path.join(full_path, 'substrate_prediction.tsv')
        cgc_file  = os.path.join(full_path, 'cgc_standard_out.tsv')
        over_file = os.path.join(full_path, 'overview.tsv')

        # ── Substrate-level hits ──────────────────────────────────────────────
        if os.path.exists(sub_file):
            sub_df = pd.read_csv(
                sub_file, sep='\t', comment='#',
                names=['cgcid', 'pulid', 'dbcan_pul_substrate',
                       'bitscore', 'sig_pairs', 'dbcansub_substrate',
                       'dbcansub_score']
            )
            sub_df['dbcan_pul_substrate'] = sub_df['dbcan_pul_substrate'].astype(str)
            sub_df['dbcansub_substrate']  = sub_df['dbcansub_substrate'].astype(str)

            mask = (
                sub_df['dbcan_pul_substrate'].str.contains(
                    term_pattern, case=False, na=False) |
                sub_df['dbcansub_substrate'].str.contains(
                    term_pattern, case=False, na=False)
            )
            hits = sub_df[mask].copy()
            if not hits.empty:
                hits['sample']              = sample
                hits['substrate_category']  = substrate
                all_substrate_hits.append(hits)
        else:
            print(f"WARNING: No substrate file for {sample}")

        # ── Load CGC gene data ────────────────────────────────────────────────
        cgc_df           = pd.DataFrame()
        cgc_gene_ids     = set()
        cgc_gene_cgc_map = {}

        if os.path.exists(cgc_file):
            cgc_df = pd.read_csv(cgc_file, sep='\t')
            cgc_df['Gene Annotation'] = cgc_df['Gene Annotation'].fillna('')
            if 'Protein ID' in cgc_df.columns and 'CGC#' in cgc_df.columns:
                cgc_gene_ids     = set(cgc_df['Protein ID'].astype(str))
                cgc_gene_cgc_map = dict(zip(
                    cgc_df['Protein ID'].astype(str),
                    cgc_df['CGC#'].astype(str)
                ))

        # ── Family-level hits from overview.tsv ───────────────────────────────
        if os.path.exists(over_file):
            over_df = pd.read_csv(over_file, sep='\t')
            over_df.columns = over_df.columns.str.strip()
            over_df['sample'] = sample

            over_df['subfamily_annotation'] = over_df['Recommend Results'].astype(str)
            over_df.loc[
                over_df['subfamily_annotation'].isin(['', 'nan', '-']),
                'subfamily_annotation'
            ] = over_df['dbCAN_hmm'].astype(str)

            over_df['all_annotations'] = (
                over_df['dbCAN_hmm'].astype(str) + '|' +
                over_df['dbCAN_sub'].astype(str) + '|' +
                over_df['DIAMOND'].astype(str)   + '|' +
                over_df['Recommend Results'].astype(str)
            )

            gene_col          = 'Gene ID'
            over_df['CGC_id'] = over_df[gene_col].astype(str).map(
                cgc_gene_cgc_map).fillna('-')
            over_df['in_any_cgc'] = over_df[gene_col].astype(str).isin(
                cgc_gene_ids)

            all_overview.append(over_df)

            # Build per-CGC classifications for this substrate
            cgc_classifications = {}
            if not cgc_df.empty:
                for cgc_id in cgc_df['CGC#'].unique():
                    cgc_classifications[cgc_id] = classify_cgc(
                        cgc_id, cgc_df, families,
                        pul_mode=pul_mode, min_cazymes=min_cazymes
                    )

            # Match genes to substrate families and assign localisation
            for family in families:
                mask = over_df['all_annotations'].str.contains(
                    rf'(?<![0-9]){family}(?![0-9])',
                    case=False, na=False, regex=True
                )
                fam_hits = over_df[mask].copy()
                if fam_hits.empty:
                    continue

                fam_hits['substrate_category'] = substrate
                fam_hits['matched_family']      = family

                # Option B: pandas map() for localisation lookup
                fam_hits['localisation'] = (
                    fam_hits['CGC_id']
                    .map(cgc_classifications)
                    .fillna('outside_CGC')
                )

                all_family_hits.append(fam_hits)
        else:
            print(f"WARNING: No overview file for {sample}")

    # ── Assemble final DataFrames ─────────────────────────────────────────────
    substrate_hits_df = (pd.concat(all_substrate_hits, ignore_index=True)
                         if all_substrate_hits else pd.DataFrame())
    family_hits_df    = (pd.concat(all_family_hits, ignore_index=True)
                         if all_family_hits else pd.DataFrame())
    overview_df       = (pd.concat(all_overview, ignore_index=True)
                         if all_overview else pd.DataFrame())

    return substrate_hits_df, family_hits_df, overview_df
