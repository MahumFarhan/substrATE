#!/usr/bin/env python3
"""
Rebuilds gh_counts_by_strain_new.tsv from current substrATE output.

For each strain, computes two kinds of count per GH/CBM family:

  {family}_raw
      Every dbCAN hit for that family, read directly from each sample's
      own overview.tsv (in cgc_output_dir), matched with the SAME
      boundary-regex logic classify_pul.py uses internally
      (rf'(?<![0-9]){family}(?![0-9])' against the concatenation of
      dbCAN_hmm | dbCAN_sub | DIAMOND | Recommend Results). This is
      NOT restricted to any substrate's FAMILY_MAP -- it reflects every
      dbCAN call for that family, full stop.

  {family}_{substrate}_filtered
      Only the subset of {family}_raw hits whose PREDICTED ACTIVITY
      (the 'activity' column in {substrate}_activity_annotated.tsv,
      derived from EC-number lookup) matches one of that substrate's
      curated activity patterns -- e.g. distinguishing a genuine
      laminarinase GH16 from an agarase GH16. This is computed by
      calling the REAL substrate.parse_substrates.load_patterns()
      function (the same one `reduced-tree` uses), not a re-derived
      approximation, so the definition here is guaranteed identical to
      what the rest of the pipeline means by "activity-filtered".

KNOWN GAP: GH30 is not in classify_pul.py's FAMILY_MAP['laminarin'],
so GH30_laminarin_filtered will be 0 for every strain -- substrATE
never tags a gene as matched_family=='GH30' for the laminarin
substrate at all. This script prints an explicit warning about this;
it is not a bug in this script, it's a gap in substrATE's own family
list for laminarin (worth deciding separately whether to add GH30 to
FAMILY_MAP['laminarin'] in classify_pul.py, if you want a real,
non-zero filtered count for it).

Run this with the same conda env / from a location where
`from substrate import ...` resolves (i.e. substrATE installed via
`pip install -e .`, or run from the substrATE repo root).
"""

import os
import re
import sys
import pandas as pd

# ─────────────────────────────────────────────────────────────────────────
# CONFIG -- edit these to match your current run
# ─────────────────────────────────────────────────────────────────────────

# Where dbCAN's own per-sample output lives (contains output_{sample}/
# directories, each with its own overview.tsv). Built by `substrate run`
# when you passed --genomes (not --dbcan_output) -- should be
# {output}/cgc_output/ under whatever --output you used.
CGC_OUTPUT_DIR = os.path.expanduser(
    "~/paper1_fasta/Kappelmann/substrate_v4_7Sep/cgc_output"
)

# Where the substrate-specific {substrate}_activity_annotated.tsv files
# live -- the same --output directory as above, one level up from
# cgc_output/ (each substrate gets its own subdirectory here).
SUBSTRATE_OUTPUT_DIR = os.path.expanduser(
    "~/paper1_fasta/Kappelmann/substrate_v4_7Sep"
)

# Path to activity_patterns.tsv -- matches _PATTERNS_FILE in cli.py.
# Adjust if your substrATE install lives somewhere else.
PATTERNS_FILE = os.path.expanduser(
    "~/substrATE/substrate/data/activity_patterns.tsv"
)

# Pattern mode: must match whatever the actual `substrate run` invocation
# used. The Kappelmann v4 rerun command did NOT pass --pattern_mode, so
# it used cli.py's default, which is 'permissive'. Change to 'strict' if
# you specifically ran with --pattern_mode strict.
PATTERN_MODE = "permissive"

# Family -> substrate mapping, matching the original comparison script's
# block_def exactly.
LAMINARIN_FAMILIES    = ["GH3", "GH16", "GH17", "GH30"]
ALPHA_GLUCAN_FAMILIES = ["GH13", "GH31", "GH65", "GH97", "CBM48"]
ALL_FAMILIES = sorted(set(LAMINARIN_FAMILIES + ALPHA_GLUCAN_FAMILIES))

OUTPUT_TSV = os.path.join(SUBSTRATE_OUTPUT_DIR, "gh_counts_by_strain_new.tsv")


# ─────────────────────────────────────────────────────────────────────────
# Make the substrATE package importable
# ─────────────────────────────────────────────────────────────────────────

try:
    from substrate import parse_substrates
except ImportError:
    print("ERROR: could not import 'substrate' package.")
    print("Run this from the substrATE repo root, or ensure it's installed")
    print("via 'pip install -e .' in the substrATE conda environment.")
    sys.exit(1)


# ─────────────────────────────────────────────────────────────────────────
# Step 1: RAW counts, from each sample's own overview.tsv
# ─────────────────────────────────────────────────────────────────────────

def count_raw_family_hits(cgc_output_dir, families):
    """
    For every sample directory in cgc_output_dir, count dbCAN hits per
    family using the same boundary-regex matching classify_pul.py uses
    internally, applied to overview.tsv directly (no FAMILY_MAP gate).
    """
    if not os.path.isdir(cgc_output_dir):
        print(f"ERROR: cgc_output_dir does not exist: {cgc_output_dir}")
        sys.exit(1)

    sample_dirs = sorted(
        d for d in os.listdir(cgc_output_dir)
        if d.startswith("output_")
        and os.path.isdir(os.path.join(cgc_output_dir, d))
    )
    print(f"Found {len(sample_dirs)} sample directories in {cgc_output_dir}")

    family_patterns = {
        fam: re.compile(rf"(?<![0-9]){re.escape(fam)}(?![0-9])")
        for fam in families
    }

    rows = []
    missing_overview = []
    for sdir in sample_dirs:
        sample = sdir[len("output_"):]
        over_path = os.path.join(cgc_output_dir, sdir, "overview.tsv")
        if not os.path.exists(over_path):
            missing_overview.append(sample)
            continue

        over_df = pd.read_csv(over_path, sep="\t")
        over_df.columns = [c.strip() for c in over_df.columns]

        for col in ["dbCAN_hmm", "dbCAN_sub", "DIAMOND", "Recommend Results"]:
            if col not in over_df.columns:
                over_df[col] = ""
            # From pandas 3, astype(str) leaves missing values as NaN, and
            # one NaN would blank the whole combined string for that gene.
            over_df[col] = over_df[col].fillna("-")

        all_annot = (
            over_df["dbCAN_hmm"].astype(str) + "|" +
            over_df["dbCAN_sub"].astype(str) + "|" +
            over_df["DIAMOND"].astype(str) + "|" +
            over_df["Recommend Results"].astype(str)
        )

        row = {"sample": sample}
        for fam in families:
            row[f"{fam}_raw"] = int(
                all_annot.str.contains(family_patterns[fam]).sum()
            )
        rows.append(row)

    if missing_overview:
        print(f"  WARNING: {len(missing_overview)} sample(s) had no "
              f"overview.tsv (excluded from raw counts):")
        for s in missing_overview:
            print(f"    - {s}")

    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────
# Step 2: FILTERED counts, from {substrate}_activity_annotated.tsv,
# using the REAL parse_substrates.load_patterns() for matching.
# ─────────────────────────────────────────────────────────────────────────

def count_filtered_family_hits(substrate_output_dir, substrate,
                                families, patterns_file, pattern_mode):
    """
    Counts, per sample, hits for `families` within `substrate`'s
    activity_annotated.tsv whose activity matches one of the substrate's
    curated activity patterns (loaded via the real pipeline function).
    """
    activity_path = os.path.join(
        substrate_output_dir, substrate, f"{substrate}_activity_annotated.tsv"
    )
    if not os.path.exists(activity_path):
        print(f"ERROR: activity file not found: {activity_path}")
        sys.exit(1)

    hits = pd.read_csv(activity_path, sep="\t")

    patterns_df = parse_substrates.load_patterns(
        substrate=substrate,
        patterns_file=patterns_file,
        pattern_mode=pattern_mode,
    )
    if patterns_df.empty:
        print(f"  WARNING: no activity patterns found for substrate "
              f"'{substrate}' (mode={pattern_mode}) -- all filtered "
              f"counts for this substrate will be 0.")
        patterns = []
    else:
        patterns = patterns_df["pattern"].str.lower().tolist()

    present_families = set(hits["matched_family"].unique())
    missing_families = sorted(set(families) - present_families)
    if missing_families:
        print(f"  NOTE: family/families with NO rows at all in "
              f"{substrate}_activity_annotated.tsv (not in substrATE's "
              f"FAMILY_MAP for '{substrate}' -- filtered counts will be "
              f"0 for every strain): {', '.join(missing_families)}")

    def activity_matches(activity_text):
        if pd.isna(activity_text) or not patterns:
            return False
        text = str(activity_text).lower()
        return any(p in text for p in patterns)

    hits = hits[hits["matched_family"].isin(families)].copy()
    hits["pattern_match"] = hits["activity"].apply(activity_matches)
    matched = hits[hits["pattern_match"]]

    # One row per (sample, Gene ID, matched_family) to avoid double-
    # counting if a gene somehow has multiple annotation rows.
    matched = matched.drop_duplicates(subset=["sample", "Gene ID", "matched_family"])

    counts = (
        matched.groupby(["sample", "matched_family"]).size()
        .unstack(fill_value=0)
        .reindex(columns=families, fill_value=0)
    )
    counts.columns = [f"{fam}_{substrate}_filtered" for fam in counts.columns]
    counts = counts.reset_index()
    return counts


# ─────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("Building gh_counts_by_strain_new.tsv")
    print("=" * 70)
    print(f"cgc_output_dir:       {CGC_OUTPUT_DIR}")
    print(f"substrate_output_dir: {SUBSTRATE_OUTPUT_DIR}")
    print(f"patterns_file:        {PATTERNS_FILE}")
    print(f"pattern_mode:         {PATTERN_MODE}")
    print()

    print("--- Step 1: raw counts (from overview.tsv) ---")
    raw_df = count_raw_family_hits(CGC_OUTPUT_DIR, ALL_FAMILIES)
    print(f"  {len(raw_df)} strains with raw counts computed.\n")

    print("--- Step 2a: laminarin filtered counts ---")
    lam_filtered_df = count_filtered_family_hits(
        SUBSTRATE_OUTPUT_DIR, "laminarin", LAMINARIN_FAMILIES,
        PATTERNS_FILE, PATTERN_MODE,
    )
    print(f"  {len(lam_filtered_df)} strains with laminarin-filtered hits.\n")

    print("--- Step 2b: glycogen filtered counts ---")
    gly_filtered_df = count_filtered_family_hits(
        SUBSTRATE_OUTPUT_DIR, "glycogen", ALPHA_GLUCAN_FAMILIES,
        PATTERNS_FILE, PATTERN_MODE,
    )
    print(f"  {len(gly_filtered_df)} strains with glycogen-filtered hits.\n")

    print("--- Step 3: merging ---")
    result = raw_df.merge(lam_filtered_df, on="sample", how="left")
    result = result.merge(gly_filtered_df, on="sample", how="left")

    # Strains with zero filtered hits won't appear in the filtered
    # dataframes at all (groupby drops them) -- fill those with 0,
    # NOT with NA, since "not in the filtered table" genuinely means
    # "zero hits passed the activity filter", not "unknown/not tested".
    filtered_cols = [c for c in result.columns if c.endswith("_filtered")]
    result[filtered_cols] = result[filtered_cols].fillna(0).astype(int)

    print(f"  Final table: {len(result)} strains x {len(result.columns)} columns\n")

    print("--- Sanity check: raw >= filtered for every family/strain ---")
    problems = 0
    for fam in LAMINARIN_FAMILIES:
        raw_col, fil_col = f"{fam}_raw", f"{fam}_laminarin_filtered"
        if raw_col in result.columns and fil_col in result.columns:
            bad = result[result[raw_col] < result[fil_col]]
            if len(bad) > 0:
                print(f"  WARNING: {fam} has {len(bad)} strain(s) where "
                      f"filtered > raw (should never happen):")
                print(bad[["sample", raw_col, fil_col]].to_string(index=False))
                problems += len(bad)
    for fam in ALPHA_GLUCAN_FAMILIES:
        raw_col, fil_col = f"{fam}_raw", f"{fam}_glycogen_filtered"
        if raw_col in result.columns and fil_col in result.columns:
            bad = result[result[raw_col] < result[fil_col]]
            if len(bad) > 0:
                print(f"  WARNING: {fam} has {len(bad)} strain(s) where "
                      f"filtered > raw (should never happen):")
                print(bad[["sample", raw_col, fil_col]].to_string(index=False))
                problems += len(bad)
    if problems == 0:
        print("  None found -- OK.\n")
    else:
        print(f"  {problems} total problem row(s) -- investigate before "
              f"trusting this TSV.\n")

    result.to_csv(OUTPUT_TSV, sep="\t", index=False)
    print(f"Written: {OUTPUT_TSV}")
    print("\nDone. Update df_new's path in the comparison script to point here.")


if __name__ == "__main__":
    main()
