"""
count.py
========
Activity-filtered GH/CBM family counting from substrATE
activity_annotated.tsv output files.

For each substrate and GH/CBM family of interest, three counts are
produced per strain:

    <FAMILY>_raw                   : all hits for that family in the
                                     activity_annotated file
    <FAMILY>_<substrate>_filtered  : hits where the assigned activity
                                     matches the substrate-specific
                                     activity whitelist

The whitelist is defined in ACTIVITY_WHITELIST below. It is keyed by
substrate name (matching the keys in classify_pul.FAMILY_MAP) and then
by family name. Each entry is a frozenset of lowercase activity strings
used as case-insensitive substring matches against the activity column
of the annotated file.

To add a new substrate or family:
    1. Add an entry to ACTIVITY_WHITELIST below.
    2. No other changes are needed — the CLI and counting logic are
       fully driven by this dictionary.

Notes on whitelist design
-------------------------
Substring matching is intentionally permissive: an activity string
"glucan endo-1,3-beta-d-glucosidase" will match whitelist term
"glucan endo-1,3-beta-d-glucosidase" exactly, but a shorter term
such as "1,3-beta-d-glucosidase" would also match it. Keep whitelist
terms specific enough to avoid false positives.

Hits with activity "-" or empty strings are never counted as filtered.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


# =============================================================================
# ACTIVITY WHITELISTS
# =============================================================================

ACTIVITY_WHITELIST: dict[str, dict[str, frozenset[str]]] = {

    "laminarin": {
        "GH16": frozenset({
            "glucan endo-1,3-beta-d-glucosidase",
            "glucan endo-1,3-beta-glucosidase",
            "endo-beta-1,3-glucanase",
            "endo-beta-1,3-glucanase / laminarinase",
            "laminarinase",
            "laminarin-degrading enzyme",
            "beta-1,3-glucanase",
            "endo-1,3-beta-glucanase",
            "endo-1,3(4)-beta-glucanase",
            "glucan 1,3-betalpha-glucosidase",
        }),
        "GH17": frozenset({
            "glucan endo-1,3-beta-d-glucosidase",
            "glucan endo-1,3-beta-glucosidase",
            "beta-1,3-glucanosyltransglycosylase",
            "endo-beta-1,3-glucanase",
            "laminarinase",
            "glucan 1,3-beta-glucosidase",
            "glucan 1,3-betalpha-glucosidase",
        }),
        "GH3": frozenset({
            "glucan 1,3-beta-glucosidase",
            "glucan 1,3-betalpha-glucosidase",
            "beta-1,3-glucosidase",
            "endo-1,3(4)-beta-glucanase",
            "glucan endo-1,3-beta-d-glucosidase",
            "glucan endo-1,3-beta-glucosidase",
        }),
        "GH30": frozenset({
            "glucan endo-1,3-beta-d-glucosidase",
            "endo-beta-1,3-glucanase",
            "laminarinase",
            "beta-1,3-glucanase",
        }),
    },

    "glycogen": {
        "GH13": frozenset({
            "alpha-amylase",
            "alpha-glucoside phosphorylase",
            "alpha glucoside phosphorylase",
            "alphalpha-glucosidase",
            "1,4-alpha-glucan branching enzyme",
            "pullulanase",
            "isoamylase",
            "cyclomaltodextrinase",
            "cyclomaltodextrin glucanotransferase",
            "oligo-1,6-glucosidase",
            "maltose alpha-d-glucosyltransferase",
            "amylosucrase",
            "sucrose phosphorylase",
            "isomaltulose synthase",
            "glucan 1,4-alpha-maltohydrolase",
            "glucan 1,4-alpha-maltotetraohydrolase",
            "glucan 1,4-alpha-maltotriohydrolase",
            "glucan 1,4-alpha-maltohexaosidase",
            "glucan 1,6-alphalpha-glucosidase",
            "glucan 1,4-alphalpha-glucosidase",
            "4-alpha-glucanotransferase",
            "neopullulanase",
            "maltogenic amylase",
            "trehalose synthase",
            "glycogen phosphorylase",
            "starch synthase (maltosyl-transferring)",
            "4-alpha-d-{(1->4)-alpha-d-glucano}trehalose trehalohydrolase",
            "(1->4)-alpha-d-glucan 1-alpha-d-glucosylmutase",
        }),
        "GH31": frozenset({
            "alpha-glucosidase",
            "alphalpha-glucosidase",
            "alpha-1,4-glucosidase",
            "sucrase-isomaltase",
            "sucrase",
            "isomaltase",
            "maltase",
            "oligo-1,6-glucosidase",
        }),
        "GH65": frozenset({
            "maltose phosphorylase",
            "trehalose phosphorylase",
            "kojibiose phosphorylase",
            "alpha,alpha-trehalose phosphorylase",
            "maltose-phosphorylase",
            "a-glucan phosphorylase",
        }),
        "GH97": frozenset({
            "alpha-glucosidase",
            "alphalpha-glucosidase",
            "glucan 1,4-alphalpha-glucosidase",
            "glucan 1,6-alphalpha-glucosidase",
        }),
        "CBM48": frozenset({
            "carbohydrate-binding module (cbm48)",
            "1,4-alpha-glucan branching enzyme",
            "pullulanase",
        }),
    },
}


# =============================================================================
# INTERNAL HELPERS
# =============================================================================

def _load_annotated(path: Path) -> pd.DataFrame:
    """
    Load an activity_annotated TSV produced by substrATE and normalise
    column names to lowercase with underscores.
    """
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]

    required = {"matched_family", "activity", "sample"}
    missing  = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Missing required columns in {path.name}: {missing}\n"
            f"Available columns: {list(df.columns)}"
        )

    df["activity_lower"] = df["activity"].fillna("").str.lower().str.strip()
    df["matched_family"] = df["matched_family"].fillna("").str.strip()
    df["sample"]         = df["sample"].fillna("").str.strip()
    return df


def _activity_passes(activity_lower: str, whitelist: frozenset) -> bool:
    """
    Return True if activity_lower contains any term from the whitelist.
    Empty or '-' activities always return False.
    """
    if not activity_lower or activity_lower == "-":
        return False
    return any(term in activity_lower for term in whitelist)


def _count_one_family(
    df: pd.DataFrame,
    family: str,
    whitelist: frozenset,
    substrate_label: str,
    all_samples,
) -> pd.DataFrame:
    """
    Count raw and filtered hits for a single family across all samples.
    """
    raw_col = f"{family}_raw"
    fil_col = f"{family}_{substrate_label}_filtered"

    fam_df = df[df["matched_family"] == family].copy()

    raw_counts = (
        fam_df.groupby("sample")
              .size()
              .reindex(all_samples, fill_value=0)
              .reset_index()
    )
    raw_counts.columns = ["sample", raw_col]

    fam_df["_passes"] = fam_df["activity_lower"].apply(
        lambda a: _activity_passes(a, whitelist)
    )
    passing = fam_df[fam_df["_passes"]]
    if passing.empty:
        fil_counts = pd.DataFrame({
            "sample": all_samples,
            fil_col:  0,
        })
    else:
        fil_counts = (
            passing.groupby("sample")
                   .size()
                   .reindex(all_samples, fill_value=0)
                   .reset_index()
        )
        fil_counts.columns = ["sample", fil_col]

    return raw_counts.merge(fil_counts, on="sample")


# =============================================================================
# PUBLIC API
# =============================================================================

def count_substrate(
    annotated_path: Path,
    substrate_label: str,
    families=None,
) -> pd.DataFrame:
    """
    Count raw and activity-filtered hits for all families of one substrate.

    Parameters
    ----------
    annotated_path  : path to <substrate>_activity_annotated.tsv
    substrate_label : substrate name e.g. 'laminarin'
    families        : optional override dict {family: whitelist_set}

    Returns
    -------
    Wide DataFrame: one row per strain, columns for each family's
    raw and filtered counts.
    """
    if families is None:
        if substrate_label not in ACTIVITY_WHITELIST:
            raise ValueError(
                f"No whitelist defined for substrate '{substrate_label}'. "
                f"Known substrates: {sorted(ACTIVITY_WHITELIST.keys())}. "
                f"Pass a custom families dict or add an entry to "
                f"ACTIVITY_WHITELIST in substrate/count.py."
            )
        families = ACTIVITY_WHITELIST[substrate_label]

    logger.info("Loading %s", annotated_path.name)
    df          = _load_annotated(annotated_path)
    all_samples = df["sample"].unique()

    result = pd.DataFrame({"sample": all_samples})

    for family, whitelist in families.items():
        logger.info("  Counting %s (%s) ...", family, substrate_label)
        counts = _count_one_family(
            df, family, whitelist, substrate_label, all_samples
        )
        result = result.merge(counts, on="sample", how="left")

    result = result.fillna(0)
    int_cols = [c for c in result.columns if c != "sample"]
    result[int_cols] = result[int_cols].astype(int)
    result = result.sort_values("sample").reset_index(drop=True)

    for family in families:
        raw_col = f"{family}_raw"
        fil_col = f"{family}_{substrate_label}_filtered"
        if raw_col in result.columns:
            total_raw = result[raw_col].sum()
            total_fil = result[fil_col].sum() if fil_col in result.columns else 0
            removed   = total_raw - total_fil
            pct       = removed / total_raw * 100 if total_raw > 0 else 0.0
            logger.info(
                "    %-8s  raw=%4d  filtered=%4d  removed=%4d  (%.1f%%)",
                family, total_raw, total_fil, removed, pct,
            )

    return result


def merge_substrates(substrate_dfs: list) -> pd.DataFrame:
    """
    Outer-join per-substrate DataFrames on the sample column.
    Strains absent from one substrate get zeros for that substrate's columns.
    """
    if not substrate_dfs:
        raise ValueError("No substrate DataFrames to merge.")

    merged = substrate_dfs[0]
    for df in substrate_dfs[1:]:
        new_cols = [c for c in df.columns
                    if c == "sample" or c not in merged.columns]
        merged = merged.merge(df[new_cols], on="sample", how="outer")

    merged = merged.fillna(0)
    int_cols = [c for c in merged.columns if c != "sample"]
    merged[int_cols] = merged[int_cols].astype(int)

    def _sort_key(col):
        parts = col.split("_")
        family = parts[0]
        is_raw = 1 if "raw" in col else 2
        return (family, is_raw, col)

    non_sample = sorted(
        [c for c in merged.columns if c != "sample"],
        key=_sort_key,
    )
    return merged[["sample"] + non_sample].sort_values("sample").reset_index(drop=True)


def audit_whitelists(annotated_paths: dict) -> pd.DataFrame:
    """
    For each substrate/family, list every observed activity with its
    count and whether it is currently whitelisted.

    Parameters
    ----------
    annotated_paths : dict mapping substrate label to annotated TSV path
                      e.g. {'laminarin': Path('laminarin_activity_annotated.tsv')}

    Returns
    -------
    DataFrame with columns:
        substrate | family | activity | count | whitelisted
    """
    rows = []
    for substrate_label, path in annotated_paths.items():
        df       = _load_annotated(path)
        wl_entry = ACTIVITY_WHITELIST.get(substrate_label, {})

        for family in sorted(df["matched_family"].unique()):
            if not family:
                continue
            fam_df    = df[df["matched_family"] == family]
            whitelist = wl_entry.get(family, frozenset())
            for activity, count in fam_df["activity"].value_counts().items():
                act_lower   = str(activity).lower().strip()
                whitelisted = _activity_passes(act_lower, whitelist)
                rows.append({
                    "substrate":   substrate_label,
                    "family":      family,
                    "activity":    activity,
                    "count":       count,
                    "whitelisted": whitelisted,
                })

    return pd.DataFrame(rows)
