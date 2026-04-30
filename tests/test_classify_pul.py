"""
Unit tests for substrate.classify_pul module.

Tests classify_cgc() across all three pul_mode values and
edge cases. Uses synthetic DataFrames — no database files needed.
"""
import pytest
import pandas as pd
from substrate.classify_pul import classify_cgc, SUSC_FAMILIES


# ── Fixtures ──────────────────────────────────────────────────────────────────

def make_cgc_df(genes):
    """
    Build a minimal cgc_standard_out.tsv DataFrame for testing.

    Args:
        genes: list of (cgc_id, gene_type, annotation) tuples

    Returns:
        DataFrame with columns CGC#, Gene Type, Gene Annotation,
        Protein ID
    """
    rows = []
    for i, (cgc_id, gene_type, annotation) in enumerate(genes):
        rows.append({
            'CGC#':            cgc_id,
            'Gene Type':       gene_type,
            'Gene Annotation': annotation,
            'Protein ID':      f'gene_{i}',
        })
    return pd.DataFrame(rows)


SUBSTRATE_FAMILIES = ['GH16', 'GH17', 'GH55']


# ── bacteroidetes mode ────────────────────────────────────────────────────────

class TestBacteroidetesMode:

    def test_canonical_pul_susc_plus_two_cazymes(self):
        """SusC + 2 substrate CAZymes -> canonical_PUL."""
        cgc_df = make_cgc_df([
            ('CGC1', 'CAZyme', 'GH16'),
            ('CGC1', 'CAZyme', 'GH17'),
            ('CGC1', 'TC',     '1.B.14'),
        ])
        result = classify_cgc(
            'CGC1', cgc_df, SUBSTRATE_FAMILIES,
            pul_mode='bacteroidetes', min_cazymes=2)
        assert result == 'canonical_PUL'

    def test_canonical_pul_susd_plus_cazymes(self):
        """SusD-type transporter (8.A.46) also gives canonical_PUL."""
        cgc_df = make_cgc_df([
            ('CGC1', 'CAZyme', 'GH16'),
            ('CGC1', 'CAZyme', 'GH55'),
            ('CGC1', 'TC',     '8.A.46'),
        ])
        result = classify_cgc(
            'CGC1', cgc_df, SUBSTRATE_FAMILIES,
            pul_mode='bacteroidetes', min_cazymes=2)
        assert result == 'canonical_PUL'

    def test_non_canonical_no_susc(self):
        """2 substrate CAZymes but no SusC/SusD -> non_canonical_CGC."""
        cgc_df = make_cgc_df([
            ('CGC1', 'CAZyme', 'GH16'),
            ('CGC1', 'CAZyme', 'GH17'),
            ('CGC1', 'TC',     '2.A.1'),  # Non-SusC transporter
        ])
        result = classify_cgc(
            'CGC1', cgc_df, SUBSTRATE_FAMILIES,
            pul_mode='bacteroidetes', min_cazymes=2)
        assert result == 'non_canonical_CGC'

    def test_outside_cgc_too_few_cazymes(self):
        """SusC present but only 1 substrate CAZyme -> outside_CGC."""
        cgc_df = make_cgc_df([
            ('CGC1', 'CAZyme', 'GH16'),
            ('CGC1', 'TC',     '1.B.14'),
        ])
        result = classify_cgc(
            'CGC1', cgc_df, SUBSTRATE_FAMILIES,
            pul_mode='bacteroidetes', min_cazymes=2)
        assert result == 'outside_CGC'

    def test_outside_cgc_no_substrate_cazymes(self):
        """SusC present but no substrate CAZymes -> outside_CGC."""
        cgc_df = make_cgc_df([
            ('CGC1', 'CAZyme', 'GH3'),   # Not in substrate families
            ('CGC1', 'TC',     '1.B.14'),
        ])
        result = classify_cgc(
            'CGC1', cgc_df, SUBSTRATE_FAMILIES,
            pul_mode='bacteroidetes', min_cazymes=2)
        assert result == 'outside_CGC'

    def test_min_cazymes_threshold_respected(self):
        """min_cazymes=3 requires 3 substrate CAZymes."""
        cgc_df = make_cgc_df([
            ('CGC1', 'CAZyme', 'GH16'),
            ('CGC1', 'CAZyme', 'GH17'),
            ('CGC1', 'TC',     '1.B.14'),
        ])
        # With min_cazymes=3, 2 CAZymes is not enough
        result = classify_cgc(
            'CGC1', cgc_df, SUBSTRATE_FAMILIES,
            pul_mode='bacteroidetes', min_cazymes=3)
        assert result == 'outside_CGC'

    def test_min_cazymes_one_accepts_single_cazyme(self):
        """min_cazymes=1 accepts a single substrate CAZyme."""
        cgc_df = make_cgc_df([
            ('CGC1', 'CAZyme', 'GH16'),
            ('CGC1', 'TC',     '1.B.14'),
        ])
        result = classify_cgc(
            'CGC1', cgc_df, SUBSTRATE_FAMILIES,
            pul_mode='bacteroidetes', min_cazymes=1)
        assert result == 'canonical_PUL'

    def test_non_substrate_cazymes_not_counted(self):
        """CAZymes not in substrate families don't count toward minimum."""
        cgc_df = make_cgc_df([
            ('CGC1', 'CAZyme', 'GH16'),  # substrate
            ('CGC1', 'CAZyme', 'GH3'),   # not substrate
            ('CGC1', 'CAZyme', 'CE1'),   # not substrate
            ('CGC1', 'TC',     '1.B.14'),
        ])
        # Only 1 substrate CAZyme despite 3 total CAZymes
        result = classify_cgc(
            'CGC1', cgc_df, SUBSTRATE_FAMILIES,
            pul_mode='bacteroidetes', min_cazymes=2)
        assert result == 'outside_CGC'

    def test_empty_cgc_returns_outside(self):
        """CGC with no genes returns outside_CGC."""
        cgc_df = make_cgc_df([
            ('CGC2', 'CAZyme', 'GH16'),  # Different CGC
        ])
        result = classify_cgc(
            'CGC1', cgc_df, SUBSTRATE_FAMILIES,
            pul_mode='bacteroidetes', min_cazymes=2)
        assert result == 'outside_CGC'


# ── generic mode ──────────────────────────────────────────────────────────────

class TestGenericMode:

    def test_canonical_any_tc_plus_cazymes(self):
        """Any TC gene + 2 substrate CAZymes -> canonical_PUL."""
        cgc_df = make_cgc_df([
            ('CGC1', 'CAZyme', 'GH16'),
            ('CGC1', 'CAZyme', 'GH17'),
            ('CGC1', 'TC',     '2.A.1'),  # Non-SusC transporter
        ])
        result = classify_cgc(
            'CGC1', cgc_df, SUBSTRATE_FAMILIES,
            pul_mode='generic', min_cazymes=2)
        assert result == 'canonical_PUL'

    def test_non_canonical_no_tc(self):
        """2 substrate CAZymes but no TC gene -> non_canonical_CGC."""
        cgc_df = make_cgc_df([
            ('CGC1', 'CAZyme', 'GH16'),
            ('CGC1', 'CAZyme', 'GH17'),
            ('CGC1', 'TF',     'LacI'),
        ])
        result = classify_cgc(
            'CGC1', cgc_df, SUBSTRATE_FAMILIES,
            pul_mode='generic', min_cazymes=2)
        assert result == 'non_canonical_CGC'

    def test_outside_cgc_tc_but_too_few_cazymes(self):
        """TC gene present but only 1 substrate CAZyme -> outside_CGC."""
        cgc_df = make_cgc_df([
            ('CGC1', 'CAZyme', 'GH16'),
            ('CGC1', 'TC',     '2.A.1'),
        ])
        result = classify_cgc(
            'CGC1', cgc_df, SUBSTRATE_FAMILIES,
            pul_mode='generic', min_cazymes=2)
        assert result == 'outside_CGC'

    def test_susc_also_canonical_in_generic_mode(self):
        """SusC counts as a TC gene in generic mode too."""
        cgc_df = make_cgc_df([
            ('CGC1', 'CAZyme', 'GH16'),
            ('CGC1', 'CAZyme', 'GH17'),
            ('CGC1', 'TC',     '1.B.14'),
        ])
        result = classify_cgc(
            'CGC1', cgc_df, SUBSTRATE_FAMILIES,
            pul_mode='generic', min_cazymes=2)
        assert result == 'canonical_PUL'


# ── cazyme_only mode ──────────────────────────────────────────────────────────

class TestCazymeOnlyMode:

    def test_canonical_two_cazymes_no_transporter(self):
        """2 substrate CAZymes -> canonical_PUL regardless of transporter."""
        cgc_df = make_cgc_df([
            ('CGC1', 'CAZyme', 'GH16'),
            ('CGC1', 'CAZyme', 'GH17'),
        ])
        result = classify_cgc(
            'CGC1', cgc_df, SUBSTRATE_FAMILIES,
            pul_mode='cazyme_only', min_cazymes=2)
        assert result == 'canonical_PUL'

    def test_outside_cgc_one_cazyme(self):
        """1 substrate CAZyme -> outside_CGC in cazyme_only mode."""
        cgc_df = make_cgc_df([
            ('CGC1', 'CAZyme', 'GH16'),
        ])
        result = classify_cgc(
            'CGC1', cgc_df, SUBSTRATE_FAMILIES,
            pul_mode='cazyme_only', min_cazymes=2)
        assert result == 'outside_CGC'

    def test_no_non_canonical_in_cazyme_only_mode(self):
        """cazyme_only mode never returns non_canonical_CGC."""
        cgc_df = make_cgc_df([
            ('CGC1', 'CAZyme', 'GH16'),
            ('CGC1', 'CAZyme', 'GH17'),
        ])
        result = classify_cgc(
            'CGC1', cgc_df, SUBSTRATE_FAMILIES,
            pul_mode='cazyme_only', min_cazymes=2)
        assert result != 'non_canonical_CGC'


# ── SUSC_FAMILIES constant ────────────────────────────────────────────────────

class TestSuscFamilies:

    def test_susc_families_defined(self):
        """SUSC_FAMILIES constant contains expected TCDB families."""
        assert '1.B.14' in SUSC_FAMILIES
        assert '8.A.46' in SUSC_FAMILIES

    def test_susc_families_is_list(self):
        assert isinstance(SUSC_FAMILIES, list)
