"""
Tests for the bundled substrate/data/activity_patterns.tsv.

These check the curated file itself, not the loading code (see
TestLoadPatterns in test_activity.py for that).
"""
import os

import pandas as pd
import pytest

import substrate
from substrate.parse_substrates import load_patterns

PATTERNS_FILE = os.path.join(
    os.path.dirname(substrate.__file__), 'data', 'activity_patterns.tsv')


def patterns(substrate_name, mode):
    df = load_patterns(substrate=substrate_name, patterns_file=PATTERNS_FILE,
                       pattern_mode=mode)
    return [p.lower() for p in df['pattern']]


def matches(activity, substrate_name, mode):
    """Same substring test that reduced-tree applies."""
    return any(p in activity.lower() for p in patterns(substrate_name, mode))


# ── File integrity ────────────────────────────────────────────────────────────

class TestPatternsFile:

    def test_expected_columns(self):
        df = pd.read_csv(PATTERNS_FILE, sep='\t', dtype=str)
        assert list(df.columns) == [
            'substrate', 'pattern', 'source', 'reviewed', 'mode']

    def test_mode_is_strict_or_permissive(self):
        df = pd.read_csv(PATTERNS_FILE, sep='\t', dtype=str)
        assert set(df['mode']) <= {'strict', 'permissive'}

    def test_no_pattern_listed_twice_for_a_substrate(self):
        """A pattern must not appear as both strict and permissive."""
        df = pd.read_csv(PATTERNS_FILE, sep='\t', dtype=str)
        dups = df[df.duplicated(['substrate', 'pattern'], keep=False)]
        assert dups.empty, dups[['substrate', 'pattern', 'mode']].to_string()

    def test_no_empty_patterns(self):
        df = pd.read_csv(PATTERNS_FILE, sep='\t', dtype=str)
        assert df['pattern'].notna().all()
        assert (df['pattern'].str.strip() != '').all()


# ── Glycogen ──────────────────────────────────────────────────────────────────

class TestGlycogenStrict:

    def test_strict_set(self):
        assert set(patterns('glycogen', 'strict')) == {
            'glycogen', 'isoamylase', 'amylo-alpha-1,6-glucosidase',
            '1,4-alpha-glucan', '4-alpha-glucanotransferase',
            'glucan phosphorylase', 'glucan 1,4-alpha'}

    @pytest.mark.parametrize('activity', [
        'glycogen phosphorylase',
        'a-glucan phosphorylase',
        'isoamylase',
        'amylo-alpha-1,6-glucosidase',
        '1,4-alpha-glucan branching enzyme',
        '4-alpha-glucanotransferase',
        'glucan 1,4-alpha-glucosidase',
        'glucan 1,4-alpha-maltohydrolase',
    ])
    def test_alpha_glucan_enzymes_kept(self, activity):
        assert matches(activity, 'glycogen', 'strict')

    @pytest.mark.parametrize('activity', [
        # beta-glucan and other unrelated enzymes
        'glucan endo-1,3-beta-D-glucosidase',
        'glucan 1,3-beta-glucosidase',
        '6-phospho-beta-glucosidase',
        'alpha-1,3-glucanase',
        'glucan 1,3-alpha-glucosidase',
        'glucan 1,6-alpha-glucosidase',
        'mannosyl-oligosaccharide glucosidase',
        'amylosucrase',
        'isomaltose a-1,6-glucosidase',
        # deliberately permissive-only: act on other alpha-glucans too
        'alpha-amylase',
        'alpha-glucosidase',
        'pullulanase',
    ])
    def test_generic_or_unrelated_enzymes_not_kept(self, activity):
        assert not matches(activity, 'glycogen', 'strict')


class TestGlycogenPermissive:

    @pytest.mark.parametrize('activity', [
        'alpha-amylase', 'alpha-glucosidase', 'isoamylase',
        'glycogen phosphorylase', '1,4-alpha-glucan branching enzyme',
    ])
    def test_broad_matches_kept(self, activity):
        assert matches(activity, 'glycogen', 'permissive')

    def test_strict_patterns_are_a_subset_of_permissive(self):
        assert set(patterns('glycogen', 'strict')) < set(
            patterns('glycogen', 'permissive'))

    def test_pullulanase_is_not_a_glycogen_pattern(self):
        assert not matches('pullulanase', 'glycogen', 'permissive')


# ── Laminarin ─────────────────────────────────────────────────────────────────

class TestLaminarinStrict:

    def test_laminarinase_names_kept(self):
        assert matches('glucan endo-1,3-beta-D-glucosidase', 'laminarin', 'strict')
        assert matches('laminarinase', 'laminarin', 'strict')

    def test_alpha_glucan_enzymes_not_kept(self):
        assert not matches('1,4-alpha-glucan branching enzyme', 'laminarin', 'strict')
        assert not matches('alpha-amylase', 'laminarin', 'strict')
