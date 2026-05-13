"""
Unit tests for substrate.activity and substrate.parse_substrates modules.

Covers:
  - normalise_activity()      — pure function
  - extract_primary_ec()      — pure function
  - get_activity_label()      — pure function (mock dicts)
  - load_ec_names()           — file I/O (minimal enzyme.dat fixture)
  - load_family_activities()  — file I/O (minimal TSV fixture)
  - annotate_references()     — file I/O (minimal TSV fixture)
  - load_patterns()           — file I/O (minimal TSV fixture),
                                including pattern_mode and back-compat
"""
import pytest
import pandas as pd

from substrate.activity import (
    normalise_activity,
    extract_primary_ec,
    get_activity_label,
    load_ec_names,
    load_family_activities,
    annotate_references,
)
from substrate.parse_substrates import load_patterns


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture
def minimal_expasy(tmp_path):
    """Write a minimal enzyme.dat with two entries and return its path."""
    content = (
        "ID   3.2.1.39\n"
        "DE   Glucan endo-1,3-beta-D-glucosidase.\n"
        "//\n"
        "ID   3.2.1.6\n"
        "DE   Endo-1,3(4)-beta-glucanase.\n"
        "//\n"
        "ID   3.2.1.91\n"
        "DE   Deleted entry.\n"
        "//\n"
        "ID   3.2.1.78\n"
        "DE   Mannan endo-1,4-beta-mannosidase.\n"
        "//\n"
    )
    path = tmp_path / "enzyme.dat"
    path.write_text(content)
    return str(path)


@pytest.fixture
def minimal_fam_sub_map(tmp_path):
    """Write a minimal fam-substrate-mapping.tsv and return its path."""
    rows = (
        "beta-glucan\tlaminarin\tGH16\tlaminarinase\t3.2.1.39\n"
        "beta-glucan\tlaminarin\tGH17\tbeta-glucosidase\t3.2.1.21\n"
        "xylan\txylan\tGH10\tendo-1,4-beta-xylanase\t3.2.1.8\n"
    )
    path = tmp_path / "fam-substrate-mapping.tsv"
    path.write_text(rows)
    return str(path)


@pytest.fixture
def minimal_ref_metadata(tmp_path):
    """Write a minimal reference_metadata.tsv (current schema) and return path."""
    content = (
        "accession\tsubstrate\tfamily\tsubfamily\tprotein_name\tec_numbers\n"
        "WP_001\tlaminarin\tGH16\tGH16_3\tLaminarinase A\t3.2.1.39\n"
        "WP_002\tlaminarin\tGH17\tGH17\tBeta-glucosidase\t3.2.1.21\n"
        "WP_003\txylan\tGH10\tGH10\tXylanase\t3.2.1.8\n"
    )
    path = tmp_path / "reference_metadata.tsv"
    path.write_text(content)
    return str(path)


@pytest.fixture
def legacy_ref_metadata(tmp_path):
    """Write a minimal reference_metadata.tsv (legacy schema) and return path."""
    content = (
        "accession\tsubstrate\tfamily\tsubfamily\tactivity\n"
        "WP_004\tlaminarin\tGH55\tGH55\texo-1,3-beta-glucosidase\n"
    )
    path = tmp_path / "reference_metadata_legacy.tsv"
    path.write_text(content)
    return str(path)


@pytest.fixture
def patterns_with_mode(tmp_path):
    """Write an activity_patterns.tsv with a mode column."""
    content = (
        "substrate\tpattern\tsource\treviewed\tmode\n"
        "laminarin\tlaminarinase\tcurated\tTrue\tstrict\n"
        "laminarin\tbeta-glucan\tcurated\tTrue\tpermissive\n"
        "laminarin\tglucanase\tcurated\tFalse\tpermissive\n"
        "xylan\txylanase\tcurated\tTrue\tstrict\n"
        "xylan\tarabinosidase\tcurated\tFalse\tpermissive\n"
    )
    path = tmp_path / "activity_patterns.tsv"
    path.write_text(content)
    return str(path)


@pytest.fixture
def patterns_without_mode(tmp_path):
    """Write a legacy activity_patterns.tsv with no mode column."""
    content = (
        "substrate\tpattern\tsource\treviewed\n"
        "laminarin\tlaminarinase\tcurated\tTrue\n"
        "laminarin\tbeta-glucan\tcurated\tTrue\n"
    )
    path = tmp_path / "activity_patterns_legacy.tsv"
    path.write_text(content)
    return str(path)


# ── normalise_activity ────────────────────────────────────────────────────────

class TestNormaliseActivity:

    def test_normal_activity_unchanged(self):
        """Short, clean activity strings are returned unchanged."""
        assert normalise_activity('laminarinase', 'GH16') == 'laminarinase'

    def test_unknown_variants_return_unknown(self):
        for val in ['unknown', '-', '', 'nan']:
            assert normalise_activity(val, 'GH16') == 'unknown'

    def test_none_returns_unknown(self):
        assert normalise_activity(None, 'GH16') == 'unknown'

    def test_verbose_cbm_replaced(self):
        """Verbose CBM description replaced with clean label."""
        long_desc = (
            'Modules of approx. 100 residues with glycogen-binding '
            'properties found in a range of enzymes'
        )
        result = normalise_activity(long_desc, 'CBM48')
        assert result == 'carbohydrate-binding module (CBM48)'

    def test_verbose_gh_replaced(self):
        """Verbose GH description replaced with clean label."""
        long_desc = (
            'Multi-domain enzyme with additional catalytic modules and '
            'substrate-binding regions linked by flexible linkers'
        )
        result = normalise_activity(long_desc, 'GH5')
        assert result == 'glycoside hydrolase (GH5)'

    def test_verbose_unknown_family_truncated(self):
        """Verbose description with unrecognised family prefix is truncated."""
        long_desc = 'X' * 65
        result = normalise_activity(long_desc, 'XX99')
        assert result.endswith('...')
        assert len(result) == 60

    def test_alpha_prefix_normalised(self):
        assert normalise_activity('a-glucosidase', 'GH31') == 'alpha-glucosidase'

    def test_beta_prefix_normalised(self):
        assert normalise_activity('b-galactosidase', 'GH2') == 'beta-galactosidase'

    def test_multiple_prefix_variants_normalised(self):
        result = normalise_activity('a-manno-b-glucosidase', 'GH130')
        assert 'alpha-manno' in result
        assert 'beta-gluco' in result

    def test_whitespace_stripped(self):
        assert normalise_activity('  laminarinase  ', 'GH16') == 'laminarinase'


# ── extract_primary_ec ────────────────────────────────────────────────────────

class TestExtractPrimaryEc:

    def test_simple_ec_returned(self):
        assert extract_primary_ec('3.2.1.39:1') == '3.2.1.39'

    def test_dash_returns_dash(self):
        assert extract_primary_ec('-') == '-'

    def test_nan_returns_dash(self):
        assert extract_primary_ec(float('nan')) == '-'

    def test_empty_returns_dash(self):
        assert extract_primary_ec('') == '-'

    def test_highest_count_wins(self):
        """EC with higher count is selected over lower count."""
        assert extract_primary_ec('3.2.1.39:1;3.2.1.6:2') == '3.2.1.6'

    def test_pipe_separated_first_segment_used(self):
        """Only the first non-empty pipe segment is used."""
        assert extract_primary_ec('3.2.1.39:1|-|-') == '3.2.1.39'

    def test_all_dash_pipe_returns_dash(self):
        assert extract_primary_ec('-|-|-') == '-'

    def test_no_count_defaults_to_one(self):
        """EC string without count still parses correctly."""
        assert extract_primary_ec('3.2.1.39') == '3.2.1.39'

    def test_partial_ec_parsed(self):
        """Partial EC numbers (e.g. 3.2.1.-) are parsed."""
        result = extract_primary_ec('3.2.1.-:1')
        assert result == '3.2.1.-'

    def test_equal_counts_first_wins(self):
        """When counts are equal, first EC encountered is returned."""
        result = extract_primary_ec('3.2.1.39:2;3.2.1.6:2')
        assert result == '3.2.1.39'


# ── get_activity_label ────────────────────────────────────────────────────────

class TestGetActivityLabel:

    EC_ACTIVITIES = {
        '3.2.1.39': 'Glucan endo-1,3-beta-D-glucosidase',
        '3.2.1.6':  'Endo-1,3(4)-beta-glucanase',
        '3.2.1.78': 'Mannan endo-1,4-beta-mannosidase',
        '3.2.1.21': 'Beta-D-glucosidase',
    }

    FAMILY_ACTIVITIES = {
        'GH16': [('3.2.1.39', 'laminarinase')],
        'GH17': [('3.2.1.6',  'beta-glucosidase')],
        'GH10': [('3.2.1.8',  'endo-1,4-beta-xylanase')],
    }

    def test_direct_ec_lookup(self):
        result = get_activity_label(
            '3.2.1.39', 'GH16',
            self.EC_ACTIVITIES, self.FAMILY_ACTIVITIES)
        assert result == 'Glucan endo-1,3-beta-D-glucosidase'

    def test_family_fallback_when_no_ec(self):
        """Falls back to family-level activity when EC is '-'."""
        result = get_activity_label(
            '-', 'GH16',
            self.EC_ACTIVITIES, self.FAMILY_ACTIVITIES)
        assert result == 'laminarinase'

    def test_unknown_when_no_ec_and_no_family(self):
        result = get_activity_label(
            '-', 'GH99',
            self.EC_ACTIVITIES, self.FAMILY_ACTIVITIES)
        assert result == 'unknown'

    def test_partial_ec_lookup(self):
        """Partial EC (3.2.1.-) resolved via prefix matching."""
        result = get_activity_label(
            '3.2.1.-', 'GH16',
            self.EC_ACTIVITIES, self.FAMILY_ACTIVITIES)
        assert result != 'unknown'

    def test_subfamily_annotation_uses_top_family(self):
        """GH16_3 should match GH16 in family_activities."""
        result = get_activity_label(
            '-', 'GH16_3',
            self.EC_ACTIVITIES, self.FAMILY_ACTIVITIES)
        assert result == 'laminarinase'

    def test_empty_family_activities_list_returns_unknown(self):
        family_activities = {'GH16': []}
        result = get_activity_label(
            '-', 'GH16',
            self.EC_ACTIVITIES, family_activities)
        assert result == 'unknown'


# ── load_ec_names ─────────────────────────────────────────────────────────────

class TestLoadEcNames:

    def test_valid_entries_loaded(self, minimal_expasy):
        ec_names = load_ec_names(minimal_expasy)
        assert '3.2.1.39' in ec_names
        assert ec_names['3.2.1.39'] == 'Glucan endo-1,3-beta-D-glucosidase'

    def test_multiple_entries_loaded(self, minimal_expasy):
        ec_names = load_ec_names(minimal_expasy)
        assert '3.2.1.6' in ec_names
        assert '3.2.1.78' in ec_names

    def test_deleted_entry_excluded(self, minimal_expasy):
        """Entries with 'Deleted' in name are skipped."""
        ec_names = load_ec_names(minimal_expasy)
        assert '3.2.1.91' not in ec_names

    def test_returns_dict(self, minimal_expasy):
        assert isinstance(load_ec_names(minimal_expasy), dict)

    def test_trailing_period_stripped(self, minimal_expasy):
        """Trailing period in DE line is stripped."""
        ec_names = load_ec_names(minimal_expasy)
        for name in ec_names.values():
            assert not name.endswith('.')


# ── load_family_activities ────────────────────────────────────────────────────

class TestLoadFamilyActivities:

    def test_families_loaded(self, minimal_fam_sub_map):
        fam_acts = load_family_activities(minimal_fam_sub_map)
        assert 'GH16' in fam_acts
        assert 'GH17' in fam_acts
        assert 'GH10' in fam_acts

    def test_ec_and_activity_tuples(self, minimal_fam_sub_map):
        fam_acts = load_family_activities(minimal_fam_sub_map)
        gh16 = fam_acts['GH16']
        assert len(gh16) == 1
        ec, activity = gh16[0]
        assert ec == '3.2.1.39'
        assert activity == 'laminarinase'

    def test_returns_dict(self, minimal_fam_sub_map):
        assert isinstance(load_family_activities(minimal_fam_sub_map), dict)

    def test_missing_ec_excluded(self, tmp_path):
        """Rows with no EC number are excluded from the list."""
        content = (
            "cat\tsub\tGH16\tlaminarinase\tnan\n"
            "cat\tsub\tGH17\tbeta-glucosidase\t3.2.1.21\n"
        )
        p = tmp_path / "fam.tsv"
        p.write_text(content)
        fam_acts = load_family_activities(str(p))
        assert fam_acts.get('GH16', []) == []
        assert len(fam_acts['GH17']) == 1


# ── annotate_references ───────────────────────────────────────────────────────

class TestAnnotateReferences:

    def test_current_schema_loaded(self, minimal_ref_metadata):
        df = annotate_references(minimal_ref_metadata, 'laminarin')
        assert len(df) == 2
        assert set(df['Gene ID']) == {'WP_001', 'WP_002'}

    def test_legacy_schema_loaded(self, legacy_ref_metadata):
        df = annotate_references(legacy_ref_metadata, 'laminarin')
        assert len(df) == 1
        assert df.iloc[0]['activity'] == 'exo-1,3-beta-glucosidase'

    def test_substrate_filter_applied(self, minimal_ref_metadata):
        """Only rows for the requested substrate are returned."""
        df = annotate_references(minimal_ref_metadata, 'xylan')
        assert len(df) == 1
        assert df.iloc[0]['Gene ID'] == 'WP_003'

    def test_empty_for_missing_substrate(self, minimal_ref_metadata):
        df = annotate_references(minimal_ref_metadata, 'starch')
        assert df.empty

    def test_missing_file_returns_empty(self, tmp_path):
        df = annotate_references(
            str(tmp_path / 'nonexistent.tsv'), 'laminarin')
        assert df.empty

    def test_required_columns_present(self, minimal_ref_metadata):
        df = annotate_references(minimal_ref_metadata, 'laminarin')
        for col in ['Gene ID', 'sample', 'substrate_category',
                    'matched_family', 'localisation', 'activity']:
            assert col in df.columns

    def test_sample_is_reference(self, minimal_ref_metadata):
        df = annotate_references(minimal_ref_metadata, 'laminarin')
        assert (df['sample'] == 'Reference').all()

    def test_localisation_is_characterised_reference(self, minimal_ref_metadata):
        df = annotate_references(minimal_ref_metadata, 'laminarin')
        assert (df['localisation'] == 'characterised_reference').all()

    def test_current_schema_activity_includes_ec(self, minimal_ref_metadata):
        """Current schema activity is formatted as 'name [ec]'."""
        df = annotate_references(minimal_ref_metadata, 'laminarin')
        wp001 = df[df['Gene ID'] == 'WP_001'].iloc[0]
        assert 'Laminarinase A' in wp001['activity']
        assert '3.2.1.39' in wp001['activity']


# ── load_patterns ─────────────────────────────────────────────────────────────

class TestLoadPatterns:

    def test_permissive_loads_all(self, patterns_with_mode):
        """Permissive mode returns all patterns regardless of mode column."""
        df = load_patterns(patterns_file=patterns_with_mode,
                           pattern_mode='permissive')
        assert len(df) == 5

    def test_strict_loads_only_strict(self, patterns_with_mode):
        """Strict mode returns only patterns with mode='strict'."""
        df = load_patterns(patterns_file=patterns_with_mode,
                           pattern_mode='strict')
        assert len(df) == 2
        assert set(df['pattern']) == {'laminarinase', 'xylanase'}

    def test_substrate_filter_applied(self, patterns_with_mode):
        df = load_patterns(substrate='laminarin',
                           patterns_file=patterns_with_mode,
                           pattern_mode='permissive')
        assert len(df) == 3
        assert (df['substrate'] == 'laminarin').all()

    def test_substrate_and_strict_combined(self, patterns_with_mode):
        df = load_patterns(substrate='laminarin',
                           patterns_file=patterns_with_mode,
                           pattern_mode='strict')
        assert len(df) == 1
        assert df.iloc[0]['pattern'] == 'laminarinase'

    def test_legacy_file_no_mode_column(self, patterns_without_mode):
        """Files without a mode column are treated as all-permissive."""
        df = load_patterns(patterns_file=patterns_without_mode,
                           pattern_mode='permissive')
        assert len(df) == 2
        assert 'mode' in df.columns
        assert (df['mode'] == 'permissive').all()

    def test_legacy_file_strict_returns_empty(self, patterns_without_mode):
        """Legacy file with no mode column returns empty under strict."""
        df = load_patterns(patterns_file=patterns_without_mode,
                           pattern_mode='strict')
        assert df.empty

    def test_invalid_pattern_mode_raises(self, patterns_with_mode):
        with pytest.raises(ValueError, match="pattern_mode must be"):
            load_patterns(patterns_file=patterns_with_mode,
                          pattern_mode='invalid')

    def test_missing_file_returns_empty_df(self, tmp_path):
        df = load_patterns(
            patterns_file=str(tmp_path / 'nonexistent.tsv'))
        assert df.empty
        assert list(df.columns) == [
            'substrate', 'pattern', 'source', 'reviewed', 'mode']

    def test_returns_dataframe(self, patterns_with_mode):
        df = load_patterns(patterns_file=patterns_with_mode)
        assert isinstance(df, pd.DataFrame)

    def test_index_reset(self, patterns_with_mode):
        """Returned DataFrame has a clean 0-based index."""
        df = load_patterns(substrate='laminarin',
                           patterns_file=patterns_with_mode)
        assert list(df.index) == list(range(len(df)))
