"""
Unit tests for colour assignment in substrate.itol.

Two different activities sharing a colour makes them indistinguishable
in a tree legend, so these tests guard against duplicates both in the
bundled palette and in the assignment logic.
"""
import os

import substrate
from substrate.itol import (
    assign_activity_colours,
    load_colour_palettes,
    _generate_hsl_palette,
)

COLOURS_FILE = os.path.join(
    os.path.dirname(substrate.__file__), 'data', 'default_colours.tsv')

PALETTE = ['#111111', '#222222', '#333333']


class TestBundledPalettes:

    def test_activity_palette_has_no_duplicates(self):
        _, activity_palette = load_colour_palettes(COLOURS_FILE)
        assert len(activity_palette) == len(set(activity_palette))

    def test_sample_palette_has_no_duplicates(self):
        sample_palette, _ = load_colour_palettes(COLOURS_FILE)
        assert len(sample_palette) == len(set(sample_palette))

    def test_palettes_not_empty(self):
        sample_palette, activity_palette = load_colour_palettes(COLOURS_FILE)
        assert sample_palette and activity_palette


class TestGenerateHslPalette:

    def test_returns_requested_number_of_distinct_colours(self):
        for n in (1, 26, 60):
            colours = _generate_hsl_palette(n)
            assert len(colours) == n
            assert len(set(colours)) == n

    def test_deterministic(self):
        assert _generate_hsl_palette(30) == _generate_hsl_palette(30)


class TestAssignActivityColours:

    def test_within_palette_uses_palette_in_alphabetical_order(self):
        result = assign_activity_colours(['zeta', 'alpha', 'mid'], PALETTE)
        assert result == {
            'alpha': '#111111', 'mid': '#222222', 'zeta': '#333333'}

    def test_repeated_activity_gets_one_entry(self):
        result = assign_activity_colours(['a', 'a', 'b'], PALETTE)
        assert set(result) == {'a', 'b'}

    def test_exactly_palette_size_has_no_duplicates(self):
        result = assign_activity_colours(['a', 'b', 'c'], PALETTE)
        assert len(set(result.values())) == 3

    def test_beyond_palette_size_has_no_duplicates(self):
        """More activities than palette colours must not reuse a colour."""
        activities = [f'activity_{i:02d}' for i in range(10)]
        result = assign_activity_colours(activities, PALETTE)
        assert len(result) == 10
        assert len(set(result.values())) == 10

    def test_same_input_gives_same_colours(self):
        activities = [f'activity_{i:02d}' for i in range(10)]
        assert (assign_activity_colours(activities, PALETTE)
                == assign_activity_colours(list(reversed(activities)), PALETTE))

    def test_references_coloured_separately_from_genomic(self):
        result = assign_activity_colours(
            ['laminarinase', 'reference: laminarin'], PALETTE)
        assert result['laminarinase'] == '#111111'
        assert result['reference: laminarin'] not in PALETTE

    def test_references_do_not_use_up_the_genomic_palette(self):
        """Reference labels are not counted against the palette size."""
        activities = ['a', 'b', 'c',
                      'reference: x', 'reference: y', 'reference: z']
        result = assign_activity_colours(activities, PALETTE)
        assert [result[k] for k in ('a', 'b', 'c')] == PALETTE
