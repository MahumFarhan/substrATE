"""
Tests for the `substrate reduced-tree` subcommand.

Each test runs the real command on a small fake output folder built in a
temporary directory: an activity table, a Newick tree and a FASTA file.
No genomes, dbCAN, MAFFT or IQ-TREE are needed.

The fake laminarin run has three genomes and two reference sequences:

    genome  gene  activity        localisation   #ofTools
    S1      g1    laminarinase    canonical_PUL  3
    S1      g2    laminarinase    outside_CGC    2
    S2      g3    laminarinase    canonical_PUL  3
    S2      g4    beta-glucanase  canonical_PUL  3
    S3      g5    xylanase        canonical_PUL  2
    S3      g1    xylanase        canonical_PUL  2   <- same gene ID as S1's g1

The test patterns file marks 'laminarinase' as strict and 'glucanase' as
permissive, so the bundled activity_patterns.tsv can change without
breaking these tests.
"""
import os

import pandas as pd
import pytest
from Bio import Phylo, SeqIO
from click.testing import CliRunner

from substrate import cli

SUB = 'laminarin'

REF1 = 'Reference__WP_1__GH16__characterised_reference'
REF2 = 'Reference__WP_2__GH16__characterised_reference'
S1G1 = 'S1__g1__GH16_3__canonical_PUL'
S1G2 = 'S1__g2__GH16_3__outside_CGC'
S2G3 = 'S2__g3__GH16_3__canonical_PUL'
S2G4 = 'S2__g4__GH16_3__canonical_PUL'
S3G5 = 'S3__g5__GH16_3__canonical_PUL'
S3G1 = 'S3__g1__GH16_3__canonical_PUL'
ALL_TIPS = [S1G1, S1G2, S2G3, S2G4, S3G5, S3G1, REF1, REF2]

GH16_TREE = (f"(({S1G1}:0.1,{S1G2}:0.1):0.1,({S2G3}:0.1,{S2G4}:0.1):0.1,"
             f"({S3G5}:0.1,{S3G1}:0.1):0.1,({REF1}:0.1,{REF2}:0.1):0.1);\n")

HITS = [
    # Gene ID, sample, activity, localisation, #ofTools
    ('g1', 'S1', 'laminarinase',   'canonical_PUL', 3),
    ('g2', 'S1', 'laminarinase',   'outside_CGC',   2),
    ('g3', 'S2', 'laminarinase',   'canonical_PUL', 3),
    ('g4', 'S2', 'beta-glucanase', 'canonical_PUL', 3),
    ('g5', 'S3', 'xylanase',       'canonical_PUL', 2),
    ('g1', 'S3', 'xylanase',       'canonical_PUL', 2),
]


# ── Fixtures and helpers ──────────────────────────────────────────────────────

@pytest.fixture
def out(tmp_path, monkeypatch):
    """Build a fake `substrate run` output folder and return its path."""
    patterns = tmp_path / 'patterns.tsv'
    patterns.write_text(
        "substrate\tpattern\tsource\treviewed\tmode\n"
        "laminarin\tlaminarinase\tcurated\tTrue\tstrict\n"
        "laminarin\tglucanase\tcurated\tTrue\tpermissive\n")
    monkeypatch.setattr(cli, '_PATTERNS_FILE', str(patterns))
    monkeypatch.setattr(cli, '_REF_METADATA', str(tmp_path / 'no_ref_metadata.tsv'))

    base = tmp_path / 'results'
    sub_dir = base / SUB
    (sub_dir / 'trees').mkdir(parents=True)
    (sub_dir / 'sequences').mkdir()

    rows = [{
        'Gene ID': gene, 'sample': sample, 'substrate_category': SUB,
        'matched_family': 'GH16', 'localisation': loc,
        'subfamily_annotation': 'GH16_3', 'EC#': '-', 'primary_ec': '-',
        'activity': act, '#ofTools': tools,
    } for gene, sample, act, loc, tools in HITS]
    pd.DataFrame(rows).to_csv(
        sub_dir / f'{SUB}_activity_annotated.tsv', sep='\t', index=False)

    (sub_dir / 'trees' / 'GH16.treefile').write_text(GH16_TREE)
    (sub_dir / 'sequences' / f'{SUB}_GH16.faa').write_text(
        ''.join(f">{tip}\nMKTAYIAKQR\n" for tip in ALL_TIPS))
    return base


def run(out, *args):
    """Run `substrate reduced-tree` on the fake folder."""
    result = CliRunner().invoke(cli.main, [
        'reduced-tree', '--substrate', SUB, '--output', str(out), *args])
    assert result.exit_code == 0, result.output
    return result


def reduced_dir(out, label='GH16'):
    return out / SUB / 'reduced_trees' / label / SUB


def tree_path(out, label='GH16', fam='GH16'):
    return reduced_dir(out, label) / 'trees' / f'{fam}.treefile'


def tips(out, label='GH16', fam='GH16'):
    tree = Phylo.read(str(tree_path(out, label, fam)), 'newick')
    return {t.name for t in tree.get_terminals()}


# ── Activity filtering ────────────────────────────────────────────────────────

class TestActivityFiltering:

    def test_default_is_strict_patterns(self, out):
        """With no --activity, only genes matching strict patterns stay."""
        result = run(out)
        assert tips(out) == {S1G1, S1G2, S2G3, REF1, REF2}
        assert 'strict patterns' in result.output

    def test_permissive_mode_adds_permissive_matches(self, out):
        run(out, '--pattern-mode', 'permissive')
        assert tips(out) == {S1G1, S1G2, S2G3, S2G4, REF1, REF2}

    def test_explicit_activity_overrides_patterns(self, out):
        """--activity selects by exact label; patterns are not applied."""
        run(out, '--activity', 'xylanase')
        assert tips(out) == {S3G5, S3G1, REF1, REF2}

    def test_pattern_mode_ignored_when_activity_given(self, out):
        run(out, '--activity', 'xylanase', '--pattern-mode', 'permissive')
        assert tips(out) == {S3G5, S3G1, REF1, REF2}

    def test_activity_with_no_matches_writes_nothing(self, out):
        result = run(out, '--family', 'GH16', '--activity', 'no such activity')
        assert not tree_path(out).exists()
        assert 'no hits after filtering' in result.output


# ── Other filters ─────────────────────────────────────────────────────────────

class TestOtherFilters:

    def test_localisation_filter(self, out):
        run(out, '--localisation', 'canonical_PUL')
        assert tips(out, 'GH16_canonical_PUL') == {S1G1, S2G3, REF1, REF2}

    def test_one_per_genome_keeps_best_supported_gene(self, out):
        """S1 has two matching genes; the one with more tools is kept."""
        run(out, '--one-per-genome')
        assert tips(out, 'GH16_1pg') == {S1G1, S2G3, REF1, REF2}

    def test_filters_combine_in_output_folder_name(self, out):
        run(out, '--localisation', 'canonical_PUL', '--one-per-genome')
        assert tree_path(out, 'GH16_canonical_PUL_1pg').exists()

    def test_exclude_sample(self, out):
        run(out, '--exclude-sample', 'S1')
        assert tips(out) == {S2G3, REF1, REF2}

    def test_exclude_samples_file(self, out, tmp_path):
        listing = tmp_path / 'exclude.txt'
        listing.write_text("S1\n\n")
        run(out, '--exclude-samples-file', str(listing))
        assert tips(out) == {S2G3, REF1, REF2}

    def test_family_option_limits_families(self, out):
        result = run(out, '--family', 'GH16')
        assert 'Processing 1 families' in result.output
        assert tree_path(out).exists()


# ── Tip matching ──────────────────────────────────────────────────────────────

class TestTipMatching:

    def test_references_always_kept(self, out):
        run(out, '--localisation', 'outside_CGC')
        assert {REF1, REF2} <= tips(out, 'GH16_outside_CGC')

    def test_same_gene_id_in_other_genome_not_kept(self, out):
        """S3 also has a gene called g1, but it is a xylanase. It must not
        ride along with S1's g1 just because the gene ID is the same."""
        run(out)
        assert S1G1 in tips(out)
        assert S3G1 not in tips(out)

    def test_shared_gene_id_not_kept_in_the_other_direction(self, out):
        """Selecting S3's xylanase g1 must not bring in S1's g1."""
        run(out, '--activity', 'xylanase')
        assert S3G1 in tips(out)
        assert S1G1 not in tips(out)

    def test_excluded_sample_stays_out_despite_shared_gene_id(self, out):
        """With S1 excluded, S1's g1 must not return because S3 is
        selected and also has a gene called g1."""
        run(out, '--activity', 'xylanase', '--exclude-sample', 'S1')
        assert tips(out) == {S3G5, S3G1, REF1, REF2}


# ── Outputs ───────────────────────────────────────────────────────────────────

class TestOutputs:

    def test_filtered_fasta_matches_tree(self, out):
        run(out)
        faa = reduced_dir(out) / 'sequences' / f'{SUB}_GH16.faa'
        ids = {r.id for r in SeqIO.parse(str(faa), 'fasta')}
        assert ids == tips(out)

    def test_filtered_activity_table_written(self, out):
        run(out)
        df = pd.read_csv(
            reduced_dir(out) / f'{SUB}_activity_annotated.tsv', sep='\t')
        assert set(zip(df['sample'], df['Gene ID'])) == {
            ('S1', 'g1'), ('S1', 'g2'), ('S2', 'g3')}

    def test_itol_annotations_generated_without_error(self, out):
        result = run(out)
        assert 'iTOL failed' not in result.output
        assert (reduced_dir(out) / 'itol_annotations').is_dir()

    def test_source_tree_not_modified(self, out):
        run(out)
        src = out / SUB / 'trees' / 'GH16.treefile'
        assert src.read_text() == GH16_TREE


# ── Finding the source tree ───────────────────────────────────────────────────

class TestSourceTree:

    def test_pruned_tree_preferred_over_full_tree(self, out):
        """If a .pruned.treefile exists it is the one that gets reduced."""
        pruned = f"({S1G1}:0.1,{S2G3}:0.1,{REF1}:0.1);\n"
        (out / SUB / 'trees' / 'GH16.pruned.treefile').write_text(pruned)
        run(out)
        assert tips(out) == {S1G1, S2G3, REF1}

    def test_tree_named_by_substrate_tree_command_is_found(self, out):
        """`substrate tree` names files {substrate}_{family}.treefile."""
        trees = out / SUB / 'trees'
        os.rename(trees / 'GH16.treefile', trees / f'{SUB}_GH16.treefile')
        run(out)
        assert tips(out) == {S1G1, S1G2, S2G3, REF1, REF2}

    def test_missing_tree_is_skipped(self, out):
        os.remove(out / SUB / 'trees' / 'GH16.treefile')
        result = run(out)
        assert 'no treefile found' in result.output
        assert not tree_path(out).exists()

    def test_too_few_tips_is_skipped(self, out):
        """Fewer than three tips left: no tree is written."""
        (out / SUB / 'trees' / 'GH16.treefile').write_text(
            f"({S1G1}:0.1,{S2G4}:0.1,{S3G5}:0.1,{S3G1}:0.1);\n")
        result = run(out)
        assert 'too few sequences' in result.output
        assert not tree_path(out).exists()

    def test_missing_activity_file_warns(self, out):
        os.remove(out / SUB / f'{SUB}_activity_annotated.tsv')
        result = run(out)
        assert 'No activity file found' in result.output


# ── --force and --force-visualise ─────────────────────────────────────────────

class TestForce:

    def test_existing_output_skipped_without_force(self, out):
        run(out)
        tree_path(out).write_text("MARKER")
        result = run(out)
        assert 'already exists' in result.output
        assert tree_path(out).read_text() == "MARKER"

    def test_force_overwrites(self, out):
        run(out)
        tree_path(out).write_text("MARKER")
        run(out, '--force')
        assert tips(out) == {S1G1, S1G2, S2G3, REF1, REF2}

    def test_force_visualise_leaves_tree_untouched(self, out):
        run(out)
        tree_path(out).write_text("MARKER")
        result = run(out, '--force-visualise')
        assert tree_path(out).read_text() == "MARKER"
        assert 'iTOL annotations written' in result.output
