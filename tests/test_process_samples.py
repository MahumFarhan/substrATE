"""
Tests for substrate.classify_pul.process_samples().

Each test runs on a small fake dbCAN output folder (overview.tsv and
cgc_standard_out.tsv for one or two samples) built in a temporary
directory. dbCAN itself is not needed.
"""
import pandas as pd
import pytest

from substrate.classify_pul import process_samples

OVERVIEW_HEADER = ("Gene ID\tEC#\tdbCAN_hmm\tdbCAN_sub\tDIAMOND"
                   "\t#ofTools\tRecommend Results\n")
CGC_HEADER = "CGC#\tGene Type\tGene Annotation\tProtein ID\n"


def write_sample(cgc_dir, sample, overview_rows, cgc_rows=()):
    d = cgc_dir / f'output_{sample}'
    d.mkdir(parents=True)
    (d / 'overview.tsv').write_text(OVERVIEW_HEADER + ''.join(overview_rows))
    if cgc_rows:
        (d / 'cgc_standard_out.tsv').write_text(CGC_HEADER + ''.join(cgc_rows))


@pytest.fixture
def cgc_dir(tmp_path):
    """One sample, S1, with a laminarin PUL (CGC1), a transporter-less
    cluster (CGC2) and one GH16 outside any cluster."""
    write_sample(
        tmp_path, 'S1',
        overview_rows=[
            "g1\t-\tGH16_3\tGH16_e1\tGH16_3\t3\tGH16_3\n",   # in CGC1
            "g2\t-\tGH17\tGH17_e2\tGH17\t3\tGH17\n",         # in CGC1
            "g4\t-\tGH16_3\tGH16_e1\tGH16_3\t3\tGH16_3\n",   # in CGC2
            "g5\t-\tGH17\tGH17_e2\tGH17\t3\tGH17\n",         # in CGC2
            "g7\t-\tGH16_3\tGH16_e1\tGH16_3\t3\tGH16_3\n",   # no cluster
            "g8\t-\tGH10\tGH10_e3\tGH10\t3\tGH10\n",         # not laminarin
        ],
        cgc_rows=[
            "CGC1\tCAZyme\tGH16_3\tg1\n",
            "CGC1\tCAZyme\tGH17\tg2\n",
            "CGC1\tTC\t1.B.14.6.1\tg3\n",
            "CGC2\tCAZyme\tGH16_3\tg4\n",
            "CGC2\tCAZyme\tGH17\tg5\n",
            "CGC2\tTF\tsome regulator\tg6\n",
        ])
    return tmp_path


def localisation(family_hits, gene):
    rows = family_hits[family_hits['Gene ID'] == gene]
    assert len(rows) >= 1, f"{gene} not in family hits"
    return set(rows['localisation'])


# ── Localisation ──────────────────────────────────────────────────────────────

class TestLocalisation:

    def test_unknown_substrate_raises(self, cgc_dir):
        with pytest.raises(ValueError):
            process_samples(str(cgc_dir), 'not_a_substrate')

    def test_gene_in_cluster_with_transporter_is_canonical(self, cgc_dir):
        _, fam, _ = process_samples(str(cgc_dir), 'laminarin')
        assert localisation(fam, 'g1') == {'canonical_PUL'}
        assert localisation(fam, 'g2') == {'canonical_PUL'}

    def test_gene_in_cluster_without_transporter_is_non_canonical(self, cgc_dir):
        _, fam, _ = process_samples(str(cgc_dir), 'laminarin')
        assert localisation(fam, 'g4') == {'non_canonical_CGC'}

    def test_gene_outside_any_cluster(self, cgc_dir):
        _, fam, _ = process_samples(str(cgc_dir), 'laminarin')
        assert localisation(fam, 'g7') == {'outside_CGC'}

    def test_other_substrates_families_not_reported(self, cgc_dir):
        _, fam, _ = process_samples(str(cgc_dir), 'laminarin')
        assert 'g8' not in set(fam['Gene ID'])

    def test_sample_name_taken_from_folder(self, cgc_dir):
        _, fam, over = process_samples(str(cgc_dir), 'laminarin')
        assert set(fam['sample']) == {'S1'}
        assert set(over['sample']) == {'S1'}

    def test_family_number_not_matched_as_prefix(self, tmp_path):
        """GH16 must not match a GH160 gene."""
        write_sample(tmp_path, 'S1', overview_rows=[
            "g1\t-\tGH160\tGH160_e1\tGH160\t3\tGH160\n"])
        _, fam, _ = process_samples(str(tmp_path), 'laminarin')
        assert fam.empty or 'GH16' not in set(fam['matched_family'])


# ── HMM transporter evidence ──────────────────────────────────────────────────

class TestHmmEvidence:

    def hmm_df(self, gene, role, source, accession):
        return pd.DataFrame([{
            'sample': 'S1', 'Protein ID': gene, 'transporter_role': role,
            'hmm_accession': accession, 'hmm_source': source}])

    def test_hmm_hit_upgrades_cluster_to_canonical(self, cgc_dir):
        hmm = self.hmm_df('g6', 'SusD', 'Pfam', 'PF07980')
        _, fam, _ = process_samples(
            str(cgc_dir), 'laminarin', transporter_hmm_df=hmm)
        assert localisation(fam, 'g4') == {'canonical_PUL'}

    def test_transporter_source_column(self, cgc_dir):
        hmm = self.hmm_df('g6', 'SusD', 'Pfam', 'PF07980')
        _, fam, _ = process_samples(
            str(cgc_dir), 'laminarin', transporter_hmm_df=hmm)
        src = fam.set_index('Gene ID')['transporter_source']
        assert src['g1'] == 'TCDB'
        assert src['g4'] == 'Pfam'
        assert src['g7'] == ''

    def test_hmm_hits_from_another_sample_ignored(self, cgc_dir):
        hmm = self.hmm_df('g6', 'SusD', 'Pfam', 'PF07980')
        hmm['sample'] = 'S2'
        _, fam, _ = process_samples(
            str(cgc_dir), 'laminarin', transporter_hmm_df=hmm)
        assert localisation(fam, 'g4') == {'non_canonical_CGC'}

    def test_empty_hmm_table_same_as_none(self, cgc_dir):
        empty = pd.DataFrame(columns=[
            'sample', 'Protein ID', 'transporter_role',
            'hmm_accession', 'hmm_source'])
        _, with_empty, _ = process_samples(
            str(cgc_dir), 'laminarin', transporter_hmm_df=empty)
        _, with_none, _ = process_samples(str(cgc_dir), 'laminarin')
        assert (list(with_empty['localisation'])
                == list(with_none['localisation']))


# ── Empty cells in overview.tsv ───────────────────────────────────────────────

class TestEmptyAnnotationCells:
    """From pandas 3, an empty cell in one annotation column used to blank
    the combined annotation string, silently dropping the gene."""

    def test_gene_with_empty_cells_still_matched(self, tmp_path):
        write_sample(tmp_path, 'S1', overview_rows=[
            "g1\t-\tGH16_3\t\tGH16_3\t2\tGH16_3\n",   # empty dbCAN_sub
            "g2\t-\t\t\tGH17\t1\tGH17\n",              # two empty cells
            "g3\t-\tGH16_3\tGH16_e1\tGH16_3\t3\tGH16_3\n",
        ])
        _, fam, _ = process_samples(str(tmp_path), 'laminarin')
        assert set(fam['Gene ID']) == {'g1', 'g2', 'g3'}

    def test_empty_recommended_result_falls_back_to_hmm(self, tmp_path):
        write_sample(tmp_path, 'S1', overview_rows=[
            "g1\t-\tGH16_3\tGH16_e1\tGH16_3\t3\t\n"])
        _, fam, _ = process_samples(str(tmp_path), 'laminarin')
        assert list(fam['subfamily_annotation']) == ['GH16_3']

    def test_dash_recommended_result_falls_back_to_hmm(self, tmp_path):
        write_sample(tmp_path, 'S1', overview_rows=[
            "g1\t-\tGH16_3\tGH16_e1\tGH16_3\t3\t-\n"])
        _, fam, _ = process_samples(str(tmp_path), 'laminarin')
        assert list(fam['subfamily_annotation']) == ['GH16_3']
