"""
Unit tests for substrate.transporter_hmm.

hmmsearch itself is never run: parse_tblout_gene_ids() is tested on small
hand-written result files, and run_hmmsearch() is replaced with a stub
for the annotate_sample()/annotate_all_samples() tests. HMMER is not
needed to run these.
"""
import pandas as pd

from substrate import transporter_hmm
from substrate.transporter_hmm import (
    parse_tblout_gene_ids,
    annotate_sample,
    annotate_all_samples,
)

COLUMNS = ['sample', 'Protein ID', 'transporter_role',
           'hmm_accession', 'hmm_source']

TBLOUT = (
    "#                                                               --- full sequence ----\n"
    "# target name        accession  query name           accession    E-value  score  bias\n"
    "#------------------- ---------- -------------------- ---------- --------- ------ -----\n"
    "gene_0001            -          TonB-Xanth-Caul      TIGR04056    1.2e-250  830.1   9.9\n"
    "gene_0420            -          TonB-Xanth-Caul      TIGR04056    3.4e-180  598.7   4.2\n"
    "\n"
    "#\n"
    "# Program:         hmmsearch\n"
    "# [ok]\n"
)


# ── parse_tblout_gene_ids ─────────────────────────────────────────────────────

class TestParseTbloutGeneIds:

    def test_gene_ids_read_from_first_column(self, tmp_path):
        p = tmp_path / "hits.tblout"
        p.write_text(TBLOUT)
        assert parse_tblout_gene_ids(str(p)) == {'gene_0001', 'gene_0420'}

    def test_comment_and_blank_lines_skipped(self, tmp_path):
        p = tmp_path / "hits.tblout"
        p.write_text("# only comments\n\n#\n")
        assert parse_tblout_gene_ids(str(p)) == set()

    def test_missing_file_returns_empty_set(self, tmp_path):
        assert parse_tblout_gene_ids(str(tmp_path / "absent.tblout")) == set()

    def test_duplicate_hits_counted_once(self, tmp_path):
        p = tmp_path / "hits.tblout"
        p.write_text("gene_1 - model ACC 1e-50 100 0\n"
                     "gene_1 - model ACC 1e-40  90 0\n")
        assert parse_tblout_gene_ids(str(p)) == {'gene_1'}


# ── annotate_sample / annotate_all_samples ────────────────────────────────────

def _fake_hmmsearch(hits_by_model):
    """Return a stand-in for run_hmmsearch that writes fixed hits,
    chosen by which HMM file was requested."""
    def _run(hmm_path, faa_path, tblout_path, cut_tc=True):
        genes = next(v for k, v in hits_by_model.items() if k in hmm_path)
        with open(tblout_path, 'w') as f:
            f.write("# header\n")
            for g in genes:
                f.write(f"{g} - model ACC 1e-50 100 0\n")
    return _run


def _make_sample(cgc_dir, sample, with_faa=True):
    d = cgc_dir / f"output_{sample}"
    d.mkdir(parents=True)
    if with_faa:
        (d / "uniInput.faa").write_text(">gene_1\nMKT\n")


class TestAnnotateSample:

    def test_missing_faa_returns_empty_with_columns(self, tmp_path):
        _make_sample(tmp_path, 'S1', with_faa=False)
        df = annotate_sample('S1', str(tmp_path), 'TIGR04056.hmm', 'PF07980.hmm')
        assert df.empty
        assert list(df.columns) == COLUMNS

    def test_hits_labelled_by_role_and_source(self, tmp_path, monkeypatch):
        _make_sample(tmp_path, 'S1')
        monkeypatch.setattr(transporter_hmm, 'run_hmmsearch', _fake_hmmsearch({
            'TIGR04056': ['gene_1'], 'PF07980': ['gene_2', 'gene_3']}))
        df = annotate_sample('S1', str(tmp_path), 'TIGR04056.hmm', 'PF07980.hmm')

        assert list(df.columns) == COLUMNS
        assert (df['sample'] == 'S1').all()
        susc = df[df['transporter_role'] == 'SusC']
        susd = df[df['transporter_role'] == 'SusD']
        assert set(susc['Protein ID']) == {'gene_1'}
        assert set(susd['Protein ID']) == {'gene_2', 'gene_3'}
        assert (susc['hmm_accession'] == 'TIGR04056').all()
        assert (susc['hmm_source'] == 'TIGRFAM').all()
        assert (susd['hmm_accession'] == 'PF07980').all()
        assert (susd['hmm_source'] == 'Pfam').all()

    def test_no_hits_returns_empty_with_columns(self, tmp_path, monkeypatch):
        _make_sample(tmp_path, 'S1')
        monkeypatch.setattr(transporter_hmm, 'run_hmmsearch', _fake_hmmsearch({
            'TIGR04056': [], 'PF07980': []}))
        df = annotate_sample('S1', str(tmp_path), 'TIGR04056.hmm', 'PF07980.hmm')
        assert df.empty
        assert list(df.columns) == COLUMNS


class TestAnnotateAllSamples:

    def test_combines_samples_and_ignores_other_entries(self, tmp_path, monkeypatch):
        _make_sample(tmp_path, 'S1')
        _make_sample(tmp_path, 'S2')
        (tmp_path / "S3.fna").write_text(">contig\nACGT\n")   # not a sample dir
        monkeypatch.setattr(transporter_hmm, 'run_hmmsearch', _fake_hmmsearch({
            'TIGR04056': ['gene_1'], 'PF07980': ['gene_2']}))
        df = annotate_all_samples(str(tmp_path), 'TIGR04056.hmm', 'PF07980.hmm')
        assert set(df['sample']) == {'S1', 'S2'}
        assert len(df) == 4

    def test_no_samples_returns_empty_with_columns(self, tmp_path):
        df = annotate_all_samples(str(tmp_path), 'TIGR04056.hmm', 'PF07980.hmm')
        assert isinstance(df, pd.DataFrame)
        assert df.empty
        assert list(df.columns) == COLUMNS
