"""
Unit tests for substrate.extract_seqs module.

Tests top_family(), _deduplicate(), and count_sequences() using
small synthetic FASTA files written to a temporary directory.
Full sequence extraction functions require dbCAN output and are
covered by the end-to-end test.
"""
import os
import pytest
import tempfile
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from substrate.extract_seqs import (
    top_family,
    is_protein_fasta,
    _deduplicate,
)
from substrate.align import count_sequences


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture
def tmp_fasta(tmp_path):
    """Write a small FASTA file to a temp directory and return its path."""
    def _write(filename, records):
        path = tmp_path / filename
        with open(path, 'w') as f:
            for seq_id, seq in records:
                f.write(f'>{seq_id}\n{seq}\n')
        return str(path)
    return _write


# ── top_family ────────────────────────────────────────────────────────────────

class TestTopFamily:

    def test_simple_family(self):
        """Simple family name returned unchanged."""
        assert top_family('GH16') == 'GH16'

    def test_subfamily_stripped(self):
        """Subfamily suffix stripped to top-level family."""
        assert top_family('GH16_3') == 'GH16'

    def test_pipe_separated_takes_first(self):
        """Pipe-separated subfamilies — first family returned."""
        assert top_family('GH16_3|GH16_5') == 'GH16'

    def test_cbm_family(self):
        assert top_family('CBM48') == 'CBM48'

    def test_pl_family(self):
        assert top_family('PL7_5') == 'PL7'

    def test_dash_returns_unknown(self):
        assert top_family('-') == 'unknown'

    def test_empty_returns_unknown(self):
        assert top_family('') == 'unknown'

    def test_nan_returns_unknown(self):
        assert top_family('nan') == 'unknown'

    def test_none_returns_unknown(self):
        assert top_family(None) == 'unknown'


# ── is_protein_fasta ──────────────────────────────────────────────────────────

class TestIsProteinFasta:

    def test_faa_accepted(self):
        assert is_protein_fasta('sequences.faa') is True

    def test_fasta_accepted(self):
        assert is_protein_fasta('sequences.fasta') is True

    def test_fa_accepted(self):
        assert is_protein_fasta('sequences.fa') is True

    def test_fna_rejected(self):
        assert is_protein_fasta('genome.fna') is False

    def test_txt_rejected(self):
        assert is_protein_fasta('file.txt') is False

    def test_no_extension_rejected(self):
        assert is_protein_fasta('sequences') is False

    def test_case_sensitive(self):
        """Extension check is case-sensitive — .FAA not accepted."""
        assert is_protein_fasta('sequences.FAA') is False


# ── _deduplicate ──────────────────────────────────────────────────────────────

class TestDeduplicate:

    def _make_records(self, ids):
        return [
            SeqRecord(Seq('ACDEFGHIKLM'), id=seq_id, description='')
            for seq_id in ids
        ]

    def test_unique_records_unchanged(self):
        """All unique IDs — all records kept."""
        records = self._make_records(['seq1', 'seq2', 'seq3'])
        result  = _deduplicate(records)
        assert len(result) == 3

    def test_duplicates_removed(self):
        """Duplicate IDs — only first occurrence kept."""
        records = self._make_records(['seq1', 'seq2', 'seq1', 'seq3'])
        result  = _deduplicate(records)
        assert len(result) == 3
        assert [r.id for r in result] == ['seq1', 'seq2', 'seq3']

    def test_all_duplicates(self):
        """All same ID — only one record kept."""
        records = self._make_records(['seq1', 'seq1', 'seq1'])
        result  = _deduplicate(records)
        assert len(result) == 1

    def test_empty_list(self):
        """Empty input returns empty list."""
        assert _deduplicate([]) == []

    def test_first_occurrence_kept(self):
        """First occurrence is always kept, not a later one."""
        records = [
            SeqRecord(Seq('AAAA'), id='seq1', description='first'),
            SeqRecord(Seq('BBBB'), id='seq1', description='second'),
        ]
        result = _deduplicate(records)
        assert len(result) == 1
        assert str(result[0].seq) == 'AAAA'


# ── count_sequences ───────────────────────────────────────────────────────────

class TestCountSequences:

    def test_counts_correct(self, tmp_fasta):
        """Correct number of sequences counted."""
        path = tmp_fasta('test.faa', [
            ('seq1', 'ACDEFGHIKLM'),
            ('seq2', 'MNPQRSTVWY'),
            ('seq3', 'ACDEFGHIKLM'),
        ])
        assert count_sequences(path) == 3

    def test_single_sequence(self, tmp_fasta):
        path = tmp_fasta('single.faa', [('seq1', 'ACDEF')])
        assert count_sequences(path) == 1

    def test_empty_file(self, tmp_fasta):
        path = tmp_fasta('empty.faa', [])
        assert count_sequences(path) == 0

    def test_multiline_sequence_counted_once(self, tmp_path):
        """Sequence split across multiple lines counted as one."""
        path = str(tmp_path / 'multiline.faa')
        with open(path, 'w') as f:
            f.write('>seq1\nACDEF\nGHIKL\nMNPQR\n')
            f.write('>seq2\nSTVWY\n')
        assert count_sequences(path) == 2

    def test_large_file(self, tmp_fasta):
        """Counting works correctly for larger files."""
        records = [(f'seq{i}', 'ACDEF') for i in range(100)]
        path = tmp_fasta('large.faa', records)
        assert count_sequences(path) == 100
