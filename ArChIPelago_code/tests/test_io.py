"""Tests for archipielago.io module."""

import io
import pytest
from archipielago.io import fasta_iter, load_fasta, save_fasta, make_train_test_beds


# ---------------------------------------------------------------------------
# fasta_iter
# ---------------------------------------------------------------------------

def test_fasta_iter_basic(simple_fasta_handle):
    records = list(fasta_iter(simple_fasta_handle))
    assert len(records) == 3
    assert records[0] == ("seq1", "ACGTACGTACGT")
    assert records[1] == ("seq2", "TTTTAAAACCCC")
    assert records[2] == ("seq3", "GGGGCCCCTTTT")


def test_fasta_iter_multiline_seq(multiline_fasta_handle):
    records = list(fasta_iter(multiline_fasta_handle))
    assert len(records) == 2
    # Multi-line sequences should be joined into one string
    assert records[0] == ("seq1", "ACGTACGT")
    assert records[1] == ("seq2", "TTTTAAAA")


def test_fasta_iter_empty():
    records = list(fasta_iter(io.StringIO("")))
    assert records == []


def test_fasta_iter_strips_whitespace():
    fasta = ">header with spaces  \nACGT  \nTTTT\n"
    records = list(fasta_iter(io.StringIO(fasta)))
    assert records[0][0] == "header with spaces"
    assert records[0][1] == "ACGTTTTT"


# ---------------------------------------------------------------------------
# load_fasta / save_fasta round-trip
# ---------------------------------------------------------------------------

def test_load_fasta(tmp_fasta):
    records = load_fasta(tmp_fasta)
    assert len(records) == 3
    assert records[0] == ("seq1", "ACGTACGTACGT")


def test_load_fasta_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_fasta(tmp_path / "nonexistent.fasta")


def test_save_fasta_roundtrip(tmp_path):
    original = [("seqA", "AAAA"), ("seqB", "CCCC")]
    out = tmp_path / "out.fasta"
    save_fasta(original, out)
    loaded = load_fasta(out)
    assert loaded == original


# ---------------------------------------------------------------------------
# make_train_test_beds
# ---------------------------------------------------------------------------

def test_make_train_test_beds_basic(peaks_df):
    train_chr = ["chr1", "chr2", "chr3"]
    test_chr = ["chr8", "chr21"]
    train, test = make_train_test_beds(peaks_df, train_chr, test_chr)

    assert set(train["chrom"].unique()).issubset(set(train_chr))
    assert set(test["chrom"].unique()).issubset(set(test_chr))
    # chr8 has 2 peaks, chr21 has 1 → 3 test rows
    assert len(test) == 3


def test_make_train_test_beds_excludes_unlisted_chroms(peaks_df):
    train, test = make_train_test_beds(peaks_df, ["chr1"], ["chr21"])
    # chr2, chr3, chr8 should not appear in either split
    assert "chr2" not in train["chrom"].values
    assert "chr2" not in test["chrom"].values


def test_make_train_test_beds_cap(peaks_df):
    """Downsampling: n_max smaller than available rows."""
    train, test = make_train_test_beds(
        peaks_df,
        train_chr=["chr1", "chr2", "chr3"],
        test_chr=["chr8", "chr21"],
        n_max=2,
        random_state=0,
    )
    assert len(train) <= 2
    assert len(test) <= 2


def test_make_train_test_beds_resets_index(peaks_df):
    train, _ = make_train_test_beds(peaks_df, ["chr1", "chr2"], ["chr21"])
    assert list(train.index) == list(range(len(train)))
