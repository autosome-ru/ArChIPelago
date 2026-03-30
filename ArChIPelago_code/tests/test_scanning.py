"""Tests for archipielago.scanning module."""

import numpy as np
import pandas as pd
import pytest
from archipielago.scanning import (
    build_feature_matrix,
    load_sarus_scores,
    run_sarus,
    select_top_features,
)


# ---------------------------------------------------------------------------
# build_feature_matrix
# ---------------------------------------------------------------------------

def test_build_feature_matrix_mono(mock_sarus_scan_dir):
    X = build_feature_matrix(mock_sarus_scan_dir, mode="mono")
    # 3 mono PWMs → 3 columns; 4 sequences → 4 rows
    assert X.shape == (4, 3)
    assert all(col.startswith("mono_") for col in X.columns)


def test_build_feature_matrix_di(mock_sarus_scan_dir):
    X = build_feature_matrix(mock_sarus_scan_dir, mode="di")
    # 2 di PWMs → 2 columns
    assert X.shape == (4, 2)
    assert all(col.startswith("di_") for col in X.columns)


def test_build_feature_matrix_mono_di(mock_sarus_scan_dir):
    X = build_feature_matrix(mock_sarus_scan_dir, mode="mono_di")
    # 3 mono + 2 di = 5 columns
    assert X.shape == (4, 5)
    mono_cols = [c for c in X.columns if c.startswith("mono_")]
    di_cols = [c for c in X.columns if c.startswith("di_")]
    assert len(mono_cols) == 3
    assert len(di_cols) == 2


def test_build_feature_matrix_no_nan(mock_sarus_scan_dir):
    X = build_feature_matrix(mock_sarus_scan_dir, mode="mono_di")
    assert not X.isnull().any().any()


def test_build_feature_matrix_values(mock_sarus_scan_dir):
    """Check that scores from conftest fixtures land in the right column."""
    X = build_feature_matrix(mock_sarus_scan_dir, mode="mono")
    # pwm_0.txt has scores [0.1, 0.2, 0.3, 0.4]
    col = "mono_pwm_0"
    assert col in X.columns
    assert pytest.approx(X[col].tolist()) == [0.1, 0.2, 0.3, 0.4]


def test_build_feature_matrix_missing_dir(tmp_path):
    with pytest.raises(FileNotFoundError):
        build_feature_matrix(tmp_path, mode="mono")


def test_build_feature_matrix_empty_dir(tmp_path):
    (tmp_path / "mono").mkdir()
    with pytest.raises(ValueError):
        build_feature_matrix(tmp_path, mode="mono")


# ---------------------------------------------------------------------------
# select_top_features
# ---------------------------------------------------------------------------

def test_select_top_features_returns_column_names(small_feature_matrix):
    X, y = small_feature_matrix
    selected = select_top_features(X, y, n=3)
    assert len(selected) == 3
    assert all(col in X.columns for col in selected)


def test_select_top_features_respects_n(small_feature_matrix):
    X, y = small_feature_matrix
    # Request more features than available — should cap at X.shape[1]
    selected = select_top_features(X, y, n=1000)
    assert len(selected) == X.shape[1]


def test_select_top_features_no_duplicates(small_feature_matrix):
    X, y = small_feature_matrix
    selected = select_top_features(X, y, n=5)
    assert len(selected) == len(set(selected))


# ---------------------------------------------------------------------------
# load_sarus_scores
# ---------------------------------------------------------------------------

def test_load_sarus_scores_basic(tmp_path):
    score_file = tmp_path / "scores.txt"
    score_file.write_text("1.5\n2.3\n0.7\n")
    s = load_sarus_scores(score_file)
    assert list(s) == pytest.approx([1.5, 2.3, 0.7])
    assert s.index.tolist() == [0, 1, 2]


def test_load_sarus_scores_empty_file(tmp_path):
    score_file = tmp_path / "empty.txt"
    score_file.write_text("")
    s = load_sarus_scores(score_file)
    assert isinstance(s, pd.Series)
    assert len(s) == 0


def test_load_sarus_scores_missing_file(tmp_path):
    s = load_sarus_scores(tmp_path / "does_not_exist.txt")
    assert isinstance(s, pd.Series)
    assert len(s) == 0


# ---------------------------------------------------------------------------
# run_sarus — validation-only (no Java invoked)
# ---------------------------------------------------------------------------

def test_run_sarus_missing_fasta(tmp_path):
    pwm_file = tmp_path / "motif.pwm"
    pwm_file.write_text("dummy")
    jar_file = tmp_path / "sarus.jar"
    jar_file.write_text("dummy")

    with pytest.raises(FileNotFoundError, match="FASTA"):
        run_sarus(
            fasta_path=tmp_path / "missing.fasta",
            pwm_path=pwm_file,
            sarus_jar=jar_file,
            output_path=tmp_path / "out.txt",
        )


def test_run_sarus_missing_jar(tmp_path):
    fasta_file = tmp_path / "seqs.fasta"
    fasta_file.write_text(">seq1\nACGT\n")
    pwm_file = tmp_path / "motif.pwm"
    pwm_file.write_text("dummy")

    with pytest.raises(FileNotFoundError, match="SARUS jar"):
        run_sarus(
            fasta_path=fasta_file,
            pwm_path=pwm_file,
            sarus_jar=tmp_path / "missing_sarus.jar",
            output_path=tmp_path / "out.txt",
        )
