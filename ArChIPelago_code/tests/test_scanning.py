"""Tests for archipielago.scanning module."""

import numpy as np
import pandas as pd
import pytest
from archipielago.scanning import build_feature_matrix, select_top_features


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
