"""Tests for archipielago.training module."""

import numpy as np
import pytest
from sklearn.ensemble import RandomForestClassifier

from archipielago.training import train_rf, evaluate_model, cross_validate_model


# ---------------------------------------------------------------------------
# train_rf
# ---------------------------------------------------------------------------

def test_train_rf_returns_fitted_model(small_feature_matrix):
    X, y = small_feature_matrix
    model = train_rf(X, y, n_estimators=10, n_jobs=1, random_state=0)
    assert isinstance(model, RandomForestClassifier)
    assert hasattr(model, "estimators_")  # fitted models have estimators_
    assert model.n_features_in_ == X.shape[1]


def test_train_rf_reproducible(small_feature_matrix):
    X, y = small_feature_matrix
    m1 = train_rf(X, y, n_estimators=10, n_jobs=1, random_state=7)
    m2 = train_rf(X, y, n_estimators=10, n_jobs=1, random_state=7)
    pred1 = m1.predict_proba(X)[:, 1]
    pred2 = m2.predict_proba(X)[:, 1]
    np.testing.assert_array_equal(pred1, pred2)


# ---------------------------------------------------------------------------
# evaluate_model
# ---------------------------------------------------------------------------

def test_evaluate_model_keys(small_feature_matrix):
    X, y = small_feature_matrix
    model = train_rf(X, y, n_estimators=10, n_jobs=1, random_state=42)
    results = evaluate_model(model, X, y)
    for key in ("roc_auc", "pr_auc", "fpr", "tpr", "precision", "recall"):
        assert key in results


def test_evaluate_model_roc_in_range(small_feature_matrix):
    X, y = small_feature_matrix
    model = train_rf(X, y, n_estimators=10, n_jobs=1, random_state=42)
    results = evaluate_model(model, X, y)
    assert 0.0 <= results["roc_auc"] <= 1.0
    assert 0.0 <= results["pr_auc"] <= 1.0


def test_evaluate_model_curve_arrays(small_feature_matrix):
    X, y = small_feature_matrix
    model = train_rf(X, y, n_estimators=10, n_jobs=1, random_state=42)
    results = evaluate_model(model, X, y)
    assert len(results["fpr"]) == len(results["tpr"])
    assert len(results["precision"]) == len(results["recall"])


def test_evaluate_model_perfect_separation():
    """Linearly separable data should give ROC-AUC == 1.0 (train == test)."""
    X = np.array([[1.0], [2.0], [3.0], [-1.0], [-2.0], [-3.0]])
    y = np.array([1, 1, 1, 0, 0, 0])
    model = train_rf(X, y, n_estimators=50, n_jobs=1, random_state=42)
    results = evaluate_model(model, X, y)
    assert results["roc_auc"] == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# cross_validate_model
# ---------------------------------------------------------------------------

def test_cross_validate_model_keys(small_feature_matrix):
    X, y = small_feature_matrix
    model = RandomForestClassifier(n_estimators=5, random_state=0)
    cv_results = cross_validate_model(model, X, y, n_splits=3)
    for key in ("roc_auc_mean", "roc_auc_std", "pr_auc_mean", "pr_auc_std",
                "roc_aucs", "pr_aucs"):
        assert key in cv_results


def test_cross_validate_model_fold_count(small_feature_matrix):
    X, y = small_feature_matrix
    model = RandomForestClassifier(n_estimators=5, random_state=0)
    cv_results = cross_validate_model(model, X, y, n_splits=3)
    assert len(cv_results["roc_aucs"]) == 3
    assert len(cv_results["pr_aucs"]) == 3


def test_cross_validate_model_mean_in_range(small_feature_matrix):
    X, y = small_feature_matrix
    model = RandomForestClassifier(n_estimators=5, random_state=0)
    cv_results = cross_validate_model(model, X, y, n_splits=3)
    assert 0.0 <= cv_results["roc_auc_mean"] <= 1.0
    assert cv_results["roc_auc_std"] >= 0.0
