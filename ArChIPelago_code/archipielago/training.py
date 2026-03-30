"""Model training and evaluation helpers for ArChIPelago pipelines."""

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve, precision_recall_curve
from sklearn.model_selection import StratifiedKFold


def train_rf(X_train, y_train, n_estimators=100, max_depth=6,
             max_samples=0.8, n_jobs=40, random_state=42):
    """Train a RandomForestClassifier on the given data.

    Parameters
    ----------
    X_train : array-like or pd.DataFrame
    y_train : array-like
        Binary labels (0 = negative, 1 = positive).
    n_estimators : int
    max_depth : int
    max_samples : float
        Fraction of samples for each tree's bootstrap.
    n_jobs : int
        Parallel jobs. Use -1 for all CPUs.
    random_state : int
        Random seed for reproducibility.

    Returns
    -------
    sklearn.ensemble.RandomForestClassifier
        Fitted model.
    """
    model = RandomForestClassifier(
        n_estimators=n_estimators,
        max_depth=max_depth,
        max_samples=max_samples,
        n_jobs=n_jobs,
        random_state=random_state,
    )
    model.fit(X_train, y_train)
    return model


def evaluate_model(model, X_test, y_test):
    """Compute ROC-AUC and PR-AUC for a fitted model on a test set.

    Parameters
    ----------
    model : fitted sklearn estimator with predict_proba()
    X_test : array-like or pd.DataFrame
    y_test : array-like
        Binary labels.

    Returns
    -------
    dict with keys:
        roc_auc : float
        pr_auc  : float
        fpr     : np.ndarray  (for ROC curve plotting)
        tpr     : np.ndarray
        recall  : np.ndarray  (for PR curve plotting)
        precision : np.ndarray
    """
    y_score = model.predict_proba(X_test)[:, 1]
    y_test = np.array(y_test)

    roc_auc = float(roc_auc_score(y_test, y_score))
    pr_auc = float(average_precision_score(y_test, y_score))
    fpr, tpr, _ = roc_curve(y_test, y_score)
    precision, recall, _ = precision_recall_curve(y_test, y_score)

    return {
        "roc_auc": roc_auc,
        "pr_auc": pr_auc,
        "fpr": fpr,
        "tpr": tpr,
        "precision": precision,
        "recall": recall,
    }


def cross_validate_model(model, X, y, n_splits=7, random_state=42):
    """Stratified k-fold cross-validation returning per-fold ROC-AUC and PR-AUC.

    Parameters
    ----------
    model : sklearn estimator (unfitted; cloned each fold)
    X : array-like or pd.DataFrame
    y : array-like
    n_splits : int
    random_state : int

    Returns
    -------
    dict with keys:
        roc_auc_mean, roc_auc_std,
        pr_auc_mean, pr_auc_std,
        roc_aucs (list), pr_aucs (list)
    """
    from sklearn.base import clone

    X = np.array(X)
    y = np.array(y)

    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    roc_aucs, pr_aucs = [], []

    for train_idx, val_idx in cv.split(X, y):
        m = clone(model)
        m.fit(X[train_idx], y[train_idx])
        y_score = m.predict_proba(X[val_idx])[:, 1]
        roc_aucs.append(float(roc_auc_score(y[val_idx], y_score)))
        pr_aucs.append(float(average_precision_score(y[val_idx], y_score)))

    return {
        "roc_auc_mean": float(np.mean(roc_aucs)),
        "roc_auc_std": float(np.std(roc_aucs)),
        "pr_auc_mean": float(np.mean(pr_aucs)),
        "pr_auc_std": float(np.std(pr_aucs)),
        "roc_aucs": roc_aucs,
        "pr_aucs": pr_aucs,
    }
