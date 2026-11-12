"""PRROC (R, through rpy2) is the primary auROC/auPRC implementation. The NumPy port below, which the v2.0.0 code used
before, is kept only as a reference: it must equal PRROC. Without rpy2 and PRROC the tests ERROR (they never skip)."""
import os
import sys

import numpy as np
from sklearn.metrics import auc, precision_recall_curve, roc_auc_score

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "1_data_preparation"))
import prroc_metrics as PM  # noqa: E402


def auprc_integral_port(y, s):
    """auPRC with Davis-Goadrich continuous interpolation between consecutive distinct-score points (verbatim port)."""
    y = np.asarray(y, float); s = np.asarray(s, float)
    o = np.argsort(-s, kind="mergesort"); y = y[o]; s = s[o]
    last = np.r_[s[1:] != s[:-1], True]                      # end of each tie block
    tp = np.r_[0.0, np.cumsum(y)[last]]; fp = np.r_[0.0, np.cumsum(1 - y)[last]]
    P = tp[-1]
    a, b = tp[:-1], tp[:-1] + fp[:-1]                        # segment start TP, TP+FP
    c, e = np.diff(tp), np.diff(tp) + np.diff(fp)            # dTP, d(TP+FP)  (e > 0 always)
    with np.errstate(divide="ignore", invalid="ignore"):
        seg = np.where(b == 0, c / e, c / e + (a - b * c / e) / e * np.log((b + e) / np.where(b == 0, 1, b)))
    seg = np.where(c == 0, 0.0, seg)
    return float(np.sum(seg * c) / P)


def _random_case(seed=0, n=5000, pos=0.1):
    rng = np.random.RandomState(seed)
    y = (rng.rand(n) < pos).astype(int)
    return y, rng.randn(n) + 1.0 * y


def _tied_case(seed=1, n=5000, pos=0.1):
    rng = np.random.RandomState(seed)
    y = (rng.rand(n) < pos).astype(int)
    return y, np.clip(np.round(rng.randn(n) + 1.0 * y + 2), 0, 5)


def _check(y, s):
    roc, prc = PM.auroc_auprc(y, s)
    assert isinstance(roc, float) and isinstance(prc, float)
    assert abs(prc - auprc_integral_port(y, s)) < 1e-9, (prc, auprc_integral_port(y, s))
    assert abs(roc - roc_auc_score(y, s)) < 1e-12, (roc, roc_auc_score(y, s))
    return roc, prc


def test_random_continuous_scores():
    for seed in range(3):
        _check(*_random_case(seed))


def test_heavily_tied_scores():
    for seed in range(3):
        _check(*_tied_case(seed))


def test_individual_functions_match_pair():
    y, s = _random_case(5)
    assert PM.auroc(y, s) == PM.auroc_auprc(y, s)[0] and PM.auprc(y, s) == PM.auroc_auprc(y, s)[1]


def test_negative_control_trapezoid_differs_on_ties():
    """The sklearn trapezoid auPRC is not the PRROC integral on tied scores, so the equality tests above can fail."""
    y, s = _tied_case()
    p, r, _ = precision_recall_curve(y, s)
    trapz = float(auc(r, p))
    prc = PM.auprc(y, s)
    assert abs(trapz - prc) > 0.01, (trapz, prc)
