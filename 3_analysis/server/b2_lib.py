"""Low-FPR metrics from saved per-sequence scores. Tie-aware: curves are built at distinct score
thresholds and joined by straight segments (= random tie-breaking in expectation). Python 3.8 compatible."""
import os
import sys
import numpy as np
from sklearn.metrics import roc_auc_score  # partial auROC (McClish) only; full auROC / auPRC come from PRROC

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "..", "1_data_preparation"))
import prroc_metrics  # noqa: E402


def curve(y, s):
    """Cumulative (tp, fp) at the end of each tie block, with (0,0) first."""
    y = np.asarray(y, float); s = np.asarray(s, float)
    o = np.argsort(-s, kind="mergesort"); y = y[o]; s = s[o]
    last = np.r_[s[1:] != s[:-1], True]
    tp = np.r_[0.0, np.cumsum(y)[last]]
    fp = np.r_[0.0, np.cumsum(1 - y)[last]]
    return tp, fp


def recall_interval_areas(tp, fp, edges):
    """Area under the Davis-Goadrich PR curve inside each recall interval [edges[i], edges[i+1]]."""
    P = tp[-1]
    a, b = tp[:-1], tp[:-1] + fp[:-1]
    c, e = np.diff(tp), np.diff(tp) + np.diff(fp)
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        t1 = np.clip(lo * P - a, 0, c); t2 = np.clip(hi * P - a, 0, c)
        dt = t2 - t1
        k = np.where(c > 0, e / np.where(c > 0, c, 1), 1.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            f = lambda t: t / k + (a - b / k) / k * np.log(b + k * t)
            seg = np.where(b == 0, dt / k, f(t2) - f(t1))
        seg = np.where(dt <= 0, 0.0, seg)
        out.append(float(np.sum(seg) / P))
    return out


def roc_points(tp, fp):
    return fp / fp[-1], tp / tp[-1]


def tpr_at_fpr(fpr, tpr, x):
    return float(np.interp(x, fpr, tpr))


def fpr_at_tpr(fpr, tpr, t):
    """smallest FPR at which the (piecewise linear) ROC reaches TPR t."""
    i = int(np.searchsorted(tpr, t - 1e-12, side="left"))
    if i == 0:
        return 0.0
    f0, f1, t0, t1 = fpr[i - 1], fpr[i], tpr[i - 1], tpr[i]
    return float(f0 + (f1 - f0) * (t - t0) / (t1 - t0))


def partial_auc(fpr, tpr, x, y, s):
    """raw = exact area under the tie-aware ROC for FPR <= x: every curve point with FPR <= x in curve order (vertical
        steps kept, so they add no area) plus the point interpolated at FPR = x on the segment that crosses x;
        std = sklearn roc_auc_score(y, s, max_fpr=x) (McClish)."""
    k = int(np.searchsorted(fpr, x, side="right"))          # fpr is non-decreasing; points 0..k-1 have FPR <= x
    xs, ys = fpr[:k], tpr[:k]
    if k < len(fpr) and xs[-1] < x:                          # segment (k-1, k) crosses x
        f0, f1, t0, t1 = fpr[k - 1], fpr[k], tpr[k - 1], tpr[k]
        xs = np.r_[xs, x]; ys = np.r_[ys, t0 + (t1 - t0) * (x - f0) / (f1 - f0)]
    raw = float(np.trapz(ys, xs))
    return raw, float(roc_auc_score(y, s, max_fpr=x))


FPR_CUTS = [0.001, 0.01, 0.05, 0.1]
TPR_AT = [0.001, 0.01, 0.05]
FPR_AT_TPR = [0.2, 0.4, 0.6, 0.8, 1.0]
EDGES = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
WIN_KB = 0.3   # negative windows are 300 bp


def all_metrics(y, s):
    tp, fp = curve(y, s)
    P, N = tp[-1], fp[-1]
    fpr, tpr = roc_points(tp, fp)
    m = {"n_pos": P, "n_neg": N}
    m["auroc"], m["auprc"] = prroc_metrics.auroc_auprc(y, s)  # PRROC roc.curve()$auc and pr.curve()$auc.integral
    for x in FPR_CUTS:
        raw, std = partial_auc(fpr, tpr, x, y, s)
        m["pauc_raw_fpr%g" % x] = raw
        m["pauc_std_fpr%g" % x] = std
    for x in TPR_AT:
        m["tpr_at_fpr%g" % x] = tpr_at_fpr(fpr, tpr, x)
    for t in FPR_AT_TPR:
        f = fpr_at_tpr(fpr, tpr, t)
        m["fpr_at_tpr%g" % t] = f
        m["fp_per_kb_at_tpr%g" % t] = f / WIN_KB
        m["prec_at_tpr%g" % t] = float(t * P / (t * P + f * N)) if (t * P + f * N) > 0 else float("nan")
    for x in TPR_AT:
        m["fp_per_kb_at_fpr%g" % x] = x / WIN_KB
    ar = recall_interval_areas(tp, fp, EDGES)
    for lo, hi, v in zip(EDGES[:-1], EDGES[1:], ar):
        m["auprc_part_r%g_%g" % (lo, hi)] = v
    return m
