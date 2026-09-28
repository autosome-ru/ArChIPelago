"""Operational relevance of the ArChIPelago gain (Sup. Table 4).

The baseline is the best monoPWM selected
by TRAINING auPRC on the human training set (score vector pwmH_byprc_test_<H|M> of the hm_mm npz), compared
with the human-trained RF (seed 0, rf_<H|M>_s0).  Ranking is tie-aware: every block of
equal scores is a unit and the expectation under a random order inside the block is taken (linear
interpolation of the PR curve through a tie block):

  * false positives to reach recall r: whole blocks above the block containing the r-th positive,
    plus the proportional share of that block;
  * precision at recall r = r*P / (r*P + FP);
  * true positives among the top N: whole blocks plus the proportional share of the block cut by N.

Reads per-TF score vectors saved by the hm_mm run (<scores_dir>/<TF>.npz; test_M = mouse chr1/8/19).
Usage: python operational_metrics.py <scores_dir> <out_csv>   (scores_dir = --scores-dir of hm_mm/run_hm_mm.py)
"""
import sys, glob, os
import numpy as np, pandas as pd


def tie_blocks(y, s):
    """Blocks of equal score in descending score order: (block sizes, positives per block)."""
    y = np.asarray(y).astype(np.int64); s = np.asarray(s, float)
    o = np.argsort(-s, kind="stable")
    s, y = s[o], y[o]
    start = np.flatnonzero(np.r_[True, s[1:] != s[:-1]])
    sizes = np.diff(np.r_[start, len(s)])
    pos = np.add.reduceat(y, start)
    return sizes, pos


def fp_at_recall(sizes, pos, r):
    """Expected false positives accepted to reach recall r under a random order inside tie blocks."""
    P = pos.sum(); target = r * P
    cum_pos = np.cumsum(pos); cum_n = np.cumsum(sizes)
    b = np.searchsorted(cum_pos, target)          # first block whose cumulative positives reach target
    n0 = cum_n[b] - sizes[b]; p0 = cum_pos[b] - pos[b]
    frac = (target - p0) / pos[b]                 # share of the block needed (0 < frac <= 1)
    fp = (n0 - p0) + frac * (sizes[b] - pos[b])
    return fp, target


def tp_top(sizes, pos, n_top):
    """Expected positives among the top n_top ranked regions (proportional share of the cut block)."""
    cum_n = np.cumsum(sizes); cum_pos = np.cumsum(pos)
    if n_top >= cum_n[-1]:
        return float(cum_pos[-1])
    b = np.searchsorted(cum_n, n_top)             # block that contains rank n_top
    n0 = cum_n[b] - sizes[b]; p0 = cum_pos[b] - pos[b]
    return p0 + (n_top - n0) / sizes[b] * pos[b]


def metrics(y, s, n_top=1000):
    y = np.asarray(y).astype(int)
    N = len(y); P = int(y.sum())
    sizes, pos = tie_blocks(y, s)
    out = {}
    for r in (0.5, 0.8):
        fp, tp = fp_at_recall(sizes, pos, r)
        out[f"fp_at_recall{r}_per10k"] = fp / N * 1e4
        out[f"precision_at_recall{r}"] = tp / (tp + fp)
    out["tp_top1000"] = tp_top(sizes, pos, n_top)
    k = max(1, int(round(0.01 * N)))
    out["tp_top1pct"] = tp_top(sizes, pos, k)
    out["n"] = N; out["n_pos"] = P
    out["frac_tied_scores"] = 1 - len(sizes) / N
    return out


def main(scores_dir, out_csv):
    if os.path.exists(out_csv):
        raise SystemExit("refusing to overwrite existing %s" % out_csv)
    keys = ["fp_at_recall0.5_per10k", "fp_at_recall0.8_per10k", "precision_at_recall0.5",
            "precision_at_recall0.8", "tp_top1000", "tp_top1pct"]
    rows = []
    for f in sorted(glob.glob(os.path.join(scores_dir, "*.npz"))):
        tf = os.path.basename(f)[:-4]
        z = np.load(f)  # our own npz of numeric arrays; no pickled objects
        for sp, yk, rk, bk in (("H", "y_test_H", "rf_H_s0", "pwmH_byprc_test_H"),
                               ("M", "y_test_M", "rf_M_s0", "pwmH_byprc_test_M")):
            y = z[yk]
            m_rf, m_mono = metrics(y, z[rk]), metrics(y, z[bk])
            row = {"tf": tf, "test": sp, "n": m_rf["n"], "n_pos": m_rf["n_pos"],
                   "baseline_pwm": str(z["pwmH_byprc_col"]),
                   "frac_tied_scores_rf": m_rf["frac_tied_scores"], "frac_tied_scores_mono": m_mono["frac_tied_scores"]}
            for k in keys:
                row[f"{k}_rf"] = m_rf[k]; row[f"{k}_mono"] = m_mono[k]
            rows.append(row)
    df = pd.DataFrame(rows)
    df["fp_reduction_recall0.5_pct_vs_mono"] = 100 * (1 - df["fp_at_recall0.5_per10k_rf"] / df["fp_at_recall0.5_per10k_mono"])
    df["fp_reduction_recall0.8_pct_vs_mono"] = 100 * (1 - df["fp_at_recall0.8_per10k_rf"] / df["fp_at_recall0.8_per10k_mono"])
    df["extra_tp_top1000_vs_mono"] = df["tp_top1000_rf"] - df["tp_top1000_mono"]
    df["precision_gain_recall0.5_vs_mono"] = df["precision_at_recall0.5_rf"] - df["precision_at_recall0.5_mono"]
    df.to_csv(out_csv, index=False)
    for sp in ("H", "M"):
        d = df[df.test == sp]
        print(f"\n=== test {sp}, baseline = best monoPWM by train auPRC (human train) (n TFs = {len(d)})")
        for c in ["fp_at_recall0.5_per10k_mono", "fp_at_recall0.5_per10k_rf", "fp_reduction_recall0.5_pct_vs_mono",
                  "fp_reduction_recall0.8_pct_vs_mono", "precision_at_recall0.5_mono", "precision_at_recall0.5_rf",
                  "precision_gain_recall0.5_vs_mono", "tp_top1000_mono", "tp_top1000_rf", "extra_tp_top1000_vs_mono"]:
            print(f"{c:40s} median {d[c].median():9.3f}   IQR {d[c].quantile(.25):8.3f}-{d[c].quantile(.75):8.3f}")
        print("TFs with fewer FP at recall 0.5:", int((d["fp_reduction_recall0.5_pct_vs_mono"] > 0).sum()), "/", len(d),
              "| at recall 0.8:", int((d["fp_reduction_recall0.8_pct_vs_mono"] > 0).sum()))


if __name__ == "__main__":
    # self-checks: perfect ranking; all-tied ranking = random expectation; untied = one block per region
    rng = np.random.default_rng(0); y = (rng.random(10000) < 0.01).astype(int); P = y.sum()
    m = metrics(y, y + rng.random(10000) * 0.1); assert m["fp_at_recall0.5_per10k"] == 0 and m["precision_at_recall0.5"] == 1
    m = metrics(y, np.zeros(10000))              # one tie block: recall r needs r*N regions
    assert abs(m["fp_at_recall0.5_per10k"] - (0.5 * 10000 - 0.5 * P) / 10000 * 1e4) < 1e-9
    assert abs(m["tp_top1000"] - 1000 * P / 10000) < 1e-9
    s = rng.random(10000)                        # no ties: block per region
    sizes, pos = tie_blocks(y, s); assert len(sizes) == 10000 and pos.sum() == P
    main(sys.argv[1], sys.argv[2])
