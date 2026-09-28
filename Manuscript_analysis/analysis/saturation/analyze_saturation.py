#!/usr/bin/env python3
"""Summary statistics of the PWM-subsampling saturation run (mouse test set = chr1/8/19).
Reads saturation_results.csv of this folder and writes into the same folder
  saturation_summary_by_k.csv   median delta (RF - best single mono PWM) per k, design, metric
  saturation_per_tf.csv         per-TF full gain, saturation k, max-over-k, "hurts" flag
  saturation_stats.md           markdown tables + headline numbers
  fig_saturation.{pdf,png}      2x2 panel figure (diagnostic preview; the shipped Fig. S5 is drawn in R)
Budget convention: a TF with P < k is used at k_eff = min(k, P) (all of its PWMs), so every
median is over all 36 TFs ("what do you get with a budget of k PWMs").
usage: python3 analyze_saturation.py   (no arguments)
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = os.path.dirname(os.path.abspath(__file__))
K_GRID = [1, 2, 4, 8, 16, 32, 64, 128]
METRICS = [("auroc_H", "auROC, human test"), ("auprc_H", "auPRC, human test"),
           ("auroc_M", "auROC, mouse test"), ("auprc_M", "auPRC, mouse test")]
HURT_TOL = 0.005
BLUE, ORANGE, GREY = "#2a78d6", "#eb6834", "#b8b8b8"

RES = "saturation_results.csv"
print("results file:", os.path.join(D, RES))
df = pd.read_csv(os.path.join(D, RES))
tfs = sorted(df.tf.unique())
if len(tfs) != 36:
    raise SystemExit("expected 36 TFs in %s, found %d" % (RES, len(tfs)))
Pmap = df.groupby("tf").P.first()
for m, _ in METRICS:
    df["d_" + m] = df[m] - df["base_" + m]
    df["dsub_" + m] = df["sub_" + m] - df["base_" + m]
dcols = ["d_" + m for m, _ in METRICS]

# mean over reps per (tf, design, k)
mean = df.groupby(["tf", "design", "k"])[dcols + [m for m, _ in METRICS]].mean().reset_index()
full = mean[(mean.design == "random") & (mean.k == mean.tf.map(Pmap))].set_index("tf")   # k = P (2 seeds)


def curve(tf, design):
    """delta at each grid k with carry-forward (k_eff = min(k, P)) plus the k = P value."""
    P = Pmap[tf]
    sub = mean[(mean.tf == tf) & (mean.design == design)].set_index("k")
    out = {}
    for k in K_GRID:
        ke = min(k, P)
        out[k] = full.loc[tf, dcols] if ke == P else sub.loc[ke, dcols]
    out["P"] = full.loc[tf, dcols]
    return pd.DataFrame(out).T.astype(float)


curves = {des: {tf: curve(tf, des) for tf in tfs} for des in ("random", "topk")}
med = {des: pd.concat(curves[des].values(), keys=tfs).groupby(level=1).median().loc[K_GRID + ["P"]]
       for des in ("random", "topk")}
q25 = {des: pd.concat(curves[des].values(), keys=tfs).groupby(level=1).quantile(0.25).loc[K_GRID + ["P"]]
       for des in ("random", "topk")}
q75 = {des: pd.concat(curves[des].values(), keys=tfs).groupby(level=1).quantile(0.75).loc[K_GRID + ["P"]]
       for des in ("random", "topk")}

rows = []
for des in ("random", "topk"):
    for k in K_GRID + ["P"]:
        for m, _ in METRICS:
            rows.append(dict(design=des, k=k, metric=m, median_delta=med[des].loc[k, "d_" + m],
                             q25=q25[des].loc[k, "d_" + m], q75=q75[des].loc[k, "d_" + m],
                             frac_of_full=med[des].loc[k, "d_" + m] / med[des].loc["P", "d_" + m],
                             n_tf_with_P_ge_k=int((Pmap >= (k if k != "P" else 0)).sum())))
byk = pd.DataFrame(rows)
byk.to_csv(os.path.join(D, "saturation_summary_by_k.csv"), index=False)


def first_k(series, thr):
    for k in K_GRID:
        if series[k] >= thr:
            return k
    return "P"


# per-TF table
pt = []
for tf in tfs:
    P = Pmap[tf]
    r = dict(tf=tf, P=P)
    c = curves["random"][tf]
    ct = curves["topk"][tf]
    for m, _ in METRICS:
        g = c.loc["P", "d_" + m]
        r["full_gain_" + m] = g
        ks = [k for k in K_GRID if k < P]
        vals = c.loc[ks, "d_" + m] if ks else pd.Series(dtype=float)
        r["max_over_k_" + m] = max(vals.max() if len(vals) else -np.inf, g)
        r["hurts_" + m] = bool(len(vals) and (vals.max() - g > HURT_TOL))
        r["sat_k_" + m] = (first_k(c["d_" + m], 0.9 * g) if g > 0 else np.nan)
        r["sat_k_topk_" + m] = (first_k(ct["d_" + m], 0.9 * g) if g > 0 else np.nan)
        # top-k minus random at k = 4 and 8 (or P if smaller)
        for k in (4, 8):
            ke = min(k, P)
            r["topk_minus_random_k%d_%s" % (k, m)] = ct.loc[k, "d_" + m] - c.loc[k, "d_" + m] if ke < P else 0.0
    pt.append(r)
pt = pd.DataFrame(pt)
pt.to_csv(os.path.join(D, "saturation_per_tf.csv"), index=False)

# ---- markdown summary
L = []
L.append("### Median delta (RF - best single mono PWM, human-train-selected) vs number of PWMs k\n")
L.append("Budget convention: TFs with P < k enter at k_eff = P; all %d TFs at every k. " % len(tfs) +
         "Values: median over TFs of the per-TF mean over replicates (random design; top-k design in brackets).\n")
L.append("| k | " + " | ".join(lbl for _, lbl in METRICS) + " |")
L.append("|---|" + "---|" * len(METRICS))
for k in K_GRID + ["P"]:
    cells = []
    for m, _ in METRICS:
        a, b = med["random"].loc[k, "d_" + m], med["topk"].loc[k, "d_" + m]
        cells.append("%+.4f [%+.4f]" % (a, b))
    L.append("| %s | " % ("all (P)" if k == "P" else k) + " | ".join(cells) + " |")
L.append("")
L.append("### k at which the median delta reaches 50 % / 90 % of the median full (k = P) delta\n")
L.append("| metric | median full delta | 50 % (random) | 90 % (random) | 50 % (top-k) | 90 % (top-k) |")
L.append("|---|---|---|---|---|---|")
for m, lbl in METRICS:
    g = med["random"].loc["P", "d_" + m]
    L.append("| %s | %+.4f | %s | %s | %s | %s |" % (
        lbl, g, first_k(med["random"]["d_" + m], 0.5 * g), first_k(med["random"]["d_" + m], 0.9 * g),
        first_k(med["topk"]["d_" + m], 0.5 * g), first_k(med["topk"]["d_" + m], 0.9 * g)))
L.append("")
L.append("### Per-TF saturation k (first grid k at which the per-TF mean delta >= 90 % of that TF's full delta; TFs with full delta <= 0 excluded)\n")
L.append("| metric | n TFs with full delta > 0 | median sat. k (random) | IQR | share saturated by k = 8 | share by k = 16 | median sat. k (top-k) | TFs where more PWMs hurt (> %.3f) |" % HURT_TOL)  # noqa
L.append("|---|---|---|---|---|---|---|---|")
for m, lbl in METRICS:
    s = pt["sat_k_" + m].dropna()
    st = pt["sat_k_topk_" + m].dropna()
    num = s.replace("P", 512).astype(float)
    numt = st.replace("P", 512).astype(float)
    med_k = lambda x: int(x.quantile(0.5, interpolation="lower")) if x.quantile(0.5, interpolation="lower") < 512 else "P"
    hurt = pt[pt["hurts_" + m]].tf.tolist()
    L.append("| %s | %d | %s | %g-%g | %.0f %% | %.0f %% | %s | %d/%d (%s) |" % (
        lbl, len(s), med_k(num), num.quantile(0.25, interpolation="lower"), num.quantile(0.75, interpolation="higher"),
        100 * (num <= 8).mean(), 100 * (num <= 16).mean(), med_k(numt),
        len(hurt), len(tfs), ", ".join(hurt) if hurt else "-"))
L.append("\n(sat. k = 512 stands for 'only at k = P')\n")
L.append("### Spearman correlation between P (number of PWMs) and the full-model gain over the best single mono PWM\n")
L.append("| metric | rho | p |")
L.append("|---|---|---|")
for m, lbl in METRICS:
    rho, p = spearmanr(pt.P, pt["full_gain_" + m])
    L.append("| %s | %.3f | %.3g |" % (lbl, rho, p))
L.append("")
L.append("### Top-k (train-auROC-ranked) minus random subsets, median over TFs of the delta difference\n")
L.append("| k | " + " | ".join(lbl for _, lbl in METRICS) + " |")
L.append("|---|" + "---|" * len(METRICS))
for k in K_GRID:
    cells = []
    for m, _ in METRICS:
        dd = [curves["topk"][tf].loc[k, "d_" + m] - curves["random"][tf].loc[k, "d_" + m] for tf in tfs if Pmap[tf] > k]
        cells.append("%+.4f (n=%d)" % (np.median(dd), len(dd)) if dd else "-")
    L.append("| %d | " % k + " | ".join(cells) + " |")
L.append("")
# seed noise at k = P
fp = df[(df.design == "random") & (df.k == df.P)].groupby("tf")[[m for m, _ in METRICS]].agg(lambda x: x.max() - x.min())
L.append("### RF seed-to-seed spread at k = P (|seed0 - seed1|, median / max over TFs)\n")
L.append("| metric | median | max |")
L.append("|---|---|---|")
for m, lbl in METRICS:
    L.append("| %s | %.4f | %.4f |" % (lbl, fp[m].median(), fp[m].max()))
L.append("")
# random-subset best single PWM vs RF on the subset (does the RF add value at small k?)
L.append("### At k = 1 and k = 2, RF on the subset vs the subset's best single PWM (median over TFs of RF - single)\n")
L.append("| k | " + " | ".join(lbl for _, lbl in METRICS) + " |")
L.append("|---|" + "---|" * len(METRICS))
for k in (1, 2, 4):
    sub = df[(df.design == "random") & (df.k == k)]
    cells = ["%+.4f" % (sub[m] - sub["sub_" + m]).groupby(sub.tf).mean().median() for m, _ in METRICS]
    L.append("| %d | " % k + " | ".join(cells) + " |")
L.append("")
n_fits = len(df)
L.append("Fits: %d rows (%d random-subset, %d top-k, %d full); total RF fit+predict time %.1f CPU-worker hours.\n" % (
    n_fits, ((df.design == "random") & (df.k < df.P)).sum(), (df.design == "topk").sum(),
    ((df.design == "random") & (df.k == df.P)).sum(), df.fit_sec.sum() / 3600))
open(os.path.join(D, "saturation_stats.md"), "w").write("\n".join(L))
print("\n".join(L))

# ---- figure
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "ps.fonttype": 42})
fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.0), sharex=True)
XP = 512  # x position of the "all P" point for the median curves
for ax, (m, lbl) in zip(axes.ravel(), METRICS):
    ax.axhline(0, color="#777777", lw=0.8, ls=":", zorder=1)
    for tf in tfs:
        P = Pmap[tf]
        c = curves["random"][tf]["d_" + m]
        xs = [k for k in K_GRID if k < P] + [P]
        ys = [c[k] for k in K_GRID if k < P] + [c["P"]]
        ax.plot(xs, ys, color=GREY, lw=0.7, alpha=0.9, zorder=2)
    for des, col, ls, lab in (("random", BLUE, "-", "median, random k-subsets"),
                              ("topk", ORANGE, "--", "median, top-k by train auROC")):
        y = med[des]["d_" + m]
        ax.plot(K_GRID + [XP], [y[k] for k in K_GRID] + [y["P"]], color=col, lw=2.2, ls=ls, zorder=4, label=lab)
        ax.plot([XP], [y["P"]], marker="o", color=col, ms=5, zorder=5)
    ax.set_xscale("log", base=2)
    ax.set_xticks(K_GRID + [XP])
    ax.set_xticklabels([str(k) for k in K_GRID] + ["all\n(P)"])
    ax.tick_params(axis="x", which="minor", bottom=False)
    ax.set_title(lbl, fontsize=9.5, loc="left")
    ax.grid(axis="y", color="#e5e5e5", lw=0.6)
for ax in axes[1]:
    ax.set_xlabel("number of PWMs used (k)")
for ax, (m, _) in zip(axes.ravel(), METRICS):
    ax.set_ylabel("$\\Delta$%s vs best single mono-PWM" % ("auROC" if m.startswith("auroc") else "auPRC"))
h, l = axes[0, 0].get_legend_handles_labels()
h.append(plt.Line2D([], [], color=GREY, lw=0.7))
l.append("single TF (mean of replicates), ends at its P")
fig.legend(h, l, loc="lower center", ncol=3, frameon=False, fontsize=8, bbox_to_anchor=(0.5, -0.01))
fig.suptitle("RF gain over the best single mono-PWM vs number of PWMs used\n(%d TFs, human-trained; P = %d-%d PWMs per TF, median %d)"
             % (len(tfs), Pmap.min(), Pmap.max(), int(Pmap.median())), fontsize=9.5, y=0.995)
fig.tight_layout(rect=(0, 0.04, 1, 0.95))
fig.savefig(os.path.join(D, "fig_saturation.pdf"))
fig.savefig(os.path.join(D, "fig_saturation.png"), dpi=200)
print("figure written")
