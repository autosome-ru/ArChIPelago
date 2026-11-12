"""B3: summary tables and the multi-panel figure. Tables are read from and written to 3_analysis/top_peaks, the figure
(Figure_S9.pdf/png) goes to this folder. usage: python3 4_visualisation/supplementary_figures/Figure_S9/Figure_S9.py
Reads the server outputs (b3_topN_results.csv, b3_peak_widths_per_dataset.csv, b3_join_test_*.csv) and b3_gain_vs_npos_*.csv.
Human (H) rows exclude p53drop.EXCLUDE_TFS (35 TFs); peak widths are restricted to the used peak sets
(human 361, mouse 276); mouse (M) top-N and join rows are the first-submission mouse data, marked stale.
The figure (Fig. S9) is redrawn from these tables; its mouse panels A-D show the stale mouse data."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import wilcoxon, spearmanr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, NullFormatter
PLAIN = FuncFormatter(lambda v, p: ("%g" % v) if v < 1e4 else format(int(v), ","))
def plain_log(ax, axis):
    a = ax.xaxis if axis == "x" else ax.yaxis
    a.set_major_formatter(PLAIN); a.set_minor_formatter(NullFormatter())

import os, sys
FIG = os.path.dirname(os.path.abspath(__file__)); R2 = os.path.dirname(os.path.dirname(os.path.dirname(FIG)))
os.chdir(os.path.join(R2, "3_analysis", "top_peaks"))  # table inputs and outputs in the analysis folder
sys.path.insert(0, os.path.join(R2, "4_visualisation", "style")); import archi_style as AS
sys.path.insert(0, os.path.join(R2, "3_analysis", "integration")); import p53drop as PD
MSTATUS = "current: rebuilt mouse data, 35 TFs"
OUTS = ["b3_topN_summary.csv", "b3_peak_widths_summary.csv", "b3_join_summary.csv", os.path.join(FIG, "Figure_S9.pdf"), os.path.join(FIG, "Figure_S9.png")]
for f in OUTS:
    PD.backup(f)
RF, PWM, GREY = AS.PALETTE["RF"], AS.PALETTE["monoPWM"], "#B0B0B0"
PFONT = "Arial Unicode MS"  # Arial lacks Unicode superscript and subscript digits
d = pd.read_csv("b3_topN_results.csv"); d["N"] = d.N.astype(str)
d = d[~((d.species == "H") & d.TF.isin(PD.EXCLUDE_TFS))]
PD.check_tfs(d[d.species == "H"].TF, "b3_topN_results H"); PD.check_tfs(d[d.species == "M"].TF, "b3_topN_results M")
NS = ["500", "1000", "2000"]

def wp(x):
    x = np.asarray(x, float); x = x[~np.isnan(x)]
    return wilcoxon(x).pvalue if len(x) > 5 and np.any(x != 0) else np.nan

rows = []
for sp in "HM":
    for var in ("allneg", "r100", "rorig"):
        base = d[(d.species == sp) & (d.variant == var)]
        allr = base[(base.selection == "top") & (base.N == "all")].set_index("TF")
        npos = allr.n_pos
        for subset in ("avail", "common_2000"):
            for N in NS + ["all"]:
                top = base[(base.selection == "top") & (base.N == N)].set_index("TF")
                if subset == "common_2000":
                    top = top[top.index.isin(npos[npos >= 2000].index)]
                if len(top) == 0: continue
                rnd = base[(base.selection == "random") & (base.N == N)].set_index("TF").reindex(top.index) if N != "all" else None
                a = allr.reindex(top.index)
                r = dict(species=sp, variant=var, subset=subset, N=N, n_tf=len(top), median_n_neg_per_pos=top.neg_per_pos.median(),
                         rf_auroc_med=top.rf_auroc.median(), pwm_auroc_med=top.pwm_auroc.median(), rf_auprc_med=top.rf_auprc.median(), pwm_auprc_med=top.pwm_auprc.median(),
                         d_auroc_med=top.d_auroc.median(), d_auprc_med=top.d_auprc.median(),
                         n_improved_auroc=int((top.d_auroc > 0).sum()), n_improved_auprc=int((top.d_auprc > 0).sum()),
                         p_wilcoxon_auroc=wp(top.d_auroc), p_wilcoxon_auprc=wp(top.d_auprc),
                         d_auroc_med_allpos_same_tfs=a.d_auroc.median(), d_auprc_med_allpos_same_tfs=a.d_auprc.median(),
                         topN_minus_allpos_auroc_med=(top.d_auroc - a.d_auroc).median(), topN_minus_allpos_auprc_med=(top.d_auprc - a.d_auprc).median(),
                         p_topN_vs_allpos_auroc=wp(top.d_auroc - a.d_auroc), p_topN_vs_allpos_auprc=wp(top.d_auprc - a.d_auprc))
                if rnd is not None:
                    r.update(d_auroc_med_random=rnd.d_auroc.median(), d_auprc_med_random=rnd.d_auprc.median(),
                             top_minus_random_auroc_med=(top.d_auroc - rnd.d_auroc).median(), top_minus_random_auprc_med=(top.d_auprc - rnd.d_auprc).median(),
                             p_top_vs_random_auroc=wp(top.d_auroc - rnd.d_auroc), p_top_vs_random_auprc=wp(top.d_auprc - rnd.d_auprc))
                rows.append(r)
summ = pd.DataFrame(rows)
summ["status"] = np.where(summ.species == "H", "current: 35 TFs, P53 excluded", MSTATUS)
summ.to_csv("b3_topN_summary.csv", index=False)

# peak widths
w = pd.read_csv("b3_peak_widths_per_dataset.csv")
used = pd.concat([PD.human_ids().assign(species="H"), PD.mouse_ids().assign(species="M")]).rename(columns={"peak_id": "peak_set"})
w = w.merge(used, on=["peak_set", "TF", "species"], how="inner")
if len(w) != len(used) or w.peak_set.duplicated().any():
    raise SystemExit("CHECK FAILED peak widths: %d rows matched of %d used peak sets" % (len(w), len(used)))
PD.check_tfs(w[w.species == "H"].TF, "peak widths H"); PD.check_tfs(w[w.species == "M"].TF, "peak widths M")
print("peak widths: %d datasets (H %d, M %d)" % (len(w), (w.species == "H").sum(), (w.species == "M").sum()))
wr = []
for sp, g in list(w.groupby("species")) + [("all", w)]:
    gg = g.dropna(subset=["width_median_top1000"])
    wr.append(dict(species=sp, n_datasets=len(g), n_peaks_total=int(g.n_peaks.sum()), width_median_of_dataset_medians=g.width_median.median(),
                   width_q25_of_dataset_medians=g.width_median.quantile(.25), width_q75_of_dataset_medians=g.width_median.quantile(.75),
                   n_datasets_median_width_le300=int((g.width_median <= 300).sum()),
                   frac_peaks_width_le300_pooled=float((g.frac_width_le300 * g.n_peaks).sum() / g.n_peaks.sum()),
                   frac_cover_ge50_median_over_datasets=g.frac_window_covers_ge50.median(), frac_cover_ge50_pooled_over_peaks=float((g.frac_window_covers_ge50 * g.n_peaks).sum() / g.n_peaks.sum()),
                   frac_cover_ge50_dataset_min=g.frac_window_covers_ge50.min(),
                   n_datasets_ge1000_peaks=len(gg), width_median_top1000_median=gg.width_median_top1000.median(), width_median_rest_median=gg.width_median_rest.median(),
                   n_datasets_top1000_wider=int((gg.width_median_top1000 > gg.width_median_rest).sum()),
                   median_ratio_top1000_to_rest=(gg.width_median_top1000 / gg.width_median_rest).median(),
                   p_wilcoxon_top1000_vs_rest=wilcoxon(gg.width_median_top1000, gg.width_median_rest).pvalue,
                   cover_ge50_top1000_median=gg.cover_ge50_top1000.median(), cover_ge50_rest_median=gg.cover_ge50_rest.median()))
pd.DataFrame(wr).to_csv("b3_peak_widths_summary.csv", index=False)

# join test
j = pd.read_csv("b3_join_test_per_tf.csv"); js = []
j = j[~((j.species == "H") & j.TF.isin(PD.EXCLUDE_TFS))]
PD.check_tfs(j[j.species == "H"].TF, "join H"); PD.check_tfs(j[j.species == "M"].TF, "join M")
for sp in "HM":
    g = j[j.species == sp]
    js.append(dict(species=sp, n_test_pos=int(g.n_test_pos.sum()), n_matched=int((g.n_matched_unique_width + g.n_matched_ambiguous_width).sum()),
                   n_ambiguous_width=int(g.n_matched_ambiguous_width.sum()), match_rate=float((g.n_matched_unique_width + g.n_matched_ambiguous_width).sum() / g.n_test_pos.sum()),
                   n_tf_full_match=int((g.match_rate == 1).sum()), min_tf_match=g.match_rate.min(), median_tf_match=g.match_rate.median()))
ps = pd.read_csv("b3_join_test_per_sequence.csv"); h = ps[(ps.species == "H") & ps.width.notna() & ~ps.TF.isin(PD.EXCLUDE_TFS)]
js[0].update(width_median_scored_positives=h.width.median(), frac_scored_width_le300=float((h.width <= 300).mean()), frac_scored_width_le600=float((h.width <= 600).mean()))
js = pd.DataFrame(js); js["status"] = np.where(js.species == "H", "current: 35 TFs, P53 excluded", MSTATUS)
js.to_csv("b3_join_summary.csv", index=False)
print(summ[(summ.variant.isin(["allneg", "rorig"]))].round(4).to_string()); print(pd.DataFrame(wr).round(3).T.to_string()); print(js.T.to_string())

# ---------------- figure
AS.apply_rc()
plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
                     "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8, "xtick.major.size": 3, "ytick.major.size": 3})
OL = dict(edgecolor="black", linewidths=0.4)
fig = plt.figure(figsize=(10.0, 8.4))
gs = fig.add_gridspec(3, 12, hspace=0.62, wspace=2.8, left=0.065, right=0.985, top=0.95, bottom=0.075)
def lab(ax, t, dx=-0.26): ax.text(dx, 1.14, t, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left")

order = ["all", "2000", "1000", "500"]
for k, (sp, met, ttl) in enumerate([("H", "auroc", "Human test set, auROC"), ("H", "auprc", "Human test set, auPRC"), ("M", "auroc", "Mouse test set, auROC"), ("M", "auprc", "Mouse test set, auPRC")]):
    ax = fig.add_subplot(gs[0, 3 * k:3 * k + 3])
    base = d[(d.species == sp) & (d.variant == "allneg") & (d.selection == "top")]
    for gi, N in enumerate(order):
        t = base[base.N == N].set_index("TF")
        xr, xp = gi + 0.17, gi - 0.17
        for tf in t.index:
            ax.plot([xp, xr], [t.loc[tf, "pwm_" + met], t.loc[tf, "rf_" + met]], color=GREY, lw=0.6, zorder=1)
        ax.scatter(np.full(len(t), xp), t["pwm_" + met], s=14, marker="^", color=PWM, zorder=3, **OL)
        ax.scatter(np.full(len(t), xr), t["rf_" + met], s=14, marker="o", color=RF, zorder=3, **OL)
    ax.set_xticks(range(4))
    ax.set_xticklabels([f"all\n({int(summ[(summ.species==sp)&(summ.variant=='allneg')&(summ['subset']=='avail')&(summ.N=='all')].n_tf.iloc[0])} TFs)"] + [f"top {n}\n({int(summ[(summ.species==sp)&(summ.variant=='allneg')&(summ['subset']=='avail')&(summ.N==n)].n_tf.iloc[0])} TFs)" for n in order[1:]], fontsize=6.5)
    ax.set_ylabel("auROC" if met == "auroc" else "auPRC"); ax.set_title(ttl, pad=6)
    lab(ax, "ABCD"[k])
    if k == 0:
        ax.scatter([], [], s=16, marker="o", color=RF, label="RF (monoPWMs+diPWMs)", **OL); ax.scatter([], [], s=16, marker="^", color=PWM, label="best monoPWM", **OL)
        ax.legend(frameon=False, loc="lower left", handletextpad=0.2, borderpad=0.1)
    lo_ = 0 if met == "auprc" else 0.6   # head room above 1.0 so that marks at 1.00 are drawn whole; no tick above 1.0
    ax.set_ylim(lo_, 1.0)
    ax.set_yticks([t_ for t_ in ax.get_yticks() if lo_ - 1e-9 <= t_ <= 1.0 + 1e-9])
    ax.set_ylim(lo_, 1.0 + 0.025 * (1.0 - lo_))
    ax.set_xlim(-0.6, 3.6)

g = pd.read_csv("b3_gain_vs_npos_per_tf.csv"); sp_ = pd.read_csv("b3_gain_vs_npos_spearman.csv")
specs = [("n_train_pos_H", "dauroc_H", "Training positives (human)", "\u0394auROC, human test set", True), ("n_train_pos_H", "dauprc_H", "Training positives (human)", "\u0394auPRC, human test set", True),
         ("n_datasets_H", "dauroc_H", "ChIP-Seq datasets (human)", "\u0394auROC, human test set", False), ("n_datasets_H", "dauprc_H", "ChIP-Seq datasets (human)", "\u0394auPRC, human test set", False)]
axB = []
for k, (x, y, xl, yl, logx) in enumerate(specs):
    ax = fig.add_subplot(gs[1, 3 * k:3 * k + 3]); axB.append(ax)
    cap = g.cap_H.values
    ax.scatter(g[x][~cap], g[y][~cap], s=16, color=RF, label="below cap", zorder=3, **OL)
    ax.scatter(g[x][cap], g[y][cap], s=20, marker="s", facecolor="white", edgecolor="black", linewidths=0.7, label="at 10,000 cap", zorder=3)
    ax.axhline(0, color="black", lw=0.6, ls=":")
    if logx: ax.set_xscale("log"); plain_log(ax, "x")
    r = sp_[(sp_.gain == y) & (sp_.predictor_col == x) & (sp_.test_set == "human test")].iloc[0]
    t = "\u03c1 = %.2f, %s" % (r.rho, AS.fmt_p(r.p))
    if x == "n_train_pos_H": t += "\n(excluding capped TFs: \u03c1 = %.2f)" % r.rho_uncapped
    ax.text(0.04, 0.96, t, transform=ax.transAxes, va="top", fontsize=7, family=PFONT)
    print("B", y, x, "rho", round(r.rho, 4), "p", r.p, "rho_uncapped", round(r.rho_uncapped, 4))
    ax.set_xlabel(xl); ax.set_ylabel(yl)
    ax.set_ylim(min(g[y].min(), 0) - 0.003, g[y].max() * 1.3)
    lab(ax, "EFGH"[k])
fig.canvas.draw()
hh, ll = axB[0].get_legend_handles_labels()
yl_ = axB[0].get_position().y1
axB[2].legend(hh, ll, loc="upper right", frameon=False, handletextpad=0.2, borderpad=0.1)

# C
H, M = w[w.species == "H"], w[w.species == "M"]
ax = fig.add_subplot(gs[2, 0:4]); lab(ax, "I", -0.16)
data = [H.width_median.values, M.width_median.values]
bp = ax.boxplot(data, widths=0.5, showfliers=False, patch_artist=True, medianprops=dict(color="black", lw=1.4))
for b in bp["boxes"]: b.set(facecolor="#CCCCCC", edgecolor="black", lw=0.8)
for wk in bp["whiskers"] + bp["caps"]: wk.set(color="black", lw=0.8)
for i, v in enumerate(data):
    ax.scatter(np.random.RandomState(0).normal(i + 1, 0.07, len(v)), v, s=6, color="#999999", edgecolors="black", linewidths=0.3, zorder=3)
ax.axhline(300, color="black", lw=0.8, ls=":"); ax.text(1.5, 330, "300 bp\nwindow", fontsize=7, ha="center", va="bottom", linespacing=1.0)   # between the two point clouds
ax.set_xticks([1, 2]); ax.set_xticklabels([f"Human\n({len(H)} datasets)", f"Mouse\n({len(M)} datasets)"]); ax.set_ylabel("Median peak width per dataset (bp)")
ax.set_title("Peak width of source datasets", pad=6)

ax = fig.add_subplot(gs[2, 4:8]); lab(ax, "J")
for sp, xo in (("H", 0), ("M", 2.2)):
    gg = w[(w.species == sp)].dropna(subset=["width_median_top1000"])
    for a, b in zip(gg.width_median_rest, gg.width_median_top1000):
        ax.plot([xo, xo + 1], [a, b], color=GREY, lw=0.4, zorder=1)
    ax.scatter(np.full(len(gg), xo), gg.width_median_rest, s=5, color="#999999", edgecolors="black", linewidths=0.2, zorder=3)
    ax.scatter(np.full(len(gg), xo + 1), gg.width_median_top1000, s=5, color="#999999", edgecolors="black", linewidths=0.2, zorder=3)
    ax.text(xo + 0.5, 1.02, "Human" if sp == "H" else "Mouse", transform=ax.get_xaxis_transform(), ha="center", fontsize=7)
ax.set_yscale("log"); plain_log(ax, "y"); ax.axhline(300, color="black", lw=0.8, ls=":")
ax.set_xticks([0, 1, 2.2, 3.2]); ax.set_xticklabels(["rest", "top 1,000", "rest", "top 1,000"], fontsize=7)
ax.set_ylabel("Median peak width (bp)")
ax.set_title("Top 1,000 by MACS significance vs rest", pad=14, family=PFONT, fontsize=8)

ax = fig.add_subplot(gs[2, 8:12]); lab(ax, "K")
data = [H.frac_window_covers_ge50.values * 100, M.frac_window_covers_ge50.values * 100]
bp = ax.boxplot(data, widths=0.5, showfliers=False, patch_artist=True, medianprops=dict(color="black", lw=1.4))
for b in bp["boxes"]: b.set(facecolor="#CCCCCC", edgecolor="black", lw=0.8)
for wk in bp["whiskers"] + bp["caps"]: wk.set(color="black", lw=0.8)
for i, v in enumerate(data):
    ax.scatter(np.random.RandomState(1).normal(i + 1, 0.07, len(v)), v, s=6, color="#999999", edgecolors="black", linewidths=0.3, zorder=3)
ax.set_xticks([1, 2]); ax.set_xticklabels(["Human", "Mouse"]); ax.set_ylabel("300 bp window covering\nat least half of the peak (%)")
ax.set_title("Coverage of the peak by the 300 bp window", pad=6)
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(FIG, "Figure_S9." + ext), dpi=300 if ext == "png" else None)
PD.write_provenance(os.path.abspath("."), "4_visualisation/supplementary_figures/Figure_S9/Figure_S9.py (copy of B3_peaks/b3_summary_figure.py)", [os.path.abspath(f) for f in OUTS],
    [os.path.abspath(f) for f in ("b3_topN_results.csv", "b3_peak_widths_per_dataset.csv", "b3_join_test_per_tf.csv", "b3_join_test_per_sequence.csv",
                                  "b3_gain_vs_npos_per_tf.csv", "b3_gain_vs_npos_spearman.csv")],
    note="W2: H rows 35 TFs; widths 361 + 276 used peak sets; M rows from the rebuilt mouse data (B11 W1G1/b3)")
