"""B4 step 6 (local): multi-panel figure from 3_analysis/feature_importance/tables/B4_*.csv -> Figure_S10.{pdf,png} in this folder. No transparency, black text, no MathText."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, NullFormatter

HERE = os.path.dirname(os.path.abspath(__file__)); R2 = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
B = os.path.join(R2, "3_analysis", "feature_importance")   # analysis folder: tables in, provenance out
T = os.path.join(B, "tables")
sys.path.insert(0, os.path.join(R2, "4_visualisation", "style"))
sys.path.insert(0, os.path.join(R2, "3_analysis", "integration"))
import p53drop as PD
OUT = [os.path.join(HERE, "Figure_S10" + e) for e in (".pdf", ".png")]
for _f in OUT:
    PD.backup(_f)
import archi_style as ST
ST.apply_rc()
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
                     "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8, "xtick.major.size": 3, "ytick.major.size": 3})
C_ = pd.read_csv(os.path.join(T, "B4_importance_concentration.csv")).sort_values(["n_features", "TF"]).reset_index(drop=True)
order = list(C_.TF); x = np.arange(len(order))
PD.check_tfs(order, "Fig. S10 concentration")
SI = pd.read_csv(os.path.join(T, "B4_subtype_importance.csv"))
SUM = pd.read_csv(os.path.join(T, "B4_subtype_summary_per_TF.csv")).set_index("TF")
AS = pd.read_csv(os.path.join(T, "B4_ablation_summary.csv"))
BP = pd.read_csv(os.path.join(T, "B4_activation_best_pwm.csv"))
GR = pd.read_csv(os.path.join(T, "B4_activation_groups.csv"))
G50, G80, GREY, DKGREEN, GREEN = "#000000", "#555555", "#999999", "#00503B", "#009E73"


def flat_ticks(ax, labels, size=6.5):
    """Horizontal tick labels, every second one a line lower (no tilted labels)."""
    ax.set_xticklabels([("\n" + l) if i % 2 else l for i, l in enumerate(labels)], rotation=0, fontsize=size, linespacing=1.1)
    for i, t in enumerate(ax.xaxis.get_major_ticks()):
        t.tick1line.set_markersize(3 if i % 2 == 0 else 10)

fig = plt.figure(figsize=(10, 11.9))   # height kept near 1.15 x width so that figure and legend fit one page of the supplement
gs = fig.add_gridspec(4, 6, height_ratios=[1, 1, 1.75, 1.65], hspace=0.74, wspace=1.25)
axA = fig.add_subplot(gs[0, :]); axB = fig.add_subplot(gs[1, :])
gsC = gs[2, :].subgridspec(2, 1, hspace=0.1)
axC1 = fig.add_subplot(gsC[0]); axC2 = fig.add_subplot(gsC[1], sharex=axC1)
axD = fig.add_subplot(gs[3, 0:3]); axE = fig.add_subplot(gs[3, 3:6])

# A: concentration
axA.bar(x, C_.n_features, color="#CCCCCC", edgecolor="black", linewidth=0.4, label="all PWMs (mono + di)")
axA.plot(x, C_.n_for_80, "s", color=G80, ms=5, markeredgecolor="black", markeredgewidth=0.4, label="PWMs carrying 80 % of Gini importance")
axA.plot(x, C_.n_for_50, "o", color=G50, ms=5, markeredgecolor="black", markeredgewidth=0.4, label="PWMs carrying 50 % of Gini importance")
axA.set_yscale("log"); axA.yaxis.set_major_formatter(FuncFormatter(lambda v, p: "%g" % v)); axA.yaxis.set_minor_formatter(NullFormatter()); axA.set_ylabel("Number of PWMs (log scale)")
axA.set_xticks(x); flat_ticks(axA, order); axA.set_xlim(-0.7, len(x) - 0.3)
axA.legend(frameon=False, loc="upper left")
axA.set_title("Gini importance concentration", loc="left", pad=6)

# B: importance per subtype at the cut of the ablation (Methods 2.10: 0.8, but 0.9 for the TFs that split only at 0.9); the cut per TF is read from the ablation summary, other TFs at 0.8
_cut = AS.groupby("TF").sim_threshold.first()
if (AS.groupby("TF").sim_threshold.nunique() != 1).any():
    raise SystemExit("CHECK FAILED ablation summary: a TF with several cuts")
CUT = {tf: float(_cut.get(tf, 0.8)) for tf in order}
S8 = pd.concat([SI[(SI.TF == tf) & (SI.sim_threshold == c)] for tf, c in CUT.items()])
if set(S8.TF) != set(order):
    raise SystemExit("CHECK FAILED Fig. S10B: TFs without subtype rows at their cut")
cols = ["#00503B", "#009E73", "#B7E3D3", "#999999"]
bot = np.zeros(len(order))
for lab, col, fn in (("subtype 1 (primary)", cols[0], lambda g: g[g.subtype == 1].gini_share.sum()),
                     ("subtype 2", cols[1], lambda g: g[g.subtype == 2].gini_share.sum()),
                     ("subtype 3", cols[2], lambda g: g[g.subtype == 3].gini_share.sum()),
                     ("subtypes 4 and higher", cols[3], lambda g: g[g.subtype >= 4].gini_share.sum())):
    v = np.array([fn(S8[S8.TF == tf]) for tf in order])
    axB.bar(x, v, bottom=bot, color=col, edgecolor="black", linewidth=0.4, label=lab); bot += v
axB.set_ylabel("Share of Gini importance"); axB.set_ylim(0, 1.0)
axB.set_xticks(x); flat_ticks(axB, [tf + "*" if CUT[tf] != 0.8 else tf for tf in order]); axB.set_xlim(-0.7, len(x) - 0.3)
_c9 = sorted(set(CUT[tf] for tf in order if CUT[tf] != 0.8))
axB.text(1.0, 1.02, "* subtypes at cut %s (as in the ablation); other TFs at 0.8" % " / ".join("%.1f" % c for c in _c9), transform=axB.transAxes,
         ha="right", va="bottom", fontsize=6.5, color="black")
axB.legend(frameon=False, ncol=4, loc="lower left", bbox_to_anchor=(0, 1.0))
n2 = int((SUM.n_ge10pct_80 >= 2).sum())
axB.set_title("Gini importance per PWM subtype", loc="left", pad=20)

# C: ablation of the secondary subtype, 29 TFs, two strips (auROC above, auPRC below) sharing the TF axis
CUTS = AS.groupby("TF").sim_threshold.first()
sec_share = {tf: SI[(SI.TF == tf) & (SI.sim_threshold == CUTS[tf]) & (SI.subtype == 2)].gini_share.sum() for tf in CUTS.index}
groups = [sorted([t for t in CUTS.index if CUTS[t] == c], key=lambda t: -sec_share[t]) for c in sorted(CUTS.unique())]
sec_order = [t for g in groups for t in g]
N8 = len(groups[0]); xi = {t: i for i, t in enumerate(sec_order)}
AF = pd.read_csv(os.path.join(T, "B4_ablation_all_fits.csv"))
full = AF[AF.cond == "full"].set_index(["TF", "seed"])
ASi = AS.set_index(["TF", "cond"])
capped = sorted(t for t in sec_order if ASi.loc[(t, "drop_random"), "n_removed"] < ASi.loc[(t, "drop_secondary"), "n_removed"])
OFF = 0.17
plotdata, pts = [], {}
for ax, m, mname in ((axC1, "auroc", "auROC"), (axC2, "auprc", "auPRC")):
    pts[ax] = []
    for tf in sec_order:
        i = xi[tf]
        r = AF[(AF.TF == tf) & (AF.cond == "drop_random")].sort_values("seed")
        dr = r[m].values - full.loc[[(tf, s_) for s_ in r.seed], m].values
        mr, sdr = ASi.loc[(tf, "drop_random"), "d_%s_mean" % m], ASi.loc[(tf, "drop_random"), "d_%s_sd" % m]
        ms_ = ASi.loc[(tf, "drop_secondary"), "d_%s_mean" % m]
        lim = mr - 3 * sdr; below = bool(ms_ < lim)
        xc, xs = i - OFF, i + OFF
        ax.plot(xc + (r.seed.values - 4.5) * 0.012, dr, "o", color="#AAAAAA", ms=2.3, mew=0, zorder=2)
        ax.plot([xc, xc], [mr, lim], "-", color="black", lw=0.9, zorder=3); ax.plot([xc - 0.07, xc + 0.07], [lim, lim], "-", color="black", lw=0.9, zorder=3)
        ax.plot([xc], [mr], "s", mfc="white", mec="black", mew=0.9, ms=4.6, zorder=4)
        ax.plot([xs], [ms_], "o", mfc=ST.RF_COL if below else "white", mec=ST.RF_COL, mew=1.1, ms=5.2, zorder=5)
        pts[ax] += [(xc, v) for v in list(dr) + [mr, lim]] + [(xs, ms_)]
        plotdata.append(dict(TF=tf, cut=CUTS[tf], metric=mname, order=i + 1, secondary_share=sec_share[tf], capped_control=tf in capped,
                             random_draws_mean=mr, random_draws_sd=sdr, limit_mean_minus_3sd=lim, secondary_removed_mean=ms_, beyond_limit=below,
                             random_draws=";".join("%.5g" % v for v in dr)))
    ax.axhline(0, color="black", lw=0.6, ls=":", zorder=1)
    ax.axvline(N8 - 0.5, color="black", lw=0.6)
    ax.set_xlim(-0.7, len(sec_order) - 0.3)
    ax.margins(y=0.06)
axC1.set_ylabel("\u0394auROC"); axC2.set_ylabel("\u0394auPRC")
for ax in (axC1, axC2):
    ax.set_xticks(range(len(sec_order)))
flat_ticks(axC2, [t + "†" if t in capped else t for t in sec_order], 6.5)   # dagger: the asterisk of panel B marks the cut 0.9
axC1.tick_params(axis="x", labelbottom=False)
for g, c, x0 in ((groups[0], 0.8, (N8 - 1) / 2), (groups[1], 0.9, N8 + (len(groups[1]) - 1) / 2)):
    axC1.annotate("cut %.1f" % c, (x0, 1.0), xycoords=("data", "axes fraction"), xytext=(0, 3), textcoords="offset points", ha="center", va="bottom", fontsize=7, color="black")
axC1.set_title("Subtype ablation (retrained minus full model)", loc="left", pad=18)
axC1.annotate("Human test set", (0, 1.0), xycoords="axes fraction", xytext=(0, 34), textcoords="offset points", ha="left", va="bottom", fontsize=8, color="black")
axC1.annotate("† fewer random PWMs removed than secondary PWMs", (1.0, 1.0), xycoords="axes fraction", xytext=(0, 18), textcoords="offset points",
              ha="right", va="bottom", fontsize=6.5, color="black")
hl = [Line2D([0], [0], marker="o", ls="", color="#AAAAAA", ms=2.8, label="10 random removals"),
      Line2D([0], [0], marker="s", ls="", mfc="white", mec="black", mew=0.9, ms=4.6, label="mean of random removals"),
      Line2D([0], [0], marker="|", ls="-", color="black", lw=0.9, ms=6, label="mean minus 3 SD"),
      Line2D([0], [0], marker="o", ls="", mfc=ST.RF_COL, mec=ST.RF_COL, ms=5.2, label="secondary removed, beyond limit"),
      Line2D([0], [0], marker="o", ls="", mfc="white", mec=ST.RF_COL, mew=1.1, ms=5.2, label="secondary removed, within limit")]
fig.canvas.draw()
def _clear(ax, leg):
    bb = leg.get_window_extent(fig.canvas.get_renderer()).expanded(1.02, 1.05)
    return not any(bb.contains(*ax.transData.transform(p)) for p in pts[ax])
leg = None
for ax in (axC1, axC2):
    for loc in ("lower right", "upper right", "lower left", "center right", "upper left", "center left"):
        if leg is not None: leg.remove()
        leg = ax.legend(handles=hl, loc=loc, ncol=1, frameon=False, fontsize=6.5, handletextpad=0.3, labelspacing=0.25, borderaxespad=0.3)
        fig.canvas.draw()
        if _clear(ax, leg):
            break
    else:
        leg.remove(); leg = None; continue
    break
if leg is None:
    raise SystemExit("CHECK FAILED Fig. S10C: no free place for the legend")
print("legend of panel C in", "upper" if ax is axC1 else "lower", "strip, loc", loc)
filled = {m: [r["TF"] for r in plotdata if r["metric"] == m and r["beyond_limit"]] for m in ("auROC", "auPRC")}
print("panel C filled:", {k: sorted(v) for k, v in filled.items()})
pd.DataFrame(plotdata).to_csv(os.path.join(HERE, "Figure_S10C_plotdata.csv"), index=False, float_format="%.6g")

# D: fraction of strong-hit sequences whose best PWM is non-primary
B8 = BP[BP.thr == 0.8]
pos = B8[B8.cls == "pos"].set_index("TF"); neg = B8[B8.cls == "neg"].set_index("TF")
tf_d = [t for t in pos.index if SUM.loc[t, "n_subtypes_80"] >= 2]
tf_d = sorted(tf_d, key=lambda t: pos.loc[t, "frac_best_nonprimary_strong"])
y = np.arange(len(tf_d))
axD.plot(pos.loc[tf_d, "frac_features_nonprimary"], y, "|", color="black", ms=7, label="share of non-primary PWMs")
axD.plot(neg.loc[tf_d, "frac_best_nonprimary_strong"], y, "o", mfc="white", mec="black", mew=0.4, ms=4.5, label="negatives")
axD.plot(pos.loc[tf_d, "frac_best_nonprimary_strong"], y, "o", color=DKGREEN, mec="black", mew=0.4, ms=4.5, label="positives")
axD.set_yticks(y); axD.set_yticklabels(tf_d, fontsize=6.5); axD.set_xlim(-0.02, 1.02)
axD.set_xlabel("Fraction with best PWM in\na non-primary subtype")
axD.set_title("Best PWM per strong-hit sequence", loc="left", pad=8)
axD.legend(frameon=False, loc="lower right", fontsize=6.5, handletextpad=0.2, borderaxespad=0.2)

# E: RF score of positives with a strong hit only in primary vs only in secondary subtypes
G8 = GR[(GR.defn == "any_nonprimary") & (GR.cls == "pos")]
po = G8[G8.group == "primary_only"].set_index("TF"); so = G8[G8.group == "secondary_only"].set_index("TF")
tf_e = [t for t in so.index if so.loc[t, "n"] >= 20 and po.loc[t, "n"] >= 20]
axE.plot([0, 1], [0, 1], ":", color="black", lw=0.8)
axE.scatter(po.loc[tf_e, "median_rf"], so.loc[tf_e, "median_rf"], s=26, color=GREEN, edgecolors="black", linewidths=0.4, zorder=3)
axE.set_xlim(0, 1); axE.set_ylim(0, 1)
axE.set_xlabel("Median RF score\n(strong hit in primary subtype only)"); axE.set_ylabel("Median RF score\n(strong hit in non-primary subtypes only)")
axE.set_title("RF score by subtype of strong hit", loc="left", pad=8)
axE.set_aspect("equal")
axE.text(0.03, 0.97, "partly expected by construction:\nthe primary subtype carries\nmost of the importance", transform=axE.transAxes,
         ha="left", va="top", fontsize=6, color="black")
for a, t in zip((axA, axB, axC1, axD, axE), "ABCDE"):
    a.text(-0.085 if a in (axA, axB, axC1) else (-0.4 if a is axE else -0.25), 1.12 if a is axA else (1.30 if a is axB else (1.55 if a is axC1 else 1.10)), t, transform=a.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left")
for a in (axD, axE):   # test-set label: panels D and E use the human test set (C has its own)
    a.text(0, 1.17, "Human test set", transform=a.transAxes, fontsize=8, ha="left", va="bottom", color="black")
fig.savefig(OUT[0], bbox_inches="tight")
fig.savefig(OUT[1], dpi=600, bbox_inches="tight")
print("saved; panel C TFs", len(sec_order), "panel D TFs", len(tf_d), "panel E TFs", len(tf_e))
PD.write_provenance(B, "4_visualisation/supplementary_figures/Figure_S10/Figure_S10.py", OUT + [os.path.join(HERE, "Figure_S10C_plotdata.csv")],
                    [os.path.join(T, n) for n in ("B4_importance_concentration.csv", "B4_subtype_importance.csv", "B4_subtype_summary_per_TF.csv",
                     "B4_ablation_summary.csv", "B4_activation_best_pwm.csv", "B4_activation_groups.csv", "B4_ablation_all_fits.csv")], note="Fig. S10")
