"""B2 figure: low-FPR performance, ArChIPelago (RF) vs best monoPWM. Main set: human-trained RF (seed 0) on human test chr1/8/21.
Reads the final-run tables of 3_analysis/low_fpr_metrics (corrected partial AUC; summaries by b2_summary_w2.py) and writes the
figure here (Figure_S8.pdf/png); its provenance block goes to 3_analysis/low_fpr_metrics/_provenance_p53drop.yaml.
usage: python3 4_visualisation/supplementary_figures/Figure_S8/Figure_S8.py"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
FIG = os.path.dirname(os.path.abspath(__file__))
R2 = os.path.dirname(os.path.dirname(os.path.dirname(FIG)))
os.chdir(os.path.join(R2, "3_analysis", "low_fpr_metrics"))   # the tables are read from the analysis folder
sys.path.insert(0, os.path.join(R2, "4_visualisation", "style"))
sys.path.insert(0, os.path.join(R2, "3_analysis", "integration"))
import archi_style as AS
import p53drop as PD
AS.apply_rc()
plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "xtick.major.width": 0.8, "ytick.major.width": 0.8, "axes.linewidth": 0.8, "legend.fontsize": 7})
PFONT = "Arial Unicode MS"   # Arial lacks the Unicode superscripts used in P values
GREEN, RED, GREY = AS.PALETTE["RF"], AS.PALETTE["monoPWM"], "#999999"
RFN = "ArChIPelago (RF)"
w = pd.read_csv("B2_metrics_wide.csv"); s = pd.read_csv("B2_summary.csv"); pc = pd.read_csv("B2_per_chromosome.csv")
OUT = [os.path.join(FIG, "Figure_S8" + e) for e in (".pdf", ".png")]
for f in OUT:
    PD.backup(f)
H = w[(w.scenario == "HH") & ~w.TF.isin(PD.EXCLUDE_TFS)]
pc = pc[~pc.TF.isin(PD.EXCLUDE_TFS)]
NTF = PD.check_tfs(H.TF, "Fig. S8 HH")
rf = H[H.model == "RF_s0"].set_index("TF"); br = H[H.model == "mono_byroc"].set_index("TF").reindex(rf.index)
bp = H[H.model == "mono_byprc"].set_index("TF").reindex(rf.index)
def pval(m):
    return s[(s.scenario == "HH") & (s.rf_seed == 0) & (s.metric == m)].wilcoxon_p.iloc[0]
def letter(ax, L):
    ax.text(-0.17, 1.10, L, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left")
def paired(ax, m, ylabel, title, L):
    for t in rf.index:
        ax.plot([0, 1], [br.loc[t, m], rf.loc[t, m]], color=GREY, lw=0.6, zorder=1)
    ax.scatter(np.zeros(len(br)), br[m], s=22, marker="^", color=RED, zorder=3, edgecolor="black", linewidth=0.4)
    ax.scatter(np.ones(len(rf)), rf[m], s=22, marker="o", color=GREEN, zorder=3, edgecolor="black", linewidth=0.4)
    for x, v in ((0, br[m]), (1, rf[m])):
        ax.hlines(np.median(v), x - 0.22, x + 0.22, color="black", lw=2, zorder=4)
    ax.set_xlim(-0.5, 1.5); ax.set_xticks([0, 1]); ax.set_xticklabels(["best monoPWM", RFN])
    ax.set_ylabel(ylabel); ax.set_title(title, loc="left", pad=6); letter(ax, L)
    lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi + 0.14 * (hi - lo))
    n = int((rf[m] > br[m]).sum())
    ax.text(0.5, 0.98, "%d of %d TFs higher, %s" % (n, NTF, AS.fmt_p(pval(m))), transform=ax.transAxes, ha="center", va="top",
            fontsize=7, family=PFONT)
    print("A/B", m, "n_higher", n, "median_best", round(br[m].median(), 8), "median_RF", round(rf[m].median(), 8), "P", pval(m))
fig, ax = plt.subplots(2, 2, figsize=(7.4, 7.2))
paired(ax[0, 0], "pauc_std_fpr0.01", "standardized partial auROC (FPR ≤ 0.01)", "Partial auROC, FPR ≤ 0.01", "A")
paired(ax[0, 1], "tpr_at_fpr0.01", "TPR at FPR = 0.01", "True positive rate at a false positive rate of 0.01", "B")
# C
parts = ["auprc_part_r0_0.2", "auprc_part_r0.2_0.4", "auprc_part_r0.4_0.6", "auprc_part_r0.6_0.8", "auprc_part_r0.8_1"]
lab = ["0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"]
x = np.arange(5); a = ax[1, 0]
mr = [rf[p].median() for p in parts]; mb = [bp[p].median() for p in parts]
a.bar(x - 0.19, mb, 0.38, color=RED, edgecolor="black", lw=0.4, label="best monoPWM")
a.bar(x + 0.19, mr, 0.38, color=GREEN, edgecolor="black", lw=0.4, label=RFN)
for i, p in enumerate(parts):
    n = int((rf[p] > bp[p]).sum())
    a.text(i, max(mr[i], mb[i]) + 0.003, "%d/%d" % (n, NTF), ha="center", fontsize=7)
    print("C", lab[i], "median_best", round(mb[i], 8), "median_RF", round(mr[i], 8), "n_higher", n)
a.set_xticks(x); a.set_xticklabels(lab); a.set_xlabel("recall interval"); a.set_ylabel("median contribution to auPRC (area)")
a.set_title("auPRC by recall interval", loc="left", pad=6); letter(a, "C")
a.legend(frameon=False, loc="upper right", handlelength=1.2)
a.set_ylim(0, max(mr + mb) * 1.15)
# D
a = ax[1, 1]; chs = ["chr1", "chr8", "chr21"]
C = pc[pc.scenario == "HH"]
def gain(ch, m, base):
    r = C[(C.chrom == ch) & (C.model == "RF_s0")].set_index("TF")[m]
    b = C[(C.chrom == ch) & (C.model == base)].set_index("TF")[m].reindex(r.index)
    return (r - b).values
# labels above the highest point of the panel, limits from the data
G = {(ch, m): gain(ch, m, base) for ch in chs for m, base in (("auroc", "mono_byroc"), ("auprc", "mono_byprc"))}
gmin, gmax = min(v.min() for v in G.values()), max(v.max() for v in G.values())
YLAB = gmax + 0.03
for i, ch in enumerate(chs):
    for j, (m, base, fill) in enumerate((("auroc", "mono_byroc", True), ("auprc", "mono_byprc", False))):
        g = G[(ch, m)]; xp = i * 3 + j
        a.boxplot([g], positions=[xp], widths=0.8, patch_artist=True, showfliers=False,
                  medianprops=dict(color="black", lw=1.2), whiskerprops=dict(color=GREEN), capprops=dict(color=GREEN),
                  boxprops=dict(edgecolor=GREEN, facecolor=GREEN if fill else "white", lw=1.0))
        jit = np.random.RandomState(2).uniform(-0.2, 0.2, len(g))
        a.scatter(xp + jit, g, s=6, color="black", zorder=3, edgecolor="none")
        a.text(xp, YLAB, "%d/%d" % (int((g > 0).sum()), NTF), ha="center", fontsize=7)
        print("D", ch, m, "median_gain", round(float(np.median(g)), 8), "n_pos", int((g > 0).sum()))
a.axhline(0, color="black", lw=0.6, ls=":")
a.set_xticks([0.5, 3.5, 6.5]); a.set_xticklabels(["%s\n(median %d positives)" % (ch, C[(C.chrom == ch) & (C.model == "RF_s0")].n_pos.median()) for ch in chs], fontsize=7)
a.set_ylabel("Δ (%s minus best monoPWM)" % RFN); a.set_ylim(gmin - 0.07, YLAB + 0.03)
a.set_title("Gain per held-out chromosome", loc="left", pad=6); letter(a, "D")
a.legend(handles=[Patch(facecolor=GREEN, edgecolor=GREEN, label="auROC (filled)"), Patch(facecolor="white", edgecolor=GREEN, label="auPRC (open)")],
         frameon=False, loc="lower left", handlelength=1.2)   # lower left: the lowest points are on chr21 (right)
fig.text(0.5, 0.995, "Human test set (human-trained model)", ha="center", va="top", fontsize=8)
fig.tight_layout(rect=(0, 0, 1, 0.97), h_pad=2.2, w_pad=2.0)
AS.save_both(fig, os.path.join(FIG, "Figure_S8"))
PD.write_provenance(os.path.abspath("."), "4_visualisation/supplementary_figures/Figure_S8/Figure_S8.py (copy of B2_metrics/b2_figure.py)", OUT,
    [os.path.abspath(f) for f in ("B2_metrics_wide.csv", "B2_summary.csv", "B2_per_chromosome.csv")], note="Fig. S8, HH only, 35 TFs")
