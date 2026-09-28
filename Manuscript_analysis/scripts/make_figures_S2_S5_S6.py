#!/usr/bin/env python3
"""Draw Fig. S2 and build the SOURCE DATA of the supplementary figures S5 and S6 (mouse test set = mouse
chr1/8/19; baseline = best single monoPWM selected on the training set).

Fig. S2  -- Random Forest vs the best single monoPWM (selected on the human training set), human test
            (left) and mouse test (right), for monoPWM / diPWM / monoPWM+diPWM models; from
            Manuscript_analysis/HUMAN_MOUSE_total_100k.csv (matplotlib; panel letters A, B = auROC row,
            C, D = auPRC row of each 2 x 2 block). Written to Figures/panels/Figure_S2.pdf.
Fig. S5  -- performance vs the number of PWMs used (subsampling experiment).
            From analysis/saturation/saturation_results.csv (reference = best single monoPWM, the
            base_* columns of that file, asserted against the results table).
Fig. S6  -- cross-species transfer and the mouse-trained control.
            From Sup. Table 5 (Sup_Tables/Sup_Table_5_*.csv).
The matplotlib versions of S5 and S6 are previews only (Figures/previews/, not tracked); Fig. S5 and S6 are
drawn in R by Figure_S5_saturation.R / Figure_S6_cross_species.R from the source-data CSVs written here to
Figures/source_data/.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

from archi_paths import BASIS, RES

OUT = os.path.join(RES, "Figures", "previews")
PREV = OUT
FIGS = os.path.join(RES, "Figures", "panels")
os.makedirs(PREV, exist_ok=True)
SRC = os.path.join(RES, "Figures", "source_data")
os.makedirs(OUT, exist_ok=True)
os.makedirs(SRC, exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"]})
COL = {"mono": "#d1478f", "di": "#5b8ff9", "mono+di": "#1f9e63"}
LBL = {"mono": "monoPWMs", "di": "diPWMs", "mono+di": "monoPWMs+diPWMs"}

TABLE = pd.read_csv(os.path.join(RES, "HUMAN_MOUSE_total_100k.csv"), sep="\t")
RF = TABLE[TABLE.Model == "RandomForestClassifier"]


# ==============================================================================================
# Figure S2
# ==============================================================================================
def panel_scatter(ax, d, xcol, ycol, xlab, ylab, by_pwm=True, label_tfs=None):
    if by_pwm:
        for pwm in ["mono", "di", "mono+di"]:
            s = d[d.PWM == pwm]
            ax.scatter(s[xcol], s[ycol], s=11, c=COL[pwm], alpha=.85, lw=0, label=LBL[pwm])
    else:
        s = d[d.PWM == "mono+di"]
        delta = s[ycol] - s[xcol] if "Delta" not in ycol else s[ycol]
        ax.scatter(s[xcol], s[ycol], s=16, c=np.where(delta < 0, "#111111", "#e6a100"), alpha=.85, lw=0)
        if label_tfs is not None:
            for _, r in s.iterrows():
                if r.TF_name in label_tfs:
                    ax.annotate(r.TF_name, (r[xcol], r[ycol]), fontsize=5.5, xytext=(2.5, 2.5),
                                textcoords="offset points", color="#333333")
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)


def figure_s2():
    fig, axes = plt.subplots(4, 4, figsize=(12.2, 10.5))
    src_rows = []
    for block in (0, 1):                                  # 0 = coloured by PWM set, 1 = mono+di labelled
        for col_block, sp, sp_lab in ((0, "H", "human"), (2, "M", "mouse")):
            for row_in, met, mlab in ((0, "roc_auc", "auROC"), (1, "pr_auc", "auPRC")):
                r = block * 2 + row_in
                x = f"{met}_test_{sp}_PWM"
                y = f"{met}_test_{sp}"
                d = RF.copy()
                d["Delta"] = d[y] - d[x]
                if block == 0:
                    panel_scatter(axes[r][col_block], d, x, y, f"{mlab} best monoPWM", mlab)
                    panel_scatter(axes[r][col_block + 1], d, x, "Delta", f"{mlab} best monoPWM", f"$\\Delta$ {mlab}")
                else:
                    s = d[d.PWM == "mono+di"]
                    lab = set(s.nlargest(6, "Delta").TF_name) | set(s.nsmallest(4, "Delta").TF_name)
                    panel_scatter(axes[r][col_block], d, x, y, f"{mlab} best monoPWM", mlab, by_pwm=False, label_tfs=lab)
                    panel_scatter(axes[r][col_block + 1], d, x, "Delta", f"{mlab} best monoPWM", f"$\\Delta$ {mlab}",
                                  by_pwm=False, label_tfs=lab)
                for c in (col_block, col_block + 1):
                    ax = axes[r][c]
                    if c % 2 == 0:
                        lo = min(d[x].min(), d[y].min()) - .03
                        hi = max(d[x].max(), d[y].max()) + .03
                        ax.plot([lo, hi], [lo, hi], ls="--", lw=.7, c="#888888", zorder=0)
                        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
                    else:
                        ax.axhline(0, ls="--", lw=.7, c="#888888", zorder=0)
                if block == 0:
                    for _, rr in d.iterrows():
                        src_rows.append(dict(TF=rr.TF_name, PWM_set=rr.PWM, test_set=sp_lab, metric=mlab,
                                             best_monoPWM=rr[x], ArChIPelago=rr[y], delta=rr["Delta"]))
    for j, t in enumerate(["Random Forest trained and tested on HUMAN data", "",
                           "Random Forest trained on HUMAN and tested on MOUSE data", ""]):
        if t:
            axes[0][j].set_title(t, fontsize=9, loc="left", pad=18, fontweight="bold")
            axes[2][j].set_title(("Human test set" if j == 0 else "Mouse test set") +
                                 ": models on monoPWMs+diPWMs, individual TFs labelled",
                                 fontsize=9, loc="left", pad=18, fontweight="bold")
    for ax, letter in zip(axes.ravel(), [c for c in "ABAB" "CDCD" "ABAB" "CDCD"]):   # A,B = auROC row; C,D = auPRC row
        ax.text(-0.30, 1.10, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))
    h = [Line2D([], [], marker="o", ls="", color=COL[p], label=LBL[p]) for p in ["mono", "di", "mono+di"]]
    h += [Line2D([], [], marker="o", ls="", color="#e6a100", label="monoPWMs+diPWMs, above the baseline"),
          Line2D([], [], marker="o", ls="", color="#111111", label="monoPWMs+diPWMs, below the baseline")]
    fig.legend(handles=h, loc="lower center", ncol=5, frameon=False, fontsize=8, bbox_to_anchor=(.5, -.004))
    fig.tight_layout(rect=(0, .025, 1, 1), h_pad=3.4, w_pad=2.4)
    fig.savefig(os.path.join(OUT, "Figure_S2.png"), dpi=160)
    fig.savefig(os.path.join(FIGS, "Figure_S2.pdf"))
    plt.close(fig)
    pd.DataFrame(src_rows).to_csv(os.path.join(SRC, "Figure_S2_source_data.csv"), index=False)


# ==============================================================================================
# Figure S5 -- saturation
# ==============================================================================================
K_GRID = [1, 2, 4, 8, 16, 32, 64, 128]
METRICS = [("auroc_H", "auROC, human test set"), ("auprc_H", "auPRC, human test set"),
           ("auroc_M", "auROC, mouse test set"), ("auprc_M", "auPRC, mouse test set")]
BLUE, ORANGE, GREY = "#2a78d6", "#eb6834", "#bdbdbd"


def figure_s5():
    df = pd.read_csv(os.path.join(BASIS, "saturation", "saturation_results.csv"))
    tfs = sorted(df.tf.unique())
    Pmap = df.groupby("tf").P.first()
    for m, _ in METRICS:
        df["d_" + m] = df[m] - df["base_" + m]
    TABLE_RF = TABLE[(TABLE.Model == "RandomForestClassifier") & (TABLE.PWM == "mono+di")].set_index("TF_name")
    for tf in tfs:      # the reference must agree with the *_PWM columns of the results table (best monoPWM)
        b = df[df.tf == tf].iloc[0]
        assert abs(b.base_auroc_H - TABLE_RF.loc[tf, "roc_auc_test_H_PWM"]) < 1e-6, tf
        assert abs(b.base_auprc_M - TABLE_RF.loc[tf, "pr_auc_test_M_PWM"]) < 1e-6, tf
    dcols = ["d_" + m for m, _ in METRICS]
    mean = df.groupby(["tf", "design", "k"])[dcols].mean().reset_index()
    full = mean[(mean.design == "random") & (mean.k == mean.tf.map(Pmap))].set_index("tf")

    def curve(tf, design):
        P = Pmap[tf]
        sub = mean[(mean.tf == tf) & (mean.design == design)].set_index("k")
        out = {k: (full.loc[tf, dcols] if min(k, P) == P else sub.loc[min(k, P), dcols]) for k in K_GRID}
        out["P"] = full.loc[tf, dcols]
        return pd.DataFrame(out).T.astype(float)

    curves = {des: {tf: curve(tf, des) for tf in tfs} for des in ("random", "topk")}
    med = {des: pd.concat(curves[des].values(), keys=tfs).groupby(level=1).median().loc[K_GRID + ["P"]]
           for des in ("random", "topk")}

    XP = 512
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.0), sharex=True)
    for ax, (m, lbl), letter in zip(axes.ravel(), METRICS, "ABCD"):
        ax.axhline(0, color="#777777", lw=.8, ls=":", zorder=1)
        for tf in tfs:
            P = Pmap[tf]
            c = curves["random"][tf]["d_" + m]
            xs = [k for k in K_GRID if k < P] + [P]
            ax.plot(xs, [c[k] for k in K_GRID if k < P] + [c["P"]], color=GREY, lw=.7, alpha=.9, zorder=2)
        for des, col, ls, lab in (("random", BLUE, "-", "median over TFs, random subsets of k PWMs"),
                                  ("topk", ORANGE, "--", "median over TFs, k best PWMs by training auROC")):
            y = med[des]["d_" + m]
            ax.plot(K_GRID + [XP], [y[k] for k in K_GRID] + [y["P"]], color=col, lw=2.2, ls=ls, zorder=4, label=lab)
            ax.plot([XP], [y["P"]], marker="o", color=col, ms=5, zorder=5)
        ax.set_xscale("log", base=2)
        ax.set_xticks(K_GRID + [XP])
        ax.set_xticklabels([str(k) for k in K_GRID] + ["all\n(P)"])
        ax.tick_params(axis="x", which="minor", bottom=False)
        ax.set_title(lbl, fontsize=9, loc="left")
        ax.text(-0.16, 1.10, letter, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")
        ax.grid(axis="y", color="#ececec", lw=.6)
        ax.set_ylabel("$\\Delta$%s vs best single monoPWM" % ("auROC" if m.startswith("auroc") else "auPRC"), fontsize=7.5)
    for ax in axes[1]:
        ax.set_xlabel("number of PWMs used (k)")
    h, l = axes[0, 0].get_legend_handles_labels()
    h.append(Line2D([], [], color=GREY, lw=.8))
    l.append("single TF (mean over replicates), ending at its own number of PWMs P")
    fig.legend(h, l, loc="lower center", ncol=1, frameon=False, fontsize=7.5, bbox_to_anchor=(.5, -.005))
    fig.tight_layout(rect=(0, .085, 1, 1))
    fig.savefig(os.path.join(PREV, "Figure_S5_saturation.pdf"))
    fig.savefig(os.path.join(PREV, "Figure_S5_saturation.png"), dpi=160)
    plt.close(fig)
    rows = []
    for des in ("random", "topk"):
        for tf in tfs:
            c = curves[des][tf]
            for k in K_GRID + ["P"]:
                rows.append(dict(design=des, TF=tf, P=Pmap[tf], k=(Pmap[tf] if k == "P" else min(k, Pmap[tf])),
                                 k_requested=k, **{f"delta_{m}": c.loc[k, "d_" + m] for m, _ in METRICS}))
    pd.DataFrame(rows).to_csv(os.path.join(SRC, "Figure_S5_source_data.csv"), index=False)
    med_out = pd.concat({des: med[des] for des in ("random", "topk")}, names=["design", "k"])
    med_out.to_csv(os.path.join(SRC, "Figure_S5_median_curves.csv"))


# ==============================================================================================
# Figure S6 -- cross-species transfer and mouse-trained control
# ==============================================================================================
SIM_COL = "Similarity of the baseline human monoPWM to the nearest mouse monoPWM (Pearson r of aligned columns)"


def figure_s6():
    t5 = pd.read_csv(os.path.join(RES, "Sup_Tables", "Sup_Table_5_cross_species_and_mouse_trained.csv")).set_index("TF")
    fail = [tf for tf in t5.index if t5.loc[tf, "H>M: below baseline on >=1 metric"]]
    fig, axes = plt.subplots(2, 2, figsize=(7.6, 6.6))

    for ax, met, letter in zip(axes[0], ["auROC", "auPRC"], "AB"):
        x = t5[f"H>M: dau{'ROC' if met == 'auROC' else 'PRC'}"]
        y = t5[f"M>M: dau{'ROC' if met == 'auROC' else 'PRC'}"]
        ax.axhline(0, ls="--", lw=.7, c="#888888"); ax.axvline(0, ls="--", lw=.7, c="#888888")
        lo, hi = min(x.min(), y.min()) - .01, max(x.max(), y.max()) + .01
        ax.plot([lo, hi], [lo, hi], ls=":", lw=.7, c="#bbbbbb")
        isf = t5.index.isin(fail)
        ax.scatter(x[~isf], y[~isf], s=18, c="#5b8ff9", lw=0, alpha=.85, label="all other TFs")
        ax.scatter(x[isf], y[isf], s=28, c="#c0392b", lw=0, label="below the baseline after transfer")
        for tf in fail:
            ax.annotate(tf, (x[tf], y[tf]), fontsize=6.5, xytext=(3, 3), textcoords="offset points", color="#c0392b")
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.set_xlabel(f"$\\Delta${met}, human-trained model on mouse\n(vs the best monoPWM selected on human data)", fontsize=7.5)
        ax.set_ylabel(f"$\\Delta${met}, mouse-trained model on mouse\n(vs best mouse monoPWM)", fontsize=7.5)
        ax.text(-0.22, 1.07, letter, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")
    axes[0][0].legend(frameon=False, fontsize=6.5, loc="lower right")

    # C: motif similarity vs cross-species gain
    ax = axes[1][0]
    sim = t5[SIM_COL]
    d = t5["H>M: dauROC"]
    isf = t5.index.isin(fail)
    ax.axhline(0, ls="--", lw=.7, c="#888888")
    ax.scatter(sim[~isf], d[~isf], s=18, c="#5b8ff9", lw=0, alpha=.85)
    ax.scatter(sim[isf], d[isf], s=28, c="#c0392b", lw=0)
    for tf in fail:
        ax.annotate(tf, (sim[tf], d[tf]), fontsize=6.5, xytext=(3, 3), textcoords="offset points", color="#c0392b")
    rho, p = spearmanr(sim, d)
    ax.set_xlabel("similarity of the baseline human monoPWM to the\nmost similar mouse monoPWM (Pearson r)", fontsize=7.5)
    ax.set_ylabel("$\\Delta$auROC, human-trained model on mouse", fontsize=7.5)
    ax.set_title(f"Spearman $\\rho$ = {rho:.2f}, P = {p:.2f}", fontsize=7.5, loc="left")
    ax.text(-0.22, 1.07, "C", transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")

    # D: distribution of the gain in the three settings
    ax = axes[1][1]
    sets = [("H>H: dauROC (reference)", "H>H: dauPRC (reference)", "human-trained,\nhuman test"),
            ("H>M: dauROC", "H>M: dauPRC", "human-trained,\nmouse test"),
            ("M>M: dauROC", "M>M: dauPRC", "mouse-trained,\nmouse test")]
    pos, data, colors = [], [], []
    for i, (rc, pc, _) in enumerate(sets):
        pos += [i * 3 + 0.6, i * 3 + 1.5]
        data += [t5[rc].values, t5[pc].values]
        colors += ["#5b8ff9", "#1f9e63"]
    bp = ax.boxplot(data, positions=pos, widths=.7, showfliers=False, patch_artist=True,
                    medianprops=dict(color="black", lw=1.2))
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c); patch.set_alpha(.45); patch.set_linewidth(.8)
    rng = np.random.default_rng(0)
    for p_, v, c in zip(pos, data, colors):
        ax.scatter(p_ + rng.uniform(-.18, .18, len(v)), v, s=7, c=c, lw=0, alpha=.75)
    ax.axhline(0, ls="--", lw=.7, c="#888888")
    ax.set_xticks([i * 3 + 1.05 for i in range(3)])
    ax.set_xticklabels([s[2] for s in sets], fontsize=7)
    ax.set_ylabel("$\\Delta$ vs the best single monoPWM", fontsize=7.5)
    ax.legend(handles=[Line2D([], [], marker="s", ls="", color="#5b8ff9", label="$\\Delta$auROC"),
                       Line2D([], [], marker="s", ls="", color="#1f9e63", label="$\\Delta$auPRC")],
              frameon=False, fontsize=7, loc="upper left")
    ax.text(-0.22, 1.07, "D", transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")
    fig.tight_layout(h_pad=2.4, w_pad=2.4)
    fig.savefig(os.path.join(PREV, "Figure_S6_cross_species.pdf"))
    fig.savefig(os.path.join(PREV, "Figure_S6_cross_species.png"), dpi=160)
    plt.close(fig)
    cols = ["H>H: dauROC (reference)", "H>H: dauPRC (reference)", "H>M: dauROC", "H>M: dauPRC",
            "M>M: dauROC", "M>M: dauPRC",
            SIM_COL, "H>M: below baseline on >=1 metric"]
    t5[cols].to_csv(os.path.join(SRC, "Figure_S6_source_data.csv"))


if __name__ == "__main__":
    figure_s2()
    figure_s5()
    figure_s6()
    print("figures written to", OUT)
