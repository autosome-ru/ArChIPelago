#!/usr/bin/env python3
"""Fig. S2: Random Forest against the best single monoPWM (selected on the human training set), human test set
(chr1, 8, 21; left) and mouse test set (chr1, 8, 19; right), for the monoPWM / diPWM / monoPWM+diPWM models;
from Manuscript_analysis/results_table.csv. Panel letters A, B = auROC row, C, D = auPRC row of each 2 x 2 block.
Writes Figures/panels/Figure_S2_RF_vs_best_monoPWM.pdf and Figures/source_data/Figure_S2_source_data.csv.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from archi_paths import RES

FIGS = os.path.join(RES, "Figures", "panels")
SRC = os.path.join(RES, "Figures", "source_data")
os.makedirs(FIGS, exist_ok=True)
os.makedirs(SRC, exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"]})
COL = {"mono": "#d1478f", "di": "#5b8ff9", "mono+di": "#1f9e63"}
LBL = {"mono": "monoPWMs", "di": "diPWMs", "mono+di": "monoPWMs+diPWMs"}

TABLE = pd.read_csv(os.path.join(RES, "results_table.csv"), sep="\t")
RF = TABLE[TABLE.Model == "RandomForestClassifier"]


# Figure S2
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
    fig.savefig(os.path.join(FIGS, "Figure_S2_RF_vs_best_monoPWM.pdf"))
    plt.close(fig)
    pd.DataFrame(src_rows).to_csv(os.path.join(SRC, "Figure_S2_source_data.csv"), index=False)


if __name__ == "__main__":
    figure_s2()
    print("written:", os.path.join(FIGS, "Figure_S2_RF_vs_best_monoPWM.pdf"))
