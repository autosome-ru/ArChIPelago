#!/usr/bin/env python3
"""Fig. S2: Random Forest against the best single monoPWM (selected on the human training set), human test set
(chr1, 8, 21; left) and mouse test set (chr1, 8, 19; right), for the monoPWM / diPWM / monoPWM+diPWM models;
from Manuscript_analysis/results_table.csv. Panel letters A, B = auROC row, C, D = auPRC row of each 2 x 2 block.
Writes Figures/panels/Figure_S2_RF_vs_best_monoPWM.pdf and Figures/source_data/Figure_S2_source_data.csv.

Restyled version of the released script; data and layout unchanged. Palette from ../archi_style.py
(PWM sets = lightness steps of the Random Forest green, black outline; above/below the baseline = dark green of monoPWMs+diPWMs / red), opaque marks,
no MathText (Unicode Delta), all text black. Usage: python3 Figure_S2_restyled.py [output_dir]; reads the released
results_table.csv (ARCHI_PKG, read only); writes Figure_S2_restyled.{pdf,png} and Figure_S2_plotdata.csv.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import archi_style as AS

RES = os.environ.get("ARCHI_PKG", os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../results/round1_package/Manuscript_analysis"))   # released package, read only
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
os.makedirs(OUT, exist_ok=True)
AS.apply_rc()
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "legend.fontsize": 9})
COL = AS.PWMSET_COL
LBL = AS.PWMSET_LABELS
ABOVE, BELOW = AS.PWMSET_COL["mono+di"], AS.PALETTE["monoPWM"]

TABLE = pd.read_csv(os.path.join(RES, "results_table.csv"), sep="\t")
RF = TABLE[TABLE.Model == "RandomForestClassifier"]


# Figure S2
MEDGREY = "#555555"


def panel_scatter(ax, d, xcol, ycol, xlab, ylab, by_pwm=True, label_tfs=None):
    if by_pwm:   # as the legend says for every panel, a model below the best monoPWM is an open black diamond
        for pwm in ["mono", "di", "mono+di"]:
            s = d[d.PWM == pwm]
            delta = (s[ycol] - s[xcol]) if "Delta" not in ycol else s[ycol]
            up, dn = s[delta >= 0], s[delta < 0]
            ax.scatter(up[xcol], up[ycol], s=16, c=COL[pwm], edgecolors="black", linewidths=0.4, label=LBL[pwm])
            ax.scatter(dn[xcol], dn[ycol], s=22, marker="D", c="white", edgecolors="black", linewidths=0.7, zorder=4)
    else:
        s = d[d.PWM == "mono+di"]
        delta = (s[ycol] - s[xcol]) if "Delta" not in ycol else s[ycol]
        up, dn = s[delta >= 0], s[delta < 0]
        ax.scatter(up[xcol], up[ycol], s=16, c=COL["mono+di"], edgecolors="black", linewidths=0.4)
        ax.scatter(dn[xcol], dn[ycol], s=22, marker="D", c="white", edgecolors="black", linewidths=0.7, zorder=4)
        ax._lab = (s, xcol, ycol, label_tfs)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)


def _cross(a, b, c, d):
    o = lambda p, q, r: np.sign((q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]))
    return o(a, b, c) * o(a, b, d) < 0 and o(c, d, a) * o(c, d, b) < 0


def _seg_hits_box(a, b, box, pad=1.0):   # does segment a-b pass through box (x0, x1, y0, y1)? sampled
    t = np.linspace(0.08, 0.92, 25)[:, None]
    q = a + t * (b - a)
    return bool(((q[:, 0] > box[0] - pad) & (q[:, 0] < box[1] + pad) & (q[:, 1] > box[2] - pad) & (q[:, 1] < box[3] + pad)).any())


def add_labels(ax, mark_r=4.5, line_r=1.5, edge=4.0, gap=2.0):
    """Deterministic greedy label placement: each label takes the nearest candidate position
    (24 directions x 12 distances around its mark) that keeps it off every mark, off the diagonal / zero line, off the other
    labels, at least `edge` px inside the axes, and whose leader line crosses no other leader line or label. The most crowded
    labels are placed first. When no position is free, the cheapest one wins: crossing a reference line costs 1, a mark,
    crossing leader or a leader through a label 3, another label 4, leaving the axes 5. Prints the labels left on a reference
    line or worse."""
    s, xcol, ycol, label_tfs = ax._lab
    lab = s[s.TF_name.isin(label_tfs)]
    old = (ax.get_xlim(), ax.get_ylim())
    for lim, setter in ((old[0], ax.set_xlim), (old[1], ax.set_ylim)):   # room for labels; data untouched
        w = lim[1] - lim[0]; setter(lim[0] - 0.10 * w, lim[1] + 0.10 * w)
    for getter, setter, lim in ((ax.get_xticks, ax.set_xticks, old[0]), (ax.get_yticks, ax.set_yticks, old[1])):
        setter([t for t in getter() if lim[0] - 1e-9 <= t <= lim[1] + 1e-9])   # no tick outside the original range
    fig = ax.figure
    fig.canvas.draw()
    T = ax.transData
    marks = T.transform(np.c_[s[xcol].values, s[ycol].values])
    ab = ax.get_window_extent()
    lines = []
    for ln in ax.lines:   # reference lines (diagonal, zero line), sampled in display coordinates
        xy = ln.get_transform().transform(ln.get_xydata())
        if len(xy) >= 2:
            q = xy[0] + np.linspace(0, 1, 200)[:, None] * (xy[-1] - xy[0])
            lines.append(q[(q[:, 0] >= ab.x0) & (q[:, 0] <= ab.x1) & (q[:, 1] >= ab.y0) & (q[:, 1] <= ab.y1)])
    lines = np.vstack(lines) if lines else np.empty((0, 2))
    texts = {r.TF_name: ax.text(r[xcol], r[ycol], r.TF_name, fontsize=7, color="black", ha="center", va="center") for _, r in lab.iterrows()}
    fig.canvas.draw()
    size = {k: (t.get_window_extent().width, t.get_window_extent().height) for k, t in texts.items()}
    tgt = {r.TF_name: T.transform((r[xcol], r[ycol])) for _, r in lab.iterrows()}
    crowd = {k: int((np.hypot(*(marks - p).T) < 40).sum()) for k, p in tgt.items()}
    px = fig.dpi / 72.0
    ang = np.deg2rad(np.arange(0, 360, 15.0))
    dist = np.array([4, 7, 10, 15, 21, 28, 36, 45, 55, 66, 80, 95]) * px
    def place(order):
      placed, segs, out, total = [], [], {}, 0
      for k in order:
          p, (w, h) = tgt[k], size[k]
          best = None
          for d in dist:
              for a in ang:
                  c = p + np.array([np.cos(a) * (d + w / 2), np.sin(a) * (d + h / 2)])
                  box = (c[0] - w / 2, c[0] + w / 2, c[1] - h / 2, c[1] + h / 2)
                  a0 = np.array([min(max(p[0], box[0]), box[1]), min(max(p[1], box[2]), box[3])])
                  bad = 0
                  parts = []
                  bad += 5 * (not (ab.x0 + edge <= box[0] and box[1] <= ab.x1 - edge and ab.y0 + edge <= box[2] and box[3] <= ab.y1 - edge))
                  parts.append(bad)
                  bad += 3 * int(((marks[:, 0] > box[0] - mark_r) & (marks[:, 0] < box[1] + mark_r) & (marks[:, 1] > box[2] - mark_r) & (marks[:, 1] < box[3] + mark_r)).sum())
                  parts.append(bad)
                  bad += bool(((lines[:, 0] > box[0] - line_r) & (lines[:, 0] < box[1] + line_r) & (lines[:, 1] > box[2] - line_r) & (lines[:, 1] < box[3] + line_r)).any())
                  parts.append(bad)
                  bad += 4 * sum(not (box[1] + gap < b[0] or b[1] + gap < box[0] or box[3] + gap < b[2] or b[3] + gap < box[2]) for b in placed)
                  parts.append(bad)
                  bad += 3 * sum(_cross(a0, p, s0, s1) for s0, s1 in segs)
                  parts.append(bad)
                  bad += 3 * (sum(_seg_hits_box(a0, p, b) for b in placed) + sum(_seg_hits_box(s0, s1, box) for s0, s1 in segs))
                  if best is None or bad < best[0] or (bad == best[0] and d < best[1]):
                      best = (bad, d, c, box, a0, parts + [bad])
                  if bad == 0:
                      break
              if best[0] == 0:
                  break
          total += best[0]
          out[k] = best
          placed.append(best[3]); segs.append((best[4], p))
      return total, out

    # greedy placement in the crowding order and in 40 seeded random orders; the order with the lowest total cost wins
    rs = np.random.RandomState(2)
    orders = [sorted(tgt, key=lambda k: -crowd[k])] + [list(rs.permutation(sorted(tgt))) for _ in range(40)]
    total, out = min((place(o) for o in orders), key=lambda r: r[0])
    n_bad = sum(b[0] > 0 for b in out.values())
    for k, b in out.items():
        texts[k].set_position(T.inverted().transform(b[2]))
        if b[0] > 0 and os.environ.get("S2_DEBUG"):
            print("  unplaced", k, "cost", b[0], "cumulative edge/mark/line/label/cross/hit", b[5])
    fig.canvas.draw()   # leader lines from the nearest edge of each label to its mark
    for k, t in texts.items():
        bb = t.get_window_extent(); px_, py_ = tgt[k]
        cx, cy = min(max(px_, bb.x0), bb.x1), min(max(py_, bb.y0), bb.y1)
        if np.hypot(px_ - cx, py_ - cy) > 7:
            ax.annotate("", xy=T.inverted().transform((px_, py_)), xycoords="data", xytext=((cx - bb.x0) / bb.width, (cy - bb.y0) / bb.height),
                        textcoords=t, arrowprops=dict(arrowstyle="-", lw=0.4, color="black", shrinkA=0, shrinkB=4, relpos=(0.5, 0.5)))
    print("labels without a conflict-free position:", n_bad, "of", len(texts))


def figure_s2():
    fig, axes = plt.subplots(4, 4, figsize=(8.4, 9.0))
    src_rows = []
    for block in (0, 1):                                  # 0 = coloured by PWM set, 1 = mono+di labelled
        for col_block, sp, sp_lab in ((0, "H", "human"), (2, "M", "mouse")):
            for row_in, met, mlab in ((0, "roc_auc", "auROC"), (1, "pr_auc", "auPRC")):
                r = block * 2 + row_in
                x = f"{met}_test_{sp}_PWM"
                y = f"{met}_test_{sp}"
                d = RF.copy()
                d["Delta"] = d[y] - d[x]
                dl = "\u0394" + mlab
                if block == 0:
                    panel_scatter(axes[r][col_block], d, x, y, f"{mlab} best monoPWM", mlab)
                    panel_scatter(axes[r][col_block + 1], d, x, "Delta", f"{mlab} best monoPWM", dl)
                else:
                    s = d[d.PWM == "mono+di"]
                    lab = set(s.nlargest(5, "Delta").TF_name) | set(s.nsmallest(3, "Delta").TF_name)
                    panel_scatter(axes[r][col_block], d, x, y, f"{mlab} best monoPWM", mlab, by_pwm=False, label_tfs=lab)
                    panel_scatter(axes[r][col_block + 1], d, x, "Delta", f"{mlab} best monoPWM", dl,
                                  by_pwm=False, label_tfs=lab)
                for c in (col_block, col_block + 1):
                    ax = axes[r][c]
                    if c % 2 == 0:
                        lo = min(d[x].min(), d[y].min()) - .03
                        hi = max(d[x].max(), d[y].max()) + .03
                        ax.plot([lo, hi], [lo, hi], ls=":", lw=0.8, c="black", zorder=0)
                        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
                    else:
                        ax.axhline(0, ls=":", lw=0.6, c="black", zorder=0)
                if block == 0:
                    for _, rr in d.iterrows():
                        src_rows.append(dict(TF=rr.TF_name, PWM_set=rr.PWM, test_set=sp_lab, metric=mlab,
                                             best_monoPWM=rr[x], ArChIPelago=rr[y], delta=rr["Delta"]))
    for ax in axes.ravel():
        ax.tick_params(labelsize=7, width=0.8, length=3)
        ax.xaxis.label.set_size(8); ax.yaxis.label.set_size(8)
        for sp_ in ax.spines.values():
            sp_.set_linewidth(0.8)
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))
    for r in range(4):
        axes[r][0].set_title("Human test set", fontsize=8, loc="left", pad=8)
        axes[r][2].set_title("Mouse test set (human-trained)", fontsize=8, loc="left", pad=8)
    fig.tight_layout(rect=(0, .045, 1, 1), h_pad=2.2, w_pad=1.4)
    fig.canvas.draw()
    for r in (2, 3):
        for c in range(4):
            add_labels(axes[r][c])
    for i, ax in enumerate(axes.ravel()):
        ax.text(-0.30, 1.12, "ABCDEFGHIJKLMNOP"[i], transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left")
    h = [Line2D([], [], marker="o", ls="", markerfacecolor=COL[p], markeredgecolor="black", markeredgewidth=0.4, markersize=6, label=LBL[p])
         for p in ["mono", "di", "mono+di"]]
    h += [Line2D([], [], marker="D", ls="", markerfacecolor="white", markeredgecolor="black", markeredgewidth=0.7, markersize=5,
                 label="below the best monoPWM")]
    fig.legend(handles=h, loc="lower center", ncol=4, frameon=False, fontsize=7.5, bbox_to_anchor=(.5, 0.0))
    AS.save_both(fig, os.path.join(OUT, "Figure_S2_restyled"))
    plt.close(fig)
    pd.DataFrame(src_rows).to_csv(os.path.join(OUT, "Figure_S2_plotdata.csv"), index=False)


if __name__ == "__main__":
    figure_s2()
    print("written:", os.path.join(OUT, "Figure_S2_restyled.pdf"))
