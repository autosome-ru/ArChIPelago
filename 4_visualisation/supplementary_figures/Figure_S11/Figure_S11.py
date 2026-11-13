#!/usr/bin/env python3
"""B5 figure (Fig. S11) from the CSVs written by 3_analysis/transfer_drivers/B5_01_drivers_and_ceiling.py; writes Figure_S11.{pdf,png} here.
Run: OMP_NUM_THREADS=1 nice -n 15 python3 4_visualisation/supplementary_figures/Figure_S11/Figure_S11.py"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats
import sys

FIG = os.path.dirname(os.path.abspath(__file__)); R2 = os.path.dirname(os.path.dirname(os.path.dirname(FIG)))
D = os.path.join(R2, "3_analysis", "transfer_drivers")   # input tables
sys.path.insert(0, os.path.join(R2, "4_visualisation", "style"))
import archi_style as ST
ST.apply_rc()
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
                     "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8, "xtick.major.size": 3, "ytick.major.size": 3})
PFONT = "Arial Unicode MS"   # Arial lacks the Unicode superscripts used in P values
GREEN, GREY = ST.PALETTE["RF"], "#999999"
C = pd.read_csv(D + "/B5_A_spearman_drivers.csv")
H = pd.read_csv(D + "/B5_B_headroom_per_TF.csv", index_col=0)
T = pd.read_csv(D + "/B5_per_TF_table.csv", index_col=0)
cm = LinearSegmentedColormap.from_list("brownpurple", ["#8C510A", "#FFFFFF", "#542788"])   # rho; green is reserved for RF
cm_sc = LinearSegmentedColormap.from_list("greys", ["#F2F2F2", "#8C8C8C", "#000000"])   # baseline ramp in greys (green is reserved for RF)

fig = plt.figure(figsize=(10.5, 7.6))
gs = fig.add_gridspec(2, 6, height_ratios=[1.25, 1], hspace=0.55, wspace=1.3)
def letter(ax, L, dx=-0.12):
    ax.text(dx, 1.10, L, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left")

# A
ax = fig.add_subplot(gs[0, 0:3])
drv = list(dict.fromkeys(C.driver)); out = list(dict.fromkeys(C.outcome))
M = C.pivot(index="driver", columns="outcome", values="rho").loc[drv, out]
P = C.pivot(index="driver", columns="outcome", values="p").loc[drv, out]
Q = C.pivot(index="driver", columns="outcome", values="q_BH").loc[drv, out]
im = ax.imshow(M.values, cmap=cm, vmin=-0.7, vmax=0.7, aspect="auto")
for i in range(len(drv)):
    for j in range(len(out)):
        s_ = "%.2f" % M.values[i, j] + ("*" if Q.values[i, j] < 0.05 else ("\u2020" if P.values[i, j] < 0.05 else ""))
        ax.text(j, i, s_, ha="center", va="center", color="black", fontsize=7, family=PFONT)
print("A n_dagger", int(((P.values < 0.05) & (Q.values >= 0.05)).sum()), "n_star", int((Q.values < 0.05).sum()), "rho_sum", round(float(M.values.sum()), 6))
ax.set_xticks(range(len(out))); ax.set_xticklabels(["\u0394auROC,\nhuman to mouse\ntransfer", "\u0394auPRC,\nhuman to mouse\ntransfer", "loss\nauROC", "loss\nauPRC"], fontsize=7)
lab = C.drop_duplicates("driver").set_index("driver").driver_label
ax.set_yticks(range(len(drv))); ax.set_yticklabels([lab[d] for d in drv], fontsize=7)
ax.set_title("Drivers of transfer, Spearman \u03c1 (n = %d TFs)" % len(T), loc="left", pad=6)
for sp_ in ax.spines.values(): sp_.set_visible(False)
cbA = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.04); cbA.ax.set_title("Spearman \u03c1", fontsize=8, pad=4); cbA.outline.set_edgecolor("black")
cbA.ax.tick_params(labelsize=7)
letter(ax, "A", -0.62)

# B: scatter
ax = fig.add_subplot(gs[0, 3:6])
sz = lambda n: 20 + 90 * np.sqrt(n / 10000.0)
sc = ax.scatter(H.dPRC, H.dROC, c=H.bROC, cmap=cm_sc, vmin=0.7, vmax=1.0, s=sz(T.n_train_pos_H), edgecolor="black", linewidth=0.4, zorder=3)
cb = fig.colorbar(sc, ax=ax, pad=0.02, fraction=0.04); cb.set_label("Baseline auROC (human test set)"); cb.outline.set_edgecolor("black")
cb.ax.tick_params(labelsize=7)
from adjustText import adjust_text
sel = H.index[(H.dROC > 0.075) | (H.dPRC > 0.15)]
txt = [ax.text(H.dPRC[tf], H.dROC[tf], tf, fontsize=7) for tf in sel]
fig.canvas.draw()
# obstacles: every mark as a centre plus a ring of its radius; try several settings, keep the one with the fewest overlaps, draw leader lines last
ang = np.linspace(0, 2 * np.pi, 7)[:-1]
disp = ax.transData.transform(np.c_[H.dPRC.values, H.dROC.values])
rad = np.sqrt(sz(T.loc[H.index, "n_train_pos_H"].values)) / 2 * fig.dpi / 72 + 2
def n_ov(texts):
    bbs = [t.get_window_extent() for t in texts]; n = 0; ab = ax.get_window_extent()
    for i, bb in enumerate(bbs):
        n += sum(bb.overlaps(bbs[j]) for j in range(i + 1, len(bbs)))
        n += int(sum((bb.x0 - r_ < px < bb.x1 + r_) and (bb.y0 - r_ < py < bb.y1 + r_) for (px, py), r_ in zip(disp, rad)))
        n += not (ab.x0 <= bb.x0 and bb.x1 <= ab.x1 and ab.y0 <= bb.y0 and bb.y1 <= ab.y1)
    return n
rng = np.random.RandomState(1); best = None
for k in range(30):
    for t, tf in zip(txt, sel):
        p0 = ax.transData.transform((H.dPRC[tf], H.dROC[tf])) + (rng.uniform(-40, 40, 2) if k else 0)
        t.set_position(ax.transData.inverted().transform(p0))
    ring = ax.transData.inverted().transform(np.vstack([disp] + [disp + (rad[:, None] * np.array([np.cos(a_), np.sin(a_)])) for a_ in ang]))
    adjust_text(txt, ax=ax, x=ring[:, 0], y=ring[:, 1], target_x=H.dPRC[sel].values, target_y=H.dROC[sel].values, expand=(1.4, 1.8),
                force_text=(1.5, 1.5), force_static=(1.5, 1.5), force_pull=(0.03, 0.03), ensure_inside_axes=True, max_move=(15, 15), iter_lim=3000)
    fig.canvas.draw(); n = n_ov(txt)
    if best is None or n < best[0]: best = (n, [t.get_position() for t in txt])
    if n == 0: break
for t, p in zip(txt, best[1]): t.set_position(p)
fig.canvas.draw(); print("B label overlaps (label-label, label-mark, label-axes):", best[0])
for t, tf in zip(txt, sel):
    bb = t.get_window_extent(); px, py = ax.transData.transform((H.dPRC[tf], H.dROC[tf]))
    cx, cy = min(max(px, bb.x0), bb.x1), min(max(py, bb.y0), bb.y1)
    if np.hypot(px - cx, py - cy) > 7:
        ax.annotate("", xy=(H.dPRC[tf], H.dROC[tf]), xycoords="data", xytext=((cx - bb.x0) / bb.width, (cy - bb.y0) / bb.height), textcoords=t,
                    arrowprops=dict(arrowstyle="-", lw=0.4, color="black", shrinkA=0, shrinkB=4, relpos=(0.5, 0.5)))
ax.set_xlabel("\u0394auPRC (RF vs best monoPWM, human test set)"); ax.set_ylabel("\u0394auROC")
ax.set_title("Gain per TF, coloured by baseline auROC", loc="left", pad=6)
hs = [ax.scatter([], [], s=sz(n), c="#999999", edgecolor="black", linewidth=0.4, label="{:,}".format(n)) for n in (1000, 5000, 10000)]
ax.legend(handles=hs, title="Training positives", title_fontsize=7, frameon=False, loc="lower right", labelspacing=0.9, borderpad=0.2, handletextpad=0.6)
letter(ax, "B", -0.2)
print("B n", len(H), "dROC_sum", round(float(H.dROC.sum()), 6), "dPRC_sum", round(float(H.dPRC.sum()), 6))

# C, D, E: headroom panels
def bpanel(pos, x, y, xl, yl, title, rho_raw=None, right=False, L=""):
    a = fig.add_subplot(gs[1, pos[0]:pos[1]])
    a.scatter(x, y, s=22, color=GREY, edgecolor="black", linewidth=0.4, zorder=3)
    r, p = stats.spearmanr(x, y)
    a.axhline(0, color="black", lw=0.6, ls=":")
    a.set_xlabel(xl); a.set_ylabel(yl); a.set_title(title, loc="left", pad=6)
    t = "\u03c1 = %.2f, %s" % (r, ST.fmt_p(p))
    if rho_raw: t += "\n" + rho_raw
    a.text(0.97 if right else 0.03, 0.97, t, transform=a.transAxes, va="top", ha="right" if right else "left", fontsize=7, family=PFONT)
    letter(a, L, -0.24)
    print(L, title, "rho", round(r, 4), "p", p)
    return a
rr = pd.read_csv(D + "/B5_B_ceiling_spearman.csv").set_index("quantity")["rho"]
r1 = rr.iloc[0]; r2 = rr.iloc[1]
bpanel((0, 2), H.bROC, H.head_roc, "Baseline auROC (best monoPWM)", "\u0394auROC / (1 - baseline auROC)", "Human test set: auROC gain per unit headroom",
       "raw \u0394auROC vs baseline: \u03c1 = %.2f" % r1, L="C")
bpanel((2, 4), H.bPRC, H.head_prc, "Baseline auPRC (best monoPWM)", "\u0394auPRC / (1 - baseline auPRC)", "Human test set: auPRC gain per unit headroom",
       "raw \u0394auPRC vs baseline: \u03c1 = %.2f" % r2, L="D")
bpanel((4, 6), H.bPRC, H.rel_prc, "Baseline auPRC (best monoPWM)", "\u0394auPRC / baseline auPRC", "Human test set: auPRC gain relative to baseline", right=True, L="E")
fig.savefig(FIG + "/Figure_S11.pdf", bbox_inches="tight"); fig.savefig(FIG + "/Figure_S11.png", dpi=600, bbox_inches="tight")
