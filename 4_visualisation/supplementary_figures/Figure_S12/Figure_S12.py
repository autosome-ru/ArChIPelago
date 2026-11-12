"""Fig. S12 and the Sup. Table 7 B9 sheets on the final data (v3 negatives, final retraining, 5 TFs).
Inputs:
  h3_inputs b9_dir/b9_random_split_per_fit.csv  -> byte copy here (md5 checked)
  h3_inputs b9_dir/b9_summary.csv               -> the summary recomputed here must equal it (max |diff| <= 1e-12)
Outputs (3_analysis/random_split/): b9_random_split_per_fit.csv, b9_summary.csv, _provenance_h3.yaml; the figure Figure_S12.pdf/png
goes to this folder.
usage: python3 4_visualisation/supplementary_figures/Figure_S12/Figure_S12.py
"""
import pandas as pd, numpy as np, matplotlib, shutil
matplotlib.use("Agg"); import matplotlib.pyplot as plt
import os, sys
FIG = os.path.dirname(os.path.abspath(__file__)); R2 = os.path.dirname(os.path.dirname(os.path.dirname(FIG)))
HERE = os.path.join(R2, "3_analysis", "random_split"); os.chdir(HERE)   # tables and provenance in the analysis folder
sys.path.insert(0, os.path.join(R2, "4_visualisation", "style")); import archi_style as AS
sys.path.insert(0, os.path.join(R2, "3_analysis/integration")); import p53drop as PD
sys.path.insert(0, os.path.join(R2, "3_analysis", "integration")); import w2; import h3_inputs as H
SRC = H.IN["b9_dir"]   # B9 rerun with the training scaler
B9_TFS = ["PRGR", "RXRA", "SPI1", "ANDR", "REST"]
FIGOUT = [os.path.join(FIG, "Figure_S12" + e) for e in (".pdf", ".png")]
for _f in ["b9_random_split_per_fit.csv", "b9_summary.csv"] + FIGOUT:
    if os.path.exists(_f):
        w2.backup(_f)
shutil.copy2(os.path.join(SRC, "b9_random_split_per_fit.csv"), "b9_random_split_per_fit.csv")
assert w2.md5("b9_random_split_per_fit.csv") == w2.md5(os.path.join(SRC, "b9_random_split_per_fit.csv"))
AS.apply_rc()
matplotlib.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
                            "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8, "xtick.major.size": 3, "ytick.major.size": 3})
d = pd.read_csv("b9_random_split_per_fit.csv")
_exp = sorted(B9_TFS)
if (d.groupby(["TF", "split"]).seed.nunique() != 3).any():
    raise SystemExit("CHECK FAILED B9: not 3 seeds per TF and split")
if sorted(d.TF.unique()) != _exp or set(d.TF) & {"P53"} or not set(d.TF) <= set(PD.reference_tfs()):
    raise SystemExit("CHECK FAILED B9 TFs %s, expected %s without P53" % (sorted(d.TF.unique()), _exp))
print("check B9: %d TFs %s, P53 absent" % (d.TF.nunique(), sorted(d.TF.unique())))
d["gain_auroc"] = d.rf_auroc - d.pwm_auroc; d["gain_auprc"] = d.rf_auprc - d.pwm_auprc
m = d.groupby(["TF", "split"])[["rf_auroc", "pwm_auroc", "gain_auroc", "rf_auprc", "pwm_auprc", "gain_auprc",
                                "frac_test_pos_overlap", "frac_test_pos_exact"]].mean().reset_index()
h3b = pd.read_csv(os.path.join(SRC, "b9_summary.csv")).set_index(["TF", "split"]).sort_index()
_dd = (m.set_index(["TF", "split"]).sort_index()[h3b.columns] - h3b).abs().max().max()
if list(h3b.columns) != [c for c in m.columns if c not in ("TF", "split")] or not _dd <= 1e-12:
    raise SystemExit("CHECK FAILED B9 summary differs from H3b b9_summary.csv (max |diff| %s)" % _dd)
print("check B9 summary = H3b b9_summary.csv: max |diff| %.1e" % _dd)
m.to_csv("b9_summary.csv", index=False)
order = m[m.split == "random"].sort_values("frac_test_pos_overlap").TF.tolist()
fig, ax = plt.subplots(1, 3, figsize=(7.4, 3.2))
for k, (col, lab) in enumerate([("gain_auroc", "ΔauROC\n(ArChIPelago (RF) minus best monoPWM)"), ("gain_auprc", "ΔauPRC\n(ArChIPelago (RF) minus best monoPWM)")]):
    a = ax[k]; x = np.arange(len(order))
    for sp, mk, c, dx in [("chromosome", "o", "#009E73", -0.13), ("random", "s", "white", 0.13)]:   # side by side: equal values stay visible
        v = m[m.split == sp].set_index("TF").loc[order, col]
        a.scatter(x + dx, v, marker=mk, s=24, facecolor=c, edgecolor="black", linewidth=0.4, zorder=3,
                  label="chromosome hold-out" if sp == "chromosome" else "random split")
    a.axhline(0, color="black", lw=0.6, ls=":"); a.set_xticks(x); a.set_xticklabels(order, rotation=45, ha="right")
    a.set_ylabel(lab); a.text(-0.26, 1.12, "AB"[k], transform=a.transAxes, fontweight="bold", fontsize=10, va="bottom"); a.set_title("Human test set", loc="left", pad=6)
    a.spines[["top", "right"]].set_visible(False)
ax[0].legend(frameon=False, loc="center right", bbox_to_anchor=(1.0, 0.52), handletextpad=0.2)
v = m[m.split == "random"].set_index("TF").loc[order]
x = np.arange(len(order))
ax[2].bar(x, 100 * v.frac_test_pos_overlap, color="#B7E3D3", edgecolor="black", lw=0.4, label="within 300 bp")
ax[2].bar(x, 100 * v.frac_test_pos_exact, color="#00503B", edgecolor="black", lw=0.4, label="same summit")
ax[2].set_xticks(x); ax[2].set_xticklabels(order, rotation=45, ha="right")
ax[2].set_ylabel("Test positives with a training\npositive at the locus (%)"); ax[2].set_ylim(0, 125); ax[2].set_yticks([0, 20, 40, 60, 80, 100]); ax[2].legend(frameon=False, loc="upper left", ncol=2, handletextpad=0.3, columnspacing=1.0)
ax[2].text(-0.26, 1.12, "C", transform=ax[2].transAxes, fontweight="bold", fontsize=10, va="bottom"); ax[2].set_title("Human test set", loc="left", pad=6); ax[2].spines[["top", "right"]].set_visible(False)
fig.tight_layout(w_pad=1.6); AS.save_both(fig, os.path.join(FIG, "Figure_S12"))
w2.provenance(HERE, "4_visualisation/supplementary_figures/Figure_S12/Figure_S12.py",
    [os.path.abspath(f) for f in ["b9_random_split_per_fit.csv", "b9_summary.csv"] + FIGOUT],
    [os.path.join(SRC, "b9_random_split_per_fit.csv"), os.path.join(SRC, "b9_summary.csv")], note="Fig. S12 (H3b, v3 negatives); TFs %s" % ", ".join(_exp))
print(m.round(4).to_string())
