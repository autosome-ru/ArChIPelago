"""Sup. Table 5 of the revision, full rebuild.

Producer: make_sup_tables_3_to_6.py of 3_analysis/package (35 TFs, the final-run results of h3_inputs.py and the cross-species
table rebuilt from the used peak sets, human 361 / mouse 276). This script copies that workbook into outputs/supplementary_tables
after checking it against Sup. Tables 1 and 2: 35 reference TFs, experiment counts per TF = Sup. Table 1 sheets, positives =
Sup. Table 2, H>M values = hm_mm_results. H>M ArChIPelago columns are labelled seed 0.
Rows of the mouse test set use the mean of seeds 0 and 1:
  per_TF has the H>M two-seed mean columns next to the seed 0 ones (ArChIPelago auROC / auPRC, dauROC / dauPRC, below baseline);
       the seed 0 gain columns are labelled '(seed 0)'.
  summary: the human-trained row is labelled '(seed 0)' and a row with the two-seed mean follows it.
  correlations: the gains named there are the seed 0 columns (labels renamed with them); a second block repeats them with the
       two-seed mean.
  Wilcoxon P of the new summary row: scipy.stats.wilcoxon defaults, as the package producer (wp).
Column 'TF family or subfamily (short label)': TFClass labels are given at family level for some TFs (AP-2, EBF, IRF ...) and at
subfamily level for others (NR2B RXR, ETS (SPI) ...); the TFClass id column names the level.
usage: python3 3_analysis/supplementary_tables/make_sup_table_5_R2.py
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

R2 = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(R2, "3_analysis/integration"))
sys.path.insert(0, os.path.join(R2, "3_analysis", "integration"))
import p53drop as PD  # noqa: E402
import w2  # noqa: E402
import h3_inputs as H  # noqa: E402

SRC = os.path.join(R2, "3_analysis", "package", "Sup_Tables", "Sup_Table_5_cross_species_and_mouse_trained.xlsx")
OUT = os.path.join(R2, "outputs/supplementary_tables", "Sup_Table_5_cross_species_and_mouse_trained.xlsx")
T1 = os.path.join(R2, "outputs/supplementary_tables", "Sup_Table_1_ChIP-Seq_experiments.xlsx")
T2 = os.path.join(R2, "outputs/supplementary_tables", "Sup_Table_2_PWMs_and_datasets.xlsx")
HM = H.IN["hm_mm"]   # run_hm_mm results of the final run

t5 = pd.read_excel(SRC, sheet_name="per_TF").set_index("TF")
NTF = PD.check_tfs(t5.index, "Sup. Table 5 per_TF")
t1 = pd.read_excel(T1, sheet_name=None)
t2 = pd.read_excel(T2, sheet_name="per_TF").set_index("TF").loc[t5.index]
hm = pd.read_csv(HM).set_index("TF").loc[t5.index]
checks = [
    (t5["Human ChIP-Seq experiments"], t1["human"].groupby("TF").size().reindex(t5.index), "human experiments vs Sup. Table 1"),
    (t5["Mouse ChIP-Seq experiments"], t1["mouse"].groupby("TF").size().reindex(t5.index), "mouse experiments vs Sup. Table 1"),
    (t5["Human training positives"], t2[H.COL["HTR"]], "human training positives vs Sup. Table 2"),
    (t5["Mouse training positives"], t2[H.COL["MTR"]], "mouse training positives vs Sup. Table 2"),
    (t5["Mouse test positives"], t2[H.COL["MTE"]], "mouse test positives vs Sup. Table 2"),
]
for a, b, lab in checks:
    if not (a.values == b.values).all():
        raise SystemExit("CHECK FAILED Sup. Table 5 %s: %s" % (lab, list(a.index[a.values != b.values])))
if not np.allclose(t5["H>M: ArChIPelago auROC (mouse test, seed 0)"], hm.HM_rf_test_M_auroc_s0, atol=5e-5):
    raise SystemExit("CHECK FAILED Sup. Table 5 H>M auROC vs hm_mm_results")
if not np.allclose(t5["M>M: ArChIPelago auPRC (mouse test, mean of 2 seeds)"], hm.MM_rf_test_M_auprc_mean, atol=5e-5):
    raise SystemExit("CHECK FAILED Sup. Table 5 M>M auPRC vs hm_mm_results")
print("Sup. Table 5: %d TFs; experiments %d human / %d mouse; checks vs Sup. Tables 1, 2 and hm_mm_results pass"
      % (NTF, t5["Human ChIP-Seq experiments"].sum(), t5["Mouse ChIP-Seq experiments"].sum()))
FAMILY_OLD, FAMILY_NEW = "TF family (TFClass)", "TF family or subfamily (short label)"
S0, MEAN = "seed 0", "mean of 2 seeds"
x5 = pd.read_excel(SRC, sheet_name=None, header=None)   # raw cells, every sheet (correlations has a text block above its table)
per = pd.read_excel(SRC, sheet_name="per_TF")
summ = pd.read_excel(SRC, sheet_name="summary", dtype=str)
if list(x5) != ["summary", "correlations", "per_TF"]:
    raise SystemExit("CHECK FAILED Sup. Table 5 sheets: %s" % list(x5))
hmo = hm.loc[per.TF]
ren = {FAMILY_OLD: FAMILY_NEW, "H>M: dauROC": "H>M: dauROC (%s)" % S0, "H>M: dauPRC": "H>M: dauPRC (%s)" % S0,
       "H>M: below baseline on >=1 metric": "H>M: below baseline on >=1 metric (%s)" % S0}
missing = [c for c in ren if c not in per.columns]
if missing:
    raise SystemExit("CHECK FAILED Sup. Table 5 per_TF: columns to rename missing %s" % missing)
per = per.rename(columns=ren)
gain = {m: hmo["HM_rf_test_M_%s_mean" % m].values - hmo["HM_%s_test_M_%s" % (b, m)].values for m, b in (("auroc", "broc"), ("auprc", "bprc"))}
new = [("H>M: ArChIPelago auROC (mouse test, %s)" % S0, "H>M: ArChIPelago auROC (mouse test, %s)" % MEAN, hmo.HM_rf_test_M_auroc_mean.values),
       ("H>M: dauROC (%s)" % S0, "H>M: dauROC (%s)" % MEAN, gain["auroc"]),
       ("H>M: ArChIPelago auPRC (mouse test, %s)" % S0, "H>M: ArChIPelago auPRC (mouse test, %s)" % MEAN, hmo.HM_rf_test_M_auprc_mean.values),
       ("H>M: dauPRC (%s)" % S0, "H>M: dauPRC (%s)" % MEAN, gain["auprc"]),
       ("H>M: below baseline on >=1 metric (%s)" % S0, "H>M: below baseline on >=1 metric (%s)" % MEAN, (gain["auroc"] < 0) | (gain["auprc"] < 0))]
for after, col, v in new:
    per.insert(per.columns.get_loc(after) + 1, col, np.round(v, 4) if v.dtype.kind == "f" else v)
# summary: seed 0 row labelled, two-seed mean row after it (gains unrounded, as add5 of the pkg producer)
lab0 = "Human-trained model on the mouse test set vs the best monoPWM (selected on human training data)"
i0 = summ.index[summ.iloc[:, 0] == lab0]
if len(i0) != 1:
    raise SystemExit("CHECK FAILED Sup. Table 5 summary: human-trained row not found once")
summ.iloc[i0[0], 0] = lab0 + ", %s" % S0
dr, dp = pd.Series(gain["auroc"], index=per.TF), pd.Series(gain["auprc"], index=per.TF)
row = [lab0 + ", %s" % MEAN, "%g" % round(dr.median(), 4), str(int((dr > 0).sum())), "%.1e" % wilcoxon(dr).pvalue, "%g" % round(dp.median(), 4),
       str(int((dp > 0).sum())), "%.1e" % wilcoxon(dp).pvalue, ", ".join(sorted(set(dr.index[dr < 0]) | set(dp.index[dp < 0]))) or "none"]
summ = pd.concat([summ.iloc[:i0[0] + 1], pd.DataFrame([row], columns=summ.columns), summ.iloc[i0[0] + 1:]], ignore_index=True)
num = [c for c in summ.columns if c.startswith(("median", "TFs with", "Wilcoxon"))]
summ[num] = summ[num].astype(float)
corr = x5["correlations"].replace({"H>M: dauROC": "H>M: dauROC (%s)" % S0, "H>M: dauPRC": "H>M: dauPRC (%s)" % S0})
corr = corr.replace(to_replace=r"'H>M: dauROC' / 'H>M: dauPRC'", value="'H>M: dauROC (%s)' / 'H>M: dauPRC (%s)'" % (S0, S0), regex=True)
if (corr.astype(str) == "H>M: dauROC").any().any() or not corr.astype(str).apply(lambda c: c.str.contains("(seed 0)' /", regex=False)).any().any():
    raise SystemExit("CHECK FAILED Sup. Table 5 correlations: seed 0 labels not applied")
# the same correlations with the two-seed mean gains, as a second block after the
# seed-0 rows (Fig. S6C shows the first pair of this block). Properties from the unrounded pkg per_TF table, as the pkg producer;
# the seed-0 rows are recomputed first and must equal the sheet (rho to 3 dp, P to 3 significant digits).
UNR = os.path.join(R2, "3_analysis", "package", "Sup_Tables", "Sup_Table_5_per_TF_unrounded.csv")
u5 = pd.read_csv(UNR).set_index("TF").loc[per.TF]
PROP = {"similarity of the baseline human monoPWM": u5["Similarity of the baseline human monoPWM to the nearest mouse monoPWM (Pearson r of aligned columns)"],
        "median, over all human monoPWMs of the TF": u5["Median similarity human monoPWMs -> nearest mouse monoPWM"],
        "number of human training positives": u5["Human training positives"], "number of mouse training positives": u5["Mouse training positives"],
        "number of human ChIP-Seq experiments": u5["Human ChIP-Seq experiments"], "number of mouse ChIP-Seq experiments": u5["Mouse ChIP-Seq experiments"],
        "number of human PWMs (monoPWMs + diPWMs)": u5["Human monoPWMs"] + u5["Human diPWMs"],
        "auROC of the best human monoPWM on the mouse test set": u5["H>M: best monoPWM auROC"]}
G0 = {"H>M: dauROC (%s)" % S0: u5["H>M: dauROC"], "H>M: dauPRC (%s)" % S0: u5["H>M: dauPRC"]}
GM = {"H>M: dauROC (%s)" % S0: ("H>M: dauROC (%s)" % MEAN, dr), "H>M: dauPRC (%s)" % S0: ("H>M: dauPRC (%s)" % MEAN, dp)}
hdr = corr.index[corr.iloc[:, 0].astype(str) == "Cross-species gain (per TF)"]
if len(hdr) != 1:
    raise SystemExit("CHECK FAILED Sup. Table 5 correlations: header row not found once")
body = corr.loc[hdr[0] + 1:].dropna(how="all")
if len(body) != 2 * len(PROP):
    raise SystemExit("CHECK FAILED Sup. Table 5 correlations: %d rows, expected %d" % (len(body), 2 * len(PROP)))
mean_rows = []
for _, r_ in body.iterrows():
    key = [k for k in PROP if str(r_[1]).startswith(k)]
    if len(key) != 1 or r_[0] not in G0:
        raise SystemExit("CHECK FAILED Sup. Table 5 correlations: row not matched: %r / %r" % (r_[0], r_[1]))
    x = PROP[key[0]].values
    rho0, p0 = spearmanr(G0[r_[0]].values, x)
    if round(rho0, 3) != float(r_[3]) or "%.3g" % p0 != str(r_[4]):
        raise SystemExit("CHECK FAILED Sup. Table 5 correlations: seed-0 row not reproduced: %s %s (%.3f %.3g vs %s %s)" % (r_[0], key[0], rho0, p0, r_[3], r_[4]))
    lab, g = GM[r_[0]]
    rho, pv = spearmanr(g.loc[per.TF].values, x)
    mean_rows.append([lab, r_[1], int(r_[2]), round(rho, 3), "%.3g" % pv])
corr = pd.concat([corr, pd.DataFrame([[np.nan] * corr.shape[1]], columns=corr.columns),
                  pd.DataFrame(mean_rows, columns=corr.columns)], ignore_index=True)
for i_ in corr.index[:hdr[0]]:   # header text: name both blocks and the pair of Fig. S6C
    v = str(corr.iloc[i_, 0])
    if v.startswith("Cross-species gain = "):
        corr.iloc[i_, 0] = v.replace("of the per_TF sheet:", "of the per_TF sheet (first block) and 'H>M: dauROC (%s)' / 'H>M: dauPRC (%s)' (second block):" % (MEAN, MEAN), 1)
    elif v == "The first pair is the relation shown in Fig. S6C.":
        corr.iloc[i_, 0] = "Fig. S6C shows the first pair of the second block (mean of two seeds)."
if not (corr.astype(str) == "Fig. S6C shows the first pair of the second block (mean of two seeds).").any().any():
    raise SystemExit("CHECK FAILED Sup. Table 5 correlations: Fig. S6C note not updated")
print("Sup. Table 5 correlations, two-seed mean: " + "; ".join("%s %s rho %.3f P %s" % (m[0][5:12], m[1][:30], m[3], m[4]) for m in mean_rows[:2]))
print("Sup. Table 5 H>M two-seed mean: auROC %d/%d up (below: %s), auPRC %d/%d up (below: %s); median gain %.4f / %.4f" % (
    (dr > 0).sum(), len(dr), ",".join(dr.index[dr < 0]), (dp > 0).sum(), len(dp), ",".join(dp.index[dp < 0]), dr.median(), dp.median()))
w2.backup(OUT)
with pd.ExcelWriter(OUT) as xw:
    summ.to_excel(xw, sheet_name="summary", index=False)
    corr.to_excel(xw, sheet_name="correlations", index=False, header=False)
    per.to_excel(xw, sheet_name="per_TF", index=False)
w2.provenance(os.path.dirname(OUT), "3_analysis/supplementary_tables/make_sup_table_5_R2.py", [OUT], [SRC, T1, T2, HM, UNR],
              note="pkg producer output (make_sup_tables_3_to_6.py) after checks; V006 H>M two-seed mean columns and summary row, V081 family header")
