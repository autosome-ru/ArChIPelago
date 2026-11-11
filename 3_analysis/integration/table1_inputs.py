"""Table 1 inputs, human and mouse rows, on the final data (35 TFs).
Convention (reproduces the first-submission Table 1 from the round 1 files, see --round1): baseline = best monoPWM selected on
the training set (by auROC for auROC, by auPRC for auPRC); two-sided Wilcoxon signed-rank on the per-TF differences.
  Human, human-trained   3_analysis/package/results_table.csv, RF mono+di (seed 0) vs *_PWM (= best monoPWM, checked)
  Mouse, human-trained   hm_mm HM_rf_test_M_*_<HM_SEED> vs HM_b*_test_M_*   (mean of the 2 seeds)
  Mouse, mouse-trained   hm_mm MM_rf_test_M_*_<MM_SEED> (mean of 2 seeds) vs MM_b*_test_M_*
  The human rows stay seed 0 (main analysis); --round1 keeps HM seed 0, as the first-submission Table 1 was built.
--round1 recomputes the old mouse rows from the round 1 hm_mm file (0.861 0.874 35/36 3e-10 ...) and writes nothing.
Output: 3_analysis/integration/table1_rows.csv
usage: ARCHI_STAGE=H5 python3 3_analysis/integration/table1_inputs.py [--round1]
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w2  # noqa: E402
import h3_inputs as H  # noqa: E402

R1 = os.path.join(w2.R2, "results/round1_package/analysis_mouse_swap/hm_mm_results.csv")
RT = os.path.join(w2.B13, "../package", "results_table.csv")
ROUND1 = "--round1" in sys.argv
HM_SEED = "s0" if ROUND1 else "mean"   # rows of the mouse test set use the mean of seeds 0 and 1
MM_SEED = "mean"
src = R1 if ROUND1 else H.IN["hm_mm"]
h = pd.read_csv(src)
n = len(h) if ROUND1 else w2.PD.check_tfs(h.TF, "hm_mm_results")
rows = []


def add(lab, r, bb, rcol, bcol, tfs):
    d = r.values - bb.values
    rows.append(dict(row=lab, baseline_median=round(float(np.median(bb)), 3), archipelago_median=round(float(np.median(r)), 3),
                     tfs_improved="%d/%d" % ((d > 0).sum(), len(d)), wilcoxon_p="%.0e" % wilcoxon(d).pvalue,
                     below=",".join(np.asarray(tfs)[d < 0]), rf_column=rcol, baseline_column=bcol))


if not ROUND1:
    t = pd.read_csv(RT, sep="\t")
    rf = t[(t.Model == "RandomForestClassifier") & (t.PWM == "mono+di")].set_index("TF_name").sort_index()
    mono = t[(t.Model == "Single best mono PWM") & (t.PWM == "mono")].drop_duplicates("TF_name").set_index("TF_name").loc[rf.index]
    for m, ml in (("roc_auc", "auROC"), ("pr_auc", "auPRC")):
        if not np.array_equal(rf[m + "_test_H_PWM"].values, mono[m + "_test_H_PWM_mono"].values):
            raise SystemExit("CHECK FAILED results_table %s_test_H_PWM is not the best monoPWM (R2-F2)" % m)
        add("Human, human-trained, " + ml, rf[m + "_test_H"], rf[m + "_test_H_PWM"], m + "_test_H", m + "_test_H_PWM", rf.index)
for sc, lab, v in (("HM", "Mouse, human-trained", HM_SEED), ("MM", "Mouse, mouse-trained", MM_SEED)):
    for m, b, ml in (("auroc", "broc", "auROC"), ("auprc", "bprc", "auPRC")):
        add("%s, %s" % (lab, ml), h["%s_rf_test_M_%s_%s" % (sc, m, v)], h["%s_%s_test_M_%s" % (sc, b, m)], "%s_rf_test_M_%s_%s" % (sc, m, v),
            "%s_%s_test_M_%s" % (sc, b, m), h.TF)
out = pd.DataFrame(rows)
print(out.to_string(index=False))
if not ROUND1:
    f = os.path.join(w2.B13, "table1_rows.csv")
    w2.backup(f)
    out.to_csv(f, index=False)
    w2.provenance(w2.B13, "3_analysis/integration/table1_inputs.py", [f], [RT, src],
                  params="{human: RF seed 0 vs best monoPWM; HM: %s; MM: %s; Wilcoxon two-sided}" % (HM_SEED, MM_SEED))
