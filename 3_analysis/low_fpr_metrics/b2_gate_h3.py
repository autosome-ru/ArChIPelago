"""B2 step-0 gate on the final data.
Recomputes full auROC / auPRC (PRROC integral) from the final-run score npz with the b2_v2 library and compares them with
  (i)   hm_mm_results.csv of the final run (the source of the npz),
  (ii)  the final results table (3_analysis/package/results_table.csv): RF mono+di seed 0 and *_PWM, which is the train-selected
        best monoPWM (asserted here),
  (iii) the b2_v2 B2_metrics_wide.csv that every B2 table and figure reads (auroc, auprc; HH, HM, MM; RF_s0, mono_byroc/byprc).
Scenarios HH (human-trained, human test), HM (human-trained, mouse test), MM (mouse-trained, mouse test), as in b2_compute.py.
Gate rule (as round 1): identity checks pass at max |diff| <= TOL_ID; the RF vs results-table check passed in round 1 only at
<= TOL_RF per TF (it failed, accepted, because the published forest was not the npz forest). AP vs integral sizes the AP difference.
Outputs (3_analysis/low_fpr_metrics/): B2_gate_h3_per_tf.csv, B2_gate_h3_summary.csv (check, max/median |diff|, tolerance, decision,
round-1 decision), _provenance_h3.yaml.
usage: python3 3_analysis/low_fpr_metrics/b2_gate_h3.py --scores DIR_WITH_35_NPZ [--scores-origin LABEL]
"""
import argparse
import glob
import importlib.util
import os
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "integration/"))
import w2  # noqa: E402
import h3_inputs as H  # noqa: E402

R2 = w2.R2
TOL_ID, TOL_RF = 1e-9, 1e-3
ap = argparse.ArgumentParser()
ap.add_argument("--scores", required=True)
ap.add_argument("--scores-origin", default="")
a = ap.parse_args()
LIB = os.path.join(H.IN["b2_dir"], "b2_lib.py")
_s = importlib.util.spec_from_file_location("b2_lib_v2", LIB)
L = importlib.util.module_from_spec(_s)
_s.loader.exec_module(L)
IN = {"hm": os.path.join(R2, H.IN["hm_mm"]), "rt": os.path.join(w2.B13, "../package", "results_table.csv"),
      "b2": H.IN["b2_wide"], "lib": LIB,
      "round1": os.path.join(R2, "results/datasets/round1_b2_gate", "B2_gate_summary.csv")}
hm = pd.read_csv(IN["hm"]).set_index("TF")
t = pd.read_csv(IN["rt"], sep="\t")
rf = t[(t.Model == "RandomForestClassifier") & (t.PWM == "mono+di")].set_index("TF_name")
mono = t[(t.Model == "Single best mono PWM") & (t.PWM == "mono")].drop_duplicates("TF_name").set_index("TF_name").loc[rf.index]
for m in ("roc_auc", "pr_auc"):
    if not np.array_equal(rf[m + "_test_H_PWM"].values, mono[m + "_test_H_PWM_mono"].values):
        raise SystemExit("CHECK FAILED results_table %s_test_H_PWM is not the best monoPWM (R2-F2)" % m)
b2 = pd.read_csv(IN["b2"]).set_index(["TF", "scenario", "model"])
files = sorted(glob.glob(os.path.join(a.scores, "*.npz")))
w2.PD.check_tfs([os.path.basename(f)[:-4] for f in files], "score npz")
SCEN = {"HH": ("y_test_H", "rf_H_s0", "pwmH_by%s_test_H", "HM", "test_H"),
        "HM": ("y_test_M", "rf_M_s0", "pwmH_by%s_test_M", "HM", "test_M"),
        "MM": ("y_test_M", "rfMM_M_s0", "pwmM_by%s_test_M", "MM", "test_M")}
rows = []
for f in files:
    tf = os.path.basename(f)[:-4]
    z = np.load(f)   # numeric arrays only (labels and scores); no pickled objects are read
    for sc, (ky, krf, kb, pre, ts) in SCEN.items():
        y = z[ky]
        r = {"TF": tf, "scenario": sc}
        for nm, key, href in (("rf_s0", krf, "%s_rf_%s_%%s_s0" % (pre, ts)), ("bROC", kb % "roc", "%s_broc_%s_%%s" % (pre, ts)),
                              ("bPRC", kb % "prc", "%s_bprc_%s_%%s" % (pre, ts))):
            s = z[key]
            r[nm + "_auroc"] = roc_auc_score(y, s)
            r[nm + "_auroc_mine"] = L.all_metrics(y, s)["auroc"] if nm == "rf_s0" else np.nan
            r[nm + "_auprc_int"] = L.auprc_integral(y, s)
            r[nm + "_auprc_ap"] = average_precision_score(y, s)
            r["ref_hm_%s_auroc" % nm] = hm.loc[tf, href % "auroc"]
            r["ref_hm_%s_auprc" % nm] = hm.loc[tf, href % "auprc"]
            bm = {"rf_s0": "RF_s0", "bROC": "mono_byroc", "bPRC": "mono_byprc"}[nm]
            r["ref_b2_%s_auroc" % nm] = b2.loc[(tf, sc, bm), "auroc"]
            r["ref_b2_%s_auprc" % nm] = b2.loc[(tf, sc, bm), "auprc"]
        if sc == "HH":
            r["rt_rf_auroc"], r["rt_rf_auprc"] = rf.loc[tf, "roc_auc_test_H"], rf.loc[tf, "pr_auc_test_H"]
            r["rt_base_auroc"], r["rt_base_auprc"] = rf.loc[tf, "roc_auc_test_H_PWM"], rf.loc[tf, "pr_auc_test_H_PWM"]
        rows.append(r)
g = pd.DataFrame(rows)
r1 = pd.read_csv(IN["round1"]).set_index("check")


def r1dec(lab):
    """Round-1 decision of the matching check (round-1 label), by the same tolerance rule."""
    if lab is None or lab not in r1.index:
        return "-"
    tol = TOL_RF if "PUBLISHED RF" in lab else (np.inf if "sklearn AP" in lab else TOL_ID)
    return "pass" if r1.loc[lab, "max_abs_diff"] <= tol else "FAIL"


CHK = [  # (label, scenario, column a, column b, tolerance, round-1 label)
    ("auROC sklearn vs own trapezoid (RF s0)", "HH", "rf_s0_auroc", "rf_s0_auroc_mine", TOL_ID, "auROC sklearn vs own trapezoid (RF s0)")]
for sc in SCEN:
    for nm in ("rf_s0", "bROC", "bPRC"):
        lab1 = {"rf_s0": "RF s0", "bROC": "bROC", "bPRC": "bPRC"}[nm]
        CHK += [("%s auROC %s vs hm_mm_results (H3a)" % (sc, lab1), sc, nm + "_auroc", "ref_hm_%s_auroc" % nm, TOL_ID,
                 "auROC %s vs hm_mm_results" % lab1 if sc == "HH" else None),
                ("%s auPRC(integral) %s vs hm_mm_results (H3a)" % (sc, lab1), sc, nm + "_auprc_int", "ref_hm_%s_auprc" % nm, TOL_ID,
                 "auPRC(integral) %s vs hm_mm_results" % lab1 if sc == "HH" else None),
                ("%s auROC %s vs b2_v2 B2_metrics_wide" % (sc, lab1), sc, nm + "_auroc", "ref_b2_%s_auroc" % nm, TOL_ID, None),
                ("%s auPRC(integral) %s vs b2_v2 B2_metrics_wide" % (sc, lab1), sc, nm + "_auprc_int", "ref_b2_%s_auprc" % nm, TOL_ID, None)]
CHK += [("HH auROC RF s0 vs results table RF (H3)", "HH", "rf_s0_auroc", "rt_rf_auroc", TOL_RF, "auROC RF s0 vs PUBLISHED RF"),
        ("HH auPRC(integral) RF s0 vs results table RF (H3)", "HH", "rf_s0_auprc_int", "rt_rf_auprc", TOL_RF, "auPRC(integral) RF s0 vs PUBLISHED RF"),
        ("HH auROC bROC vs results table best monoPWM (roc_auc_test_H_PWM)", "HH", "bROC_auroc", "rt_base_auroc", TOL_ID,
         "auROC bROC vs PUBLISHED single best mono (test_H_PWM)"),
        ("HH auPRC(integral) bPRC vs results table best monoPWM (pr_auc_test_H_PWM)", "HH", "bPRC_auprc_int", "rt_base_auprc", TOL_ID,
         "auPRC(integral) bPRC vs PUBLISHED single best mono (pr_auc_test_H_PWM)"),
        ("HH auPRC integral vs sklearn AP, RF s0 (size only)", "HH", "rf_s0_auprc_int", "rf_s0_auprc_ap", np.inf, "auPRC integral vs sklearn AP, RF s0"),
        ("HH auPRC integral vs sklearn AP, bPRC (size only)", "HH", "bPRC_auprc_int", "bPRC_auprc_ap", np.inf, "auPRC integral vs sklearn AP, bPRC")]
out = []
for lab, sc, ca, cb, tol, lab_r1 in CHK:
    d = (g.loc[g.scenario == sc, ca] - g.loc[g.scenario == sc, cb]).abs()
    if d.isna().any():
        raise SystemExit("CHECK FAILED %s: %d missing values" % (lab, d.isna().sum()))
    out.append(dict(check=lab, n_TF=len(d), max_abs_diff=d.max(), median_abs_diff=d.median(), tolerance=tol,
                    decision=("pass" if d.max() <= tol else "FAIL") if np.isfinite(tol) else "size only",
                    round1_check=lab_r1 or "-", round1_max_abs_diff=r1.loc[lab_r1, "max_abs_diff"] if lab_r1 in r1.index else np.nan,
                    round1_decision=r1dec(lab_r1) if np.isfinite(tol) else ("size only" if lab_r1 else "-")))
c = pd.DataFrame(out)
fo, fs = os.path.join(HERE, "B2_gate_h3_per_tf.csv"), os.path.join(HERE, "B2_gate_h3_summary.csv")
for f in (fo, fs):
    w2.backup(f)
g.to_csv(fo, index=False)
c.to_csv(fs, index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 80)
print(c[["check", "max_abs_diff", "decision", "round1_decision"]].to_string(index=False))
hh = g[g.scenario == "HH"]
print("HH medians: RF auROC %.4f auPRC %.4f | bROC auROC %.4f | bPRC auPRC %.4f" % (hh.rf_s0_auroc.median(), hh.rf_s0_auprc_int.median(),
      hh.bROC_auroc.median(), hh.bPRC_auprc_int.median()))
w2.provenance(HERE, "3_analysis/low_fpr_metrics/b2_gate_h3.py", [fo, fs], list(IN.values()) + files,
              params="{TOL_ID: %g, TOL_RF: %g, scores: '%s', scores_origin: '%s'}" % (TOL_ID, TOL_RF, a.scores, a.scores_origin))
