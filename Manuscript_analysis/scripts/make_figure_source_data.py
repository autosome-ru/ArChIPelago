#!/usr/bin/env python3
"""Source-data CSVs for the figures that are drawn by the R scripts Figure_2_H_H.R, Figure_S1_H_M.R, Figure_3_H_H.R, Figure_S3_H_M.R.
Values come from Manuscript_analysis/HUMAN_MOUSE_total_100k.csv (human test set = chr1/8/21; mouse test set =
mouse chr1/8/19; baseline = best single monoPWM). Each CSV holds exactly the values that are plotted, so that
every panel can be checked without re-running R. Written to Manuscript_analysis/Figures/source_data/.
The two single-best-PWM reference series of Fig. 2 / S1 carry one row per TF (the results table repeats them
once per algorithm; the plotters draw one dot per TF for both series).
"""
import os
import pandas as pd

from archi_paths import RES

SRC = os.path.join(RES, "Figures", "source_data")
os.makedirs(SRC, exist_ok=True)

F = pd.read_csv(os.path.join(RES, "HUMAN_MOUSE_total_100k.csv"), sep="\t")

# ---- Fig. 2 (human test) and Fig. S1 (mouse test): all models x PWM sets, plus the baseline lines
for sp, name in (("H", "Figure_2"), ("M", "Figure_S1")):
    d = F[F.Seq_count > 100][["TF_name", "Model", "PWM", f"roc_auc_test_{sp}", f"pr_auc_test_{sp}", "Count", "Seq_count"]].copy()
    d = d.rename(columns={f"roc_auc_test_{sp}": "auROC", f"pr_auc_test_{sp}": "auPRC", "PWM": "PWM_set",
                          "Count": "n_monoPWMs", "Seq_count": "n_training_positives"})
    d["test_set"] = "human" if sp == "H" else "mouse"
    d = d.drop_duplicates(["TF_name", "Model", "PWM_set"])      # one dot per TF for the two reference series
    assert len(d) == 36 * (5 * 3 + 2), len(d)
    d.to_csv(os.path.join(SRC, f"{name}_source_data.csv"), index=False)
    ref = d[d.Model.isin(["Single best mono PWM", "Single best di PWM"])].groupby("Model")[["auROC", "auPRC"]].median()
    ref.to_csv(os.path.join(SRC, f"{name}_reference_lines.csv"))

# ---- Fig. 3 (human) and Fig. S3 (mouse): RF per PWM set vs its baseline + the delta panel
for sp, name in (("H", "Figure_3"), ("M", "Figure_S3")):
    r = F[F.Model == "RandomForestClassifier"].copy()
    n_pwm = F[F.Model == "RandomForestClassifier"].groupby("TF_name").Count.apply(lambda s: sum(set(s)))
    out = r[["TF_name", "PWM", f"roc_auc_test_{sp}_PWM", f"roc_auc_test_{sp}", f"pr_auc_test_{sp}_PWM", f"pr_auc_test_{sp}", "Seq_count"]].copy()
    out.columns = ["TF", "PWM_set", "auROC_best_PWM", "auROC_ArChIPelago", "auPRC_best_PWM", "auPRC_ArChIPelago", "n_training_positives"]
    out["delta_auROC"] = out.auROC_ArChIPelago - out.auROC_best_PWM
    out["delta_auPRC"] = out.auPRC_ArChIPelago - out.auPRC_best_PWM
    out["n_PWMs_total_mono_plus_di"] = out.TF.map(n_pwm)
    out["test_set"] = "human" if sp == "H" else "mouse"
    out.sort_values(["TF", "PWM_set"]).to_csv(os.path.join(SRC, f"{name}_source_data.csv"), index=False)

# ---- Fig. 4 / S4: written by analysis/fig4/assemble_fig4.py, not here.
print("source data written to", SRC)
