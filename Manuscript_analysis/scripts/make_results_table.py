#!/usr/bin/env python3
"""Results table of the manuscript (Manuscript_analysis/results_table.csv, 900 x 47).

Two sources with the same rows and columns:
  notebook 4        analysis/inputs/notebook4_results_table.csv (the table written by notebook 4)
  model evaluation  analysis/model_evaluation/evaluation_table.csv (analysis/model_evaluation/)
From the notebook 4 table: every HUMAN value of every ArChIPelago model, the human best-monoPWM baselines,
  Count, Seq_count, Names, Unnamed: 0, row order.
From the model evaluation table:
  1. every MOUSE column of every row
     (mouse test set = chromosomes 1, 8 and 19; models trained on the human training set, evaluated on
      the aligned mouse feature matrix; mouse baselines selected on the human training set by PWM identity);
  2. the 'Single best di PWM' rows and all *_PWM_di columns
     (the best diPWM of each TF, selected on the human training set).
The baseline columns used by the figures, *_PWM, are the best single monoPWM (*_PWM_mono) of the TF in every
row: the common reference for models built on monoPWMs, on diPWMs and on both, as the figure legends state.

Output: Manuscript_analysis/results_table.csv (+ the headline medians on stdout).
"""
import os
import numpy as np
import pandas as pd

from archi_paths import INPUTS, BASIS, RES

EVAL = os.path.join(BASIS, "model_evaluation", "evaluation_table.csv")
NB4 = os.path.join(INPUTS, "notebook4_results_table.csv")
OUT = os.path.join(RES, "results_table.csv")

KEY = ["TF_name", "Model", "PWM"]
nb4 = pd.read_csv(NB4, sep="\t")
ev = pd.read_csv(EVAL, sep="\t")
assert nb4.shape == ev.shape == (900, 47)
assert list(nb4.columns) == list(ev.columns)
assert (nb4[KEY].values == ev[KEY].values).all(), "row order must be identical"

out = nb4.copy()

# ---------------------------------------------------------------- 1. all mouse columns
for c in [c for c in nb4.columns if "test_M" in c]:
    out[c] = ev[c].values

# ---------------------------------------------------------------- 2. diPWM baseline (human side)
for c in [c for c in nb4.columns if c.endswith("_PWM_di") and "test_M" not in c]:
    out[c] = ev[c].values
is_di_row = out.Model == "Single best di PWM"
for c in ["roc_auc_test_H", "pr_auc_test_H", "roc_auc_train_H", "pr_auc_train_H"]:
    out.loc[is_di_row, c] = ev.loc[is_di_row, c].values
for c in [c for c in nb4.columns if c.startswith(("mean_", "median_", "std_"))]:
    out.loc[is_di_row, c] = ev.loc[is_di_row, c].values
# the mean_/median_/std_ columns are not read by any figure script; on the mouse side they carry the
# point value of the mouse basis

# ---------------------------------------------------------------- 3. baselines, one convention for
# every row: *_PWM_mono = best single monoPWM of the TF, *_PWM_di = best single diPWM of the TF,
# *_PWM (the column the figures compare against) = best single monoPWM, as the figure legends state.
mono_rows = out[(out.Model == "Single best mono PWM") & (out.PWM == "mono")].drop_duplicates("TF_name").set_index("TF_name")
di_rows = ev[(ev.Model == "Single best di PWM") & (ev.PWM == "di")].drop_duplicates("TF_name").set_index("TF_name")
for met in ["roc_auc", "pr_auc"]:
    for sp in ["train_H", "test_H", "test_M"]:
        out[f"{met}_{sp}_PWM_mono"] = out.TF_name.map(mono_rows[f"{met}_{sp}_PWM_mono"]).values
        out[f"{met}_{sp}_PWM_di"] = out.TF_name.map(di_rows[f"{met}_{sp}_PWM_di"]).values
    for sp in ["test_H", "test_M"]:
        out[f"{met}_{sp}_PWM"] = out.TF_name.map(mono_rows[f"{met}_{sp}_PWM_mono"]).values

# ---------------------------------------------------------------- checks of the merge
rf_nb4 = nb4[(nb4.Model == "RandomForestClassifier") & (nb4.PWM == "mono+di")].set_index("TF_name")
rf_out = out[(out.Model == "RandomForestClassifier") & (out.PWM == "mono+di")].set_index("TF_name")
assert np.allclose(rf_out.roc_auc_test_H, rf_nb4.roc_auc_test_H), "human model values come from notebook 4"
assert np.allclose(rf_out.pr_auc_test_H, rf_nb4.pr_auc_test_H)
mono_nb4 = nb4[(nb4.Model == "Single best mono PWM") & (nb4.PWM == "mono")].drop_duplicates("TF_name").set_index("TF_name")
mono_out = out[(out.Model == "Single best mono PWM") & (out.PWM == "mono")].drop_duplicates("TF_name").set_index("TF_name")
assert np.allclose(mono_out.roc_auc_test_H, mono_nb4.roc_auc_test_H), "human monoPWM baseline comes from notebook 4"
assert np.allclose(rf_out.roc_auc_test_H_PWM, mono_nb4.loc[rf_out.index, "roc_auc_test_H_PWM_mono"]), "Fig. 3 baseline must be the best monoPWM"
assert (out[KEY].values == nb4[KEY].values).all() and out.shape == (900, 47)
assert np.allclose(out.Count, nb4.Count) and np.allclose(out.Seq_count, nb4.Seq_count)

os.makedirs(RES, exist_ok=True)
out.to_csv(OUT, sep="\t", index=False)
print("written:", OUT)

# ---------------------------------------------------------------- headline medians
rf_out = out[(out.Model == "RandomForestClassifier") & (out.PWM == "mono+di")].set_index("TF_name").sort_index()
print("\nheadline medians (36 TFs, Random Forest on monoPWMs+diPWMs, baseline = best single monoPWM)")
for sp, lab in (("H", "human test set"), ("M", "mouse test set, chr1/8/19")):
    b_roc, b_prc = rf_out[f"roc_auc_test_{sp}_PWM"], rf_out[f"pr_auc_test_{sp}_PWM"]
    dr = rf_out[f"roc_auc_test_{sp}"] - b_roc
    dp = rf_out[f"pr_auc_test_{sp}"] - b_prc
    below = sorted(set(dr.index[dr < 0]) | set(dp.index[dp < 0]))
    print(f"  {lab}: RF {rf_out[f'roc_auc_test_{sp}'].median():.4f}/{rf_out[f'pr_auc_test_{sp}'].median():.4f}"
          f"  baseline {b_roc.median():.4f}/{b_prc.median():.4f}"
          f"  median delta {dr.median():+.4f}/{dp.median():+.4f}  improved {(dr > 0).sum()}/36, {(dp > 0).sum()}/36  below: {below or 'none'}")
