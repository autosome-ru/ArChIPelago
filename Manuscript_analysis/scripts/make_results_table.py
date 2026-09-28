#!/usr/bin/env python3
"""Results table of the manuscript (Manuscript_analysis/HUMAN_MOUSE_total_100k.csv, 900 x 47).

Base            : the table written by notebook 4 (analysis/inputs/HUMAN_MOUSE_total_100k_notebook4.csv).
Taken from the recomputed table (analysis/results_table/HUMAN_MOUSE_total_100k_recomputed.csv):
  1. every MOUSE column of every row
     (mouse test set = chromosomes 1, 8 and 19; models trained on the human training set, evaluated on
      the aligned mouse feature matrix; mouse baselines selected on the human training set by PWM identity);
  2. the 'Single best di PWM' rows and all *_PWM_di columns  <-  the same table
     (the best diPWM of each TF, selected on the human training set);
  3. the baseline columns used by the figures, *_PWM, are the best single monoPWM (*_PWM_mono) of the TF
     in every row -- the common reference for models built on monoPWMs, on diPWMs and on both, as the
     figure legends state.
Kept            : every HUMAN value of every ArChIPelago model (notebook 4 run), the human best-monoPWM
                  baselines, Count, Seq_count, Names, Unnamed: 0, row order.

Output: Manuscript_analysis/HUMAN_MOUSE_total_100k.csv (+ a provenance report on stdout).
"""
import os
import numpy as np
import pandas as pd

from archi_paths import INPUTS, BASIS, RES

FIX = os.path.join(BASIS, "results_table", "HUMAN_MOUSE_total_100k_recomputed.csv")
NB4 = os.path.join(INPUTS, "HUMAN_MOUSE_total_100k_notebook4.csv")
OUT = os.path.join(RES, "HUMAN_MOUSE_total_100k.csv")

KEY = ["TF_name", "Model", "PWM"]
pub = pd.read_csv(NB4, sep="\t")
fix = pd.read_csv(FIX, sep="\t")
print("notebook 4 table:", NB4)
print("recomputed table:", FIX)
assert pub.shape == fix.shape == (900, 47)
assert list(pub.columns) == list(fix.columns)
assert (pub[KEY].values == fix[KEY].values).all(), "row order must be identical"

out = pub.copy()
prov = {c: "notebook 4" for c in pub.columns}

# ---------------------------------------------------------------- 1. all mouse columns
mouse_cols = [c for c in pub.columns if "test_M" in c]
for c in mouse_cols:
    out[c] = fix[c].values
    prov[c] = "recomputed table (mouse test set = chr1/8/19; baseline by PWM identity)"

# ---------------------------------------------------------------- 2. diPWM baseline (human side)
di_cols_H = [c for c in pub.columns if c.endswith("_PWM_di") and "test_M" not in c]
for c in di_cols_H:
    out[c] = fix[c].values
    prov[c] = "best diPWM selected on the training set"
is_di_row = out.Model == "Single best di PWM"
for c in ["roc_auc_test_H", "pr_auc_test_H", "roc_auc_train_H", "pr_auc_train_H"]:
    out.loc[is_di_row, c] = fix.loc[is_di_row, c].values
    prov[c] = "notebook 4 for the ArChIPelago models; 'Single best di PWM' rows from the recomputed table"
for c in [c for c in pub.columns if c.startswith(("mean_", "median_", "std_"))]:
    out.loc[is_di_row, c] = fix.loc[is_di_row, c].values
    if "test_M" not in c:
        prov[c] = "notebook 4 for the ArChIPelago models; 'Single best di PWM' rows from the recomputed table"
# the mean_/median_/std_ columns are not read by any figure script; on the mouse side they carry the
# point value of the mouse basis

# ---------------------------------------------------------------- 3. baselines, one convention for
# every row: *_PWM_mono = best single monoPWM of the TF, *_PWM_di = best single diPWM of the TF,
# *_PWM (the column the figures compare against) = best single monoPWM, as the figure legends state.
mono_rows = out[(out.Model == "Single best mono PWM") & (out.PWM == "mono")].drop_duplicates("TF_name").set_index("TF_name")
di_rows = fix[(fix.Model == "Single best di PWM") & (fix.PWM == "di")].drop_duplicates("TF_name").set_index("TF_name")
for met in ["roc_auc", "pr_auc"]:
    for sp in ["train_H", "test_H", "test_M"]:
        base_mono = out.TF_name.map(mono_rows[f"{met}_{sp}_PWM_mono"])
        base_di = out.TF_name.map(di_rows[f"{met}_{sp}_PWM_di"])
        out[f"{met}_{sp}_PWM_mono"] = base_mono.values
        out[f"{met}_{sp}_PWM_di"] = base_di.values
        prov[f"{met}_{sp}_PWM_mono"] = ("best single monoPWM of the TF in every row "
                                        + ("(notebook 4 values)" if sp != "test_M" else "(recomputed table, chr1/8/19 test set)"))
        prov[f"{met}_{sp}_PWM_di"] = "best single diPWM of the TF in every row (selected on the training set)"
    for sp in ["test_H", "test_M"]:
        out[f"{met}_{sp}_PWM"] = out.TF_name.map(mono_rows[f"{met}_{sp}_PWM_mono"]).values
        prov[f"{met}_{sp}_PWM"] = "best single monoPWM (the baseline of the figure legends; common reference for all model types)"

# ---------------------------------------------------------------- checks
rf_pub = pub[(pub.Model == "RandomForestClassifier") & (pub.PWM == "mono+di")].set_index("TF_name")
rf_out = out[(out.Model == "RandomForestClassifier") & (out.PWM == "mono+di")].set_index("TF_name")
assert np.allclose(rf_out.roc_auc_test_H, rf_pub.roc_auc_test_H), "human model values must stay as in notebook 4"
assert np.allclose(rf_out.pr_auc_test_H, rf_pub.pr_auc_test_H)
mono_pub = pub[(pub.Model == "Single best mono PWM") & (pub.PWM == "mono")].drop_duplicates("TF_name").set_index("TF_name")
mono_chk = out[(out.Model == "Single best mono PWM") & (out.PWM == "mono")].drop_duplicates("TF_name").set_index("TF_name")
assert np.allclose(mono_chk.roc_auc_test_H, mono_pub.roc_auc_test_H), "human monoPWM baseline must stay as in notebook 4"
assert np.allclose(rf_out.roc_auc_test_H_PWM, mono_pub.loc[rf_out.index, "roc_auc_test_H_PWM_mono"]), "Fig. 3 baseline must be the best monoPWM"
assert (out[KEY].values == pub[KEY].values).all() and out.shape == (900, 47)
assert np.allclose(out.Count, pub.Count) and np.allclose(out.Seq_count, pub.Seq_count)

os.makedirs(RES, exist_ok=True)
out.to_csv(OUT, sep="\t", index=False)

# ---------------------------------------------------------------- report
print("written:", OUT)
print("\ncolumn provenance")
for c in pub.columns:
    changed = not np.allclose(pd.to_numeric(out[c], errors="coerce"), pd.to_numeric(pub[c], errors="coerce"), equal_nan=True) \
        if pd.api.types.is_numeric_dtype(pub[c]) else not (out[c].astype(str) == pub[c].astype(str)).all()
    print(f"  {c:28s} {'CHANGED' if changed else 'notebook 4':14s} {prov[c]}")

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
