#!/usr/bin/env python3
"""Source-data CSVs of the figures drawn by the R scripts (Figures/source_data/); each CSV holds exactly the
plotted values, so every panel can be checked without re-running R.
  Fig. 2 / S1, 3 / S3   from Manuscript_analysis/results_table.csv (human test set = chr1/8/21; mouse test set =
                        mouse chr1/8/19; baseline = best single monoPWM). The two single-best-PWM reference
                        series of Fig. 2 / S1 carry one row per TF (the results table repeats them once per
                        algorithm; the plotters draw one dot per TF for both series).
  Fig. S5               from analysis/saturation/saturation_results.csv (reference = best single monoPWM, the
                        base_* columns of that file, asserted against the results table).
  Fig. S6               from Sup. Table 5 (make_sup_tables_3_to_6.py).
  Fig. S7               the four TFs drawn in Fig. S7, from analysis/motif_subtypes/Figure_S7_source_data.csv.
Fig. 4 / S4: written by analysis/slim_dichipmunk/assemble_Figure_4_S4.py; Fig. S2: by Figure_S2_RF_vs_best_monoPWM.py.
"""
import os
import pandas as pd

from archi_paths import BASIS, RES

SRC = os.path.join(RES, "Figures", "source_data")
os.makedirs(SRC, exist_ok=True)

F = pd.read_csv(os.path.join(RES, "results_table.csv"), sep="\t")

# Fig. 2 (human test) and Fig. S1 (mouse test): all models x PWM sets, plus the baseline lines
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

# Fig. 3 (human) and Fig. S3 (mouse): RF per PWM set vs its baseline + the delta panel
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

# Fig. S5: RF gain over the best single monoPWM vs the number of PWMs used
K_GRID = [1, 2, 4, 8, 16, 32, 64, 128]
METRICS = ["auroc_H", "auprc_H", "auroc_M", "auprc_M"]
df = pd.read_csv(os.path.join(BASIS, "saturation", "saturation_results.csv"))
tfs = sorted(df.tf.unique())
Pmap = df.groupby("tf").P.first()
for m in METRICS:
    df["d_" + m] = df[m] - df["base_" + m]
F_RF = F[(F.Model == "RandomForestClassifier") & (F.PWM == "mono+di")].set_index("TF_name")
for tf in tfs:      # the reference must agree with the *_PWM columns of the results table (best monoPWM)
    b = df[df.tf == tf].iloc[0]
    assert abs(b.base_auroc_H - F_RF.loc[tf, "roc_auc_test_H_PWM"]) < 1e-6, tf
    assert abs(b.base_auprc_M - F_RF.loc[tf, "pr_auc_test_M_PWM"]) < 1e-6, tf
dcols = ["d_" + m for m in METRICS]
mean = df.groupby(["tf", "design", "k"])[dcols].mean().reset_index()
full = mean[(mean.design == "random") & (mean.k == mean.tf.map(Pmap))].set_index("tf")


def curve(tf, design):
    P = Pmap[tf]
    sub = mean[(mean.tf == tf) & (mean.design == design)].set_index("k")
    out = {k: (full.loc[tf, dcols] if min(k, P) == P else sub.loc[min(k, P), dcols]) for k in K_GRID}
    out["P"] = full.loc[tf, dcols]
    return pd.DataFrame(out).T.astype(float)


curves = {des: {tf: curve(tf, des) for tf in tfs} for des in ("random", "topk")}
med = {des: pd.concat(curves[des].values(), keys=tfs).groupby(level=1).median().loc[K_GRID + ["P"]]
       for des in ("random", "topk")}
rows = []
for des in ("random", "topk"):
    for tf in tfs:
        c = curves[des][tf]
        for k in K_GRID + ["P"]:
            rows.append(dict(design=des, TF=tf, P=Pmap[tf], k=(Pmap[tf] if k == "P" else min(k, Pmap[tf])),
                             k_requested=k, **{f"delta_{m}": c.loc[k, "d_" + m] for m in METRICS}))
pd.DataFrame(rows).to_csv(os.path.join(SRC, "Figure_S5_source_data.csv"), index=False)
med_out = pd.concat({des: med[des] for des in ("random", "topk")}, names=["design", "k"])
med_out.to_csv(os.path.join(SRC, "Figure_S5_median_curves.csv"))

# Fig. S6: cross-species transfer and the mouse-trained control (extract of Sup. Table 5)
t5 = pd.read_csv(os.path.join(RES, "Sup_Tables", "Sup_Table_5_cross_species_and_mouse_trained.csv")).set_index("TF")
t5[["H>H: dauROC (reference)", "H>H: dauPRC (reference)", "H>M: dauROC", "H>M: dauPRC", "M>M: dauROC", "M>M: dauPRC",
    "Similarity of the baseline human monoPWM to the nearest mouse monoPWM (Pearson r of aligned columns)",
    "H>M: below baseline on >=1 metric"]].to_csv(os.path.join(SRC, "Figure_S6_source_data.csv"))

# Fig. S7: the four TFs of the figure
s7 = pd.read_csv(os.path.join(BASIS, "motif_subtypes", "Figure_S7_source_data.csv"))
s7[s7.TF.isin(["E2F4", "RXRA", "TAL1", "TFE2"])].to_csv(os.path.join(SRC, "Figure_S7_source_data.csv"), index=False)

print("source data written to", SRC)
