"""B3 analysis 2: per-TF RF gain vs number of training
positives / datasets / PWMs (Spearman), 35 TFs.
Inputs (all explicit, recorded in provenance):
  human gain  3_analysis/package/results_table.csv, RandomForestClassifier mono+di, test_H minus *_PWM, where *_PWM is the
              train-selected best monoPWM (make_results_table resets it; checked here against the 'Single best mono PWM' rows)
  mouse gain  hm_mm_results (final run): human-trained RF on the mouse test set, mean of seeds 0 and 1, minus the
              best monoPWM (selected on human train); the human-test columns of the retrained model stay seed 0
  positives   per_tf_positive_counts.csv of results/server/human_data and results/server/mouse_data (file "train" = human training, mouse file "control" = training)
  PWMs        hm_mm P_mono_H + P_di_H (monoPWMs + diPWMs), checked against the PWM counts of the final run
  datasets    used peak-set lists (human 361, mouse 276)
usage: python3 3_analysis/top_peaks/b3_gain_vs_npos_h3.py OUTDIR
"""
import os, sys, pandas as pd, numpy as np
from scipy.stats import spearmanr
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "integration/"))
import w2  # noqa: E402
import h3_inputs as H  # noqa: E402
PD = w2.PD
OUT = sys.argv[1]
HM = H.IN["hm_mm"]
RT = os.path.join(w2.B13, "../package", "results_table.csv")
OUTS = [f"{OUT}/b3_gain_vs_npos_per_tf.csv", f"{OUT}/b3_gain_vs_npos_spearman.csv", f"{OUT}/b3_gain_capped_vs_uncapped.csv"]
for f in OUTS:
    PD.backup(f)
res_t = pd.read_csv(RT, sep="\t")
rf = res_t[(res_t.Model == "RandomForestClassifier") & (res_t.PWM == "mono+di")].set_index("TF_name").sort_index()
mono = res_t[(res_t.Model == "Single best mono PWM") & (res_t.PWM == "mono")].drop_duplicates("TF_name").set_index("TF_name").sort_index()
PD.check_tfs(rf.index, "b3_gain_vs_npos RF rows")
for c in ("roc_auc_test_H", "pr_auc_test_H"):          # the baseline must be the train-selected best monoPWM
    if not np.allclose(rf[c + "_PWM"], mono.loc[rf.index, c + "_PWM_mono"], rtol=0, atol=0):
        raise SystemExit("CHECK FAILED b3_gain_vs_npos: %s_PWM is not the best monoPWM" % c)
hm = pd.read_csv(HM).set_index("TF").loc[rf.index]
ph = pd.read_csv(H.IN["pos_H"]).set_index("TF").loc[rf.index]
pm = pd.read_csv(H.IN["pos_M"]).set_index("TF").loc[rf.index]
ch = pd.read_csv(H.IN["counts_H"]).set_index("TF").loc[rf.index]
df = pd.DataFrame(index=rf.index)
df["n_train_pos_H"] = ph.train_file_pos
df["n_train_pos_M"] = pm.control_file_pos
if not ((df.n_train_pos_H == hm.n_train_pos_H).all() and (df.n_train_pos_M == hm.n_train_pos_M).all()):
    raise SystemExit("CHECK FAILED training positives: counts tables vs hm_mm_results")
# dataset counts from the used peak-set lists (human 361, mouse 276)
hid, mid = PD.human_ids(), PD.mouse_ids()
df["n_datasets_H"] = hid.groupby("TF").size().reindex(df.index)
df["n_datasets_M"] = mid.groupby("TF").size().reindex(df.index)
if df.n_datasets_H.isna().any() or df.n_datasets_M.isna().any() or df.n_datasets_H.sum() != len(hid) or df.n_datasets_M.sum() != len(mid):
    raise SystemExit("CHECK FAILED dataset counts: human %s of %d, mouse %s of %d" % (df.n_datasets_H.sum(), len(hid), df.n_datasets_M.sum(), len(mid)))
print("dataset counts: human %d, mouse %d" % (df.n_datasets_H.sum(), df.n_datasets_M.sum()))
if not ((hm.P_mono_H == ch.n_HUMAN_pwm_mono).all() and (hm.P_di_H == ch.n_HUMAN_pwm_di).all()):
    raise SystemExit("CHECK FAILED PWM counts: hm_mm vs B14 counts")
df["n_pwms_H"] = hm.P_mono_H + hm.P_di_H
df["cap_H"] = df.n_train_pos_H == 10000
df["cap_M"] = df.n_train_pos_M == 10000
df["dauroc_H"] = rf.roc_auc_test_H - rf.roc_auc_test_H_PWM
df["dauprc_H"] = rf.pr_auc_test_H - rf.pr_auc_test_H_PWM
# mouse test set: model trained on human; training positives of that model = human
df["dauroc_M"] = hm.HM_rf_test_M_auroc_mean - hm.HM_broc_test_M_auroc   # mean of two seeds
df["dauprc_M"] = hm.HM_rf_test_M_auprc_mean - hm.HM_bprc_test_M_auprc
df["dauroc_H_retrained"] = hm.HM_rf_test_H_auroc_s0 - hm.HM_broc_test_H_auroc
df["dauprc_H_retrained"] = hm.HM_rf_test_H_auprc_s0 - hm.HM_bprc_test_H_auprc
df.index.name = "TF"
PD.check_tfs(df.index, "b3_gain_vs_npos per TF")
if df[["dauroc_H", "dauprc_H"]].isna().any().any():
    raise SystemExit("CHECK FAILED human gains missing")
df.to_csv(f"{OUT}/b3_gain_vs_npos_per_tf.csv")
rows = []
for test, g, x, lab in [("human test", "dauroc_H", "n_train_pos_H", "train positives (human)"), ("human test", "dauprc_H", "n_train_pos_H", "train positives (human)"),
                        ("human test", "dauroc_H", "n_datasets_H", "datasets (human)"), ("human test", "dauprc_H", "n_datasets_H", "datasets (human)"),
                        ("human test", "dauroc_H", "n_pwms_H", "PWMs (human)"), ("human test", "dauprc_H", "n_pwms_H", "PWMs (human)"),
                        ("mouse test", "dauroc_M", "n_train_pos_H", "train positives (human; model trained on human)"), ("mouse test", "dauprc_M", "n_train_pos_H", "train positives (human; model trained on human)"),
                        ("mouse test", "dauroc_M", "n_datasets_H", "datasets (human)"), ("mouse test", "dauprc_M", "n_datasets_H", "datasets (human)"),
                        ("mouse test", "dauroc_M", "n_datasets_M", "datasets (mouse)"), ("mouse test", "dauprc_M", "n_datasets_M", "datasets (mouse)")]:
    r, p = spearmanr(df[x], df[g])
    sub = df[~df["cap_H"]] if x == "n_train_pos_H" else df
    r2, p2 = spearmanr(sub[x], sub[g]) if len(sub) > 3 else (np.nan, np.nan)
    rows.append(dict(test_set=test, gain=g, predictor=lab, predictor_col=x, n_tf=len(df), rho=r, p=p, n_tf_uncapped=len(sub), rho_uncapped=r2, p_uncapped=p2))
res = pd.DataFrame(rows)
res.to_csv(f"{OUT}/b3_gain_vs_npos_spearman.csv", index=False)
# capped vs uncapped medians
cs = []
for g in ["dauroc_H", "dauprc_H", "dauroc_M", "dauprc_M"]:
    c = df.cap_H
    cs.append(dict(gain=g, n_capped=int(c.sum()), median_capped=df.loc[c, g].median(), n_uncapped=int((~c).sum()), median_uncapped=df.loc[~c, g].median()))
pd.DataFrame(cs).to_csv(f"{OUT}/b3_gain_capped_vs_uncapped.csv", index=False)
w2.provenance(os.path.abspath(OUT), "3_analysis/top_peaks/b3_gain_vs_npos_h3.py", [os.path.abspath(f) for f in OUTS],
              [RT, HM, H.IN["pos_H"], H.IN["pos_M"], H.IN["counts_H"], PD.HUMAN_IDS, PD.MOUSE_IDS],
              params="{gain: RF mono+di minus best monoPWM (train-selected); mouse test set: mean of RF seeds 0 and 1; Spearman; uncapped = human training positives < 10,000}")
print(res.round(3).to_string()); print(pd.DataFrame(cs).round(4).to_string())
print(df[["dauroc_H","dauprc_H","dauroc_M","dauprc_M"]].median())
