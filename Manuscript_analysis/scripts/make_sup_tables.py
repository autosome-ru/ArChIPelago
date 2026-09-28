#!/usr/bin/env python3
"""Supplementary Tables 3-6 of the manuscript.

Basis: human values of every ArChIPelago model = the notebook 4 run; mouse test set = mouse chromosomes 1, 8
and 19 (mouse training set = mouse chromosomes 2-7, 9, 10 and 13-18); baseline everywhere = the best single monoPWM
selected on the training set (human training set for the human-trained models, mouse training set for the
mouse-trained control). No other reference PWM is reported.

Inputs (under Manuscript_analysis/analysis/, see archi_paths.py):
  ../HUMAN_MOUSE_total_100k.csv                             results table (make_results_table.py)
  operational/operational_metrics.csv                       tie-aware operational metrics
  hm_mm/hm_mm_results.csv                                   cross-species + mouse-trained control
  crossspecies/crossspecies_table.csv                       TF metadata + motif similarity (crossspecies_analysis.py)
  saturation/saturation_summary_by_k.csv, saturation_per_tf.csv, saturation_results.csv   (analyze_saturation.py)

Sup. Table 3  -- the performance table (+ README, summary)
Sup. Table 4  -- operational metrics (+ README, summary)
Sup. Table 5  -- cross-species diagnostics + mouse-trained control (+ README, summary, correlations)
Sup. Table 6  -- runtime / memory benchmark (values measured with the raw timings in analysis/runtime/)
headline_numbers.json -- every headline number quoted in the manuscript (written to Sup_Tables/)
"""
import json
import os
import numpy as np
import pandas as pd
from scipy import stats

from archi_paths import BASIS, RES

OUT = os.path.join(RES, "Sup_Tables")
os.makedirs(OUT, exist_ok=True)
MOUSE_TEST = "chr1, 8, 19"
HUMAN_TEST = "chr1, 8, 21"

T = pd.read_csv(os.path.join(RES, "HUMAN_MOUSE_total_100k.csv"), sep="\t")
assert T.shape == (900, 47)
TFS = sorted(T.TF_name.unique())
assert len(TFS) == 36

rf = T[(T.Model == "RandomForestClassifier") & (T.PWM == "mono+di")].set_index("TF_name").sort_index()
mono = T[(T.Model == "Single best mono PWM") & (T.PWM == "mono")].drop_duplicates("TF_name").set_index("TF_name").sort_index()
# the baseline: best single monoPWM (*_PWM columns == *_PWM_mono)
BASE = rf[["roc_auc_test_H_PWM", "roc_auc_test_M_PWM", "pr_auc_test_H_PWM", "pr_auc_test_M_PWM"]].copy()
for sp in ["H", "M"]:
    assert np.allclose(BASE[f"roc_auc_test_{sp}_PWM"], mono[f"roc_auc_test_{sp}_PWM_mono"]) and np.allclose(BASE[f"pr_auc_test_{sp}_PWM"], mono[f"pr_auc_test_{sp}_PWM_mono"])
headline = {"basis": f"human test set {HUMAN_TEST}; mouse test set {MOUSE_TEST}; baseline = best single monoPWM selected on the training set"}


def wp(x):
    return float(stats.wilcoxon(x).pvalue)


def desc(d):
    """TF names in descending order of the value."""
    return list(d.sort_values(ascending=False).index)


# ----------------------------------------------------------------------------------------------
# Sup. Table 3
# ----------------------------------------------------------------------------------------------
t3 = T.drop(columns=["Unnamed: 0", "Names"] + [c for c in T.columns if c.startswith(("mean_", "median_", "std_"))])
t3 = t3.rename(columns={"Count": "Number of PWMs", "PWM": "PWM type"})
readme3 = pd.DataFrame({"Sup. Table 3 -- ArChIPelago performance": [
    "One row per TF x model x PWM set (900 rows, 36 TFs). Models: RandomForestClassifier, XGBClassifier, LogisticRegression, "
    "BaggingClassifier_XGBClassifier, BaggingClassifier_LogisticRegression; the 'Single best mono PWM' / 'Single best di PWM' rows give the "
    "single-PWM values (repeated once per model block). All ArChIPelago models of this table are trained on human data only.",
    f"roc_auc_* / pr_auc_* = auROC / auPRC (PRROC integral) on the human training set (train_H), the human test set (test_H, {HUMAN_TEST}) "
    f"and the mouse test set (test_M, mouse {MOUSE_TEST}; the mouse training set of the mouse-trained control of Sup. Table 5 is the remaining "
    "mouse autosomes).",
    "*_PWM_mono / *_PWM_di = the best single monoPWM / diPWM of the TF, selected on the human training set (independently by auROC and by "
    "auPRC) and evaluated by PWM identity on each test set; the same value is repeated in every row of that TF. *_PWM = the baseline the "
    "figures and the text compare against = the best single monoPWM, the common reference for models built on monoPWMs, on diPWMs and on "
    "both, as stated in the figure legends. Sheet 'summary': the Random Forest on monoPWMs + diPWMs against this baseline on both test sets.",
    "The human values of the ArChIPelago models are those of the published run. The mouse values are computed with the same models and "
    "hyperparameters (random_state = 0) trained on the human training set and evaluated on the aligned mouse feature matrix of the mouse "
    "test set; the mouse baselines are the same human-train-selected PWMs evaluated on the mouse test set.",
    "'Number of PWMs' = number of monoPWMs of the TF; Seq_count = number of human training positives.",
]})
sum3 = []
for sp, lab in (("H", f"human test set ({HUMAN_TEST})"), ("M", f"mouse test set ({MOUSE_TEST})")):
    b_r, b_p = mono[f"roc_auc_test_{sp}_PWM_mono"], mono[f"pr_auc_test_{sp}_PWM_mono"]
    d_r, d_p = rf[f"roc_auc_test_{sp}"] - b_r, rf[f"pr_auc_test_{sp}"] - b_p
    sum3.append([lab, "best single monoPWM (selected on the human training set)", round(rf[f"roc_auc_test_{sp}"].median(), 4), round(b_r.median(), 4),
                 round(d_r.median(), 4), int((d_r > 0).sum()), f"{wp(d_r):.1e}",
                 round(rf[f"pr_auc_test_{sp}"].median(), 4), round(b_p.median(), 4), round(d_p.median(), 4), int((d_p > 0).sum()), f"{wp(d_p):.1e}",
                 ", ".join(sorted(set(d_r.index[d_r < 0]) | set(d_p.index[d_p < 0]))) or "none"])
sum3 = pd.DataFrame(sum3, columns=["Test set", "Baseline", "ArChIPelago auROC (median)", "baseline auROC (median)", "median per-TF dauROC", "TFs with dauROC > 0",
                                   "Wilcoxon P (auROC)", "ArChIPelago auPRC (median)", "baseline auPRC (median)", "median per-TF dauPRC", "TFs with dauPRC > 0",
                                   "Wilcoxon P (auPRC)", "TFs below the baseline on >= 1 metric"])
with pd.ExcelWriter(os.path.join(OUT, "Sup_Table_3_ArChIPelago_performance.xlsx")) as xw:
    sum3.to_excel(xw, sheet_name="summary", index=False)
    t3.to_excel(xw, sheet_name="performance", index=False)
sum3.to_csv(os.path.join(OUT, "Sup_Table_3_summary.csv"), index=False)
t3.to_csv(os.path.join(OUT, "Sup_Table_3_ArChIPelago_performance.csv"), index=False)

# headline numbers, baseline = best single monoPWM
for sp in ["H", "M"]:
    b_roc, b_prc = BASE[f"roc_auc_test_{sp}_PWM"], BASE[f"pr_auc_test_{sp}_PWM"]
    dr = rf[f"roc_auc_test_{sp}"] - b_roc
    dp = rf[f"pr_auc_test_{sp}"] - b_prc
    headline[f"RF_{sp}"] = dict(
        rf_auroc_median=round(rf[f"roc_auc_test_{sp}"].median(), 4), rf_auprc_median=round(rf[f"pr_auc_test_{sp}"].median(), 4),
        base_auroc_median=round(b_roc.median(), 4), base_auprc_median=round(b_prc.median(), 4),
        diff_of_medians_auroc=round(rf[f"roc_auc_test_{sp}"].median() - b_roc.median(), 4),
        diff_of_medians_auprc=round(rf[f"pr_auc_test_{sp}"].median() - b_prc.median(), 4),
        median_per_tf_auroc=round(dr.median(), 4), median_per_tf_auprc=round(dp.median(), 4),
        n_improved_auroc=int((dr > 0).sum()), n_improved_auprc=int((dp > 0).sum()),
        below_baseline={tf: [round(dr[tf], 3), round(dp[tf], 3)] for tf in sorted(set(dr.index[dr < 0]) | set(dp.index[dp < 0]))},
        wilcoxon_p=[wp(dr), wp(dp)],
        gt005_both=sorted(dr.index[(dr > 0.05) & (dp > 0.05)]),
        gt005_auroc=desc(dr[dr > 0.05]),            # TFs with dauROC > 0.05, descending
        gt008_auprc=desc(dp[dp >= 0.08]),           # TFs with dauPRC >= 0.08, descending
        gt01_auprc=sorted(dp.index[dp > 0.1]), gt02_auprc=sorted(dp.index[dp > 0.2]),
        smallest_auprc={tf: round(dp[tf], 3) for tf in dp.nsmallest(3).index},
        per_tf={tf: [round(dr[tf], 3), round(dp[tf], 3)] for tf in ["ANDR", "E2F4", "IRF4", "SRF", "RXRA", "P53", "MAFK", "SPI1", "SOX2", "STAT1", "GATA3", "FLI1", "PRGR", "STAT3"]},
        spearman_gain_vs_baseline=[round(v, 4) for v in stats.spearmanr(dr, b_roc)],
        spearman_gain_vs_baseline_P=float(stats.spearmanr(dr, b_roc)[1]),
        spearman_auprc_gain_vs_baseline=[round(stats.spearmanr(dp, b_prc)[0], 4), float(stats.spearmanr(dp, b_prc)[1])],
        baseline_auroc_ge_097=sorted(b_roc.index[b_roc >= 0.97]),
    )
for pwm in ["mono", "di", "mono+di"]:
    r = T[(T.Model == "RandomForestClassifier") & (T.PWM == pwm)].set_index("TF_name").sort_index()
    headline[f"RF_{pwm}_H"] = dict(auroc=round(r.roc_auc_test_H.median(), 4), auprc=round(r.pr_auc_test_H.median(), 4),
                                   n_impr_auroc=int(((r.roc_auc_test_H - BASE.roc_auc_test_H_PWM) > 0).sum()),
                                   n_impr_auprc=int(((r.pr_auc_test_H - BASE.pr_auc_test_H_PWM) > 0).sum()))
allm = T[(T.PWM == "mono+di") & (~T.Model.str.startswith("Single"))]
headline["models_monodi_medians_H"] = allm.groupby("Model")[["roc_auc_test_H", "pr_auc_test_H"]].median().round(4).to_dict()

# ----------------------------------------------------------------------------------------------
# Sup. Table 4 -- operational metrics (baseline: best single monoPWM, selected by training auPRC)
# ----------------------------------------------------------------------------------------------
op = pd.read_csv(os.path.join(BASIS, "operational", "operational_metrics.csv"))
assert op.shape[0] == 72 and sorted(op.tf.unique()) == TFS
assert not any((c.endswith(("_di", "_pwm")) and c != "baseline_pwm") or c.startswith(("baseline_class", "train_auprc_best")) for c in op.columns), \
    "operational_metrics.csv must carry the best-monoPWM baseline only"
ren = {"tf": "TF", "test": "Test set", "n": "Candidate regions in test set", "n_pos": "True binding sites (positives)",
       "frac_tied_scores_rf": "Fraction of regions sharing their score with another region: ArChIPelago RF",
       "frac_tied_scores_mono": "Fraction of regions sharing their score with another region: best monoPWM"}
cols4 = ["TF", "Test set", "Candidate regions in test set", "True binding sites (positives)", ren["frac_tied_scores_rf"], ren["frac_tied_scores_mono"]]
BASES = [("mono", "best monoPWM"), ("rf", "ArChIPelago RF")]
for rec, lab in [("0.5", "50 %"), ("0.8", "80 %")]:
    for b, bl in BASES:
        ren[f"fp_at_recall{rec}_per10k_{b}"] = f"FP per 10,000 regions at {lab} recall: {bl}"
        cols4.append(ren[f"fp_at_recall{rec}_per10k_{b}"])
    for b, bl in BASES:
        ren[f"precision_at_recall{rec}_{b}"] = f"Precision at {lab} recall: {bl}"
        cols4.append(ren[f"precision_at_recall{rec}_{b}"])
    if rec == "0.5":
        ren["precision_gain_recall0.5_vs_mono"] = "Precision gain at 50 % recall, ArChIPelago minus best monoPWM"
        cols4.append(ren["precision_gain_recall0.5_vs_mono"])
    ren[f"fp_reduction_recall{rec}_pct_vs_mono"] = f"FP reduction at {lab} recall, ArChIPelago vs best monoPWM (%)"
    cols4.append(ren[f"fp_reduction_recall{rec}_pct_vs_mono"])
for b, bl in BASES:
    ren[f"tp_top1000_{b}"] = f"True sites among top 1,000 predictions: {bl}"
    cols4.append(ren[f"tp_top1000_{b}"])
ren["extra_tp_top1000_vs_mono"] = "Additional true sites in top 1,000, ArChIPelago vs best monoPWM"
cols4.append(ren["extra_tp_top1000_vs_mono"])
for b, bl in BASES:
    ren[f"tp_top1pct_{b}"] = f"True sites among top 1 % of predictions: {bl}"
    cols4.append(ren[f"tp_top1pct_{b}"])
t4 = op.rename(columns=ren)[cols4].copy()
t4["Test set"] = t4["Test set"].map({"H": f"human ({HUMAN_TEST})", "M": f"mouse ({MOUSE_TEST})"})
t4 = t4.sort_values(["Test set", "TF"]).reset_index(drop=True)
t4[t4.select_dtypes("number").columns] = t4.select_dtypes("number").round(3)
rows = []
for sp, lab in [("H", f"human test set ({HUMAN_TEST})"), ("M", f"mouse test set ({MOUSE_TEST})")]:
    o = op[op.test == sp]
    bl = "best monoPWM (selected on the human training set by auPRC)"
    for rec in ["0.5", "0.8"]:
        red = o[f"fp_reduction_recall{rec}_pct_vs_mono"]
        rows.append([lab, bl, f"FP reduction at {float(rec) * 100:.0f} % recall (%)", round(red.median(), 1),
                     f"{red.quantile(.25):.1f} to {red.quantile(.75):.1f}", int((red > 0).sum()),
                     round(o[f"fp_at_recall{rec}_per10k_mono"].median()), round(o[f"fp_at_recall{rec}_per10k_rf"].median()),
                     round(o[f"precision_at_recall{rec}_mono"].median(), 3), round(o[f"precision_at_recall{rec}_rf"].median(), 3)])
    ex = o["extra_tp_top1000_vs_mono"]
    rows.append([lab, bl, "additional true sites in top 1,000", round(ex.median(), 1),
                 f"{ex.quantile(.25):.1f} to {ex.quantile(.75):.1f}", int((ex > 0).sum()),
                 round(o["tp_top1000_mono"].median()), round(o["tp_top1000_rf"].median()), np.nan, np.nan])
    pg = o["precision_gain_recall0.5_vs_mono"]
    rows.append([lab, bl, "precision gain at 50 % recall (paired, per TF)", round(pg.median(), 3),
                 f"{pg.quantile(.25):.3f} to {pg.quantile(.75):.3f}", int((pg > 0).sum()),
                 round(o["precision_at_recall0.5_mono"].median(), 3), round(o["precision_at_recall0.5_rf"].median(), 3), np.nan, np.nan])
sum4 = pd.DataFrame(rows, columns=["Test set", "Baseline", "Quantity", "Median over 36 TFs", "IQR", "TFs with improvement (of 36)",
                                   "Median baseline value", "Median ArChIPelago value", "Median baseline precision", "Median ArChIPelago precision"])
oh, om = op[op.test == "H"], op[op.test == "M"]
b = "mono"
headline["operational"] = {b: dict(
    H_fp50=[round(oh[f"fp_reduction_recall0.5_pct_vs_{b}"].median(), 1), round(oh[f"fp_reduction_recall0.5_pct_vs_{b}"].quantile(.25), 1),
            round(oh[f"fp_reduction_recall0.5_pct_vs_{b}"].quantile(.75), 1), int((oh[f"fp_reduction_recall0.5_pct_vs_{b}"] > 0).sum())],
    H_fp80=[round(oh[f"fp_reduction_recall0.8_pct_vs_{b}"].median(), 1), int((oh[f"fp_reduction_recall0.8_pct_vs_{b}"] > 0).sum())],
    H_precision50=[round(oh[f"precision_at_recall0.5_{b}"].median(), 3), round(oh["precision_at_recall0.5_rf"].median(), 3)],
    H_precision50_gain=[round(oh[f"precision_gain_recall0.5_vs_{b}"].median(), 3), round(oh[f"precision_gain_recall0.5_vs_{b}"].quantile(.25), 3),
                        round(oh[f"precision_gain_recall0.5_vs_{b}"].quantile(.75), 3), int((oh[f"precision_gain_recall0.5_vs_{b}"] > 0).sum())],
    M_precision50=[round(om[f"precision_at_recall0.5_{b}"].median(), 3), round(om["precision_at_recall0.5_rf"].median(), 3)],
    M_precision50_gain=[round(om[f"precision_gain_recall0.5_vs_{b}"].median(), 3), int((om[f"precision_gain_recall0.5_vs_{b}"] > 0).sum())],
    M_extra_tp1000=[round(om[f"extra_tp_top1000_vs_{b}"].median(), 1), int((om[f"extra_tp_top1000_vs_{b}"] > 0).sum())],
    H_extra_tp1000=[round(oh[f"extra_tp_top1000_vs_{b}"].median(), 1), int((oh[f"extra_tp_top1000_vs_{b}"] > 0).sum())],
    H_fp50_per10k=[round(oh[f"fp_at_recall0.5_per10k_{b}"].median()), round(oh["fp_at_recall0.5_per10k_rf"].median())],
    M_fp50=[round(om[f"fp_reduction_recall0.5_pct_vs_{b}"].median(), 1), int((om[f"fp_reduction_recall0.5_pct_vs_{b}"] > 0).sum())],
    M_fp80=[round(om[f"fp_reduction_recall0.8_pct_vs_{b}"].median(), 1), int((om[f"fp_reduction_recall0.8_pct_vs_{b}"] > 0).sum())],
)}
readme4 = pd.DataFrame({"Sup. Table 4 -- operational relevance of the ArChIPelago gain": [
    "For each TF and test set, the predictions of the human-trained ArChIPelago Random Forest (monoPWMs + diPWMs) and of the best single "
    "monoPWM are ranked by score over all candidate regions of the held-out test chromosomes (positives = ChIP-Seq peak summits of the TF; "
    f"negatives = GC-matched peaks of unrelated TFs, about 1:100). Test sets: human {HUMAN_TEST}; mouse {MOUSE_TEST}.",
    "Ties: single-PWM scores are heavily tied (a large share of the candidate regions share their best-hit score with other regions; the fraction is "
    "given per TF). Every block of equal scores is therefore treated as a unit and counts are the expectation under a random order inside "
    "the block (the linear interpolation of the precision-recall curve through the tie block), identically for the Random Forest and for the "
    "single PWM.",
    "Reported: false positives (per 10,000 candidate regions) that must be accepted to recover 50 % and 80 % of the true sites; the precision at "
    "these recall levels; true sites among the top 1,000 and the top 1 % of the ranked regions; the relative reduction of false positives and "
    "the additional true sites of ArChIPelago over the baseline; and the paired per-TF precision gain at 50 % recall.",
    "Baseline: the best single monoPWM, selected on the human training set by auPRC (PRROC integral).",
    f"Top-1,000 counts on the mouse test set are close to saturation ({int(om.n.min()):,}-{int(om.n.max()):,} regions with up to "
    f"{int(om.n_pos.max()):,} positives), so the fixed-recall statistics are the informative ones there. Sheet 'summary' gives medians and "
    "interquartile ranges over the 36 TFs.",
]})
with pd.ExcelWriter(os.path.join(OUT, "Sup_Table_4_operational_metrics.xlsx")) as xw:
    sum4.to_excel(xw, sheet_name="summary", index=False)
    t4.to_excel(xw, sheet_name="per_TF", index=False)
t4.to_csv(os.path.join(OUT, "Sup_Table_4_operational_metrics.csv"), index=False)
sum4.to_csv(os.path.join(OUT, "Sup_Table_4_operational_metrics_summary.csv"), index=False)

# ----------------------------------------------------------------------------------------------
# Sup. Table 5 -- cross-species + mouse-trained control (baseline: best single monoPWM)
# ----------------------------------------------------------------------------------------------
hm = pd.read_csv(os.path.join(BASIS, "hm_mm", "hm_mm_results.csv")).set_index("TF").sort_index()
cs = pd.read_csv(os.path.join(BASIS, "crossspecies", "crossspecies_table.csv")).set_index("TF").sort_index()
assert list(hm.index) == TFS == list(cs.index) == list(rf.index) and hm.MM_available.all()
assert np.allclose(rf.roc_auc_test_M, hm.HM_rf_test_M_auroc_s0) and np.allclose(rf.pr_auc_test_M, hm.HM_rf_test_M_auprc_s0)
assert np.allclose(mono.roc_auc_test_M_PWM_mono, hm.HM_broc_test_M_auroc) and np.allclose(mono.pr_auc_test_M_PWM_mono, hm.HM_bprc_test_M_auprc)

t5 = pd.DataFrame(index=TFS)
t5.index.name = "TF"
t5["TF family (TFClass)"] = cs.family
t5["TFClass id"] = cs.TFclass_id
t5["Human monoPWMs"] = hm.P_mono_H
t5["Human diPWMs"] = hm.P_di_H
t5["Mouse monoPWMs"] = hm.P_mono_M
t5["Mouse diPWMs"] = hm.P_di_M
t5["Human ChIP-Seq experiments"] = cs.n_exp_H
t5["Mouse ChIP-Seq experiments"] = cs.n_exp_M
t5["Human training positives"] = hm.n_train_pos_H
t5["Mouse training positives"] = hm.n_train_pos_M
t5["Mouse test positives"] = hm.n_test_pos_M
t5["Human cell types (experiments)"] = cs.cells_H
t5["Mouse cell types (experiments)"] = cs.cells_M
t5["Cell-type overlap human/mouse"] = cs.celltype_overlap
t5["H>M: ArChIPelago auROC (mouse test)"] = rf.roc_auc_test_M
t5["H>M: best monoPWM auROC"] = mono.roc_auc_test_M_PWM_mono
t5["H>M: dauROC"] = rf.roc_auc_test_M - mono.roc_auc_test_M_PWM_mono
t5["H>M: ArChIPelago auPRC (mouse test)"] = rf.pr_auc_test_M
t5["H>M: best monoPWM auPRC"] = mono.pr_auc_test_M_PWM_mono
t5["H>M: dauPRC"] = rf.pr_auc_test_M - mono.pr_auc_test_M_PWM_mono
t5["H>M: below baseline on >=1 metric"] = (t5["H>M: dauROC"] < 0) | (t5["H>M: dauPRC"] < 0)
t5["H>H: dauROC (reference)"] = rf.roc_auc_test_H - mono.roc_auc_test_H_PWM_mono
t5["H>H: dauPRC (reference)"] = rf.pr_auc_test_H - mono.pr_auc_test_H_PWM_mono
t5["M>M: ArChIPelago auROC (mouse test, mean of 2 seeds)"] = hm.MM_rf_test_M_auroc_mean
t5["M>M: best mouse monoPWM auROC"] = hm.MM_broc_test_M_auroc
t5["M>M: dauROC"] = hm.MM_rf_test_M_auroc_mean - hm.MM_broc_test_M_auroc
t5["M>M: ArChIPelago auPRC (mouse test, mean of 2 seeds)"] = hm.MM_rf_test_M_auprc_mean
t5["M>M: best mouse monoPWM auPRC"] = hm.MM_bprc_test_M_auprc
t5["M>M: dauPRC"] = hm.MM_rf_test_M_auprc_mean - hm.MM_bprc_test_M_auprc
t5["M>M: below baseline on >=1 metric"] = (t5["M>M: dauROC"] < 0) | (t5["M>M: dauPRC"] < 0)
t5["M>M: auROC seed 0"] = hm.MM_rf_test_M_auroc_s0
t5["M>M: auROC seed 1"] = hm.MM_rf_test_M_auroc_s1
t5["M>M: auPRC seed 0"] = hm.MM_rf_test_M_auprc_s0
t5["M>M: auPRC seed 1"] = hm.MM_rf_test_M_auprc_s1
t5["M>M minus H>M: dauROC (same mouse test set, mean of 2 seeds for both)"] = hm.MM_rf_test_M_auroc_mean - hm.HM_rf_test_M_auroc_mean
t5["M>M minus H>M: dauPRC (same mouse test set, mean of 2 seeds for both)"] = hm.MM_rf_test_M_auprc_mean - hm.HM_rf_test_M_auprc_mean
t5["Baseline human monoPWM (best training auROC)"] = cs.topH_PWM
t5["Similarity of the baseline human monoPWM to the nearest mouse monoPWM (Pearson r of aligned columns)"] = cs.topH_vs_nearestM_pcc
t5["Median similarity human monoPWMs -> nearest mouse monoPWM"] = cs.median_H_vs_nearestM_pcc
unrounded5 = t5.select_dtypes("number").astype(float).copy()
t5[t5.select_dtypes("number").columns] = t5.select_dtypes("number").astype(float).round(4)
for c in ["Human monoPWMs", "Human diPWMs", "Mouse monoPWMs", "Mouse diPWMs", "Human ChIP-Seq experiments", "Mouse ChIP-Seq experiments",
          "Human training positives", "Mouse training positives", "Mouse test positives"]:
    t5[c] = t5[c].astype(int)

dHM_roc, dHM_prc = unrounded5["H>M: dauROC"], unrounded5["H>M: dauPRC"]
dMM_roc, dMM_prc = unrounded5["M>M: dauROC"], unrounded5["M>M: dauPRC"]
dD_roc = unrounded5["M>M minus H>M: dauROC (same mouse test set, mean of 2 seeds for both)"]
dD_prc = unrounded5["M>M minus H>M: dauPRC (same mouse test set, mean of 2 seeds for both)"]
seed_spread_roc = (hm.MM_rf_test_M_auroc_s0 - hm.MM_rf_test_M_auroc_s1).abs().max()
seed_spread_prc = (hm.MM_rf_test_M_auprc_s0 - hm.MM_rf_test_M_auprc_s1).abs().max()
summ = []


def add5(setting, dr, dp, extra=""):
    summ.append([setting, round(dr.median(), 4), int((dr > 0).sum()), f"{wp(dr):.1e}", round(dp.median(), 4), int((dp > 0).sum()), f"{wp(dp):.1e}",
                 ", ".join(sorted(set(dr.index[dr < 0]) | set(dp.index[dp < 0]))) or "none"])


add5("Human-trained model on the mouse test set vs the best monoPWM (selected on human training data)", dHM_roc, dHM_prc)
add5("Mouse-trained model on the mouse test set vs the best mouse monoPWM (selected on mouse training data)", dMM_roc, dMM_prc)
add5("Mouse-trained minus human-trained model on the same mouse test set (positive = the mouse-trained model is better)", dD_roc, dD_prc)
add5("Human-trained model on the human test set vs the best monoPWM [reference]", unrounded5["H>H: dauROC (reference)"], unrounded5["H>H: dauPRC (reference)"])
sum5 = pd.DataFrame(summ, columns=["Comparison (36 TFs)", "median dauROC", "TFs with dauROC > 0", "Wilcoxon P (auROC)", "median dauPRC",
                                   "TFs with dauPRC > 0", "Wilcoxon P (auPRC)", "TFs below baseline on >= 1 metric"])
sim = unrounded5["Similarity of the baseline human monoPWM to the nearest mouse monoPWM (Pearson r of aligned columns)"]
corr = []
CORR_VARS = [("similarity of the baseline human monoPWM to its most similar mouse monoPWM (per_TF column 'Similarity of the baseline human monoPWM ...')", sim),
             ("median, over all human monoPWMs of the TF, of the similarity to the most similar mouse monoPWM (per_TF column 'Median similarity ...')",
              unrounded5["Median similarity human monoPWMs -> nearest mouse monoPWM"]),
             ("number of human training positives (per_TF column 'Human training positives')", t5["Human training positives"]),
             ("number of mouse training positives (per_TF column 'Mouse training positives')", t5["Mouse training positives"]),
             ("number of human ChIP-Seq experiments", t5["Human ChIP-Seq experiments"]), ("number of mouse ChIP-Seq experiments", t5["Mouse ChIP-Seq experiments"]),
             ("number of human PWMs (monoPWMs + diPWMs)", t5["Human monoPWMs"] + t5["Human diPWMs"]),
             ("auROC of the best human monoPWM on the mouse test set (quality of the baseline; per_TF column 'H>M: best monoPWM auROC')", unrounded5["H>M: best monoPWM auROC"])]
for lab, x in CORR_VARS:
    for dl, dd in [("H>M: dauROC", dHM_roc), ("H>M: dauPRC", dHM_prc)]:
        r, p = stats.spearmanr(dd, x)
        corr.append([dl, lab, 36, round(r, 3), f"{p:.3g}"])
CORR_HEADER = [
    "Spearman rank correlation, across the 36 TFs, between the cross-species gain of a TF and one TF-level property (one row per pair; two-sided P; no multiple-testing correction).",
    "Cross-species gain = 'H>M: dauROC' / 'H>M: dauPRC' of the per_TF sheet: the human-trained ArChIPelago Random Forest (monoPWMs + diPWMs) minus the best single human monoPWM "
    "(selected on the human training set), both evaluated on the mouse test set (mouse chromosomes 1, 8 and 19).",
    "Motif similarity = maximum over all alignments (both strands, at least 5 overlapping columns) of the mean per-column Pearson correlation between the log-odds matrices of a "
    "human monoPWM and a mouse monoPWM of the same TF; for each human monoPWM the most similar mouse monoPWM of the TF is taken. The baseline human monoPWM is the one with the "
    "best training auROC (the *_PWM baseline of Sup. Table 3; per_TF column 'Baseline human monoPWM').",
    "The first pair is the relation shown in Fig. S6C.",
]
corr5 = pd.DataFrame(corr, columns=["Cross-species gain (per TF)", "TF property", "n (TFs)", "Spearman rho", "P (two-sided)"])
fail = sorted(t5.index[t5["H>M: below baseline on >=1 metric"]])
headline["crossspecies"] = dict(
    HM=dict(rf_auroc=round(rf.roc_auc_test_M.median(), 4), base_auroc=round(mono.roc_auc_test_M_PWM_mono.median(), 4),
            rf_auprc=round(rf.pr_auc_test_M.median(), 4), base_auprc=round(mono.pr_auc_test_M_PWM_mono.median(), 4),
            d_auroc=round(dHM_roc.median(), 4), d_auprc=round(dHM_prc.median(), 4),
            n_improved=[int((dHM_roc > 0).sum()), int((dHM_prc > 0).sum())], p=[wp(dHM_roc), wp(dHM_prc)],
            fail={tf: [round(dHM_roc[tf], 3), round(dHM_prc[tf], 3)] for tf in fail}),
    MM=dict(rf_auroc=round(hm.MM_rf_test_M_auroc_mean.median(), 4), base_auroc=round(hm.MM_broc_test_M_auroc.median(), 4),
            rf_auprc=round(hm.MM_rf_test_M_auprc_mean.median(), 4), base_auprc=round(hm.MM_bprc_test_M_auprc.median(), 4),
            d_auroc=round(dMM_roc.median(), 4), d_auprc=round(dMM_prc.median(), 4),
            n_improved=[int((dMM_roc > 0).sum()), int((dMM_prc > 0).sum())], p=[wp(dMM_roc), wp(dMM_prc)],
            fail={tf: [round(dMM_roc[tf], 3), round(dMM_prc[tf], 3)] for tf in sorted(set(dMM_roc.index[dMM_roc < 0]) | set(dMM_prc.index[dMM_prc < 0]))}),
    MM_minus_HM=dict(d_auroc=round(dD_roc.median(), 4), d_auprc=round(dD_prc.median(), 4),
                     n=[int((dD_roc > 0).sum()), int((dD_prc > 0).sum())]),
    MM_seed_spread_max=[round(float(seed_spread_roc), 4), round(float(seed_spread_prc), 4)],
    fail_recovered_MM={tf: [round(dMM_roc[tf], 3), round(dMM_prc[tf], 3)] for tf in fail},
    fail_train_pos={tf: [int(t5.loc[tf, "Human training positives"]), int(t5.loc[tf, "Mouse training positives"])] for tf in fail},
    fail_similarity={tf: round(sim[tf], 3) for tf in fail}, similarity_range=[round(sim.min(), 2), round(sim.max(), 2)],
    similarity_median=round(sim.median(), 3),
    similarity_mwu_fail_vs_rest=[round(v, 3) for v in stats.mannwhitneyu(sim[fail], sim.drop(fail))] if fail else None,
    spearman_sim=[[round(v, 3) for v in stats.spearmanr(dHM_roc, sim)], [round(v, 3) for v in stats.spearmanr(dHM_prc, sim)]],
    top_gain_auroc={tf: round(dHM_roc[tf], 3) for tf in dHM_roc.nlargest(4).index},
    top_gain_auprc={tf: round(dHM_prc[tf], 3) for tf in dHM_prc.nlargest(4).index},
    gt005_auroc=desc(dHM_roc[dHM_roc > 0.05]), gt006_auprc=desc(dHM_prc[dHM_prc > 0.06]),
    mouse_train_pos_median=int(hm.n_train_pos_M.median()),
)
readme5 = pd.DataFrame({"Sup. Table 5 -- cross-species transfer (human -> mouse) and the mouse-trained control": [
    "H>M = the ArChIPelago Random Forest trained on human ChIP-Seq data with human monoPWMs + diPWMs (the model of the paper), evaluated on the "
    f"mouse test set (mouse {MOUSE_TEST}). Baseline: the best single monoPWM selected on the HUMAN training set (by PWM identity, "
    "independently by auROC and by auPRC) and evaluated on the mouse test set.",
    "M>M = the ArChIPelago Random Forest trained on MOUSE ChIP-Seq data (mouse training set = mouse chromosomes 2-7, 9, 10 and 13-18; mouse monoPWMs + "
    "diPWMs of the same TF) and evaluated on the same mouse test set; baseline = the best single mouse monoPWM selected on the mouse training set. "
    "'M>M minus H>M' compares the two models on identical test rows (positive = the mouse-trained model is better).",
    "H>M values are the mouse numbers of Sup. Table 3 (random_state = 0); M>M values are the mean of two fits (random_state 0 and 1, given "
    f"separately); the 'M>M minus H>M' columns use the two-seed mean for both models. Between the two seeds a per-TF auROC changes by at most "
    f"{seed_spread_roc:.3f} and an auPRC by at most {seed_spread_prc:.3f}, so differences of that size should be read as ties.",
    "Metadata: TF family and TFClass id, PWM counts, numbers of ChIP-Seq experiments and cell types (from Sup. Tables 1-2 and the GTRD metadata), "
    "training and test positives. Motif similarity = maximum over alignments (both strands, at least 5 overlapping columns) of the mean per-column "
    "Pearson correlation between the log-odds matrices of the baseline human monoPWM (the best by training auROC, i.e. the *_PWM baseline of Sup. Table 3) "
    "and its most similar mouse monoPWM of the same TF (maximum over all mouse monoPWMs of the TF).",
    "Sheets: 'summary' (medians, counts and Wilcoxon signed-rank tests over the 36 TFs), 'correlations' (Spearman correlations of the cross-species "
    "gain with TF properties), 'per_TF'.",
]})
with pd.ExcelWriter(os.path.join(OUT, "Sup_Table_5_cross_species_and_mouse_trained.xlsx")) as xw:
    sum5.to_excel(xw, sheet_name="summary", index=False)
    pd.DataFrame({"Spearman correlations of the cross-species gain with TF properties": CORR_HEADER}).to_excel(xw, sheet_name="correlations", index=False)
    corr5.to_excel(xw, sheet_name="correlations", index=False, startrow=len(CORR_HEADER) + 2)
    t5.reset_index().to_excel(xw, sheet_name="per_TF", index=False)
t5.reset_index().to_csv(os.path.join(OUT, "Sup_Table_5_cross_species_and_mouse_trained.csv"), index=False)
sum5.to_csv(os.path.join(OUT, "Sup_Table_5_summary.csv"), index=False)

# ----------------------------------------------------------------------------------------------
# Sup. Table 6 -- runtime / memory (values measured with the raw timings in analysis/runtime/)
# ----------------------------------------------------------------------------------------------
hw = pd.DataFrame({"Item": ["CPU", "RAM", "OS", "Java / SARUS", "Python", "Timing", "Benchmark input", "Training matrices", "Thread policy"],
                   "Value": ["Intel Xeon E5-4607 v2 @ 2.60 GHz (4 sockets x 6 cores, 2 threads per core; shared server, load 3-8 of 48 threads during the measurements)",
                             "503 GB total", "Ubuntu 20.04.6 LTS, kernel 5.4.0",
                             "OpenJDK 1.8.0_232; SPRY-SARUS 2.0.2, command line identical to the pipeline (--skipn --show-non-matching --output-scoring-mode score besthit), JVM pinned to one CPU (-XX:ActiveProcessorCount=1, -Xmx2G)",
                             "Python 3.8; scikit-learn 1.3.0, pandas 2.0.3, numpy 1.24.3, joblib 1.2.0",
                             "GNU time (/usr/bin/time -f '%e %U %S %M'): wall-clock s, user s, system s, maximum resident set size; in-process timers for load / scaler / fit / predict",
                             "50,000 sequences x 300 bp = 15.0 Mb (10,000 CTCF human test positives + 40,000 GC-matched negatives); PWM sets of STAT3 (5 mono / 5 di), HNF4A (26 / 13) and CTCF (198 / 99)",
                             "the pipeline's human training matrices: STAT3 229,482 rows (3,000 positives), HNF4A 201,495 (10,000), CTCF 231,291 (10,000); features = best-hit log-odds of the monoPWMs",
                             "SARUS with the JVM limited to one active processor (JIT/GC helper threads still run: CPU time = 1.4-1.8 x wall-clock); Random Forest with n_jobs = 1 "
                             "(genuinely single-threaded, CPU/wall = 0.96-0.99) and n_jobs = 4 (OMP_NUM_THREADS = 4); all jobs at nice 10. Repeats: single-PWM scans 3 "
                             "(medians reported), PWM loops and RF runs 1 (CTCF single-call inference 2). Feature matrices of the RF benchmarks contain the "
                             "monoPWMs only (5 / 26 / 198 features); the mono+di models of the paper have 10 / 39 / 297 features and scale about linearly."]})
scan = pd.DataFrame([
    ["JVM start-up (10 sequences, 3 kb)", "any", "mono", 1, 0.003, 0.18, 0.18, 0.06],
    ["1 PWM", "STAT3", "mono", 1, 15.0, 2.36, 2.50, 0.31], ["1 PWM", "HNF4A", "mono", 1, 15.0, 2.19, 2.56, 0.29], ["1 PWM", "CTCF", "mono", 1, 15.0, 2.25, 2.87, 0.30],
    ["1 PWM, pipeline JVM flags (no CPU cap)", "CTCF", "mono", 1, 15.0, 2.34, 3.40, 0.23],
    ["1 diPWM", "STAT3", "di", 1, 15.0, 2.54, 2.75, 0.53], ["1 diPWM", "HNF4A", "di", 1, 15.0, 2.94, 3.17, 0.63], ["1 diPWM", "CTCF", "di", 1, 15.0, 2.95, 3.25, 0.64],
    ["all monoPWMs of the TF (sequential loop)", "STAT3", "mono", 5, 15.0, 11.5, 14.1, 0.31], ["all monoPWMs of the TF (sequential loop)", "HNF4A", "mono", 26, 15.0, 58.2, 72.9, 0.31],
    ["all monoPWMs of the TF (sequential loop)", "CTCF", "mono", 198, 15.0, 478.4, 614.7, 0.35],
    ["all diPWMs of the TF (sequential loop)", "STAT3", "di", 5, 15.0, 14.0, 15.8, 0.63], ["all diPWMs of the TF (sequential loop)", "HNF4A", "di", 13, 15.0, 36.2, 40.6, 0.64],
    ["all diPWMs of the TF (sequential loop)", "CTCF", "di", 99, 15.0, 306.3, 354.1, 0.64],
], columns=["Task (SPRY-SARUS, JVM limited to one active processor)", "TF", "PWM type", "Number of PWMs", "Scanned sequence (Mb)", "Wall-clock (s)", "User CPU (s)", "Peak memory (GB)"])
scan_notes = pd.DataFrame({"Derived quantities": [
    "Per-PWM cost is independent of the TF: 0.146-0.157 s WALL-CLOCK per monoPWM per Mb (44-47 us per 300-bp sequence); 0.17-0.20 s per diPWM per Mb; JVM start-up + class loading = 0.18 s per invocation (included in these figures; the marginal rate is 0.13-0.15 s/Mb).",
    "CPU time is higher than wall-clock because the JVM keeps JIT/GC helper threads even with -XX:ActiveProcessorCount=1: 0.22-0.24 s CPU per monoPWM per Mb for single scans and 0.27-0.29 s in the loops (user + system time; CPU/wall = 1.4-1.8). On a genuinely single core the wall-clock cost would approach the CPU figure.",
    "Scan time is linear in the number of PWMs k: T_scan ~ k x (0.18 s + 0.15 s/Mb x Mb) for monoPWMs (0.19 s/Mb for diPWMs).",
    "Relative to one PWM (2.3 s per 15 Mb) the full ArChIPelago feature step costs k times more: 11 s (STAT3, 5 monoPWMs), 58 s (HNF4A, 26) and 8.0 min (CTCF, 198); with the diPWMs of the mono+di models 25 s, 94 s and 13.1 min.",
    "Memory is flat (0.23-0.35 GB per mono scan, 0.53-0.64 GB per di scan) because PWMs are scanned one at a time; the scans are independent, so wall time divides by the number of cores used.",
    "Extrapolation (not measured): one monoPWM over the human genome (3.1 Gb) ~ 8 min wall-clock (~12 min CPU); the 198 CTCF monoPWMs ~ 26 h wall-clock-equivalent / ~49 CPU-h at 0.29 s CPU per PWM per Mb (~1-2 h on 24 cores). Single-PWM rows are medians of 3 runs, loops 1 run.",
]})
train = pd.DataFrame([
    ["STAT3", 5, "229,482 (3,000 + 226,482)", 1, 28.4, 0.03, 32.4, 32.0, 0.22, 0.98], ["STAT3", 5, "229,482", 4, 7.2, 0.03, 9.8, 31.7, 0.30, 0.98],
    ["HNF4A", 26, "201,495 (10,000 + 191,495)", 1, 38.2, 0.30, 43.0, 42.2, 0.36, 1.02], ["HNF4A", 26, "201,495", 4, 10.7, 0.14, 14.5, 45.3, 0.41, 1.02],
    ["CTCF", 198, "231,291 (10,000 + 221,291)", 1, 104.1, 1.62, 120.0, 115.1, 1.80, 0.97], ["CTCF", 198, "231,291", 4, 32.3, 0.98, 43.2, 138.3, 1.87, 0.97],
], columns=["TF", "Features (monoPWMs)", "Training rows (positives + negatives)", "Threads (n_jobs)", "RandomForest fit (s, in-process)", "StandardScaler (s)",
            "Whole process wall-clock incl. CSV reading (s)", "Whole process user CPU (s)", "Peak memory (GB)", "Model size on disk (MB)"])
train_notes = pd.DataFrame({"Notes": [
    "RandomForestClassifier(max_depth = 6, max_samples = 0.8, n_estimators = 100) on StandardScaler-transformed features, as in the paper. Training auROC (sanity check): STAT3 0.84, HNF4A 0.96, CTCF 0.98.",
    "The pickled model is about 1 MB for every TF, the fitted scaler 0.7-5.4 kB. Peak memory is dominated by the feature matrix, not by the forest. Speed-up with 4 threads: 3.2-3.9x. Each configuration was run once; the matrices hold the monoPWM features only (the mono+di models have about twice as many features).",
    "Reading the feature matrix takes 0.6 s (STAT3, 28 MB), 1.0-2.4 s (HNF4A, 45 MB) and 7-10 s warm / 26 s cold cache (CTCF, 257 MB).",
]})
infer = pd.DataFrame([
    ["chunked, 10 x 100,000 rows", "STAT3", 5, 1, 0.07, 6.7, 6.9, 8.8, 0.17, "146k"], ["chunked, 10 x 100,000 rows", "STAT3", 5, 4, 0.06, 1.9, 2.1, 4.0, 0.19, "482k"],
    ["chunked, 10 x 100,000 rows", "HNF4A", 26, 1, 0.33, 7.4, 8.2, 10.7, 0.29, "122k"], ["chunked, 10 x 100,000 rows", "HNF4A", 26, 4, 0.35, 2.7, 3.6, 5.9, 0.30, "278k"],
    ["chunked, 10 x 100,000 rows", "CTCF", 198, 1, 3.1, 19.6, 30.8, 41.3, 1.65, "32k"], ["chunked, 10 x 100,000 rows", "CTCF", 198, 4, 3.4, 6.8, 18.3, 29.0, 1.65, "55k"],
    ["single call on the full 1e6 x k matrix", "STAT3", 5, 1, 0.07, 7.7, np.nan, 10.0, 0.34, "130k"], ["single call on the full 1e6 x k matrix", "HNF4A", 26, 1, 2.1, 8.3, np.nan, 18.4, 0.81, "120k"],
    ["single call on the full 1e6 x k matrix (2 runs)", "CTCF", 198, 1, "13.8 / 21.3", "21.5 / 40.2", np.nan, "50.9 / 79.6", 5.13, "25-46k"],
    ["single call on the full 1e6 x k matrix (2 runs)", "CTCF", 198, 4, "31.5 / 49.4", "16.2 / 17.1", np.nan, "66.1 / 85.4", 5.21, "15-20k"],
], columns=["Variant (1,000,000 sequences)", "TF", "Features (monoPWMs)", "Threads (n_jobs)", "scaler.transform (s)", "predict_proba (s)", "Total in-process (s)",
            "Whole process wall-clock (s)", "Peak memory (GB)", "Sequences per second"])
infer_notes = pd.DataFrame({"Notes": [
    "predict_proba itself takes 7-20 s per 10^6 sequences single-threaded (2-7 s on 4 threads); with scaler.transform included 7-23 s, the whole chunked loop 7-31 s (the extra time for CTCF is the gathering of rows from the resident pool, a benchmark artefact). 'Sequences per second' = 10^6 / total in-process time for the chunked rows and 10^6 / predict_proba for the single-call rows. In the chunked runs peak memory is dominated by the resident matrix used as the sampling pool (not measured separately); the inference working set proper is the ~1 MB model plus one 100,000 x k chunk (arithmetic, a few hundred MB for CTCF).",
    "Feeding one monolithic 10^6 x 198 matrix raises memory to 5.1-5.2 GB and the timings become system-time dominated on the shared node, so the chunked numbers are the representative deployment figures (identical results).",
    "Companion scanner (ArChIPelago-TFBS-finder): cost ~ (1 + n_null) x k x (0.18 s + 0.15 s/Mb x Mb) + negligible RF time; the empirical-null calibration (default 50 shuffles) dominates and can be reduced with --n_null or switched off with --no-null.",
]})
with pd.ExcelWriter(os.path.join(OUT, "Sup_Table_6_runtime_memory.xlsx")) as xw:
    hw.to_excel(xw, sheet_name="setup", index=False)
    scan.to_excel(xw, sheet_name="PWM_scanning", index=False)
    scan_notes.to_excel(xw, sheet_name="PWM_scanning", index=False, startrow=len(scan) + 3)
    train.to_excel(xw, sheet_name="RF_training", index=False)
    train_notes.to_excel(xw, sheet_name="RF_training", index=False, startrow=len(train) + 3)
    infer.to_excel(xw, sheet_name="RF_inference", index=False)
    infer_notes.to_excel(xw, sheet_name="RF_inference", index=False, startrow=len(infer) + 3)
scan.to_csv(os.path.join(OUT, "Sup_Table_6_runtime_PWM_scanning.csv"), index=False)
train.to_csv(os.path.join(OUT, "Sup_Table_6_runtime_RF_training.csv"), index=False)
infer.to_csv(os.path.join(OUT, "Sup_Table_6_runtime_RF_inference.csv"), index=False)


# ----------------------------------------------------------------------------------------------
# saturation (Fig. S5) headline numbers -- against the best single monoPWM
# ($ARCHI_BASIS_DIR/saturation/, summaries written by analyze_saturation.py)
# ----------------------------------------------------------------------------------------------
SAT = os.path.join(BASIS, "saturation")
sat = pd.read_csv(os.path.join(SAT, "saturation_summary_by_k.csv"))
satk = {}
for metric in ["auroc_H", "auprc_H", "auroc_M", "auprc_M"]:
    s = sat[(sat.design == "random") & (sat.metric == metric)].set_index("k")
    full = s.loc["all", "median_delta"] if "all" in s.index else s.median_delta.iloc[-1]
    satk[metric] = {str(k): [round(v, 4), round(v / full, 2)] for k, v in s.median_delta.items()}
headline["saturation_random_median_delta_and_frac"] = satk
per_tf = pd.read_csv(os.path.join(SAT, "saturation_per_tf.csv"))
sat_per_tf = {}
for m in ["auroc_H", "auprc_H", "auroc_M", "auprc_M"]:
    col = per_tf[f"sat_k_{m}"]
    k = pd.to_numeric(col.where(col.astype(str) != "P", per_tf.P), errors="coerce")
    k = k[per_tf[f"full_gain_{m}"] > 0]
    # strict convention: the 90 % level is reached with a GENUINE subset (k < P); a TF whose curve reaches it only
    # at the first grid point >= P (carry-forward of the full model) is not counted
    P_ = per_tf.P.loc[k.index]
    strict = k.where(k < P_)
    sat_per_tf[m] = dict(n_tf=int(k.notna().sum()), median_k=float(k.median()), q25=float(k.quantile(.25)), q75=float(k.quantile(.75)),
                         median_k_lower=float(k.quantile(.5, interpolation="lower")),
                         frac_sat_by_8=round(float((k <= 8).mean()), 2), frac_sat_by_16=round(float((k <= 16).mean()), 2),
                         frac_sat_by_32=round(float((k <= 32).mean()), 2),
                         n_only_at_full_set_le16=int(((k <= 16) & (k >= P_)).sum()),
                         n_only_at_full_set_le32=int(((k <= 32) & (k >= P_)).sum()),
                         frac_sat_by_16_strict=round(float((strict <= 16).mean()), 2), frac_sat_by_32_strict=round(float((strict <= 32).mean()), 2),
                         n_more_pwms_hurt=int(per_tf[f"hurts_{m}"].sum()),
                         tfs_more_pwms_hurt=sorted(per_tf.tf[per_tf[f"hurts_{m}"]]))
headline["saturation_per_tf"] = sat_per_tf
headline["saturation_frac_at_32_all_four"] = {m: satk[m]["32"][1] for m in ["auroc_H", "auprc_H", "auroc_M", "auprc_M"]}
headline["saturation_first_k_ge_90pct"] = {m: int(min(int(k) for k, (v, f) in satk[m].items() if k not in ("all", "P") and f >= 0.9))
                                           for m in ["auroc_H", "auprc_H", "auroc_M", "auprc_M"]}
headline["saturation_n_fits"] = int(len(pd.read_csv(os.path.join(SAT, "saturation_results.csv"))))
headline["saturation_spearman_P_vs_gain"] = {m: [round(v, 3) for v in stats.spearmanr(per_tf.P, per_tf[f"full_gain_{m}"])]
                                             for m in ["auroc_H", "auprc_H", "auroc_M", "auprc_M"]}
headline["saturation_k_beats_baseline"] = {}
for metric in ["auroc_H", "auprc_H", "auroc_M", "auprc_M"]:
    for des in ["random", "topk"]:
        s_ = sat[(sat.design == des) & (sat.metric == metric)]
        s_ = s_[~s_.k.astype(str).isin(["all", "P"])]
        pos = s_[s_.median_delta > 0]
        headline["saturation_k_beats_baseline"][f"{metric}_{des}"] = int(pd.to_numeric(pos.k).min()) if len(pos) else None
sat_topk_minus_random = {}
for metric in ["auroc_H", "auprc_H"]:
    r_ = sat[(sat.design == "random") & (sat.metric == metric)].set_index("k").median_delta
    t_ = sat[(sat.design == "topk") & (sat.metric == metric)].set_index("k").median_delta
    sat_topk_minus_random[metric] = {str(k): round(float(t_[k] - r_[k]), 4) for k in r_.index if str(k) not in ("all", "P")}
headline["saturation_topk_minus_random_median_curve"] = sat_topk_minus_random

with open(os.path.join(OUT, "headline_numbers.json"), "w") as fh:
    json.dump(headline, fh, indent=1, default=str)
print("written to", OUT)
print(json.dumps({k: headline[k] for k in ["RF_H", "RF_M", "operational", "crossspecies"]}, indent=1, default=str))
