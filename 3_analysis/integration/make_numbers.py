"""NUMBERS_human_35TF.csv: every human number of the inventory (and the rows that quote them): old = first submission
(36 TFs), new = final data (35 TFs).

Old values are read from the frozen copies under results/previous_values/<git ref>/ (3_analysis/integration/h3_old.py), from
the released package results/round1_package/Manuscript_analysis or from the first-submission Sup. Tables
(results/round1_package/Resubmission), never from live tables. New values come from the final outputs:
  results table   3_analysis/integration/inputs/results_table.csv (= 3_analysis/package/results_table.csv; *_PWM = best monoPWM)
  Fig. 2          4_visualisation/main_figures/Figure_2/Figure_2_plotdata.csv   operational   outputs/supplementary_tables/Sup_Table_4 per_TF
  low-FPR metrics 3_analysis/low_fpr_metrics                              top peaks     3_analysis/top_peaks
  importance      3_analysis/feature_importance/tables                    drivers       3_analysis/transfer_drivers
  Fig. S5         3_analysis/package/Figures/source_data                  Sup. Table 2  outputs/supplementary_tables
  Fig. 4          4_visualisation/main_figures/Figure_4/Figure_4_stats.csv; random split (Fig. S12) from 3_analysis/random_split
Each old value that the manuscript quotes is asserted against the quoted text value (the old computation must reproduce
the first submission). Only rows whose value changes are written (unchanged checks are printed).
usage: python3 3_analysis/integration/make_numbers.py
"""
import hashlib
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

import p53drop as PD
sys.path.insert(0, os.path.join(PD.R2, "3_analysis", "integration"))
import h3_old as O  # noqa: E402

R2 = PD.R2
R1 = os.path.join(R2, "results", "round1_package", "Resubmission")
PKG = os.path.join(R2, "results/round1_package/Manuscript_analysis")
RR = "manuscript/"          # repository path prefix of the old-value tables
B13 = "3_analysis/"
PENDING = "pending (H3b)"
EX = list(PD.EXCLUDE_TFS)
rows, same = [], []


def p(path):
    return os.path.join(R2, path)


def old(path, **kw):
    """old (first-submission) copy of a table, from results/previous_values"""
    return O.csv(RR + path, **kw)


def src2(new, oldpath):
    return "%s (new); %s" % (new, O.src(RR + oldpath))


def fmt(v):
    if isinstance(v, str):
        return v
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    v = float(v)
    if v != 0 and abs(v) < 1e-3:
        return "%.2e" % v
    return "%.4g" % v


def add(item, old, new, src, inv, quoted=None):
    o, n = fmt(old), fmt(new)
    if quoted is not None:   # the old value must reproduce the number quoted in the first-submission text
        dec = len(quoted.split(".")[1]) if "." in quoted else 0
        if abs(float(old) - float(quoted)) > 0.5 * 10 ** -dec + 1e-12:
            raise SystemExit("CHECK FAILED %s: old %s does not reproduce the quoted %s" % (item, o, quoted))
    (rows if o != n else same).append(dict(item=item, old_value_36TF=o, new_value_35TF=n, source_table=src, inventory_row=inv))


# ---------------------------------------------------------------- headline (released-package results_table vs P53-free copy)
def rf_block(path):
    t = pd.read_csv(path, sep="\t")
    return t[(t.Model == "RandomForestClassifier") & (t.PWM == "mono+di")].set_index("TF_name")
RF_old, RF_new = rf_block(os.path.join(PKG, "results_table.csv")), rf_block(p("3_analysis/integration/inputs/results_table.csv"))
_t = pd.read_csv(p("3_analysis/integration/inputs/results_table.csv"), sep="\t")
_mono = _t[(_t.Model == "Single best mono PWM") & (_t.PWM == "mono")].drop_duplicates("TF_name").set_index("TF_name").loc[RF_new.index]
for _c in ("roc_auc_test_H", "pr_auc_test_H"):     # the baseline column must be the train-selected best monoPWM
    if not np.array_equal(RF_new[_c + "_PWM"].values, _mono[_c + "_PWM_mono"].values):
        raise SystemExit("CHECK FAILED results table: %s_PWM is not the best monoPWM" % _c)
# the new values are the final results table; the input copy is byte-identical (checked here)
if hashlib.md5(open(p("3_analysis/integration/inputs/results_table.csv"), "rb").read()).digest() != hashlib.md5(open(p(B13 + "package/results_table.csv"), "rb").read()).digest():
    raise SystemExit("CHECK FAILED 3_analysis/integration/inputs/results_table.csv is not the H3a pkg results table")
src = "PKG/results_table.csv (old); 3_analysis/package/results_table.csv (new, H3a results table; 3_analysis/integration/inputs/results_table.csv is a byte-identical copy); RF mono+di vs *_PWM, human test"
inv = "H01, M18, M25, M32, F01, R11"
for met, b, qb, qr, qg in (("roc_auc_test_H", "roc_auc_test_H_PWM", "0.861", "0.891", "0.024"),
                           ("pr_auc_test_H", "pr_auc_test_H_PWM", "0.255", "0.298", "0.041")):
    o, n = RF_old, RF_new
    nm = "auROC" if met.startswith("roc") else "auPRC"
    add("Human RF mono+di median %s, best monoPWM" % nm, o[b].median(), n[b].median(), src, inv, qb)
    add("Human RF mono+di median %s, ArChIPelago RF" % nm, o[met].median(), n[met].median(), src, inv, qr)
    add("Human median per-TF gain in %s" % nm, (o[met] - o[b]).median(), (n[met] - n[b]).median(), src, inv, qg)
    add("Human TFs improved, %s (n of N)" % nm, "%d/%d" % (((o[met] - o[b]) > 0).sum(), len(o)), "%d/%d" % (((n[met] - n[b]) > 0).sum(), len(n)), src, inv)
    add("Human Wilcoxon signed-rank P (two-sided), %s gain" % nm, wilcoxon(o[met] - o[b]).pvalue, wilcoxon(n[met] - n[b]).pvalue, src, inv)

# ---------------------------------------------------------------- Fig. 2 boxes
src = src2("4_visualisation/main_figures/Figure_2/Figure_2_plotdata.csv", "B6_figures/Fig2_S1/Figure_2_plotdata.csv")
f2o, f2n = old("B6_figures/Fig2_S1/Figure_2_plotdata.csv"), pd.read_csv(p("4_visualisation/main_figures/Figure_2/Figure_2_plotdata.csv"))
for (mdl, pw), g in f2o.groupby(["Model", "PWM"], sort=False):
    h = f2n[(f2n.Model == mdl) & (f2n.PWM == pw)]
    for m in ("auROC", "auPRC"):
        if mdl != "best monoPWM":
            add("Fig. 2 %s, %s on %s: TFs above best monoPWM" % (m, mdl, pw), "%d/%d" % ((g[m] > g["base_" + m]).sum(), len(g)),
                "%d/%d" % ((h[m] > h["base_" + m]).sum(), len(h)), src, "H02, M17, M19, F02")
        add("Fig. 2 %s, %s on %s: median" % (m, mdl, pw), g[m].median(), h[m].median(), src, "H02, F02")

# ---------------------------------------------------------------- operational metrics (first-submission Sup. Table 4 per-TF rows)
op = pd.read_csv(os.path.join(R1, "Sup_Tables", "Sup_Table_4_operational_metrics.csv"))
op = op[op["Test set"].str.startswith("human")]
on = pd.read_excel(p("outputs/supplementary_tables/Sup_Table_4_operational_metrics.xlsx"), sheet_name="per_TF")
on = on[on["Test set"].str.startswith("human")]
PD.check_tfs(on.TF, "Sup. Table 4 human rows")
src = "outputs/supplementary_tables/Sup_Table_4_operational_metrics.xlsx per_TF, human rows (new); R1/Sup_Tables/Sup_Table_4_operational_metrics.csv (old)"
inv = "H03, M20, M25, M32, F01, R02"
c50 = "FP reduction at 50 % recall, ArChIPelago vs best monoPWM (%)"
c80 = "FP reduction at 80 % recall, ArChIPelago vs best monoPWM (%)"
add("Human FP reduction at 50 % recall, median (%)", op[c50].median(), on[c50].median(), src, inv, "23.9")
add("Human FP reduction at 50 % recall, IQR low (%)", op[c50].quantile(.25), on[c50].quantile(.25), src, inv, "18.6")
add("Human FP reduction at 50 % recall, IQR high (%)", op[c50].quantile(.75), on[c50].quantile(.75), src, inv, "43.1")
add("Human TFs with fewer FP at 50 % recall (n of N)", "%d/%d" % ((op[c50] > 0).sum(), len(op)), "%d/%d" % ((on[c50] > 0).sum(), len(on)), src, inv)
add("Human FP reduction at 80 % recall, median (%)", op[c80].median(), on[c80].median(), src, inv, "24.1")
for c, q in (("Precision gain at 50 % recall, ArChIPelago minus best monoPWM", "0.042"), ("Precision at 50 % recall: best monoPWM", "0.176"),
             ("Precision at 50 % recall: ArChIPelago RF", "0.276"), ("Additional true sites in top 1,000, ArChIPelago vs best monoPWM", None)):
    add("Human %s, median" % c, op[c].median(), on[c].median(), src, inv, q)

# ---------------------------------------------------------------- headroom and Fig. 3C groups (transfer drivers)
cso, csn = old("B5_drivers/B5_B_ceiling_spearman.csv"), pd.read_csv(p("3_analysis/transfer_drivers/B5_B_ceiling_spearman.csv"))
src = src2("3_analysis/transfer_drivers/B5_B_ceiling_spearman.csv", "B5_drivers/B5_B_ceiling_spearman.csv")
inv = "H04, M21, R13, L15"
h = 1.96 / np.sqrt(np.array(cso.n, float) - 3)
cso["lo"], cso["hi"] = np.tanh(np.arctanh(cso.rho) - h), np.tanh(np.arctanh(cso.rho) + h)
for i, lab, qr, qp in ((0, "dauROC vs baseline auROC", "-0.65", None), (1, "dauPRC vs baseline auPRC", "0.52", None),
                       (2, "dauROC/(1 - baseline auROC) vs baseline auROC", "0.07", "0.68"), (4, "dauPRC/baseline auPRC vs baseline auPRC", "-0.38", "0.02")):
    add("Headroom Spearman rho, %s" % lab, cso.rho[i], csn.rho[i], src, inv, qr)
    add("Headroom Spearman P, %s" % lab, cso.p[i], csn.p[i], src, inv, qp)
add("Headroom rho 95 % CI (Fisher), dauROC/(1 - baseline auROC)", "%.2f to %.2f" % (cso.lo[2], cso.hi[2]),
    "%.2f to %.2f" % (csn.rho_ci95_lo_fisher[2], csn.rho_ci95_hi_fisher[2]), src, inv)
cmo, cmn = old("B5_drivers/B5_B_clouds_per_TF.csv"), pd.read_csv(p("3_analysis/transfer_drivers/B5_B_clouds_per_TF.csv"))
src = src2("3_analysis/transfer_drivers/B5_B_clouds_per_TF.csv", "B5_drivers/B5_B_clouds_per_TF.csv")
for col, lab in (("cloud_prc_high", "dauPRC > 0.15"), ("cloud_roc_high", "dauROC > 0.075"), ("cloud_ward2", "Ward outlying")):
    sel = lambda d: sorted(d.TF[d[col].isin(["dauPRC>0.15", "dauROC>0.075", "outlying"])])
    add("Fig. 3C group %s: TFs" % lab, "%d (%s)" % (len(sel(cmo)), ", ".join(sel(cmo))), "%d (%s)" % (len(sel(cmn)), ", ".join(sel(cmn))), src, inv)
sepo, sepn = old("B5_drivers/B5_B_cloud_separation.csv"), pd.read_csv(p("3_analysis/transfer_drivers/B5_B_cloud_separation.csv"))
src = src2("3_analysis/transfer_drivers/B5_B_cloud_separation.csv", "B5_drivers/B5_B_cloud_separation.csv") + " (q_BH over all Mann-Whitney tests)"
for sp_ in ("cloud_roc_high", "cloud_prc_high"):
    for v in ("base_auROC_H", "base_auPRC_H"):
        a, b = sepo[(sepo.split == sp_) & (sepo.variable == v)].iloc[0], sepn[(sepn.split == sp_) & (sepn.variable == v)].iloc[0]
        add("Fig. 3C %s vs rest, %s: median group / rest" % (sp_, v), "%.3f / %.3f" % (a.median_group, a.median_rest), "%.3f / %.3f" % (b.median_group, b.median_rest), src, inv)
        add("Fig. 3C %s vs rest, %s: Mann-Whitney P (q)" % (sp_, v), "%.2g (%.2g)" % (a.p_MWU, a.q_BH), "%.2g (%.2g)" % (b.p_MWU, b.q_BH), src, inv)
add("Fig. 3C Mann-Whitney tests in the BH family", int(sepo.p_MWU.notna().sum()), int(sepn.p_MWU.notna().sum()), src, inv)

# ---------------------------------------------------------------- low-FPR metrics
so, sn = old("B2_metrics/B2_summary.csv"), pd.read_csv(p(B13 + "low_fpr_metrics/B2_summary.csv"))
src = src2(B13 + "low_fpr_metrics/B2_summary.csv", "B2_metrics/B2_summary.csv") + ", HH, RF seed 0"
def b2(d, m):
    return d[(d.scenario == "HH") & (d.rf_seed == 0) & (d.metric == m)].iloc[0]
Q = {"tpr_at_fpr0.01": ("0.199", "0.273"), "pauc_std_fpr0.01": ("0.553", "0.574"), "prec_at_tpr0.4": ("0.255", "0.331"),
     "auroc": (None, "0.891"), "auprc": (None, "0.303"), "fp_per_kb_at_tpr0.8": ("0.72", "0.58")}
for m in ["auroc", "auprc", "pauc_std_fpr0.01", "tpr_at_fpr0.01", "tpr_at_fpr0.001", "prec_at_tpr0.4",
          "fpr_at_tpr0.2", "fpr_at_tpr0.4", "fpr_at_tpr0.6", "fpr_at_tpr0.8",
          "fp_per_kb_at_tpr0.2", "fp_per_kb_at_tpr0.4", "fp_per_kb_at_tpr0.6", "fp_per_kb_at_tpr0.8",
          "auprc_part_r0_0.2", "auprc_part_r0.2_0.4", "auprc_part_r0.4_0.6", "auprc_part_r0.6_0.8", "auprc_part_r0.8_1"]:
    a, b = b2(so, m), b2(sn, m)
    inv = "H05, H07, M13, M24" if m in ("auroc", "auprc") else "H05, M01, M22, M25, R01, R05, R14, F11, L13"
    qb, qr = Q.get(m, (None, None))
    add("B2 HH %s, best monoPWM median" % m, a.base_median, b.base_median, src, inv, qb)
    add("B2 HH %s, RF seed 0 median" % m, a.rf_median, b.rf_median, src, inv, qr)
    add("B2 HH %s, TFs improved (n of N)" % m, "%d/%d" % (a.n_improved, a.n_TF), "%d/%d" % (b.n_improved, b.n_TF), src, inv)
    add("B2 HH %s, Wilcoxon P" % m, a.wilcoxon_p, b.wilcoxon_p, src, inv)
parts = [m for m in so.metric.unique() if m.startswith("auprc_part")]
rng = lambda d: "%d to %d of %d" % (min(b2(d, m).n_improved for m in parts), max(b2(d, m).n_improved for m in parts), b2(d, parts[0]).n_TF)
add("B2 HH TFs improved per recall interval (range)", rng(so), rng(sn), src, "H05, M22, R05")
pco, pcn = old("B2_metrics/B2_per_chromosome_summary.csv"), pd.read_csv(p(B13 + "low_fpr_metrics/B2_per_chromosome_summary.csv"))
src = src2(B13 + "low_fpr_metrics/B2_per_chromosome_summary.csv", "B2_metrics/B2_per_chromosome_summary.csv") + ", HH"
for ch, q in (("chr1", "0.025"), ("chr8", "0.025"), ("chr21", "0.022")):
    for m in ("auroc", "auprc"):
        a = pco[(pco.scenario == "HH") & (pco.chrom == ch) & (pco.metric == m)].iloc[0]
        b = pcn[(pcn.scenario == "HH") & (pcn.chrom == ch) & (pcn.metric == m)].iloc[0]
        add("B2 HH %s %s gain, median" % (ch, m), a.diff_median, b.diff_median, src, "H06, M22, R07, R10, F11", q if m == "auroc" else None)
        add("B2 HH %s %s, TFs improved (n of N)" % (ch, m), "%d/%d" % (a.n_improved, a.n_TF), "%d/%d" % (b.n_improved, b.n_TF), src, "H06, M22, R07, R10, F11")
        add("B2 HH %s %s, Wilcoxon P" % (ch, m), a.wilcoxon_p, b.wilcoxon_p, src, "H06, F11")
wo, wn = old("B2_metrics/B2_metrics_wide.csv"), pd.read_csv(p(B13 + "low_fpr_metrics/B2_metrics_wide.csv"))
wo, wn = wo[(wo.scenario == "HH") & (wo.model == "RF_s0")], wn[(wn.scenario == "HH") & (wn.model == "RF_s0")]
PD.check_tfs(wn.TF, "B2 wide HH RF_s0")
r_o, r_n = wo.n_neg / wo.n_pos, wn.n_neg / wn.n_pos
src = src2(B13 + "low_fpr_metrics/B2_metrics_wide.csv", "B2_metrics/B2_metrics_wide.csv") + ", HH RF_s0, n_neg / n_pos"
add("Human test class ratio, median negatives per positive", r_o.median(), r_n.median(), src, "H08, M09, M17, L13", "42.3")
add("Human test class ratio, range", "%.1f to %.1f" % (r_o.min(), r_o.max()), "%.1f to %.1f" % (r_n.min(), r_n.max()), src, "H08, M09, L13")

# ---------------------------------------------------------------- strongest peaks (top-N)
to, tn = old("B3_peaks/b3_topN_summary.csv"), pd.read_csv(p(B13 + "top_peaks/b3_topN_summary.csv"))
src = src2(B13 + "top_peaks/b3_topN_summary.csv", "B3_peaks/b3_topN_summary.csv") + ", H, allneg, avail"
for N in ("500", "1000", "2000", "all"):
    a = to[(to.species == "H") & (to.variant == "allneg") & (to.subset == "avail") & (to.N.astype(str) == N)].iloc[0]
    b = tn[(tn.species == "H") & (tn.variant == "allneg") & (tn.subset == "avail") & (tn.N.astype(str) == N)].iloc[0]
    for c in ("n_tf", "n_improved_auroc", "n_improved_auprc", "pwm_auroc_med", "rf_auroc_med", "d_auroc_med", "d_auprc_med", "p_wilcoxon_auroc",
              "d_auroc_med_allpos_same_tfs", "topN_minus_allpos_auroc_med", "p_topN_vs_allpos_auroc"):
        if c in a.index and not (pd.isna(a[c]) and pd.isna(b[c])):
            add("Top-N human N = %s: %s" % (N, c), a[c], b[c], src, "H09, M22, M25, R02")

# ---------------------------------------------------------------- gain vs positives / datasets (top peaks)
go, gn = old("B3_peaks/b3_gain_vs_npos_spearman.csv"), pd.read_csv(p(B13 + "top_peaks/b3_gain_vs_npos_spearman.csv"))
src = src2(B13 + "top_peaks/b3_gain_vs_npos_spearman.csv", "B3_peaks/b3_gain_vs_npos_spearman.csv") + ", human test rows"
Q10 = {("dauroc_H", "n_train_pos_H"): ("-0.33", "0.053"), ("dauprc_H", "n_train_pos_H"): ("0.25", "0.14"),
       ("dauprc_H", "n_datasets_H"): ("0.43", "0.009"), ("dauroc_H", "n_datasets_H"): ("-0.06", "0.71")}
for (gg, x), g in go[go.test_set == "human test"].groupby(["gain", "predictor_col"], sort=False):
    a, b = g.iloc[0], gn[(gn.test_set == "human test") & (gn.gain == gg) & (gn.predictor_col == x)].iloc[0]
    qr, qp = Q10.get((gg, x), (None, None))
    add("Spearman %s vs %s, rho" % (gg, x), a.rho, b.rho, src, "H10, R12, F12", qr)
    add("Spearman %s vs %s, P" % (gg, x), a.p, b.p, src, "H10, R12, F12", qp)

# ---------------------------------------------------------------- feature importance
T4 = "3_analysis/feature_importance/tables/"
co, cn = old("B4_importance/tables/B4_importance_concentration.csv"), pd.read_csv(p(T4 + "B4_importance_concentration.csv"))
src = src2(T4 + "B4_importance_concentration.csv", "B4_importance/tables/B4_importance_concentration.csv")
inv = "H11, M27, R03, R15"
for c, q in (("n_for_50", "6"), ("n_for_80", "15")):
    add("PWMs carrying %s %% of the Gini importance, median" % c[-2:], co[c].median(), cn[c].median(), src, inv, q)
    add("PWMs carrying %s %% of the Gini importance, range" % c[-2:], "%d-%d" % (co[c].min(), co[c].max()), "%d-%d" % (cn[c].min(), cn[c].max()), src, inv)
add("Top PWM is a diPWM (n of N)", "%d/%d" % ((co.top1_kind == "di").sum(), len(co)), "%d/%d" % ((cn.top1_kind == "di").sum(), len(cn)), src, "H11")
add("Feature PWMs (gini_all_PWMs rows)", int(co.n_features.sum()), int(cn.n_features.sum()), src, "H11, L07, R03, T15")
sso, ssn = old("B4_importance/tables/B4_similarity_summary.csv"), pd.read_csv(p(T4 + "B4_similarity_summary.csv"))
add("Median pairwise similarity of the monoPWMs of a TF (median over TFs)", sso.median_pairwise_mono.median(), ssn.median_pairwise_mono.median(),
    src2(T4 + "B4_similarity_summary.csv", "B4_importance/tables/B4_similarity_summary.csv"), inv, "0.91")
suo, sun = old("B4_importance/tables/B4_subtype_summary_per_TF.csv"), pd.read_csv(p(T4 + "B4_subtype_summary_per_TF.csv"))
src = src2(T4 + "B4_subtype_summary_per_TF.csv", "B4_importance/tables/B4_subtype_summary_per_TF.csv") + ", cut 0.8"
sec = lambda d: sorted(d.TF[d.n_ge10pct_80 >= 2])
add("TFs with a secondary subtype carrying >= 10 % of the importance", "%d (%s)" % (len(sec(suo)), ", ".join(sec(suo))), "%d (%s)" % (len(sec(sun)), ", ".join(sec(sun))), src, "H11, M14, M27, R03, L14")
add("TFs with a non-primary subtype (>= 2 subtypes at cut 0.8; Fig. S10D)", int((suo.n_subtypes_80 >= 2).sum()), int((sun.n_subtypes_80 >= 2).sum()), src, "H11, L14")
abo, abn = old("B4_importance/tables/B4_ablation_summary.csv"), pd.read_csv(p(T4 + "B4_ablation_summary.csv"))
def ablation_pass(d):
    out = []
    for tf, g in d.groupby("TF"):
        s, r = g[g.cond == "drop_secondary"].iloc[0], g[g.cond == "drop_random"].iloc[0]
        if s.d_auroc_mean < r.d_auroc_mean - 3 * r.d_auroc_sd or s.d_auprc_mean < r.d_auprc_mean - 3 * r.d_auprc_sd:
            out.append(tf)
    return out
pa, pn = ablation_pass(abo), ablation_pass(abn)
if len(pa) != 8:
    raise SystemExit("CHECK FAILED ablation rule reproduces %d TFs, not the quoted 8" % len(pa))
src = src2(T4 + "B4_ablation_summary.csv", "B4_importance/tables/B4_ablation_summary.csv") + "; rule: secondary removal below random removal by > 3 sd (auROC or auPRC)"
add("Ablated TFs", abo.TF.nunique(), abn.TF.nunique(), src, "H11, M14, M27, R03, R15, L07, L14")
add("Ablated TFs where the loss exceeds 3 sd of the random removals", "%d of %d (%s)" % (len(pa), abo.TF.nunique(), ", ".join(pa)),
    "%d of %d (%s)" % (len(pn), abn.TF.nunique(), ", ".join(pn)), src, "H11, M27, R03, R15")
gro, grn = old("B4_importance/tables/B4_activation_groups.csv"), pd.read_csv(p(T4 + "B4_activation_groups.csv"))
def act(d):
    g = d[(d.defn == "any_nonprimary") & (d.cls == "pos")]
    po, so_ = g[g.group == "primary_only"].set_index("TF"), g[g.group == "secondary_only"].set_index("TF")
    t = [x for x in so_.index if so_.loc[x, "n"] >= 20 and po.loc[x, "n"] >= 20]
    return sum(so_.loc[x, "median_rf"] < po.loc[x, "median_rf"] for x in t), len(t)
add("Secondary-only strong hit gives the lower median RF score (n of TFs with >= 20 positives in both groups)", "%d of %d" % act(gro), "%d of %d" % act(grn),
    src2(T4 + "B4_activation_groups.csv", "B4_importance/tables/B4_activation_groups.csv"), "H11, M27, R15, L14")

# ---------------------------------------------------------------- random split: retraining on negatives v3 (5 TFs)
# old values from the frozen copy (first submission as quoted, P53 included); new from 3_analysis/random_split/b9_summary.csv
bo = old("B9_random_split/b9_summary.csv")
bn = pd.read_csv(p(B13 + "random_split/b9_summary.csv"))
src = B13 + "random_split/b9_summary.csv (new, H3b v3 negatives); " + O.src(RR + "B9_random_split/b9_summary.csv") + " (old); means of 3 seeds"
inv = "H12, M11, M13, R07, F15, L16"


def b9(b):
    r_, c_ = b[b.split == "random"].set_index("TF"), b[b.split == "chromosome"].set_index("TF")
    return (r_, r_.gain_auprc - c_.gain_auprc.reindex(r_.index), r_.rf_auroc - c_.rf_auroc.reindex(r_.index),
            r_.pwm_auroc - c_.pwm_auroc.reindex(r_.index), r_.gain_auroc - c_.gain_auroc.reindex(r_.index))


(ro, dgo, dro, dpo, gao), (rn, dgn, drn, dpn, gan) = b9(bo), b9(bn)
if not set(rn.index) <= set(PD.reference_tfs()) or len(rn) != 5:
    raise SystemExit("CHECK FAILED B9 (H3b): TFs %s, expected 5 reference TFs without P53" % sorted(rn.index))
add("Random-split retraining: TFs", len(ro), len(rn), src, inv)
add("Random split: median share of test positives with a training positive within 300 bp (%)", 100 * ro.frac_test_pos_overlap.median(), 100 * rn.frac_test_pos_overlap.median(), src, inv, "81")
add("Random split: median share of test positives at the same summit (%)", 100 * ro.frac_test_pos_exact.median(), 100 * rn.frac_test_pos_exact.median(), src, inv, "52")
add("Random split raised the auPRC gain (n of N TFs)", "%d of %d" % ((dgo > 0).sum(), len(dgo)), "%d of %d" % ((dgn > 0).sum(), len(dgn)), src, inv)
add("Random split: largest rise of the auPRC gain", dgo.max(), dgn.max(), src, inv, "0.045")
add("Random minus chromosome split, RF auROC (range)", "%+.3f to %+.3f" % (dro.min(), dro.max()), "%+.3f to %+.3f" % (drn.min(), drn.max()), src, inv)
add("Random minus chromosome split, best monoPWM auROC (range)", "%+.3f to %+.3f" % (dpo.min(), dpo.max()), "%+.3f to %+.3f" % (dpn.min(), dpn.max()), src, inv)
add("Random minus chromosome split, auROC gain (range)", "%+.3f to %+.3f" % (gao.min(), gao.max()), "%+.3f to %+.3f" % (gan.min(), gan.max()), src, inv)

fo = old("B6_figures/Fig4_S4/Figure_4_stats.csv")
s4o = pd.read_csv(os.path.join(PKG, "Figures/source_data/Figure_4_source_data.csv"))
# new values from the Fig. 4 stats table that the figure producer writes from the Fig. 4 source data (via the package copy);
# medians, TFs above 0 and the one-sided P equal Figure_4_and_S4_numbers.csv of the server run
F4N = "4_visualisation/main_figures/Figure_4/Figure_4_stats.csv"
fn = pd.read_csv(p(F4N)).set_index(["metric", "model"])
if set(fn.index) != set(zip(fo.metric, fo.model)) or (fn.n != 35).any():
    raise SystemExit("CHECK FAILED Fig. 4 (H3b): stats rows differ from round 1 or n != 35")
src = F4N + " (new, H3b: Fig. 4 models on negatives v3); " + O.src(RR + "B6_figures/Fig4_S4/Figure_4_stats.csv") + " (old); old two-sided P from PKG/Figures/source_data/Figure_4_source_data.csv"
Q4 = {("auROC", "ArChIPelago RF2f"): "0.009", ("auPRC", "ArChIPelago RF2f"): "0.004", ("auROC", "ArChIPelago RF2f + diChIPMunk"): "0.015",
      ("auROC", "ArChIPelago RF2f + Slim m=1"): "0.020", ("auROC", "ArChIPelago RF2f + LSlim m=-5"): "0.018",
      ("auROC", "ArChIPelago RF2f + Slim m=1, LSlim m=-5 + diChIPMunk"): "0.023", ("auPRC", "ArChIPelago RF2f + Slim m=1, LSlim m=-5 + diChIPMunk"): "0.030",
      ("auROC", "ArChIPelago RF on all PWMs"): "0.024", ("auPRC", "ArChIPelago RF on all PWMs"): "0.041"}
for _, a in fo.iterrows():
    nw = fn.loc[(a.metric, a.model)]
    add("Fig. 4 %s, %s: median delta" % (a.metric, a.model), a.med, nw.med, src, "H14, M29, F07", Q4.get((a.metric, a.model)))
    add("Fig. 4 %s, %s: TFs above 0" % (a.metric, a.model), "%d/%d" % (a.n_pos, a.n), "%d/%d" % (nw.n_pos, nw.n), src, "H14, M29, M31, F07")
    add("Fig. 4 %s, %s: one-sided Wilcoxon P" % (a.metric, a.model), a.p, nw.p, src, "H14, M29, F07")
    add("Fig. 4 %s, %s: two-sided Wilcoxon P" % (a.metric, a.model), wilcoxon(s4o[s4o.model == a.model]["d_" + a.metric]).pvalue, nw.p_two_sided, src, "H14, M29")

# ---------------------------------------------------------------- Fig. S5 A, B median curves (human columns)
S5N = p(B13 + "package/Figures/source_data/Figure_S5_source_data.csv")
s5o = pd.read_csv(os.path.join(PKG, "Figures/source_data/Figure_S5_source_data.csv")); s5n = pd.read_csv(S5N)
mc = pd.read_csv(os.path.join(PKG, "Figures/source_data/Figure_S5_median_curves.csv"))
src = B13 + "package/Figures/source_data/Figure_S5_source_data.csv (new, H3a saturation); PKG/Figures/source_data/Figure_S5_source_data.csv (old; checked against the released median curves)"
PD.check_tfs(s5n.TF, "Fig. S5 source data (H3)")
def s5_curve(s5, design, k, c):   # analyze_saturation.py: carry-forward k_eff = min(k, P); k = P from the random design
    mean5 = s5.groupby(["TF", "design", "k"])[["delta_auroc_H", "delta_auprc_H"]].mean()
    Pmap = s5.groupby("TF").P.first()
    v = []
    for tf in sorted(Pmap.index):
        P = Pmap[tf]
        ke = P if k == "P" else min(int(k), P)
        v.append(mean5.loc[(tf, "random" if ke == P else design, ke), c])
    return float(np.median(v))
for design in ("random", "topk"):
    for k in mc[mc.design == design].k:
        for c in ("delta_auroc_H", "delta_auprc_H"):
            o = s5_curve(s5o, design, k, c)
            ref = float(mc[(mc.design == design) & (mc.k == k)]["d_" + c[6:]].iloc[0])
            if abs(ref - o) > 1e-9:
                raise SystemExit("CHECK FAILED Fig. S5 median curve %s k=%s %s: %s vs released %s" % (design, k, c, o, ref))
            add("Fig. S5 %s k = %s, median %s" % (design, k, c), o, s5_curve(s5n, design, k, c), src, "H15, M32, F09, L11")

# ---------------------------------------------------------------- Sup. Table 2 human (and mouse PWM counts), Sup. Table 1 counts
T2R1 = os.path.join(R1, "Sup_Tables", "Sup_Table_2_PWMs_and_datasets.xlsx")
T1R1 = os.path.join(R1, "Sup_Tables", "Sup_Table_1_ChIP-Seq_experiments.xlsx")
t2o = pd.read_excel(T2R1, sheet_name="per_TF")
t2n = pd.read_excel(p("outputs/supplementary_tables/Sup_Table_2_PWMs_and_datasets.xlsx"), sheet_name="per_TF")
PD.check_tfs(t2n.TF, "Sup. Table 2 (H3)")
src = "outputs/supplementary_tables/Sup_Table_2_PWMs_and_datasets.xlsx per_TF (new); R1/Sup_Tables/Sup_Table_2_PWMs_and_datasets.xlsx (old, 36 TFs)"
inv = "H16, M06, M15, M16, F01, R04"
Q16 = {"Human monoPWMs": "1495", "Human diPWMs": "780", "Human PWMs (mono + di)": "2275", "Mouse monoPWMs": "1192", "Mouse diPWMs": "624", "Mouse PWMs (mono + di)": "1816"}
# the training column of the new table carries the chromosome label of the final run (h3_inputs COL); the old column keeps its name
import h3_inputs as HI  # noqa: E402
NEWCOL = {"Human training positives (chr2-7, 9, 10, 13-20)": HI.COL["HTR"], "Human test positives (chr1, 8, 21)": HI.COL["HTE"]}
for c in ("Human monoPWMs", "Human diPWMs", "Human PWMs (mono + di)", "Mouse monoPWMs", "Mouse diPWMs", "Mouse PWMs (mono + di)",
          "Human training positives (chr2-7, 9, 10, 13-20)", "Human test positives (chr1, 8, 21)", "Human pooled peak set (peaks, all chromosomes)"):
    cn = NEWCOL.get(c, c)
    nv = (int(t2n[cn].sum()), t2n[cn].median()) if cn in t2n.columns else ("removed", "removed")
    # item label: the column name of the final run (h3_inputs COL), so that the label names the chromosomes the value counts
    add("Sup. Table 2 total, %s" % cn, int(t2o[c].sum()), nv[0], src, inv, Q16.get(c))
    add("Sup. Table 2 median per TF, %s" % cn, t2o[c].median(), nv[1], src, inv, "34.5" if c == "Human PWMs (mono + di)" else ("37" if c == "Mouse PWMs (mono + di)" else None))
for c in [c for c in t2n.columns if c.startswith("Human") and "GC" in c]:      # the first submission gave mean GC; the new table gives medians (negatives v3)
    add("Sup. Table 2 median per TF, %s" % c, "-", t2n[c].median(), src, "H16")
t1o = pd.read_excel(T1R1, sheet_name=None)
_miss = set(pd.read_csv(p("results/datasets/top_peaks/b3_missing_datasets.csv")).peak_set)   # 19 human sets without a peak file, removed from the earlier copy
t1o = {k: v[~v["GTRD peak set ID"].isin(_miss)] for k, v in t1o.items()}
hnew = PD.human_ids()
src1 = "outputs/supplementary_tables/Sup_Table_1_ChIP-Seq_experiments.xlsx (new: sets used); R1/Sup_Tables/Sup_Table_1_ChIP-Seq_experiments.xlsx (old: rows listed, minus the 19 sets of results/datasets/top_peaks/b3_missing_datasets.csv)"
add("Human ChIP-Seq datasets used (Sup. Table 1)", len(t1o["human"]), len(hnew), src1, "M01, M03, M16, T01, L01, F01", "395")
add("Human datasets per TF, median (used sets)", t2o["Human ChIP-Seq experiments (Sup. Table 1)"].median(), hnew.groupby("TF").size().median(),
    "results/datasets/peak_ids/human_used_peak_ids.csv without P53 (new); R1/Sup_Tables/Sup_Table_2_PWMs_and_datasets.xlsx (old)", "H16, M03")
add("Human datasets per TF, minimum (used sets)", int(t2o["Human ChIP-Seq experiments (Sup. Table 1)"].min()), int(hnew.groupby("TF").size().min()),
    "results/datasets/peak_ids/human_used_peak_ids.csv without P53 (new); R1/Sup_Tables/Sup_Table_2_PWMs_and_datasets.xlsx (old)", "M03, G9")
add("Datasets of both species (manuscript total)", len(t1o["human"]) + len(t1o["mouse"]), len(hnew) + len(PD.mouse_ids()), src1, "M01, M02, M15, M16, R01, L06, F01", "685")

# ---------------------------------------------------------------- peak widths (top peaks)
wo, wn = old("B3_peaks/b3_peak_widths_summary.csv"), pd.read_csv(p(B13 + "top_peaks/b3_peak_widths_summary.csv"))
src = src2(B13 + "top_peaks/b3_peak_widths_summary.csv", "B3_peaks/b3_peak_widths_summary.csv") + " (new: 361 human + 276 mouse used peak sets)"
Q5 = {("H", "width_median_of_dataset_medians"): "575", ("M", "width_median_of_dataset_medians"): "551", ("H", "width_median_top1000_median"): "898",
      ("M", "width_median_top1000_median"): "812.75", ("all", "width_median_top1000_median"): "857.5", ("all", "width_median_rest_median"): "553",
      ("all", "n_datasets_ge1000_peaks"): "683"}
for sp_ in ("H", "M", "all"):
    a, b = wo[wo.species == sp_].iloc[0], wn[wn.species == sp_].iloc[0]
    for c in ("n_datasets", "width_median_of_dataset_medians", "width_median_top1000_median", "width_median_rest_median", "n_datasets_ge1000_peaks",
              "n_datasets_top1000_wider", "frac_cover_ge50_median_over_datasets"):
        add("Peak widths %s: %s" % (sp_, c), a[c], b[c], src, "M05, M32, R02, R09, T14, L06", Q5.get((sp_, c)))

out = pd.DataFrame(rows).rename(columns={"old_value_36TF": "old value (36 TFs)", "new_value_35TF": "new value (35 TFs)",
                                          "source_table": "source table", "inventory_row": "inventory row"})
dst = os.path.join(PD.B12, "../../outputs/numbers/NUMBERS_human_35TF.csv")
out.to_csv(dst, index=False)
pd.DataFrame(same).to_csv(os.path.join(PD.B12, "../../outputs/numbers/NUMBERS_human_35TF_unchanged.csv"), index=False)
print("%d changed numbers -> %s; %d unchanged -> NUMBERS_human_35TF_unchanged.csv" % (len(out), dst, len(same)))
srcs = [os.path.join(PKG, "results_table.csv"), p("3_analysis/integration/inputs/results_table.csv"), os.path.join(PKG, "Figures/source_data/Figure_4_source_data.csv"),
        os.path.join(R1, "Sup_Tables", "Sup_Table_4_operational_metrics.csv"), p("outputs/supplementary_tables/Sup_Table_4_operational_metrics.xlsx"),
        os.path.join(PKG, "Figures/source_data/Figure_S5_source_data.csv"), S5N, T2R1, T1R1, p("results/datasets/top_peaks/b3_missing_datasets.csv"), p("outputs/supplementary_tables/Sup_Table_2_PWMs_and_datasets.xlsx"),
        p(B13 + "low_fpr_metrics/B2_summary.csv"), p(B13 + "low_fpr_metrics/B2_per_chromosome_summary.csv"), p(B13 + "low_fpr_metrics/B2_metrics_wide.csv"), p(B13 + "top_peaks/b3_topN_summary.csv"),
        p(B13 + "top_peaks/b3_gain_vs_npos_spearman.csv"), p(B13 + "top_peaks/b3_peak_widths_summary.csv"), p("3_analysis/transfer_drivers/B5_B_ceiling_spearman.csv"),
        p("3_analysis/transfer_drivers/B5_B_clouds_per_TF.csv"), p("3_analysis/transfer_drivers/B5_B_cloud_separation.csv"), p("4_visualisation/main_figures/Figure_2/Figure_2_plotdata.csv")]
PD.write_provenance(PD.B12, "3_analysis/integration/make_numbers.py", [dst, os.path.join(PD.B12, "../../outputs/numbers/NUMBERS_human_35TF_unchanged.csv")], srcs,
                    note="H3: old values from git tag %s, PKG or R1 (asserted against the quoted text numbers); new values from the H3 outputs" % O.TAG)
