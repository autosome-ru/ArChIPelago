#!/usr/bin/env python3
"""Drivers of human->mouse transfer, and the ceiling effect with the two clouds of Fig. 3C.
Read-only on inputs; writes B5_* files next to this script. Run: OMP_NUM_THREADS=1 nice -n 15 python3 3_analysis/transfer_drivers/B5_01_drivers_and_ceiling.py
Python 3.9 compatible.
TFs in p53drop.EXCLUDE_TFS are removed (35 TFs); human dataset counts come from the 361 used peak sets.
All three tables are required arguments, without defaults, and are recorded in _provenance.yaml:
  --results-table   results table with *_PWM = the train-selected best monoPWM (3_analysis/package/results_table.csv)
  --mouse-transfer  run_hm_mm results table (results/server/human/hm_mm_results.csv)
  --crossspecies    cross-species table with celltype_overlap, topH_vs_nearestM_pcc, family, TFclass_id
The H->M gains are computed here from the unrounded hm_mm values and checked against the cross-species table;
n_PWMs_H = monoPWMs + diPWMs of hm_mm (P_mono_H + P_di_H).
"""
import os, hashlib, platform, subprocess, sys, datetime, argparse
import numpy as np, pandas as pd, scipy
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster

OUT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(OUT, "..", "integration/"))
import p53drop as PD
MS = os.path.join(OUT, "..", "package")
ap = argparse.ArgumentParser()
ap.add_argument("--results-table", required=True)
ap.add_argument("--mouse-transfer", required=True)
ap.add_argument("--crossspecies", required=True)
ARGS = ap.parse_args()
MOUSE_STATUS = "current"
INP = {
    "crossspecies_table": ARGS.crossspecies,
    "results_table": ARGS.results_table,
    "mouse_transfer": ARGS.mouse_transfer,
    "single_PWM_features": MS + "/analysis/inputs/single_PWM_features.csv",
    "human_used_peak_ids": PD.HUMAN_IDS,
}
ROC_GAP = 0.075    # Fig. 3C: empty band of dauROC above the core (0.060 -> 0.079 in the first submission)
PRC_GAP = 0.15     # Fig. 3C: auPRC gain separating the right-hand group
MOUSE_VARS = {"base_auROC_M", "base_auPRC_M", "celltype_overlap_ord"}   # depend on the mouse test set / mouse datasets
OUTFILES = ["B5_per_TF_table.csv", "B5_A_spearman_drivers.csv", "B5_A_sensitivity_auPRC_baseline_IC.csv", "B5_A_multivariable_rank_ols.csv",
            "B5_A_multivariable_summary.csv", "B5_B_ceiling_spearman.csv", "B5_B_headroom_per_TF.csv", "B5_B_clouds_per_TF.csv",
            "B5_B_cloud_separation.csv", "B5_B_cloud_members.txt", "B5_checks.csv", "_provenance.yaml"]
for _f in OUTFILES:
    PD.backup(os.path.join(OUT, _f))
PWM_DIR = os.path.join(os.environ.get("ARCHI_ZENODO_DIR", os.path.join(OUT, "../../results/datasets/zenodo_v1.0.0")), "PWMs_mono_HUMAN")
BG = 0.25          # background probability per base
SEED = 1
NBOOT = 2000


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); prev = 1.0
    for rank in range(n, 0, -1):
        i = o[rank - 1]; prev = min(prev, p[i] * n / rank); q[i] = prev
    return q


def md5(f):
    h = hashlib.md5()
    with open(f, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def pwm_ic(path):
    """PWM files hold log-odds w = ln(p / 0.25) (natural log; checked: sum_j 0.25*exp(w_ij) = 1 in every column).
    PFM p_ij = BG*exp(w_ij); IC_i = 2 + sum_j p_ij log2 p_ij (bits, background 0.25); total IC = sum_i IC_i."""
    m = np.loadtxt(path, skiprows=1).reshape(-1, 4)
    p = BG * np.exp(m)
    assert np.allclose(p.sum(1), 1, atol=1e-6), path
    ic = 2 + (p * np.log2(np.clip(p, 1e-12, None))).sum(1)
    return ic.sum(), len(m), ic.mean()


# ---------------------------------------------------------------- per-TF table
cs = pd.read_csv(INP["crossspecies_table"]).set_index("TF").sort_index()
res = pd.read_csv(INP["results_table"], sep="\t")
rf = res[(res.Model == "RandomForestClassifier") & (res.PWM == "mono+di")].set_index("TF_name").sort_index()
hm = pd.read_csv(INP["mouse_transfer"]).set_index("TF").sort_index()
cs, rf, hm = [x[~x.index.isin(PD.EXCLUDE_TFS)] for x in (cs, rf, hm)]
feat = pd.read_csv(INP["single_PWM_features"]).set_index(["TF", "feature"])
for _n, _x in (("crossspecies", cs), ("results_table RF", rf), ("mouse_transfer", hm)):
    PD.check_tfs(_x.index, "B5 " + _n)
assert list(cs.index) == list(rf.index) == list(hm.index)
NTF = len(cs)

T = pd.DataFrame(index=cs.index)
T["family"] = cs.family
T["tfclass_superclass"] = cs.TFclass_id.astype(str).str.split(".").str[0]
for kind, colname in (("auROC", "HM_broc_col"), ("auPRC", "HM_bprc_col")):
    names, tot, ln, mean = [], [], [], []
    for tf in T.index:
        col = hm.loc[tf, colname]
        nm = feat.loc[(tf, col), "pwm_name"]
        path = "%s/%s/%s.pwm" % (PWM_DIR, tf, col.split("_")[1])
        assert open(path).readline().strip().lstrip(">") == nm
        a, b, c = pwm_ic(path); assert b == feat.loc[(tf, col), "motif_length"]
        names.append(nm); tot.append(a); ln.append(b); mean.append(c)
    T["base_pwm_" + kind] = names; T["IC_total_" + kind] = tot; T["motif_len_" + kind] = ln; T["IC_mean_per_pos_" + kind] = mean
T["same_baseline_pwm"] = T.base_pwm_auROC == T.base_pwm_auPRC
T["n_train_pos_H"] = rf.Seq_count
if not (rf.Count == hm.P_mono_H).all():   # the Count of the mono+di row is the monoPWM count
    raise SystemExit("CHECK FAILED results_table Count of the mono+di rows vs hm_mm P_mono_H")
T["n_PWMs_H"] = hm.P_mono_H + hm.P_di_H       # monoPWMs + diPWMs, as colour in Fig. 3C
# human dataset counts: the used peak sets (human_used_peak_ids.csv without the excluded TFs, 361)
_hid = PD.human_ids()
T["n_datasets_H"] = _hid.groupby("TF").size().reindex(T.index)
if T.n_datasets_H.isna().any() or T.n_datasets_H.sum() != len(_hid):
    raise SystemExit("CHECK FAILED n_datasets_H: %s of %d" % (T.n_datasets_H.sum(), len(_hid)))
T["base_auROC_H"] = rf.roc_auc_test_H_PWM
T["base_auPRC_H"] = rf.pr_auc_test_H_PWM
T["base_auROC_M"] = hm.HM_broc_test_M_auroc
T["base_auPRC_M"] = hm.HM_bprc_test_M_auprc
T["motif_sim_H_M"] = cs.topH_vs_nearestM_pcc
T["celltype_overlap"] = cs.celltype_overlap.str.replace("?", "", regex=False)
T["celltype_overlap_ord"] = T.celltype_overlap.map({"yes": 1.0, "no": 0.0})   # the written rule (3_analysis/integration/celltype_overlap.py) gives yes/no only
assert T.celltype_overlap_ord.notna().all(), "CHECK FAILED celltype_overlap label outside yes/no: %s" % sorted(set(T.celltype_overlap))
T["prevalence_H_test"] = hm.n_test_pos_H / (hm.n_test_pos_H + hm.n_test_neg_H)
# outcomes
# human-trained RF on the mouse test set vs the best human monoPWM, unrounded from hm_mm (the table rounds to 4 dp).
# the outcome is the mean of the two RF seeds, as for the mouse-trained models; the seed-0 gains are
# kept only to check the input against the cross-species table (which holds seed 0) and for the seed spread in B5_checks.csv
s0_roc = hm.HM_rf_test_M_auroc_s0 - hm.HM_broc_test_M_auroc; s0_prc = hm.HM_rf_test_M_auprc_s0 - hm.HM_bprc_test_M_auprc
if not (np.allclose(s0_roc, cs.dROC_M, rtol=0, atol=5e-4) and np.allclose(s0_prc, cs.dPRC_M, rtol=0, atol=5e-4)):
    raise SystemExit("CHECK FAILED H->M gains: hm_mm seed 0 vs cross-species table")
T["dROC_HM"] = hm.HM_rf_test_M_auroc_mean - hm.HM_broc_test_M_auroc; T["dPRC_HM"] = hm.HM_rf_test_M_auprc_mean - hm.HM_bprc_test_M_auprc
hm_roc = hm.HM_rf_test_M_auroc_mean - hm.HM_broc_test_M_auroc
hm_prc = hm.HM_rf_test_M_auprc_mean - hm.HM_bprc_test_M_auprc
mm_roc = hm.MM_rf_test_M_auroc_mean - hm.MM_broc_test_M_auroc
mm_prc = hm.MM_rf_test_M_auprc_mean - hm.MM_bprc_test_M_auprc
T["gain_HM_roc_2seed"] = hm_roc; T["gain_HM_prc_2seed"] = hm_prc
T["gain_MM_roc_2seed"] = mm_roc; T["gain_MM_prc_2seed"] = mm_prc
T["loss_roc"] = mm_roc - hm_roc        # transfer loss = mouse-trained gain minus human-trained gain, mouse test set
T["loss_prc"] = mm_prc - hm_prc
T["loss_rf_only_roc"] = hm.MM_rf_test_M_auroc_mean - hm.HM_rf_test_M_auroc_mean   # Sup. Table 5 definition (RF vs RF)
T["loss_rf_only_prc"] = hm.MM_rf_test_M_auprc_mean - hm.HM_rf_test_M_auprc_mean
T["dROC_HH"] = rf.roc_auc_test_H - rf.roc_auc_test_H_PWM; T["dPRC_HH"] = rf.pr_auc_test_H - rf.pr_auc_test_H_PWM
T["mouse_columns_status"] = MOUSE_STATUS
T.to_csv(OUT + "/B5_per_TF_table.csv", float_format="%.6g")
if not (np.array_equal(T.dROC_HM, T.gain_HM_roc_2seed) and np.array_equal(T.dPRC_HM, T.gain_HM_prc_2seed)):
    raise SystemExit("CHECK FAILED V006: driver outcomes dROC_HM / dPRC_HM are not the two-seed mean gains")
chk = {"max_abs_dROC_HM_seed0_vs_2seed": float((s0_roc - T.gain_HM_roc_2seed).abs().max()),
       "max_abs_dPRC_HM_seed0_vs_2seed": float((s0_prc - T.gain_HM_prc_2seed).abs().max())}

# ---------------------------------------------------------------- Analysis 1: Spearman drivers x outcomes
DRV = [("IC_total_auROC", "Baseline monoPWM total IC (bits)"), ("motif_len_auROC", "Baseline motif length (positions)"),
       ("IC_mean_per_pos_auROC", "Baseline mean IC per position (bits)"), ("n_train_pos_H", "Human training positives"),
       ("n_PWMs_H", "Number of PWMs (mono + di)"), ("n_datasets_H", "Number of human datasets"),
       ("base_auROC_H", "Baseline auROC, human test"), ("base_auPRC_H", "Baseline auPRC, human test"),
       ("base_auROC_M", "Baseline auROC, mouse test"), ("base_auPRC_M", "Baseline auPRC, mouse test"),
       ("motif_sim_H_M", "Motif similarity human-mouse"), ("celltype_overlap_ord", "Cell-type overlap (shared tissue or lineage: yes=1, no=0)")]
OUTC = [("dROC_HM", "ΔauROC, human-trained on mouse test (mean of 2 seeds)"), ("dPRC_HM", "ΔauPRC, human-trained on mouse test (mean of 2 seeds)"),
        ("loss_roc", "Transfer loss auROC (M>M gain - H>M gain)"), ("loss_prc", "Transfer loss auPRC (M>M gain - H>M gain)")]
rows = []
for d, dl in DRV:
    for o, ol in OUTC:
        r, p = stats.spearmanr(T[d], T[o])
        rows.append(dict(driver=d, driver_label=dl, outcome=o, outcome_label=ol, n=int(T[[d, o]].dropna().shape[0]), rho=r, p=p))
C = pd.DataFrame(rows); C["q_BH"] = bh(C.p.values)   # BH over all 48 tests of the table
C["status"] = MOUSE_STATUS                           # every outcome is measured on the mouse test set
C.to_csv(OUT + "/B5_A_spearman_drivers.csv", index=False, float_format="%.6g")
# sensitivity: IC of the auPRC-selected baseline PWM against the auPRC outcomes
sens = []
for d in ["IC_total_auPRC", "motif_len_auPRC", "IC_mean_per_pos_auPRC"]:
    for o in ["dPRC_HM", "loss_prc"]:
        r, p = stats.spearmanr(T[d], T[o]); sens.append(dict(driver=d, outcome=o, rho=r, p=p))
pd.DataFrame(sens).assign(status=MOUSE_STATUS).to_csv(OUT + "/B5_A_sensitivity_auPRC_baseline_IC.csv", index=False, float_format="%.6g")


# multivariable: rank-transformed OLS, 4 predeclared predictors, standardised ranks
def rank_ols(y, X, names, rng):
    Y = stats.rankdata(y); Z = np.column_stack([stats.rankdata(X[:, j]) for j in range(X.shape[1])])
    Y = (Y - Y.mean()) / Y.std(); Z = (Z - Z.mean(0)) / Z.std(0)
    n, k = Z.shape; A = np.column_stack([np.ones(n), Z])
    beta = np.linalg.lstsq(A, Y, rcond=None)[0]; resid = Y - A @ beta
    s2 = resid @ resid / (n - k - 1); se = np.sqrt(np.diag(s2 * np.linalg.inv(A.T @ A)))
    tt = beta / se; pv = 2 * stats.t.sf(abs(tt), n - k - 1)
    r2 = 1 - resid @ resid / (Y @ Y); adj = 1 - (1 - r2) * (n - 1) / (n - k - 1)
    boots = []
    for _ in range(NBOOT):
        i = rng.integers(0, n, n)
        yy = stats.rankdata(Y[i]); zz = np.column_stack([stats.rankdata(Z[i, j]) for j in range(k)])
        yy = (yy - yy.mean()) / yy.std(); zz = (zz - zz.mean(0)) / (zz.std(0) + 1e-12)
        boots.append(np.linalg.lstsq(np.column_stack([np.ones(n), zz]), yy, rcond=None)[0][1:])
    lo, hi = np.percentile(boots, [2.5, 97.5], axis=0)
    out = [dict(predictor=names[j], beta_std_rank=beta[j + 1], se=se[j + 1], p=pv[j + 1], boot_ci_lo=lo[j], boot_ci_hi=hi[j]) for j in range(k)]
    return out, r2, adj


mv, mvsum = [], []
rng = np.random.default_rng(SEED)
for o, pred_roc in [("dROC_HM", "base_auROC_M"), ("loss_roc", "base_auROC_M"), ("dPRC_HM", "base_auPRC_M"), ("loss_prc", "base_auPRC_M")]:
    names = ["motif_sim_H_M", "n_train_pos_H", pred_roc, "IC_total_auROC"]
    out, r2, adj = rank_ols(T[o].values, T[names].values, names, rng)
    for r in out: r["outcome"] = o; mv.append(r)
    mvsum.append(dict(outcome=o, n=NTF, k=4, R2=r2, adj_R2=adj, status=MOUSE_STATUS))
    Xc = pd.DataFrame(T[names].rank()).corr()
    # variance inflation of each predictor
    Z = T[names].rank().values; Z = (Z - Z.mean(0)) / Z.std(0)
    for j, nm in enumerate(names):
        o2 = np.delete(Z, j, axis=1); A = np.column_stack([np.ones(len(Z)), o2]); b = np.linalg.lstsq(A, Z[:, j], rcond=None)[0]
        rr = 1 - ((Z[:, j] - A @ b) ** 2).sum() / (Z[:, j] ** 2).sum()
        for r in mv:
            if r["outcome"] == o and r["predictor"] == nm: r["VIF"] = 1 / (1 - rr)
MV = pd.DataFrame(mv)
MV["q_BH"] = bh(MV.p.values)
MV["status"] = MOUSE_STATUS
MV.to_csv(OUT + "/B5_A_multivariable_rank_ols.csv", index=False, float_format="%.6g")
pd.DataFrame(mvsum).to_csv(OUT + "/B5_A_multivariable_summary.csv", index=False, float_format="%.6g")

# ---------------------------------------------------------------- Analysis 2
H = pd.DataFrame(index=T.index)
H["dROC"] = T.dROC_HH; H["dPRC"] = T.dPRC_HH; H["bROC"] = T.base_auROC_H; H["bPRC"] = T.base_auPRC_H
H["head_roc"] = H.dROC / (1 - H.bROC); H["head_prc"] = H.dPRC / (1 - H.bPRC); H["rel_prc"] = H.dPRC / H.bPRC
H["rel_roc"] = H.dROC / H.bROC
# TFs with the CTCF-like ceiling: report
num = None
rep = []
def sp(name, x, y):
    r, p = stats.spearmanr(x, y); h = 1.96 / np.sqrt(len(x) - 3)   # 95 % CI by Fisher transformation (quoted in the text)
    rep.append(dict(quantity=name, rho=r, p=p, n=len(x), rho_ci95_lo_fisher=np.tanh(np.arctanh(r) - h), rho_ci95_hi_fisher=np.tanh(np.arctanh(r) + h)))
sp("dauROC vs baseline auROC", H.dROC, H.bROC)
sp("dauPRC vs baseline auPRC", H.dPRC, H.bPRC)
sp("dauROC/(1-auROC_base) vs baseline auROC", H.head_roc, H.bROC)
sp("dauPRC/(1-auPRC_base) vs baseline auPRC", H.head_prc, H.bPRC)
sp("dauPRC/auPRC_base (relative) vs baseline auPRC", H.rel_prc, H.bPRC)
sp("dauROC/auROC_base (relative) vs baseline auROC", H.rel_roc, H.bROC)
sp("dauPRC vs baseline auROC", H.dPRC, H.bROC)
sp("dauROC vs baseline auPRC", H.dROC, H.bPRC)
sp("baseline auROC vs baseline auPRC", H.bROC, H.bPRC)
sp("dauROC vs dauPRC", H.dROC, H.dPRC)
_out = sorted(H.index[H.dPRC > PRC_GAP])   # with the first-submission 36 TFs: ANDR, MAFK, P53
sp("dauPRC vs baseline auPRC, excluding the TFs with dauPRC > %g (%s)" % (PRC_GAP, ", ".join(_out)), H.drop(_out).dPRC, H.drop(_out).bPRC)
sp("dauPRC vs prevalence of positives (human test)", H.dPRC, T.prevalence_H_test)
sp("baseline auPRC vs prevalence of positives (human test)", H.bPRC, T.prevalence_H_test)
sp("auPRC gain / baseline auPRC vs prevalence", H.rel_prc, T.prevalence_H_test)
# lift over chance: baseline auPRC / prevalence, headroom against the chance-normalised scale
H["prevalence"] = T.prevalence_H_test
H["pr_norm_gain"] = H.dPRC / (1 - H.bPRC)   # same as head_prc
R = pd.DataFrame(rep); R["q_BH"] = bh(R.p.values)
R["status"] = "current: human test, %d TFs" % NTF
R.to_csv(OUT + "/B5_B_ceiling_spearman.csv", index=False, float_format="%.6g")
H.to_csv(OUT + "/B5_B_headroom_per_TF.csv", float_format="%.6g")
s_roc, s_prc = [stats.spearmanr(H.dROC, H.bROC), stats.spearmanr(H.dPRC, H.bPRC)]
if not PD.EXCLUDE_TFS:   # regression check: with all 36 TFs the first-submission values come back
    assert round(s_roc[0], 4) == -0.6515 and round(s_prc[0], 4) == 0.5243 and abs(s_roc[1] - 1.68e-5) < 1e-7 and abs(s_prc[1] - 0.00103) < 1e-5
print("headroom: rho dROC vs bROC %.4f (P %.3g), dPRC vs bPRC %.4f (P %.3g), n %d" % (s_roc[0], s_roc[1], s_prc[0], s_prc[1], NTF))
# the two cut-offs must still fall in empty bands of the gains (else the group definition would need revisiting)
_r = np.sort(H.dROC.values); _p = np.sort(H.dPRC.values)
print("dauROC just below / above %.3f: %.4f / %.4f; dauPRC just below / above %.2f: %.4f / %.4f" % (
    ROC_GAP, _r[_r <= ROC_GAP].max(), _r[_r > ROC_GAP].min(), PRC_GAP, _p[_p <= PRC_GAP].max(), _p[_p > PRC_GAP].min()))

# clouds in Fig. 3C (x = dauPRC, y = dauROC, human test set, RF mono+di vs best monoPWM)
Zs = np.column_stack([stats.zscore(H.dPRC), stats.zscore(H.dROC)])
lab = fcluster(linkage(Zs, "ward"), 2, "maxclust")
big = pd.Series(lab).value_counts().idxmax()
H["cloud_ward2"] = np.where(lab == big, "core", "outlying")
# simple gap-based definition for the vertical split: dauROC > 0.075 (largest empty band 0.060 -> 0.079 above the core)
H["cloud_roc_high"] = np.where(H.dROC > ROC_GAP, "dauROC>%g" % ROC_GAP, "rest")
H["cloud_prc_high"] = np.where(H.dPRC > PRC_GAP, "dauPRC>%g" % PRC_GAP, "rest")
H["quadrant"] = np.where(H.dROC > ROC_GAP, "top (dauROC>%g)" % ROC_GAP, np.where(H.dPRC > PRC_GAP, "right (dauPRC>%g)" % PRC_GAP, "core"))
H.to_csv(OUT + "/B5_B_clouds_per_TF.csv", float_format="%.6g")
VARS = [d for d, _ in DRV] + ["prevalence_H_test"]
V = T.join(H[["prevalence"]], rsuffix="_x")
sep = []
rng = np.random.default_rng(SEED)
for grp, col in [("ward2 outlying vs core", "cloud_ward2"), ("top: dauROC>%g vs rest" % ROC_GAP, "cloud_roc_high"), ("right: dauPRC>%g vs rest" % PRC_GAP, "cloud_prc_high")]:
    g = H[col]; pos = g.isin(["outlying", "dauROC>%g" % ROC_GAP, "dauPRC>%g" % PRC_GAP]).values
    for v in VARS:
        x = T[v].values.astype(float)
        u = stats.mannwhitneyu(x[pos], x[~pos], alternative="two-sided")
        auc = u.statistic / (pos.sum() * (~pos).sum())
        sep.append(dict(split=col, n_group=int(pos.sum()), n_rest=int((~pos).sum()), variable=v, median_group=np.median(x[pos]), median_rest=np.median(x[~pos]),
                        AUC_group_higher=auc, p_MWU=u.pvalue))
    # family: superclass
    for sc, n in T.loc[pos, "tfclass_superclass"].value_counts().items():
        sep.append(dict(split=col, n_group=int(pos.sum()), n_rest=int((~pos).sum()), variable="TFClass superclass %s (n in group / in rest)" % sc,
                        median_group=n, median_rest=int((T.loc[~pos, "tfclass_superclass"] == sc).sum()), AUC_group_higher=np.nan, p_MWU=np.nan))
S = pd.DataFrame(sep)
S["q_BH"] = np.nan
m = S.p_MWU.notna(); S.loc[m, "q_BH"] = bh(S.loc[m, "p_MWU"].values)
S["status"] = np.where(S.variable.isin(MOUSE_VARS), MOUSE_STATUS, "current: human")
S.to_csv(OUT + "/B5_B_cloud_separation.csv", index=False, float_format="%.6g")
# membership list
with open(OUT + "/B5_B_cloud_members.txt", "w") as f:
    for g, d in H.groupby("quadrant"): f.write("%s (n=%d): %s\n" % (g, len(d), ", ".join(d.sort_values("dROC", ascending=False).index)))
    for g, d in H.groupby("cloud_ward2"): f.write("ward2 %s (n=%d): %s\n" % (g, len(d), ", ".join(d.index)))
pd.Series(chk).to_csv(OUT + "/B5_checks.csv", header=["value"])

# ---------------------------------------------------------------- provenance
try: git = subprocess.run(["git", "-C", MS, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip() or "n/a"
except Exception: git = "n/a"
with open(OUT + "/_provenance.yaml", "w") as f:
    f.write("p53drop: {exclude_tfs: [%s], n_TF: %d, mouse_inputs: '%s'}\n" % (", ".join(PD.EXCLUDE_TFS), NTF, MOUSE_STATUS))
    f.write("script: %s\nrun: %s\npython: %s\nnumpy: %s\npandas: %s\nscipy: %s\nrandom_seed: %d\nbootstrap_n: %d\nbackground_per_base: %.2f\n" %
            (os.path.abspath(__file__), datetime.datetime.now().isoformat(timespec="seconds"), platform.python_version(), np.__version__, pd.__version__, scipy.__version__, SEED, NBOOT, BG))
    f.write("ArChIPelago_repo_git_head: %s\ninputs:\n" % git)
    for k, v in INP.items(): f.write("  %s: {path: %s, md5: %s}\n" % (k, v, md5(v)))
    f.write("  PWM_dir: {path: %s, note: 'monoPWM files read for the 72 baseline PWMs (%d TFs x auROC/auPRC selection)'}\n" % (PWM_DIR, NTF))
    f.write("outputs:\n")
    for _f in OUTFILES[:-1]: f.write("  - {file: %s, md5: %s}\n" % (_f, md5(OUT + "/" + _f)))
    f.write("pwm_to_pfm: 'files hold ln(p/0.25); p = 0.25*exp(w); columns sum to 1 (asserted); IC_i = 2 + sum p log2 p'\n")
print("done", chk)
