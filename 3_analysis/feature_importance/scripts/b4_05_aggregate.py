"""Step 5 of the feature-importance analysis (local, miniconda python3): aggregate the per-TF server outputs into the result tables.
Input : {gini,perm,metrics,subtypes,usage,ablation}/<TF>_*.csv under the folder B4_IN (and <TF>_sim.npy under B4_SIM)
Output: 3_analysis/feature_importance/tables/B4_*.csv
TFs in p53drop.EXCLUDE_TFS are not read (35 TFs; ablation 11 TFs).
usage: B4_IN=<folder> B4_SIM=<folder> python3 3_analysis/feature_importance/scripts/b4_05_aggregate.py
"""
import os, sys, glob, csv
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
B = os.path.normpath(os.path.join(HERE, ".."))
T = os.path.join(B, "tables"); os.makedirs(T, exist_ok=True)
# the per-TF inputs are read from B4_IN (the mirror of the server results); no default, so a
# run cannot silently aggregate the first-submission per-TF files
B_IN = os.environ.get("B4_IN") or sys.exit("CHECK FAILED b4_05_aggregate: set B4_IN to the per-TF input folder")
# the PWM similarity matrices (<TF>_sim.npy) are not in the mirror (no .npy); they depend on the PWMs and the
# Gini order only, so they are read from B4_SIM after checking that its sibling gini files equal the B4_IN ones
B_SIM = os.environ.get("B4_SIM") or sys.exit("CHECK FAILED b4_05_aggregate: set B4_SIM to the folder of the <TF>_sim.npy files")
_gs = sorted(os.listdir(os.path.join(B_IN, "gini")))
_bad = [f for f in _gs if open(os.path.join(B_IN, "gini", f), "rb").read() != open(os.path.join(B_SIM, "..", "gini", f), "rb").read()]
if _bad or sorted(os.listdir(os.path.join(B_SIM, "..", "gini"))) != _gs:
    sys.exit("CHECK FAILED b4_05_aggregate: gini files of B4_IN and of the B4_SIM tree differ (%s): the similarity cache is stale" % _bad[:5])
print("sim cache: %d gini files byte-identical between B4_IN and B4_SIM" % len(_gs))
META = os.path.join(B, "../package/analysis/inputs/GTRD_metadata.csv")
sys.path.insert(0, os.path.join(B, "..", "integration/"))
import p53drop as PD
for f in sorted(glob.glob(os.path.join(T, "B4_*.csv"))):
    if "_backup_" not in f:
        PD.backup(f)
INPUTS = []


def cat(sub, pat):
    fs = [f for f in sorted(glob.glob(os.path.join(B_IN, sub, pat))) if os.path.getsize(f) > 2]   # TFs without a secondary subtype write empty usage files
    fs = [f for f in fs if os.path.basename(f).split("_")[0] not in PD.EXCLUDE_TFS]
    INPUTS.extend(fs)
    return pd.concat([pd.read_csv(f) for f in fs], ignore_index=True), len(fs)


G, n = cat("gini", "*_gini.csv"); print("gini files", n)
PD.check_tfs(G.TF, "B4 gini")
M, _ = cat("metrics", "*_metrics.csv")
PD.check_tfs(M.TF, "B4 metrics")
# ---------- 1. concentration
rows = []
for tf, g in G.groupby("TF"):
    g = g.sort_values("rank"); p = g.gini.values / g.gini.sum()
    c = np.cumsum(p)
    rows.append(dict(TF=tf, n_features=len(g), n_mono=int((g.kind == "mono").sum()), n_di=int((g.kind == "di").sum()),
                     n_for_50=int(np.argmax(c >= 0.5) + 1), n_for_80=int(np.argmax(c >= 0.8) + 1),
                     top1_share=p[0], top3_share=p[:3].sum(), top1_feature=g.feature.iloc[0], top1_kind=g.kind.iloc[0],
                     di_share=float(p[(g.kind == "di").values].sum()),
                     effective_n=float(1 / (p ** 2).sum())))
C = pd.DataFrame(rows).merge(M[["TF", "model_source", "n_test", "n_pos", "auroc_sav", "auprc_sav"]], on="TF")
C["frac_features_for_50"] = C.n_for_50 / C.n_features; C["frac_features_for_80"] = C.n_for_80 / C.n_features
C.to_csv(os.path.join(T, "B4_importance_concentration.csv"), index=False, float_format="%.5g")
G.to_csv(os.path.join(T, "B4_gini_all_features.csv"), index=False, float_format="%.6g")
# ---------- 2. PWM origin
O = G.groupby("TF").agg(n_features=("feature", "size"), n_header_tf_is_target=("header_tf_is_target", "sum"),
                        n_file_matches_header=("file_matches_header", "sum"), n_datasets=("dataset", "nunique"),
                        n_variant_M=("variant", lambda v: int((v == "M").sum())), n_variant_S=("variant", lambda v: int((v == "S").sum()))).reset_index()
O.to_csv(os.path.join(T, "B4_pwm_origin_check.csv"), index=False)
# ---------- 3. permutation
P, _ = cat("perm", "*_perm.csv")
PS = P.groupby(["TF", "feature"]).agg(d_auroc_mean=("d_auroc", "mean"), d_auroc_sd=("d_auroc", "std"),
                                      d_auprc_mean=("d_auprc", "mean"), d_auprc_sd=("d_auprc", "std"), n_repeats=("repeat", "size")).reset_index()
PS = PS.merge(G[["TF", "feature", "rank", "gini", "kind"]], on=["TF", "feature"]).sort_values(["TF", "rank"])
PS.to_csv(os.path.join(T, "B4_permutation_top20.csv"), index=False, float_format="%.6g")
rows = []
for tf, g in PS.groupby("TF"):
    rho = spearmanr(g.gini, g.d_auroc_mean)[0] if len(g) > 3 else np.nan
    rows.append(dict(TF=tf, n_tested=len(g), max_d_auroc=g.d_auroc_mean.max(), max_d_auprc=g.d_auprc_mean.max(),
                     feature_max_d_auroc=g.loc[g.d_auroc_mean.idxmax(), "feature"], n_d_auroc_gt_0p001=int((g.d_auroc_mean > 0.001).sum()),
                     n_d_auroc_gt_0p005=int((g.d_auroc_mean > 0.005).sum()), n_d_auprc_gt_0p005=int((g.d_auprc_mean > 0.005).sum()),
                     sum_d_auroc=g.d_auroc_mean.sum(), spearman_gini_vs_perm_auroc=rho,
                     perm_rank_of_top_gini=int(g.d_auroc_mean.rank(ascending=False)[g["rank"].idxmin()])))
PT = pd.DataFrame(rows); PT.to_csv(os.path.join(T, "B4_permutation_summary_per_TF.csv"), index=False, float_format="%.5g")
# ---------- 4. subtypes
S, _ = cat("subtypes", "*_subtypes.csv"); SS, _ = cat("subtypes", "*_sim_summary.csv")
meta = {}
with open(META, newline="") as fh:
    for r in csv.reader(fh):
        if len(r) > 6 and r[6].startswith("PEAKS"):
            meta[r[6]] = r[5]
S["cell_line"] = S.dataset.map(meta).fillna("not in metadata")
S.to_csv(os.path.join(T, "B4_feature_subtypes.csv"), index=False, float_format="%.6g")
rows = []
for tf, s in S.groupby("TF"):
    for thr in (70, 80, 90):
        c = "subtype_%d" % thr
        for st, g in s.groupby(c):
            rows.append(dict(TF=tf, sim_threshold=thr / 100, subtype=st, n_features=len(g), n_mono=int((g.kind == "mono").sum()),
                             n_di=int((g.kind == "di").sum()), gini_share=g.gini.sum(), n_datasets=g.dataset.nunique(),
                             n_cell_lines=g[g.cell_line != "not in metadata"].cell_line.nunique(),
                             cell_lines=";".join(sorted(set(g.cell_line)))))
SI = pd.DataFrame(rows); SI.to_csv(os.path.join(T, "B4_subtype_importance.csv"), index=False, float_format="%.5g")
SS.to_csv(os.path.join(T, "B4_similarity_summary.csv"), index=False, float_format="%.5g")
rows = []
for tf, g in SI.groupby("TF"):
    r = dict(TF=tf)
    for thr in (0.7, 0.8, 0.9):
        h = g[g.sim_threshold == thr].sort_values("subtype")
        r["n_subtypes_%d" % round(thr * 100)] = len(h)
        r["n_ge10pct_%d" % round(thr * 100)] = int((h.gini_share >= 0.10).sum())
        r["primary_share_%d" % round(thr * 100)] = h.gini_share.iloc[0]
        r["secondary_share_%d" % round(thr * 100)] = h.gini_share.iloc[1] if len(h) > 1 else 0.0
    h = g[g.sim_threshold == 0.8].sort_values("subtype")
    cl1 = set(h.cell_lines.iloc[0].split(";")); cl2 = set(h.cell_lines.iloc[1].split(";")) if len(h) > 1 else set()
    r["cell_lines_sec_not_in_primary_80"] = ";".join(sorted(cl2 - cl1))
    rows.append(r)
SUM = pd.DataFrame(rows).merge(SS[["TF", "median_pairwise_mono", "median_pairwise_all", "min_pairwise_all", "median_sim_to_top_feature"]], on="TF")
SUM.to_csv(os.path.join(T, "B4_subtype_summary_per_TF.csv"), index=False, float_format="%.5g")
# ---------- 5. usage
BP, _ = cat("usage", "*_best_pwm.csv"); GR, _ = cat("usage", "*_groups.csv"); GP, _ = cat("usage", "*_grouped_perm.csv")
BP.to_csv(os.path.join(T, "B4_activation_best_pwm.csv"), index=False, float_format="%.5g")
GR.to_csv(os.path.join(T, "B4_activation_groups.csv"), index=False, float_format="%.5g")
GPS = GP.groupby(["TF", "subtype", "n_features", "gini_share"]).agg(d_auroc_mean=("d_auroc", "mean"), d_auroc_sd=("d_auroc", "std"),
                                                                   d_auprc_mean=("d_auprc", "mean"), d_auprc_sd=("d_auprc", "std")).reset_index()
GPS.to_csv(os.path.join(T, "B4_grouped_permutation.csv"), index=False, float_format="%.5g")
# ---------- 6. ablation
A, _ = cat("ablation", "*_ablation.csv")
if set(A.TF) & set(PD.EXCLUDE_TFS):
    raise SystemExit("CHECK FAILED ablation contains %s" % sorted(set(A.TF) & set(PD.EXCLUDE_TFS)))
print("ablation TFs:", A.TF.nunique())
rows = []
for tf, a in A.groupby("TF"):
    full = a[a.cond == "full"].set_index("seed")
    for cond, c in a.groupby("cond"):
        c = c.set_index("seed")
        dr = (c.auroc - full.auroc.reindex(c.index)); dp = (c.auprc - full.auprc.reindex(c.index))
        rows.append(dict(TF=tf, sim_threshold=c.sim_threshold.iloc[0], cond=cond, n_fits=len(c), n_removed=c.n_removed.iloc[0],
                         n_prim=c.n_prim.iloc[0], n_sec=c.n_sec.iloc[0], auroc_mean=c.auroc.mean(), auprc_mean=c.auprc.mean(),
                         d_auroc_mean=dr.mean(), d_auroc_sd=dr.std(), d_auprc_mean=dp.mean(), d_auprc_sd=dp.std(),
                         d_auroc_min=dr.min(), d_auroc_max=dr.max(), d_auprc_min=dp.min(), d_auprc_max=dp.max()))
AS = pd.DataFrame(rows); AS.to_csv(os.path.join(T, "B4_ablation_summary.csv"), index=False, float_format="%.5g")
A.to_csv(os.path.join(T, "B4_ablation_all_fits.csv"), index=False, float_format="%.6g")
print("tables:", sorted(os.listdir(T)))

# ---------- 7. similarity of the top secondary subtype to the primary subtype (similarity cut 0.8; needs cache_sim/<TF>_sim.npy)
rows = []
for tf, s in S.groupby("TF"):
    sim = np.load(os.path.join(B_SIM, "%s_sim.npy" % tf))   # features in gini-rank order == row order of subtypes csv
    s = s.reset_index(drop=True)
    for thr in (80, 90):
        st = s["subtype_%d" % thr].values
        if (st == 2).sum() == 0:
            continue
        i1, i2 = np.flatnonzero(st == 1), np.flatnonzero(st == 2)
        sub = sim[np.ix_(i2, i1)]
        g2 = s[st == 2]
        rows.append(dict(TF=tf, sim_threshold=thr / 100, n_primary=len(i1), n_secondary=len(i2), gini_share_secondary=g2.gini.sum(),
                         max_sim_secondary_to_primary=sub.max(), median_of_best_match=float(np.median(sub.max(axis=1))),
                         n_secondary_mono=int((g2.kind == "mono").sum()), n_secondary_di=int((g2.kind == "di").sum()),
                         secondary_datasets=";".join(sorted(set(g2.dataset))), secondary_cell_lines=";".join(sorted(set(g2.cell_line))),
                         primary_cell_lines=";".join(sorted(set(s[st == 1].cell_line))),
                         secondary_cell_lines_not_in_primary=";".join(sorted(set(g2.cell_line) - set(s[st == 1].cell_line)))))
pd.DataFrame(rows).to_csv(os.path.join(T, "B4_secondary_vs_primary_subtype.csv"), index=False, float_format="%.4g")
print("secondary table written")
OUTS = sorted(f for f in glob.glob(os.path.join(T, "B4_*.csv")) if "_backup_" not in f)
for f in OUTS:
    d = pd.read_csv(f)
    if "TF" in d.columns and set(d.TF) & set(PD.EXCLUDE_TFS):
        raise SystemExit("CHECK FAILED %s contains %s" % (f, sorted(set(d.TF) & set(PD.EXCLUDE_TFS))))
PD.write_provenance(B, "3_analysis/feature_importance/scripts/b4_05_aggregate.py", OUTS,
                    sorted(set(INPUTS)) + [os.path.join(B_SIM, "%s_sim.npy" % t) for t in sorted(S.TF.unique())],
                    note="per-TF files of excluded TFs skipped; %d output tables" % len(OUTS))
