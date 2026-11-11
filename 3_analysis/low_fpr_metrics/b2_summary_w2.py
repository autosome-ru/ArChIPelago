"""Summary tables from B2_metrics_wide.csv / B2_per_chromosome.csv (runs locally, py3.9)."""
import os, sys
import numpy as np, pandas as pd, warnings
from scipy.stats import wilcoxon
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "3_analysis/integration"))
os.chdir(os.path.dirname(os.path.abspath(__file__)))  # inputs and outputs in this folder (3_analysis/low_fpr_metrics)
import p53drop as PD
# inputs are the final-run tables copied by 3_analysis/integration/make_b2_h3.py (HH, HM, MM, 35 TFs each).
STATUS = {"HH": "current: H3a", "HM": "current: H3a", "MM": "current: H3a"}
for f in ("B2_summary.csv", "B2_per_chromosome_summary.csv"):
    PD.backup(f)
w = pd.read_csv("B2_metrics_wide.csv")
w = w[~((w.scenario == "HH") & w.TF.isin(PD.EXCLUDE_TFS))]
for _sc in ("HH", "HM", "MM"):
    PD.check_tfs(w[w.scenario == _sc].TF, "B2_metrics_wide " + _sc)
LOWER_BETTER = lambda m: m.startswith("fpr_at_tpr") or m.startswith("fp_per_kb_at_tpr")
PRTYPE = lambda m: m == "auprc" or m.startswith("prec_at_tpr") or m.startswith("auprc_part")
skip = {"n_pos", "n_neg"} | {c for c in w.columns if c.startswith("fp_per_kb_at_fpr")}
metrics = [c for c in w.columns if c not in ("TF", "scenario", "model") and c not in skip]

def wp(d):
    d = d[~np.isnan(d)]
    if len(d) == 0 or np.all(d == 0): return np.nan
    return wilcoxon(d, zero_method="wilcox").pvalue
def bh(p):
    p = np.asarray(p, float); q = np.full(len(p), np.nan); ok = ~np.isnan(p)
    pp = p[ok]; o = np.argsort(pp); r = pp[o] * len(pp) / (np.arange(len(pp)) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]; out = np.empty(len(pp)); out[o] = np.minimum(r, 1); q[ok] = out
    return q
rows = []
for sc in ("HH", "HM", "MM"):
    for seed in ("RF_s0", "RF_s1"):
        for m in metrics:
            base = "mono_byprc" if PRTYPE(m) else "mono_byroc"
            a = w[(w.scenario == sc) & (w.model == seed)].set_index("TF")[m]
            b = w[(w.scenario == sc) & (w.model == base)].set_index("TF")[m].reindex(a.index)
            d = a - b
            better = (d < 0) if LOWER_BETTER(m) else (d > 0)
            q = lambda x, p: float(np.nanpercentile(x, p))
            rows.append(dict(scenario=sc, rf_seed=seed[-1], metric=m, baseline=base, n_TF=int(a.notna().sum()),
                rf_median=a.median(), rf_q25=q(a, 25), rf_q75=q(a, 75), base_median=b.median(), base_q25=q(b, 25), base_q75=q(b, 75),
                diff_median=d.median(), diff_q25=q(d, 25), diff_q75=q(d, 75),
                n_improved=int(better.sum()), n_equal=int((d == 0).sum()), n_worse=int((~better & (d != 0)).sum()),
                wilcoxon_p=wp(d.values)))
s = pd.DataFrame(rows)
s["wilcoxon_q_BH"] = s.groupby(["scenario", "rf_seed"])["wilcoxon_p"].transform(lambda x: bh(x.values))
s["status"] = s.scenario.map(STATUS)
if not (s.n_TF == len(PD.reference_tfs())).all():
    raise SystemExit("CHECK FAILED B2_summary HH n_TF %s" % sorted(s[s.scenario == "HH"].n_TF.unique()))
s.to_csv("B2_summary.csv", index=False)
# per chromosome
c = pd.read_csv("B2_per_chromosome.csv")
c = c[~((c.scenario == "HH") & c.TF.isin(PD.EXCLUDE_TFS))]
for _sc in ("HH", "HM", "MM"):
    PD.check_tfs(c[c.scenario == _sc].TF, "B2_per_chromosome " + _sc)
rows = []
for sc in ("HH", "HM", "MM"):
    for ch in sorted(c[c.scenario == sc].chrom.unique(), key=lambda x: int(x[3:])):
        for m, base in [x for x in (("auroc", "mono_byroc"), ("auprc", "mono_byprc"), ("tpr_at_fpr0.01", "mono_byroc"), ("pauc_std_fpr0.01", "mono_byroc"))
                        if x[0] in c.columns]:   # pauc_* columns are skipped when absent
            a = c[(c.scenario == sc) & (c.chrom == ch) & (c.model == "RF_s0")].set_index("TF")[m]
            b = c[(c.scenario == sc) & (c.chrom == ch) & (c.model == base)].set_index("TF")[m].reindex(a.index)
            d = a - b
            rows.append(dict(scenario=sc, chrom=ch, metric=m, baseline=base, n_TF=len(a), rf_median=a.median(), base_median=b.median(),
                             diff_median=d.median(), diff_q25=np.percentile(d, 25), diff_q75=np.percentile(d, 75), n_improved=int((d > 0).sum()),
                             wilcoxon_p=wp(d.values), median_n_pos=c[(c.scenario == sc) & (c.chrom == ch) & (c.model == "RF_s0")].n_pos.median()))
pc = pd.DataFrame(rows); pc["wilcoxon_q_BH"] = bh(pc.wilcoxon_p.values)   # BH over all rows as before
pc["status"] = pc.scenario.map(STATUS)
pc.to_csv("B2_per_chromosome_summary.csv", index=False)
PD.write_provenance(os.path.dirname(os.path.abspath(__file__)), "3_analysis/low_fpr_metrics/b2_summary_w2.py (copy of B2_metrics/b2_summary.py)",
    [os.path.abspath(f) for f in ("B2_summary.csv", "B2_per_chromosome_summary.csv")],
    [os.path.abspath(f) for f in ("B2_metrics_wide.csv", "B2_per_chromosome.csv")],
    params="{baselines: mono_byroc for ROC-type, mono_byprc for PR-type; Wilcoxon two-sided, zero_method wilcox; BH per scenario and seed; per chromosome BH over all rows}")
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
cols = ["metric", "n_TF", "rf_median", "base_median", "diff_median", "n_improved", "n_equal", "wilcoxon_p"]
for sc in ("HH", "HM", "MM"):
    print(sc); print(s[(s.scenario == sc) & (s.rf_seed == "0")][cols].round(4).to_string(index=False))
print(pc.round(4).to_string(index=False))
