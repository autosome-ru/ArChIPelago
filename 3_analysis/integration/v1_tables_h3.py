"""Source tables for the mono vs mono+di comparison, the GC comparison and the positive-count ties, on the final data.
  mono vs mono+di   mono_vs_monodi_h3.csv: RF on monoPWMs + diPWMs vs RF on monoPWMs only (seed 0, human-trained), human and mouse test sets,
         auROC and auPRC: median difference, TFs improved, two-sided Wilcoxon signed-rank P (pkg results_table.csv)
  GC     gc_per_tf_h3.csv: median GC of positives and negatives per TF, species and split (v3 negatives, gc_check_v3.csv),
         difference in percentage points; gc_vs_gain_h3.csv: Spearman of |dGC| (test split) with the gain of the RF over the
         best monoPWM on the same test set (human: pkg results table; mouse: mouse-trained, hm_mm MM_*, mean of 2 seeds)
  ties   positive_count_ties_h3.csv: TFs that share a positive count (per species and split; Sup. Table 2), with their pool sizes
         (per_tf_positive_counts.csv of the human and mouse data builds) and whether the count is the 10,000 cap
usage: ARCHI_STAGE=H5 python3 3_analysis/integration/v1_tables_h3.py
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w2  # noqa: E402
import h3_inputs as H  # noqa: E402

PD, R2, B13 = w2.PD, w2.R2, w2.B13
CAP = 10000
IN = {"rt": os.path.join(B13, "../package", "results_table.csv"), "gc": os.path.join(R2, H.IN["gc_v3"]), "hm": os.path.join(R2, H.IN["hm_mm"]),
      "t2": os.path.join(R2, "outputs/supplementary_tables", "Sup_Table_2_PWMs_and_datasets.xlsx"),
      "cH": os.path.join(R2, H.IN["pos_H"]), "cM": os.path.join(R2, H.IN["pos_M"])}
OUT = {k: os.path.join(B13, f) for k, f in (("mvd", "mono_vs_monodi_h3.csv"), ("gc", "gc_per_tf_h3.csv"), ("gcs", "gc_vs_gain_h3.csv"),
                                            ("ties", "positive_count_ties_h3.csv"))}
for f in OUT.values():
    w2.backup(f)

# mono vs mono+di
t = pd.read_csv(IN["rt"], sep="\t")
rf = t[t.Model == "RandomForestClassifier"]
if rf.duplicated(["TF_name", "PWM"]).any():
    raise SystemExit("CHECK FAILED results_table: duplicate RF rows")
md, mo = rf[rf.PWM == "mono+di"].set_index("TF_name").sort_index(), rf[rf.PWM == "mono"].set_index("TF_name")
PD.check_tfs(md.index, "results_table RF mono+di")
mo = mo.loc[md.index]
rows = []
for ts, lab in (("H", "human"), ("M", "mouse")):
    for m, ml in (("roc_auc", "auROC"), ("pr_auc", "auPRC")):
        c = "%s_test_%s" % (m, ts)
        d = md[c] - mo[c]
        rows.append(dict(test_set=lab, metric=ml, n_TF=len(d), median_monodi=md[c].median(), median_mono=mo[c].median(), median_diff=d.median(),
                         n_improved=int((d > 0).sum()), n_worse=int((d < 0).sum()), wilcoxon_p_two_sided=wilcoxon(d).pvalue,
                         worse_TFs=",".join(d.index[d < 0])))
pd.DataFrame(rows).to_csv(OUT["mvd"], index=False)

# GC
g = pd.read_csv(IN["gc"])
g = g[g.pool.isin([H.GC_POOL["H"], H.GC_POOL["M"]])].copy()   # pool names of the stage
# logical split: the mouse files are swapped on disk, the chromosome list names the split
g["split"] = np.where(g.chromosomes.isin(["chr1/8/21", "chr1/8/19"]), "test", "training")
g["species"] = g.species.str.lower()
g["gc_pos_pct"], g["gc_neg_pct"], g["gc_diff_points"] = 100 * g.median_gc_pos, 100 * g.median_gc_neg, 100 * g.median_gc_diff
for (sp, sl), x in g.groupby(["species", "split"]):
    PD.check_tfs(x.TF, "GC v3 %s %s" % (sp, sl))
    if not np.allclose(x.median_gc_pos - x.median_gc_neg, x.median_gc_diff, atol=1e-4):
        raise SystemExit("CHECK FAILED GC %s %s: diff != pos - neg" % (sp, sl))
gt = g[["TF", "species", "split", "chromosomes", "n_pos", "n_neg", "gc_pos_pct", "gc_neg_pct", "gc_diff_points", "ks_D", "ks_p"]].sort_values(["species", "split", "TF"])
gt.to_csv(OUT["gc"], index=False)
hm = pd.read_csv(IN["hm"]).set_index("TF")
gain = {"human": {"auROC": md.roc_auc_test_H - md.roc_auc_test_H_PWM, "auPRC": md.pr_auc_test_H - md.pr_auc_test_H_PWM},
        "mouse": {"auROC": hm.MM_rf_test_M_auroc_mean - hm.MM_broc_test_M_auroc, "auPRC": hm.MM_rf_test_M_auprc_mean - hm.MM_bprc_test_M_auprc}}
srows = []
for sp in ("human", "mouse"):
    x = gt[(gt.species == sp) & (gt.split == "test")].set_index("TF")
    for ml, gv in gain[sp].items():
        r = spearmanr(x.gc_diff_points.abs(), gv.reindex(x.index))
        srows.append(dict(species=sp, test_set="%s test" % sp, model="human-trained RF" if sp == "human" else "mouse-trained RF (mean of 2 seeds)",
                          metric=ml, n_TF=len(x), spearman_rho=r.correlation, p=r.pvalue))
pd.DataFrame(srows).to_csv(OUT["gcs"], index=False)

# positive-count ties
t2 = pd.read_excel(IN["t2"], sheet_name="per_TF").set_index("TF")
pools = {"Human": pd.read_csv(IN["cH"]).set_index("TF"), "Mouse": pd.read_csv(IN["cM"]).set_index("TF")}
trows = []
for c in [c for c in t2.columns if "positives (" in c and "GC" not in c]:
    sp = c.split()[0]
    v = t2[c]
    for n, grp in v[v.duplicated(keep=False)].groupby(v):
        trows.append(dict(species=sp.lower(), column=c, count=int(n), TFs=",".join(sorted(grp.index)), at_cap=int(n) == CAP,
                          pool_after_mask=",".join("%s %d" % (tf, pools[sp].loc[tf, "after_mask"]) for tf in sorted(grp.index))))
pd.DataFrame(trows).to_csv(OUT["ties"], index=False)
for k, f in OUT.items():
    print(k, os.path.basename(f)); print(pd.read_csv(f).to_string(index=False) if k != "gc" else pd.read_csv(f).groupby(["species", "split"]).gc_diff_points.describe().round(2))
w2.provenance(B13, "3_analysis/integration/v1_tables_h3.py", list(OUT.values()), list(IN.values()),
              params="{cap: %d, gc_threshold_points: 5, wilcoxon: two-sided signed-rank}" % CAP)
