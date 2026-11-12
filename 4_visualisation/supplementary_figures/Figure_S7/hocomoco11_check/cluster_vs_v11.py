"""Agreement of the clustering subtypes used for the ablation (Methods 2.10: cut 0.8, but 0.9 for the TFs ablated at 0.9, i.e. MAFK
and SPI1; the cut per TF is read from B4_ablation_summary.csv) with the nearest HOCOMOCO v11 subtype, for the
11 ablated TFs. Every feature PWM of the TF's final model (diPWMs as their mononucleotide marginal, as
for Fig. S7) is matched to the TF's own v11 HUMAN mono models with pcc_max (subtype_vs_hocomoco11.py); its v11 subtype is the
index of the nearest own model. Agreement per TF = Rand index of the two partitions (share of PWM pairs on which both agree:
same cluster and same v11 subtype, or different and different); also the adjusted Rand index, the contingency table and the
number of PWMs whose two nearest own subtypes are within TIE in r, and 'informative' (both partitions have at least 2 groups, so
that the ARI is defined). Summary rows: median Rand over the 11 TFs, median ARI over the informative TFs, pooled pairs.
Output: cluster_vs_v11.csv (one row per TF, then "all"), cluster_vs_v11_per_PWM.csv, _provenance_h3.yaml (this folder).
usage: ARCHI_STAGE=H5 python3 4_visualisation/figures/FigS7/hocomoco11_check/cluster_vs_v11.py
"""
import importlib.util
import os
import sys
from itertools import combinations

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

HERE = os.path.dirname(os.path.abspath(__file__))
S7D = os.path.dirname(HERE)
R2 = os.path.dirname(os.path.dirname(os.path.dirname(S7D)))
_s = importlib.util.spec_from_file_location("s7h3", os.path.join(S7D, "build_Figure_S7_source_data_h3.py"))
S7 = importlib.util.module_from_spec(_s)
_s.loader.exec_module(S7)
P, V11, w2 = S7.P, S7.V11, S7.w2
TIE = S7.TIE
IN = {"sub": os.path.join(R2, "3_analysis/feature_importance", "tables", "B4_feature_subtypes.csv"),
      "abl": os.path.join(R2, "3_analysis/feature_importance", "tables", "B4_ablation_summary.csv"), "v11": V11.PCM}
OUT = {"tf": os.path.join(HERE, "cluster_vs_v11.csv"), "pwm": os.path.join(HERE, "cluster_vs_v11_per_PWM.csv")}


def rand(a, b):
    pairs = list(combinations(range(len(a)), 2))
    return sum((a[i] == a[j]) == (b[i] == b[j]) for i, j in pairs), len(pairs)


def main():
    for f in OUT.values():
        w2.backup(f)
    ab = pd.read_csv(IN["abl"])
    cut = ab.groupby("TF").sim_threshold.unique()
    if len(cut) != 11 or (cut.str.len() != 1).any():
        raise SystemExit("CHECK FAILED %d ablated TFs (expected 11), or a TF with several cuts" % len(cut))
    cut = cut.str[0]
    tfs = sorted(cut.index)
    sub = pd.read_csv(IN["sub"])
    z11 = {k: P.zrows(m) for k, m in V11.read_pcms(IN["v11"]).items()}
    per, rows = [], []
    for tf in tfs:
        x = sub[sub.TF == tf]
        col = "subtype_%d" % round(100 * cut[tf])
        g = S7.importances(tf)
        if set(x.feature) != set(g.feature):
            raise SystemExit("CHECK FAILED %s: subtype table features differ from the H3a Gini table" % tf)
        own = {k: z for k, z in z11.items() if k.startswith(tf + "_HUMAN.")}
        if not own:
            raise SystemExit("CHECK FAILED no HOCOMOCO v11 model for %s" % tf)
        for _, r in x.iterrows():
            h, w = P.read_pwm(P.pwm_path(tf, r.feature))
            za = P.zrows(P.to_prob(r.feature.split("_")[0], w))
            rr = {k: V11.pcc_max(za, z, V11.MIN_OVERLAP)[0] for k, z in own.items()}
            s = sorted(rr.values(), reverse=True)
            bo = max(rr, key=rr.get)
            per.append(dict(TF=tf, feature=r.feature, kind=r.kind, dataset=r.dataset, gini=r.gini, cut=cut[tf], cluster=int(r[col]), nearest_own_v11=bo,
                            v11_subtype=int(bo.split(".")[-2]), r_nearest=round(rr[bo], 3),
                            r_second=round(s[1], 3) if len(s) > 1 else np.nan, tie=len(s) > 1 and s[0] - s[1] < TIE))
        d = pd.DataFrame([p for p in per if p["TF"] == tf])
        a, n = rand(d.cluster.values, d.v11_subtype.values)
        ct = pd.crosstab(d.cluster, d.v11_subtype)
        inf = d.cluster.nunique() >= 2 and d.v11_subtype.nunique() >= 2
        rows.append(dict(TF=tf, cut=cut[tf], n_PWMs=len(d), n_clusters=d.cluster.nunique(), n_v11_models=len(own),
                         n_v11_subtypes_nearest=d.v11_subtype.nunique(), n_tie=int(d.tie.sum()), agreeing_pairs=a, pairs=n,
                         rand_index=a / n, adjusted_rand=adjusted_rand_score(d.cluster, d.v11_subtype), informative=inf,
                         gini_share_majority_v11_of_primary=float(d[d.cluster == 1].groupby("v11_subtype").gini.sum().max() / d[d.cluster == 1].gini.sum()),
                         contingency=" ; ".join("cluster %d: %s" % (c, ", ".join("v11 %d x%d" % (k, v) for k, v in ct.loc[c].items() if v))
                                                for c in ct.index)))
    t = pd.DataFrame(rows)
    inf = t[t.informative]
    t = pd.concat([t, pd.DataFrame([
        dict(TF="median, all TFs", n_PWMs=t.n_PWMs.sum(), rand_index=t.rand_index.median()),   # ARI only over informative TFs
        dict(TF="median, informative TFs (%d)" % len(inf), n_PWMs=inf.n_PWMs.sum(), rand_index=inf.rand_index.median(),
             adjusted_rand=inf.adjusted_rand.median()),
        dict(TF="pooled pairs, all TFs", n_PWMs=t.n_PWMs.sum(), n_tie=t.n_tie.sum(), agreeing_pairs=t.agreeing_pairs.sum(),
             pairs=t.pairs.sum(), rand_index=t.agreeing_pairs.sum() / t.pairs.sum())])], ignore_index=True)
    # summary rows: label in 'summary', TF left empty (the TF column of a Sup. Table holds TF names only)
    t.insert(0, "summary", np.where(t.TF.isin(tfs), "", t.TF)); t.loc[~t.TF.isin(tfs), "TF"] = np.nan
    t.to_csv(OUT["tf"], index=False)
    pd.DataFrame(per).to_csv(OUT["pwm"], index=False)
    pd.set_option("display.width", 250)
    print(t.drop(columns=["contingency"]).round(3).to_string(index=False))
    w2.provenance(HERE, "4_visualisation/figures/FigS7/hocomoco11_check/cluster_vs_v11.py", list(OUT.values()),
                  [IN["sub"], IN["abl"], IN["v11"]] + [os.path.join(S7.IN["gini"], "%s_gini.csv" % t_) for t_ in tfs],
                  params="{clustering: ablation cut per TF (%s), tie: %.2f, v11_min_overlap: %d, agreement: Rand index, ARI}" % (
                      ", ".join("%s %.1f" % kv for kv in cut.items() if kv[1] != 0.8) or "all 0.8", TIE, V11.MIN_OVERLAP))


if __name__ == "__main__":
    main()
