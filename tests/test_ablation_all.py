"""Ablation of the secondary motif subtype: 29 TFs, ten seeds (0-9) in every condition.
1. the tables hold 29 TFs x 4 conditions x 10 fits (116 summary rows, 1,160 fits), 20 TFs at cut 0.8 and 9 at 0.9, and the rule
   (mean change after removing the secondary subtype below the mean of the random removals by more than 3 SD, auROC or auPRC)
   flags the 9 TFs of the text: 6 in auROC, 8 in auPRC, 9 of the 20 at cut 0.8 and none of the 9 at cut 0.9.
   Red: T1_MUTATE=ablation_rule uses 2 SD instead of 3, which flags more TFs.
2. Sup. Table 9 sheet 'ablation' is the summary (116 rows).
3. Sup. Table 9 sheet 'ablation_fits' holds the 1,160 fits and gives the chance rate of the rule.
   Red: T1_MUTATE=ablation_rule.
"""
import os

import numpy as np
import pandas as pd
import pytest

import t1_helpers as T

pytestmark = pytest.mark.r2
TAB = os.path.join(T.R2, "3_analysis/feature_importance/tables")


def test_summary_has_29_tfs_ten_seeds_and_the_rule_flags_the_9_tfs():
    s = pd.read_csv(os.path.join(TAB, "B4_ablation_summary.csv")); f = pd.read_csv(os.path.join(TAB, "B4_ablation_all_fits.csv"))
    assert len(s) == 116 and s.TF.nunique() == 29 and len(f) == 1160 and (s.n_fits == 10).all()
    assert (f.groupby(["TF", "cond"]).seed.apply(lambda x: sorted(x) == list(range(10)))).all()
    cut = s.groupby("TF").sim_threshold.first()
    assert (cut == 0.8).sum() == 20 and (cut == 0.9).sum() == 9
    k = 2 if T.mutate("ablation_rule") else 3
    sec = s[s.cond == "drop_secondary"].set_index("TF"); rnd = s[s.cond == "drop_random"].set_index("TF")
    lost = {m: sorted(sec.index[sec["d_%s_mean" % m] < rnd["d_%s_mean" % m] - k * rnd["d_%s_sd" % m]]) for m in ("auroc", "auprc")}
    assert lost["auroc"] == ["GATA3", "IRF4", "JUND", "STAT1", "TAL1", "TFE2"], lost["auroc"]
    assert lost["auprc"] == ["GATA3", "IRF4", "PRGR", "SOX2", "STA5A", "STAT1", "TAL1", "TFE2"], lost["auprc"]
    either = sorted(set(lost["auroc"]) | set(lost["auprc"]))
    assert len(either) == 9 and sum(cut[t] == 0.8 for t in either) == 9 and sum(cut[t] == 0.9 for t in either) == 0
    print("flagged", either)


def test_sup_table_9_ablation_sheet_is_the_summary():
    s = pd.read_csv(os.path.join(TAB, "B4_ablation_summary.csv"))
    x = pd.read_excel(os.path.join(T.R2, "outputs/supplementary_tables/Sup_Table_9_feature_importance.xlsx"), sheet_name="ablation")
    assert len(x) == len(s) == 116 and list(x.columns) == list(s.columns)
    assert x.TF.nunique() == 29 and set(x.cond) == {"full", "drop_secondary", "drop_random", "drop_primary"}
    assert np.allclose(x.d_auroc_mean.astype(float), s.d_auroc_mean.astype(float), rtol=1e-9, atol=0)


def test_sup_table_9_ablation_fits_sheet_gives_the_chance_rate():
    """Sheet 'ablation_fits' holds every fit; the chance rate quoted in Methods 2.10 (each random draw against the other nine
    of its TF: 12 of 290 draws flagged in either metric) is recomputed from the sheet alone.
    Red: T1_MUTATE=ablation_rule (2 SD) flags more draws."""
    x = pd.read_excel(os.path.join(T.R2, "outputs/supplementary_tables/Sup_Table_9_feature_importance.xlsx"), sheet_name="ablation_fits")
    assert len(x) == 1160 and x.TF.nunique() == 29 and sorted(x.seed.unique()) == list(range(10))
    k = 2 if T.mutate("ablation_rule") else 3
    full = x[x.cond == "full"].set_index(["TF", "seed"])
    r = x[x.cond == "drop_random"].copy()
    key = pd.MultiIndex.from_frame(r[["TF", "seed"]])
    flagged = 0
    for m in ("auroc", "auprc"):
        r["d_" + m] = r[m].values - full[m].reindex(key).values
    for tf, g in r.groupby("TF"):
        g = g.sort_values("seed")
        for d in range(10):
            either = False
            for m in ("auroc", "auprc"):
                v = g["d_" + m].values; rest = np.delete(v, d)
                either |= bool(v[d] < rest.mean() - k * rest.std(ddof=1))
            flagged += either
    assert flagged == 12, flagged
