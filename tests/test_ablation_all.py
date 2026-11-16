"""Ablation of the secondary motif subtype: 29 TFs, ten seeds (0-9) in every condition.
1. the new summary restricted to seeds 0-2 reproduces the three-seed summary of the 11 TFs of the earlier release (fixture
   ablation_summary_3seeds_11TF.csv, 44 rows), recomputed here from the per-fit files of the server mirror
   Red: T1_MUTATE=ablation_seeds keeps seeds 0-3, so the means and SDs differ from the fixture.
2. the tables hold 29 TFs x 4 conditions x 10 fits (116 summary rows, 1,160 fits), 20 TFs at cut 0.8 and 9 at 0.9, and the rule
   (mean change after removing the secondary subtype below the mean of the random removals by more than 3 SD, auROC or auPRC)
   flags the 9 TFs of the text: 6 in auROC, 8 in auPRC, 9 of the 20 at cut 0.8 and none of the 9 at cut 0.9.
   Red: T1_MUTATE=ablation_rule uses 2 SD instead of 3, which flags more TFs.
3. Sup. Table 9 sheet 'ablation' is the summary (116 rows).
"""
import glob
import os

import numpy as np
import pandas as pd
import pytest

import t1_helpers as T

pytestmark = pytest.mark.r2
HERE = os.path.dirname(os.path.abspath(__file__))
TAB = os.path.join(T.R2, "3_analysis/feature_importance/tables")
MIRROR = os.path.join(T.R2, "results/server/human/B4_importance/ablation")
EARLIER = ["E2F4", "IRF4", "MAFK", "PRGR", "RXRA", "SOX2", "SPI1", "SRF", "STA5A", "TAL1", "TFE2"]


def summarise(a):
    """Per TF and condition: change against the full model of the same seed; mean, SD (ddof 1), min, max."""
    rows = []
    for tf, g in a.groupby("TF"):
        full = g[g.cond == "full"].set_index("seed")
        for cond, c in g.groupby("cond"):
            c = c.set_index("seed")
            dr = c.auroc - full.auroc.reindex(c.index); dp = c.auprc - full.auprc.reindex(c.index)
            rows.append(dict(TF=tf, cond=cond, n_fits=len(c), d_auroc_mean=dr.mean(), d_auroc_sd=dr.std(), d_auprc_mean=dp.mean(), d_auprc_sd=dp.std(),
                             d_auroc_min=dr.min(), d_auroc_max=dr.max(), d_auprc_min=dp.min(), d_auprc_max=dp.max()))
    return pd.DataFrame(rows).set_index(["TF", "cond"])


def test_seeds_0_to_2_of_the_new_files_reproduce_the_three_seed_summary_of_the_11_earlier_tfs():
    top = 3 if T.mutate("ablation_seeds") else 2
    a = pd.concat([pd.read_csv(os.path.join(MIRROR, "%s_ablation.csv" % tf)) for tf in EARLIER], ignore_index=True)
    s = summarise(a[a.seed <= top])
    fx = pd.read_csv(os.path.join(HERE, "ablation_summary_3seeds_11TF.csv")).set_index(["TF", "cond"])
    assert len(fx) == 44 and len(s) == 44
    s = s.loc[fx.index]
    worst = 0.0
    for c in ("n_fits", "d_auroc_mean", "d_auroc_sd", "d_auprc_mean", "d_auprc_sd", "d_auroc_min", "d_auroc_max", "d_auprc_min", "d_auprc_max"):
        rel = np.abs(s[c].values - fx[c].values) / np.maximum(np.abs(fx[c].values), 1e-12)
        worst = max(worst, rel.max())
    assert worst < 6e-5, "largest relative difference to the three-seed summary: %.2e" % worst   # the fixture carries 5 significant digits
    print("largest relative difference %.2e over 44 rows" % worst)


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
