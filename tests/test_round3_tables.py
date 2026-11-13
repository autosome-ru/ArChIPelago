"""Fig. S11B label thresholds (= B5 group thresholds; Fig. 3C and S11B draw no lines), Sup. Table 3 sheets with the per-TF values of Fig. 4,
Fig. S4 and Fig. S5 and the paired RF vs Slim / diChIPMunk counts, and the negatives of Sup. Table 2.
Red: T1_MUTATE=f3c_thr reads the Fig. S11B ROC label threshold + 0.005; T1_MUTATE=t3_paired adds one TF to the human
Slim m=1 auPRC count of the Sup. Table 3 sheet; T1_MUTATE=t2_neg adds one negative to the first TF's human test negatives."""
import os
import re

import numpy as np
import pandas as pd
import pytest
from scipy.stats import wilcoxon

import t1_helpers as T

pytestmark = pytest.mark.r2


def J(*p):
    return os.path.join(T.R2, *p)


def test_fig3c_thresholds_equal_b5():
    r = open(J("4_visualisation", "main_figures", "Figure_3", "Figure_3_and_S3.R")).read()
    b5 = open(J("3_analysis/transfer_drivers", "B5_01_drivers_and_ceiling.py")).read()
    f11 = open(J("4_visualisation", "supplementary_figures", "Figure_S11", "Figure_S11.py")).read()
    roc_s, prc_s = re.search(r"H\.dROC > ([0-9.]+)\) \| \(H\.dPRC > ([0-9.]+)\)", f11).groups()   # Fig. S11B labels
    roc_r, prc_r = float(roc_s), float(prc_s)
    if T.mutate("f3c_thr"):
        roc_r += 0.005
    roc_b = float(re.search(r"^ROC_GAP = ([0-9.]+)", b5, re.M).group(1))
    prc_b = float(re.search(r"^PRC_GAP = ([0-9.]+)", b5, re.M).group(1))
    assert (roc_r, prc_r) == (roc_b, prc_b), ((roc_r, prc_r), (roc_b, prc_b))
    assert "group thresholds" not in r and "F3C_" not in r and "dashed" not in r, "Fig. 3C draws no group-threshold lines"
    assert "axhline(0.075" not in f11 and "group thresholds" not in f11, "Fig. S11B draws no group-threshold lines"


def test_sup_table3_fig4_sheets_and_paired_counts():
    t3 = pd.read_excel(J("outputs/supplementary_tables", "Sup_Table_3_ArChIPelago_performance.xlsx"), sheet_name=None)
    IN = T.h3in().IN
    src = {"Fig4_human_test": IN["fig4_src"], "FigS4_mouse_test": IN["figS4_src"],
           "FigS5_per_TF": J("3_analysis", "package", "Figures", "source_data", "Figure_S5_source_data.csv")}
    for sh, f in src.items():
        a, b = t3[sh], pd.read_csv(f)
        assert list(a.columns) == list(b.columns) and len(a) == len(b), sh
        num = a.select_dtypes("number").columns
        assert np.allclose(a[num].values, b[num].values, equal_nan=True), sh
    p = t3["Fig4_S4_paired_counts"].set_index(["test_set", "comparator"])
    if T.mutate("t3_paired"):
        p.loc[("human", "Slim m=1"), "n_TF_RF_above_auPRC"] += 1
    assert len(p) == 8
    for ts, f in (("human", IN["fig4_src"]), ("mouse", IN["figS4_src"])):
        d = pd.read_csv(f)
        rf = d[d.model == "ArChIPelago RF on all PWMs"].set_index("TF")
        for c in ("Slim m=0", "Slim m=1", "LSlim m=-5", "diChIPMunk"):
            for m in ("auROC", "auPRC"):
                dd = rf[m] - d[d.model == c].set_index("TF").loc[rf.index, m]
                assert p.loc[(ts, c), "n_TF_RF_above_" + m] == int((dd > 0).sum()), (ts, c, m)
                assert abs(p.loc[(ts, c), "wilcoxon_P_two_sided_" + m] / wilcoxon(dd).pvalue - 1) < 1e-9, (ts, c, m)
    # the human Slim counts are the ones NUMBERS and the text quote
    hp = p.loc["human"]
    if T.h3in().STAGE != "H4":   # the values the text quotes for the previous run; the flip report lists the changed counts of the final run
        return
    assert [int(hp.loc[c, "n_TF_RF_above_auROC"]) for c in ("Slim m=0", "Slim m=1", "LSlim m=-5")] == [35, 33, 34]
    assert [int(hp.loc[c, "n_TF_RF_above_auPRC"]) for c in ("Slim m=0", "Slim m=1", "LSlim m=-5")] == [33, 27, 33]


def test_sup_table2_negatives_equal_the_models_and_ratios_recompute():
    t2 = pd.read_excel(J("outputs/supplementary_tables", "Sup_Table_2_PWMs_and_datasets.xlsx"), sheet_name="per_TF").set_index("TF")
    hm = pd.read_csv(T.h3in().IN["hm_mm"]).set_index("TF").loc[t2.index]
    C = T.h3in().COL   # split labels of the stage
    cols = {C[k].replace(" positives ", " negatives "): (h, C[k]) for k, h in (("HTR", "n_train_neg_H"), ("HTE", "n_test_neg_H"),
                                                                              ("MTR", "n_train_neg_M"), ("MTE", "n_test_neg_M"))}
    if T.mutate("t2_neg"):
        t2.iloc[0, t2.columns.get_loc(C["HTE"].replace(" positives ", " negatives "))] += 1
    for c, (h, pc) in cols.items():
        assert (t2[c].values == hm[h].values).all(), c
        r = t2[c.replace(" negatives ", " negatives per positive ")]
        assert np.allclose(r.values, (t2[c] / t2[pc]).round(1).values), c
    # the realized ratios of Methods 2.2 (NUMBERS 'realised negatives per positive used by the models'): human 1:4.5-1:100, mouse 1:4.2-1:82
    rh = pd.concat([t2[c] / t2[cols[c][1]] for c in list(cols)[:2]])
    rm = pd.concat([t2[c] / t2[cols[c][1]] for c in list(cols)[2:]])
    if T.h3in().STAGE == "H4":   # the values the text quotes for the previous run; the flip report and NUMBERS carry the new range of the final run
        assert ("%.1f" % rh.min(), "%.1f" % rh.max(), "%.1f" % rm.min(), "%.1f" % rm.max()) == ("4.5", "100.0", "4.2", "82.0")
