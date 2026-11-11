"""Content checks of Sup. Tables 2, 4, 5, 6, 8, 9 and 10. Each test recomputes the checked cells from the source table of the
producer and prints what it saw.
Red: T1_MUTATE=<name> (named in each test) feeds an earlier variant of the table.
"""
import os

import numpy as np
import pandas as pd
import pytest

import t1_helpers as T

pytestmark = pytest.mark.r2


def J(*p):
    return os.path.join(T.R2, *p)


def xl(n, sheet):
    f = {2: "Sup_Table_2_PWMs_and_datasets", 4: "Sup_Table_4_operational_metrics", 5: "Sup_Table_5_cross_species_and_mouse_trained",
         6: "Sup_Table_6_runtime_memory", 8: "Sup_Table_8_top_peaks_and_peak_widths", 9: "Sup_Table_9_feature_importance",
         10: "Sup_Table_10_transfer_drivers_and_headroom"}[n]
    return pd.read_excel(J("outputs/supplementary_tables", f + ".xlsx"), sheet_name=sheet)


def test_sup5_hm_two_seed_mean_columns_and_family_header():
    """The H>M two-seed mean columns equal hm_mm (mean of seeds 0 and 1), the summary row recomputes, and the family column header is the short label.
    Red: T1_MUTATE=sup5_hm_s0 expects the seed 0 values in the mean columns."""
    from scipy.stats import wilcoxon
    p = xl(5, "per_TF").set_index("TF")
    h = pd.read_csv(T.h3in().IN["hm_mm"]).set_index("TF").loc[p.index]
    v = "s0" if T.mutate("sup5_hm_s0") else "mean"
    d = {}
    for m, M, b in (("auroc", "auROC", "broc"), ("auprc", "auPRC", "bprc")):
        rf = h["HM_rf_test_M_%s_%s" % (m, v)]
        d[M] = rf - h["HM_%s_test_M_%s" % (b, m)]
        got = p["H>M: ArChIPelago %s (mouse test, mean of 2 seeds)" % M]
        print("%s: max |table - hm_mm %s| %.2g; TFs below %s" % (M, v, (got - rf).abs().max(), list(d[M].index[d[M] < 0])))
        assert np.allclose(got, rf, atol=5e-5), M
        assert np.allclose(p["H>M: d%s (mean of 2 seeds)" % M], d[M], atol=1e-4), M
        assert "H>M: d%s (seed 0)" % M in p.columns and "H>M: d%s" % M not in p.columns
    s = xl(5, "summary")
    r = s[s.iloc[:, 0].str.endswith("mean of 2 seeds") & s.iloc[:, 0].str.startswith("Human-trained")].iloc[0]
    print("summary mean row: %s up %d / %d, P %s / %s, below %s" % (r.iloc[0][-15:], r["TFs with dauROC > 0"], r["TFs with dauPRC > 0"],
                                                                    r["Wilcoxon P (auROC)"], r["Wilcoxon P (auPRC)"], r["TFs below baseline on >= 1 metric"]))
    for M in ("auROC", "auPRC"):
        assert r["TFs with d%s > 0" % M] == (d[M] > 0).sum()
        assert abs(r["median d%s" % M] - d[M].median()) <= 5e-5
        assert float("%.1e" % wilcoxon(d[M]).pvalue) == float(r["Wilcoxon P (%s)" % M])
    assert r["TFs below baseline on >= 1 metric"] == ", ".join(sorted(set(d["auROC"].index[d["auROC"] < 0]) | set(d["auPRC"].index[d["auPRC"] < 0])))
    assert "TF family or subfamily (short label)" in p.columns and "TF family (TFClass)" not in p.columns


def test_sup2_gc_difference_reproduces_the_methods_correlations():
    """|GC difference| of the test sets in Sup. Table 2 gives the Spearman rho of 3_analysis/integration/gc_vs_gain_h3.csv.
    Red: T1_MUTATE=sup2_gc_rounded uses the difference of the rounded GC medians (an earlier version of the table)."""
    from scipy.stats import spearmanr
    t2 = xl(2, "per_TF").set_index("TF")
    ref = pd.read_csv(J("3_analysis/integration", "gc_vs_gain_h3.csv")).set_index(["species", "metric"])
    rt = pd.read_csv(J("3_analysis", "package", "results_table.csv"), sep="\t")
    rf = rt[(rt.Model == "RandomForestClassifier") & (rt.PWM == "mono+di")].set_index("TF_name").loc[t2.index]
    h = pd.read_csv(T.h3in().IN["hm_mm"]).set_index("TF").loc[t2.index]
    gain = {("human", "auROC"): rf.roc_auc_test_H - rf.roc_auc_test_H_PWM, ("human", "auPRC"): rf.pr_auc_test_H - rf.pr_auc_test_H_PWM,
            ("mouse", "auROC"): h.MM_rf_test_M_auroc_mean - h.MM_broc_test_M_auroc, ("mouse", "auPRC"): h.MM_rf_test_M_auprc_mean - h.MM_bprc_test_M_auprc}
    for (sp, m), g in gain.items():
        lab = sp.capitalize() + " test"
        if T.mutate("sup2_gc_rounded"):
            pos = [c for c in t2.columns if c.startswith(lab + " positives, median GC")][0]
            neg = [c for c in t2.columns if c.startswith(lab + " negatives, median GC")][0]
            x = t2[pos] - t2[neg]
        else:
            x = t2[[c for c in t2.columns if c.startswith(lab + ", median GC difference")][0]]
        r = spearmanr(x.abs(), g)
        print("%s %s: rho %.4f (table %.4f), P %.3f" % (sp, m, r[0], ref.loc[(sp, m), "spearman_rho"], r[1]))
        assert abs(r[0] - ref.loc[(sp, m), "spearman_rho"]) < 1e-9


def test_sup4_tfs_in_the_median_and_precision_gain_decimals():
    """'TFs in the median' = TFs with a defined value (34 for the false positives at 50 % recall, REST undefined);
    the per-TF precision gain has 4 decimals, so that its count of TFs improved equals the summary.
    Red: T1_MUTATE=sup4_round3 rounds the per-TF precision gain to 3 decimals (an earlier version of the sheet)."""
    s, p = xl(4, "summary"), xl(4, "per_TF")
    col = {"FP reduction at 50 % recall (%)": "FP reduction at 50 % recall, ArChIPelago vs best monoPWM (%)",
           "FP reduction at 80 % recall (%)": "FP reduction at 80 % recall, ArChIPelago vs best monoPWM (%)",
           "additional true sites in top 1,000": "Additional true sites in top 1,000, ArChIPelago vs best monoPWM",
           "precision gain at 50 % recall (paired, per TF)": "Precision gain at 50 % recall, ArChIPelago minus best monoPWM"}
    pg = col["precision gain at 50 % recall (paired, per TF)"]
    if T.mutate("sup4_round3"):
        p = p.assign(**{pg: p[pg].round(3)})
    for _, x in s.iterrows():
        v = p[p["Test set"].str[:5] == x["Test set"][:5]][col[x.Quantity]]
        n_up = int((v > 0).sum())
        print("%s %-46s TFs in the median %d (defined %d); improved %d (sheet %d)" % (x["Test set"][:5], x.Quantity, x["TFs in the median"],
                                                                                      v.notna().sum(), n_up, x["TFs with improvement (of 35)"]))
        assert x["TFs in the median"] == v.notna().sum()
        assert n_up == x["TFs with improvement (of 35)"]
    assert (s[s.Quantity == "FP reduction at 50 % recall (%)"]["TFs in the median"] == 34).all()
    assert "Median over TFs" in s.columns


def test_sup6_scanner_rate_from_the_scanning_sheet():
    """The companion-scanner formula quotes the per-Mb rates of sheet PWM_scanning (single-scan mean, lowest loop rate).
    Red: T1_MUTATE=sup6_old_rate reads the first-submission note (0.15 s/Mb)."""
    import re
    note = xl(6, "RF_inference").iloc[:, 0].dropna().iloc[-1]
    if T.mutate("sup6_old_rate"):
        note = note.replace(re.search(r"r x Mb, with r = [^)]*\)", note).group(0), "0.15 s/Mb x Mb")
    ps = " ".join(xl(6, "PWM_scanning").iloc[:, 0].dropna().astype(str))
    single = float(re.search(r"per monoPWM per Mb \(mean ([\d.]+)\)", ps).group(1))
    loop_lo = float(re.search(r"\): ([\d.]+)-[\d.]+ s per monoPWM per Mb", ps).group(1))
    want = "r = %.2f-%.2f s/Mb" % (single, loop_lo)
    print("PWM_scanning: single-scan mean %.3f, loops from %.3f -> %r; note: %s" % (single, loop_lo, want, note[:120]))
    assert want in note and "0.15 s/Mb x Mb" not in note


def test_sup8_uncapped_columns_only_for_the_positives_predictor():
    """n_tf_uncapped / rho_uncapped / p_uncapped are filled only for the training positives (n_train_pos_H), and there
    they leave out the capped TFs (n < 35). Red: T1_MUTATE=sup8_unblank fills them from the full columns (an earlier version of the sheet)."""
    d = xl(8, "gain_vs_positives_rho")
    if T.mutate("sup8_unblank"):
        d = d.assign(n_tf_uncapped=d.n_tf_uncapped.fillna(d.n_tf), rho_uncapped=d.rho_uncapped.fillna(d.rho))
    pos = d.predictor_col == "n_train_pos_H"
    print("uncapped filled: positives rows %d of %d (n_tf_uncapped %s), other rows %d of %d" % (
        d[pos].rho_uncapped.notna().sum(), pos.sum(), sorted(set(d[pos].n_tf_uncapped)), d[~pos].rho_uncapped.notna().sum(), (~pos).sum()))
    assert d[pos][["n_tf_uncapped", "rho_uncapped", "p_uncapped"]].notna().all().all()
    assert (d[pos].n_tf_uncapped < d[pos].n_tf).all()
    assert d[~pos][["n_tf_uncapped", "rho_uncapped", "p_uncapped"]].isna().all().all()


def test_sup9_activation_best_pwm_sheet_is_the_fig_s10d_source():
    """Sheet 'activation_best_pwm' holds 3_analysis/feature_importance/tables/B4_activation_best_pwm.csv (read by b4_06_figure.py for
    Fig. S10D). Red: T1_MUTATE=sup9_no_sheet drops the sheet."""
    sheets = pd.read_excel(J("outputs/supplementary_tables", "Sup_Table_9_feature_importance.xlsx"), sheet_name=None)
    if T.mutate("sup9_no_sheet"):
        sheets.pop("activation_best_pwm", None)
    src = pd.read_csv(J("3_analysis/feature_importance", "tables", "B4_activation_best_pwm.csv"))
    print("Sup. Table 9 sheets: %s" % list(sheets))
    assert "activation_best_pwm" in sheets
    got = sheets["activation_best_pwm"]
    print("activation_best_pwm %d rows, source %d rows" % (len(got), len(src)))
    pd.testing.assert_frame_equal(got.reset_index(drop=True), src[got.columns].reset_index(drop=True), check_dtype=False, atol=1e-9)
    assert "Fig. S10D" not in "".join(got.columns)
    assert "B4_activation_best_pwm.csv" in open(J("3_analysis/feature_importance", "scripts", "b4_06_figure.py")).read()


def test_sup10_multivariable_sheet_has_coefficients_and_bootstrap_intervals():
    """Sheet 'multivariable' holds per outcome and predictor the standardized coefficient, its bootstrap 95 % interval, VIF
    and BH q (B5_A_multivariable_rank_ols.csv), with n, k, R2 and adjusted R2 of the outcome (B5_A_multivariable_summary.csv).
    Red: T1_MUTATE=sup10_summary_only reads an earlier sheet layout (summary rows only)."""
    m = xl(10, "multivariable")
    if T.mutate("sup10_summary_only"):
        m = pd.read_csv(J("3_analysis/transfer_drivers", "B5_A_multivariable_summary.csv")).drop(columns="status")
    co = pd.read_csv(J("3_analysis/transfer_drivers", "B5_A_multivariable_rank_ols.csv"))
    su = pd.read_csv(J("3_analysis/transfer_drivers", "B5_A_multivariable_summary.csv")).set_index("outcome")
    print("multivariable: %d rows, columns %s" % (len(m), list(m.columns)))
    need = ["outcome", "predictor", "beta_std_rank", "boot_ci_lo", "boot_ci_hi", "VIF", "q_BH", "R2", "adj_R2"]
    assert all(c in m.columns for c in need) and len(m) == len(co)
    k = m.set_index(["outcome", "predictor"]).loc[list(zip(co.outcome, co.predictor))]
    for c in ("beta_std_rank", "boot_ci_lo", "boot_ci_hi", "q_BH"):
        assert np.allclose(k[c].values, co[c].values), c
    assert np.allclose(k.R2.values, su.loc[co.outcome, "R2"].values)
