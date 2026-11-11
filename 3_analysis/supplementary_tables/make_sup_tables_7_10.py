"""Sup. Tables 7-10 of the revision: one xlsx per table, sheets copied unchanged from the result CSVs of the low-FPR, top-peak,
feature-importance and transfer-driver analyses (35 TFs, final data).

Output: outputs/supplementary_tables/Sup_Table_<n>_<name>.xlsx
usage: python3 3_analysis/supplementary_tables/make_sup_tables_7_10.py [7] [8] [9] [10]   (default: all four)
  7  B2 summaries of 3_analysis/low_fpr_metrics (HH, HM, MM; corrected partial AUC), random_split_overlap (human 361, mouse 276
     sets, repeat mask as the pipeline; shares only, no window counts), top1000_check from the check logs of the rebuilds
     (3_analysis/integration/top1000_check_h3.py; columns that count unique windows or hold the shifted-TF control are not
     shown), random-split retraining of 3_analysis/random_split (v3 negatives, 5 TFs)
  8  3_analysis/top_peaks (B3 summaries); gain_vs_positives_rho: the uncapped columns (n_tf_uncapped, rho_uncapped, p_uncapped)
     are blanked for predictors other than the training positives, because the cap (10,000 positives) concerns only that predictor
  9  3_analysis/feature_importance/tables (b4_05_aggregate.py); sheet 'activation_best_pwm' (B4_activation_best_pwm.csv, the
     source of Fig. S10D); 'multivariable' of Sup. Table 10 is described below
  10 3_analysis/transfer_drivers (run with the package results table and the final hm_mm); 'multivariable': the standardised
     rank-regression coefficients with their 2,000-resample bootstrap 95 % intervals, VIF and BH q per predictor and outcome
     (B5_A_multivariable_rank_ols.csv), joined with n, k, R2 and adjusted R2 of each outcome (B5_A_multivariable_summary.csv)
  Sup. Table 8 sheet datasets_without_peaks is not written: the 19 human sets without a peak file are not in Sup. Table 1.
Columns status / mouse_columns_status must all read "current" and are then dropped.
A table is only written when no sheet holds a TF of p53drop.EXCLUDE_TFS.
A sheet entry may carry a third element: "drop" removes the rows of EXCLUDE_TFS, a list keeps only those columns, a function
transforms the table. A source may be a tuple of CSVs, merged on their shared key columns (status columns checked first).
"""
import os
import sys
import pandas as pd

R2 = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(R2, "outputs/supplementary_tables")
sys.path.insert(0, os.path.join(R2, "3_analysis/integration"))
import p53drop as PD
sys.path.insert(0, os.path.join(R2, "3_analysis", "integration"))
import w2
BUILD = set(sys.argv[1:]) or {"7", "8", "9", "10"}
os.makedirs(OUT, exist_ok=True)
UNCAPPED_PREDICTOR = "n_train_pos_H"   # the only predictor whose uncapped variant is defined (cap of the training positives)


def blank_uncapped(df):
    """n_tf_uncapped / rho_uncapped / p_uncapped only for the positives predictor."""
    cols = ["n_tf_uncapped", "rho_uncapped", "p_uncapped"]
    other = df.predictor_col != UNCAPPED_PREDICTOR
    if not other.any() or other.all():
        raise SystemExit("CHECK FAILED gain_vs_positives_rho: predictor %s present in %d of %d rows" % (UNCAPPED_PREDICTOR, (~other).sum(), len(df)))
    df = df.copy()
    df[cols] = df[cols].astype(float)
    df.loc[other, cols] = float("nan")
    print("gain_vs_positives_rho: uncapped columns kept for %d rows (%s), blanked for %d" % ((~other).sum(), UNCAPPED_PREDICTOR, other.sum()))
    return df

TABLES = {
    "Sup_Table_7_low_FPR_metrics": [
        ("summary", "3_analysis/low_fpr_metrics/B2_summary.csv"),
        ("per_TF", "3_analysis/low_fpr_metrics/B2_metrics_wide.csv"),
        ("per_chromosome_summary", "3_analysis/low_fpr_metrics/B2_per_chromosome_summary.csv"),
        ("per_chromosome", "3_analysis/low_fpr_metrics/B2_per_chromosome.csv"),
        ("random_split_overlap", "results/datasets/split_leakage/b8w1g2_leakage_per_TF_masked.csv",
         ["TF", "species", "split", "seed", "n_datasets", "frac_overlap", "frac_half", "frac_other_ds"]),   # window counts not shown
        ("top1000_check", "3_analysis/integration/top1000_check_h3.csv"),   # check C of the two rebuilds (the earlier check counted the published human positives)
        ("random_split_retraining", "3_analysis/random_split/b9_random_split_per_fit.csv"),   # final retraining, v3 negatives, 5 TFs
        ("random_split_retraining_summary", "3_analysis/random_split/b9_summary.csv"),
    ],
    "Sup_Table_8_top_peaks_and_peak_widths": [
        ("topN_summary", "3_analysis/top_peaks/b3_topN_summary.csv"),
        ("topN_per_TF", "3_analysis/top_peaks/b3_topN_results.csv"),
        ("gain_vs_positives", "3_analysis/top_peaks/b3_gain_vs_npos_per_tf.csv"),
        ("gain_vs_positives_rho", "3_analysis/top_peaks/b3_gain_vs_npos_spearman.csv", blank_uncapped),
        ("peak_width_summary", "3_analysis/top_peaks/b3_peak_widths_summary.csv"),
        # n_peaks_table (peak count of the first-submission Sup. Table 1) is empty for every dataset because Sup. Table 1
        # carries no peak count; dropped, n_peaks is the count of the peak file itself
        ("peak_width_per_dataset", "3_analysis/top_peaks/b3_peak_widths_per_dataset.csv",
         ["TF", "species", "peak_set", "n_peaks", "width_median", "width_q25", "width_q75", "frac_width_le300", "frac_window_covers_ge50",
          "median_cover", "width_median_top1000", "width_median_rest", "cover_ge50_top1000", "cover_ge50_rest"]),
    ],
    "Sup_Table_9_feature_importance": [
        ("concentration", "3_analysis/feature_importance/tables/B4_importance_concentration.csv"),
        ("gini_all_PWMs", "3_analysis/feature_importance/tables/B4_gini_all_features.csv"),
        ("permutation_top20", "3_analysis/feature_importance/tables/B4_permutation_top20.csv"),
        ("subtypes_per_TF", "3_analysis/feature_importance/tables/B4_subtype_summary_per_TF.csv"),
        ("subtype_importance", "3_analysis/feature_importance/tables/B4_subtype_importance.csv"),
        ("activation", "3_analysis/feature_importance/tables/B4_activation_groups.csv"),
        ("activation_best_pwm", "3_analysis/feature_importance/tables/B4_activation_best_pwm.csv"),   # Fig. S10D
        ("ablation", "3_analysis/feature_importance/tables/B4_ablation_summary.csv"),
        ("PWM_origin", "3_analysis/feature_importance/tables/B4_pwm_origin_check.csv"),
        # ablation subtypes vs nearest HOCOMOCO v11 subtype, per TF and per PWM (265 PWMs of the 11 ablated TFs)
        ("subtypes_vs_HOCOMOCO_v11", "4_visualisation/figures/FigS7/hocomoco11_check/cluster_vs_v11.csv"),
        ("HOCOMOCO_v11_per_PWM", "4_visualisation/figures/FigS7/hocomoco11_check/cluster_vs_v11_per_PWM.csv"),
    ],
    "Sup_Table_10_transfer_drivers_and_headroom": [
        ("per_TF", "3_analysis/transfer_drivers/B5_per_TF_table.csv"),
        ("drivers_spearman", "3_analysis/transfer_drivers/B5_A_spearman_drivers.csv"),
        ("multivariable", ("3_analysis/transfer_drivers/B5_A_multivariable_summary.csv", "3_analysis/transfer_drivers/B5_A_multivariable_rank_ols.csv")),
        ("headroom_spearman", "3_analysis/transfer_drivers/B5_B_ceiling_spearman.csv"),
        ("headroom_per_TF", "3_analysis/transfer_drivers/B5_B_headroom_per_TF.csv"),
        ("Fig3C_groups", "3_analysis/transfer_drivers/B5_B_cloud_separation.csv"),
    ],
}

FULL_TF_SHEETS = {"per_TF", "per_chromosome", "random_split_overlap", "top1000_check", "topN_per_TF", "gain_vs_positives",
                  "peak_width_per_dataset", "headroom_per_TF"}
for name, sheets in TABLES.items():
    if name.split("_")[2] not in BUILD:
        continue
    path = os.path.join(OUT, name + ".xlsx")
    dfs = []
    for entry in sheets:
        sheet, src = entry[:2]
        if isinstance(src, tuple):   # tables merged on their shared columns after the status check
            parts = []
            for s_ in src:
                p_ = pd.read_csv(os.path.join(R2, s_))
                if "status" in p_.columns:
                    if not p_.status.astype(str).str.startswith("current").all():
                        raise SystemExit("CHECK FAILED %s / %s: status not current in %s" % (name, sheet, s_))
                    p_ = p_.drop(columns="status")
                parts.append(p_)
            df = parts[0]
            for p_ in parts[1:]:
                key = [c for c in df.columns if c in p_.columns]
                n_ = len(p_)
                df = df.merge(p_, on=key, how="right", validate="one_to_many")
                if len(df) != n_ or df.isna().any().any():
                    raise SystemExit("CHECK FAILED %s / %s: merge on %s gave %d rows for %d" % (name, sheet, key, len(df), n_))
            src = " + ".join(src)
        else:
            df = pd.read_csv(os.path.join(R2, src))
        opt = entry[2] if len(entry) > 2 else None
        if callable(opt):
            df = opt(df)
        elif opt == "drop":
            n0 = len(df); df = PD.drop(df)
            print("%s / %s: %d rows of %s dropped" % (name, sheet, n0 - len(df), ",".join(PD.EXCLUDE_TFS)))
        elif opt is not None:
            df = df[opt]
        _empty = [c for c in df.columns if df[c].isna().all()]
        if _empty:   # no published column may be empty
            raise SystemExit("CHECK FAILED %s / %s: empty column(s) %s" % (name, sheet, _empty))
        for _st in [c for c in df.columns if c in ("status", "mouse_columns_status")]:   # work-status columns are not published
            if not df[_st].astype(str).str.startswith("current").all():
                raise SystemExit("CHECK FAILED %s / %s: %s not current: %s" % (name, sheet, _st, sorted(set(df[_st]))))
            df = df.drop(columns=_st)
        if sheet in FULL_TF_SHEETS:   # these sheets must hold the 35 reference TFs in every scenario / species
            grp = next((c for c in ("scenario", "species") if c in df.columns), None)
            for _sc, _g in (df.groupby(grp) if grp else [("all", df)]):
                PD.check_tfs(_g.TF, "%s / %s %s" % (name, sheet, _sc))
        bad = sorted(set(df.TF) & set(PD.EXCLUDE_TFS)) if "TF" in df.columns else []
        if bad:
            raise SystemExit("CHECK FAILED %s / %s contains %s (from %s)" % (name, sheet, bad, src))
        dfs.append((sheet, src, df))
    w2.backup(path)
    with pd.ExcelWriter(path) as xw:
        for sheet, src, df in dfs:
            df.to_excel(xw, sheet_name=sheet, index=False)
            print("%s / %s: %d rows from %s" % (name, sheet, len(df), src))
    w2.provenance(OUT, "3_analysis/supplementary_tables/make_sup_tables_7_10.py (%s)" % name, [path], [os.path.join(R2, s_) for _, src, _ in dfs for s_ in src.split(" + ")])
