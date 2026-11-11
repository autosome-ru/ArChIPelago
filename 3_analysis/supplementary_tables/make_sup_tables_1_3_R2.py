"""Supplementary Tables 1-4 and 6 of the revision (35 reference TFs, final data).
Sup. Table 1 lists exactly the peak sets used, human 361 (results/datasets/peak_ids/human_used_peak_ids.csv without P53) and
  mouse 276 (mouse_sym_peak_ids.csv). The 19 human datasets without a peak file (ESR1 11, GCR 8; 3_analysis/top_peaks/
  b3_missing_datasets.csv) are therefore absent: none of their peaks is among the positives (results/datasets/split_leakage:
  all ESR1 and GCR training and test positives match the other datasets). A column with the number of MACS peaks in the GTRD
  interval file is added (b3_peak_widths_per_dataset.csv); the GEM column is kept as it was.
Sup. Table 2 per_TF (P53 row dropped):
  PWM counts                  first-submission Sup. Table 2 (unchanged; mouse counts asserted equal to counts_per_tf.csv)
  experiments                 the Sup. Table 1 sheets written above (361 human, 276 mouse)
  pooled peak sets            not shown: the columns of the first submission are removed for both species
  human positives             orig_pos_* of human_counts_per_tf.csv = the published human training tables actually used
                              (asserted equal to the first-submission Sup. Table 2 and to hm_mm_results n_*_pos_H)
  mouse positives             counts_per_tf.csv pos_chr1_8_19_file_train (test) / pos_other_file_control (training), asserted
                              equal to hm_mm_results.csv n_test_pos_M / n_train_pos_M
  GC content                  per species and split the median GC difference, positives minus negatives (percentage points,
                              unrounded), from median_gc_diff of gc_check_v3 (negatives v3); the Spearman correlations of
                              Methods 2.2 (3_analysis/integration/gc_vs_gain_h3.csv) are computed from these values
Sup. Table 3: from 3_analysis/package/Sup_Tables (make_sup_tables_3_to_6.py on the 35-TF results table) plus the 'Number of PWMs'
  of the monoPWMs + diPWMs rows set to the total of monoPWMs and diPWMs (Sup. Table 2).
Sup. Table 4: copy of 3_analysis/package/Sup_Tables/Sup_Table_4_operational_metrics.xlsx. Summary: 'Median over TFs' with a
  column 'TFs in the median' (the false-positive reduction at 50 % recall is undefined for REST, which has no false positives
  for either model, so that median is over 34 TFs). per_TF: the precision gain at 50 % recall with PREC_GAIN_DECIMALS (4)
  decimals from the operational metrics table, so that small gains (PRGR +0.0002 on the mouse test set) are visible.
Every current input comes from 3_analysis/integration/h3_inputs.py and is recorded in provenance:
  Sup. Table 1 MACS peak counts  b3_widths_H / b3_widths_M (n_peaks of the interval files)
  Sup. Table 2 human positives   counts_H (each peak set read once: pos_train = file "train", pos_test = chr1/8/21),
                                 asserted equal to hm_mm n_train_pos_H / n_test_pos_H; human PWM counts asserted vs counts_H
               mouse positives   counts_M, asserted vs hm_mm n_*_pos_M
Output: outputs/supplementary_tables/Sup_Table_{1,2,3,4,6}_*.xlsx.
usage: python3 3_analysis/supplementary_tables/make_sup_tables_1_3_R2.py [1] [2] [3] [4] [6]  (default: all five)
"""
import os
import sys
import pandas as pd
R2 = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(R2, "3_analysis/integration"))
import p53drop as PD
sys.path.insert(0, os.path.join(R2, "3_analysis", "integration"))
import w2
import h3_inputs as H
BUILD = set(sys.argv[1:]) or {"1", "2", "3", "4", "6"}
SRC = os.path.join(R2, "results", "round1_package", "Resubmission", "Sup_Tables")
OUT = os.path.join(R2, "outputs/supplementary_tables")

miss = pd.read_csv(os.path.join(R2, "results/datasets/top_peaks", "b3_missing_datasets.csv"))
widths = pd.concat([pd.read_csv(H.IN["b3_widths_H"]), pd.read_csv(H.IN["b3_widths_M"])], ignore_index=True)
nmacs = dict(zip(widths.peak_set, widths.n_peaks))

t1 = pd.read_excel(os.path.join(SRC, "Sup_Table_1_ChIP-Seq_experiments.xlsx"), sheet_name=None)
USED = {"human": PD.human_ids(), "mouse": PD.mouse_ids()}
for sh in t1:
    d = t1[sh]
    n0 = len(d)
    d = d[~d["GTRD peak set ID"].isin(set(miss.peak_set))].copy()
    u = USED[sh]
    d = d[d["GTRD peak set ID"].isin(set(u.peak_id))].copy()
    chk = d[["TF", "GTRD peak set ID"]].merge(u, left_on=["TF", "GTRD peak set ID"], right_on=["TF", "peak_id"])
    if len(d) != len(u) or len(chk) != len(u) or d["GTRD peak set ID"].duplicated().any():
        raise SystemExit("CHECK FAILED Sup. Table 1 %s: %d rows, %d used sets, %d TF-matched" % (sh, len(d), len(u), len(chk)))
    PD.check_tfs(d.TF, "Sup. Table 1 " + sh)
    d.insert(d.columns.get_loc("GEM peaks (GTRD peak calls)") + 1, "MACS peaks (GTRD interval file, used in this study)",
             d["GTRD peak set ID"].map(nmacs))
    assert d["MACS peaks (GTRD interval file, used in this study)"].notna().all(), sh
    t1[sh] = d
    print("Sup. Table 1", sh, n0, "->", len(d))
if "1" in BUILD:
    p1 = os.path.join(OUT, "Sup_Table_1_ChIP-Seq_experiments.xlsx")
    PD.backup(p1)
    with pd.ExcelWriter(p1) as xw:
        for sh, d in t1.items():
            d.to_excel(xw, sheet_name=sh, index=False)
    PD.write_provenance(OUT, "3_analysis/supplementary_tables/make_sup_tables_1_3_R2.py (Sup. Table 1)", [p1],
                        [os.path.join(SRC, "Sup_Table_1_ChIP-Seq_experiments.xlsx"), H.IN["b3_widths_H"], H.IN["b3_widths_M"],
                         os.path.join(R2, "results/datasets/top_peaks", "b3_missing_datasets.csv")],
                        note="human %d, mouse %d used peak sets" % (len(t1["human"]), len(t1["mouse"])))
if not BUILD & {"2", "3", "4", "6"}:
    sys.exit(0)
import numpy as np
B11 = os.path.join(R2, "B11_mouse_sym")
PKGC = os.path.join(R2, "3_analysis", "package", "Sup_Tables")
IN2 = {"t2_round1": os.path.join(SRC, "Sup_Table_2_PWMs_and_datasets.xlsx"), "counts_mouse": H.IN["counts_M"],
       "counts_human": H.IN["counts_H"], "gc": H.IN["gc_v3"], "hm_mm": H.IN["hm_mm"]}
t2 = pd.read_excel(IN2["t2_round1"], sheet_name=None)
per = PD.drop(t2["per_TF"]).set_index("TF").sort_index()
NTF = PD.check_tfs(per.index, "Sup. Table 2 per_TF")
TFS = list(per.index)
cm = pd.read_csv(IN2["counts_mouse"]).set_index("TF").drop("POOLED").loc[TFS]
ch = pd.read_csv(IN2["counts_human"]).set_index("TF").loc[TFS]
hm = pd.read_csv(IN2["hm_mm"]).set_index("TF").loc[TFS]
gc = pd.read_csv(IN2["gc"])
HTR, HTE, MTR, MTE = H.COL["HTR"], H.COL["HTE"], H.COL["MTR"], H.COL["MTE"]   # split labels of the final run (h3_inputs)
# checks of the sources
for a, b, lab in ((per["Mouse monoPWMs"], cm.n_MOUSE_pwm_mono, "mouse monoPWMs"), (per["Mouse diPWMs"], cm.n_MOUSE_pwm_di, "mouse diPWMs"),
                  (per["Human monoPWMs"], cm.n_HUMAN_pwm_mono, "human monoPWMs"), (per["Human diPWMs"], cm.n_HUMAN_pwm_di, "human diPWMs"),
                  (ch.pos_train, hm.n_train_pos_H, "human training positives vs hm_mm"), (ch.pos_test, hm.n_test_pos_H, "human test positives vs hm_mm"),
                  (cm.pos_other_file_control, hm.n_train_pos_M, "mouse training positives vs hm_mm"), (cm.pos_chr1_8_19_file_train, hm.n_test_pos_M, "mouse test positives vs hm_mm")):
    if not (a.values == b.values).all():
        raise SystemExit("CHECK FAILED Sup. Table 2 source: %s differ for %s" % (lab, list(a.index[a.values != b.values])))
print("Sup. Table 2 sources agree (PWM counts, rebuilt human and mouse positives vs hm_mm)")
col = "Human ChIP-Seq experiments (Sup. Table 1)"
old_exp = per[col].copy()
per[col] = t1["human"].groupby("TF").size().reindex(TFS)
per["Mouse ChIP-Seq experiments (Sup. Table 1)"] = t1["mouse"].groupby("TF").size().reindex(TFS)
if per[[col, "Mouse ChIP-Seq experiments (Sup. Table 1)"]].isna().any().any():
    raise SystemExit("CHECK FAILED Sup. Table 2: a TF without experiments")
# no pooled peak count (Methods describe the pooling without a count)
per = per.drop(columns=["Human pooled peak set (peaks, all chromosomes)", "Mouse pooled peak set (peaks, all chromosomes)"])
# first-submission count columns carry the old chromosome labels; they would survive next to the new ones
R1_POS = ["Human training positives (chr2-7, 9, 10, 13-20)", "Human test positives (chr1, 8, 21)",
          "Mouse training positives (chr2-7, 9, 10, 13-18)", "Mouse test positives (chr1, 8, 19)"]
per = per.drop(columns=R1_POS)
per[HTR], per[HTE] = ch.pos_train, ch.pos_test
per[MTR], per[MTE] = cm.pos_other_file_control, cm.pos_chr1_8_19_file_train
# negatives used by the models (hm_mm) and the realised negatives per positive
NEG = {c.replace(" positives ", " negatives "): (h, c) for c, h in ((HTR, "n_train_neg_H"), (HTE, "n_test_neg_H"), (MTR, "n_train_neg_M"), (MTE, "n_test_neg_M"))}
RATIO_COLS = []
for _c, (_h, _pc) in NEG.items():
    per[_c] = hm[_h]
for _c, (_h, _pc) in NEG.items():
    _r = _c.replace(" negatives ", " negatives per positive ")
    per[_r] = (hm[_h] / per[_pc]).round(1)
    RATIO_COLS.append(_r)
# GC: median GC (%) from gc_check.csv; the split is identified by the chromosome column, not by the file name
GCMAP = {(H.GC_POOL[sp], H.SPLIT_GC["%s_%s" % (sp, part)]): (lab, "training" if part == "train" else "test", H.SPLIT["%s_%s" % (sp, part)])
         for sp, sl, lab in (("H", "human", "Human"), ("M", "mouse", "Mouse")) for part in ("train", "test")}
gc = gc[gc.pool.isin([H.GC_POOL["H"], H.GC_POOL["M"]])]
per = per.drop(columns=[c for c in per.columns if "GC content" in c])
GC_COLS = []
for (pool, chrom), (lab, slab, chs) in GCMAP.items():
    g = gc[(gc.pool == pool) & (gc.chromosomes == chrom)].set_index("TF")
    PD.check_tfs(g.index, "gc_check %s %s" % (pool, chrom))
    g = g.loc[TFS]
    n_ref = {("Human", "training"): per[HTR], ("Human", "test"): per[HTE], ("Mouse", "training"): per[MTR], ("Mouse", "test"): per[MTE]}[(lab, slab)]
    if not (g.n_pos.values == n_ref.values).all():
        raise SystemExit("CHECK FAILED gc_check n_pos %s %s differs from the positives" % (pool, chrom))
    for cls, c in (("positives", "median_gc_pos"), ("negatives", "median_gc_neg")):
        name = "%s %s %s, median GC content (%%; %s)" % (lab, slab, cls, chs)
        per[name] = (100 * g[c]).round(2)
        GC_COLS.append(name)
    name = "%s %s, median GC difference positives minus negatives (percentage points, unrounded; %s)" % (lab, slab, chs)
    per[name] = 100 * g.median_gc_diff
    GC_COLS.append(name)
cnt = [c for c in per.columns if c not in GC_COLS + RATIO_COLS]
per[cnt] = per[cnt].astype(int)
print("Sup. Table 2 human experiments changed:", {tf: (int(old_exp[tf]), int(per.loc[tf, col])) for tf in TFS if int(old_exp[tf]) != int(per.loc[tf, col])})
summ = []
for st, fn in (("Total over %d TFs" % NTF, lambda v: int(v.sum())), ("Median per TF", lambda v: float(v.median())), ("Mean per TF", lambda v: round(float(v.mean()), 2))):
    row = {"Statistic": st}
    for c in per.columns:
        row[c] = np.nan if (c in GC_COLS + RATIO_COLS and st.startswith("Total")) else fn(per[c])
    summ.append(row)
summ = pd.DataFrame(summ)
per = per.reset_index()
print(summ.set_index("Statistic").T.to_string())
OUT2 = os.path.join(OUT, "Sup_Table_2_PWMs_and_datasets.xlsx")
if "2" in BUILD:
    w2.backup(OUT2)
    with pd.ExcelWriter(OUT2) as xw:
        summ.to_excel(xw, sheet_name="summary", index=False)
        per.to_excel(xw, sheet_name="per_TF", index=False)
    w2.provenance(OUT, "3_analysis/supplementary_tables/make_sup_tables_1_3_R2.py (Sup. Table 2)", [OUT2], list(IN2.values()) + [os.path.join(OUT, "Sup_Table_1_ChIP-Seq_experiments.xlsx")],
                  params="{pooled peak set: not shown; human positives: B14 rebuild (pos_train, pos_test); GC: median of 300 bp windows, negatives v3}")

t3 = pd.read_excel(os.path.join(PKGC, "Sup_Table_3_ArChIPelago_performance.xlsx"), sheet_name=None)
perf = t3["performance"].copy()
PD.check_tfs(perf.TF_name, "Sup. Table 3 performance")
tot = per.set_index("TF")["Human PWMs (mono + di)"]
m = perf["PWM type"] == "mono+di"
perf.loc[m, "Number of PWMs"] = perf.loc[m, "TF_name"].map(tot).astype(int)
print("Sup. Table 3 rows corrected:", int(m.sum()), "ANDR mono+di:", perf.loc[m & (perf.TF_name == "ANDR"), "Number of PWMs"].unique())
OUT3 = os.path.join(OUT, "Sup_Table_3_ArChIPelago_performance.xlsx")
# per-TF values of Fig. 4 (human test), Fig. S4 (mouse test) and Fig. S5 (number of PWMs), and the
# paired counts RF on all PWMs vs Slim / diChIPMunk (3_analysis/integration/fig4_paired_h4.py)
SD = os.path.join(R2, "3_analysis", "package", "Figures", "source_data")
IN3X = {"Fig4_human_test": os.path.join(SD, "Figure_4_source_data.csv"), "FigS4_mouse_test": os.path.join(SD, "Figure_S4_source_data.csv"),
        "Fig4_S4_paired_counts": os.path.join(R2, "4_visualisation/figures/output", "Fig4", "Figure_4_S4_paired_RF_vs_Slim.csv"),
        "FigS5_per_TF": os.path.join(SD, "Figure_S5_source_data.csv"), "FigS5_median_curves": os.path.join(SD, "Figure_S5_median_curves.csv")}
t3x = {sh: pd.read_csv(f) for sh, f in IN3X.items()}
for sh in ("Fig4_human_test", "FigS4_mouse_test", "FigS5_per_TF"):
    PD.check_tfs(t3x[sh].TF.unique(), "Sup. Table 3 %s" % sh)
for sh, f in (("Fig4_human_test", H.IN["fig4_src"]), ("FigS4_mouse_test", H.IN["figS4_src"])):
    if w2.md5(IN3X[sh]) != w2.md5(f):
        raise SystemExit("CHECK FAILED Sup. Table 3 %s: pkg source data differ from %s" % (sh, f))
if "3" in BUILD:
    w2.backup(OUT3)
    with pd.ExcelWriter(OUT3) as xw:
        for sh, d in t3.items():
            (perf if sh == "performance" else d).to_excel(xw, sheet_name=sh, index=False)
        for sh, d in t3x.items():
            d.to_excel(xw, sheet_name=sh, index=False)
    w2.provenance(OUT, "3_analysis/supplementary_tables/make_sup_tables_1_3_R2.py (Sup. Table 3)", [OUT3],
                  [os.path.join(PKGC, "Sup_Table_3_ArChIPelago_performance.xlsx"), OUT2] + list(IN3X.values()),
                  params="{mono+di rows: Number of PWMs = monoPWMs + diPWMs; Fig. 4 / S4 / S5 sheets: source data as plotted}")
if "4" in BUILD:
    t4 = pd.read_excel(os.path.join(PKGC, "Sup_Table_4_operational_metrics.xlsx"), sheet_name=None)
    PD.check_tfs(t4["per_TF"].TF, "Sup. Table 4 per_TF")
    if len(t4["per_TF"]) != 2 * NTF:
        raise SystemExit("CHECK FAILED Sup. Table 4 per_TF rows %d" % len(t4["per_TF"]))
    PREC_GAIN_DECIMALS = 4
    op4 = pd.read_csv(H.IN["operational"])
    pt = t4["per_TF"]
    key = pt.TF + "|" + pt["Test set"].str[0].str.upper()          # 'human (chr...)' -> H, 'mouse (chr...)' -> M
    opk = op4.set_index(op4.tf + "|" + op4.test)
    if sorted(key) != sorted(opk.index):
        raise SystemExit("CHECK FAILED Sup. Table 4 per_TF rows vs operational_metrics.csv")
    PG = "Precision gain at 50 % recall, ArChIPelago minus best monoPWM"
    pg = opk.loc[key, "precision_gain_recall0.5_vs_mono"].values
    if not np.allclose(pt[PG].values, pg, atol=5e-4 + 1e-9):   # the pkg column is the same value rounded to 3 decimals
        raise SystemExit("CHECK FAILED Sup. Table 4 %s differs from operational_metrics.csv by %.2g" % (PG, np.abs(pt[PG].values - pg).max()))
    pt[PG] = np.round(pg, PREC_GAIN_DECIMALS)                    
    # number of TFs behind each summary median (the quantity is undefined where both models have no false positive)
    QCOL = {"FP reduction at 50 % recall (%)": "fp_reduction_recall0.5_pct_vs_mono", "FP reduction at 80 % recall (%)": "fp_reduction_recall0.8_pct_vs_mono",
            "additional true sites in top 1,000": "extra_tp_top1000_vs_mono", "precision gain at 50 % recall (paired, per TF)": "precision_gain_recall0.5_vs_mono"}
    sm = t4["summary"]
    MOLD = "Median over %d TFs" % NTF
    nmed, rmed = [], []
    for _, x in sm.iterrows():
        v = op4[op4.test == x["Test set"][0].upper()][QCOL[x.Quantity]].dropna()
        nmed.append(len(v)); rmed.append(v.median())
    tol = np.where(sm.Quantity.str.startswith("precision"), 5e-4, 5e-2) + 1e-9   # summary rounds precision to 3, the rest to 1 decimal
    if not (np.abs(sm[MOLD].values - np.array(rmed)) <= tol).all():
        raise SystemExit("CHECK FAILED Sup. Table 4 summary medians differ from operational_metrics.csv: %s vs %s" % (list(sm[MOLD]), rmed))
    sm = sm.rename(columns={MOLD: "Median over TFs"})
    sm.insert(sm.columns.get_loc("Median over TFs") + 1, "TFs in the median", nmed)
    # exact two-sided Wilcoxon signed-rank P of the paired per-TF values (TFs with equal values left out),
    # so the Table 1 false-positive row has its test in this sheet (same values and test as make_numbers_w2.py)
    from scipy.stats import wilcoxon
    PAIR = {"FP reduction at 50 % recall (%)": "FP per 10,000 regions at 50 % recall: {}",
            "FP reduction at 80 % recall (%)": "FP per 10,000 regions at 80 % recall: {}",
            "additional true sites in top 1,000": "True sites among top 1,000 predictions: {}"}
    wp, wn = [], []
    for _, x in sm.iterrows():
        sub = pt[pt["Test set"].str.split().str[0] == x["Test set"].split()[0]]
        if len(sub) != NTF:
            raise SystemExit("CHECK FAILED Sup. Table 4 per_TF rows for %s: %d" % (x["Test set"], len(sub)))
        d = (sub[PG] if x.Quantity.startswith("precision") else
             sub[PAIR[x.Quantity].format("ArChIPelago RF")] - sub[PAIR[x.Quantity].format("best monoPWM")]).values
        nz = d[d != 0]
        wp.append(float("%.2g" % wilcoxon(nz, method="exact").pvalue)); wn.append(len(nz))
    sm["Wilcoxon P (exact, two-sided; paired per TF, equal values left out)"] = wp
    sm["TFs in the test (non-zero differences)"] = wn
    t4["summary"] = sm
    print("Sup. Table 4: TFs in the median %s; precision gain per TF with %d decimals (mouse PRGR %+.4f)" % (
        dict(zip(sm["Test set"].str[:5] + " " + sm.Quantity, nmed)), PREC_GAIN_DECIMALS, opk.loc["PRGR|M", "precision_gain_recall0.5_vs_mono"]))
    OUT4 = os.path.join(OUT, "Sup_Table_4_operational_metrics.xlsx")
    w2.backup(OUT4)
    with pd.ExcelWriter(OUT4) as xw:
        for sh, d in t4.items():
            d.to_excel(xw, sheet_name=sh, index=False)
    w2.provenance(OUT, "3_analysis/supplementary_tables/make_sup_tables_1_3_R2.py (Sup. Table 4)", [OUT4], [os.path.join(PKGC, "Sup_Table_4_operational_metrics.xlsx"), H.IN["operational"]],
                  params="{precision gain decimals: %d; summary: TFs in the median (V188); exact two-sided Wilcoxon P on non-zero paired differences (G7-3)}" % PREC_GAIN_DECIMALS)
if "6" in BUILD:   # RF rows and single-PWM scans rerun on the final training matrices
    import sup_table_6_h3 as S6
    S6.build(os.path.join(OUT, "Sup_Table_6_runtime_memory.xlsx"))
