"""Builds outputs/numbers/NUMBERS_R2_35TF.csv = the human rows of outputs/numbers/NUMBERS_human_35TF.csv (unchanged except the
columns mapped and the rows that the rows below supersede, marked in statement_note) + the mouse rebuild, dataset counts, split
leakage, GC check, Fig. S4, operational, Table 1 mouse, cross-species, drivers and Fig. 1 rows. Every new value is read from a
rebuilt table here; old values are the quoted text of the inventory. Columns: item, old text value, new value, source file,
inventory row, statement_note.
Full regeneration on the final data. Inputs come from h3_inputs.py; old values of computed rows come from the frozen copies under
results/previous_values/<git ref>/ or from the released package (h3_old.py); notes are computed from the tables, not typed; rows
of Fig. 4, S4 and S12 read 'pending' with a statement_note starting 'pending' while their source tables are missing.
usage: python3 3_analysis/integration/make_numbers_w2.py"""
import os, re, ast, sys
import numpy as np, pandas as pd
from scipy.stats import wilcoxon
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w2
import h3_inputs as H
import h3_old as O
PD, R2, B13 = w2.PD, w2.R2, w2.B13
J = lambda *p: os.path.join(R2, *p)
SRC = {k: J(*v.split("/")) for k, v in {
    "b12": "outputs/numbers/NUMBERS_human_35TF.csv", "t1": "outputs/supplementary_tables/Sup_Table_1_ChIP-Seq_experiments.xlsx",
    "t2": "outputs/supplementary_tables/Sup_Table_2_PWMs_and_datasets.xlsx", "t3": "outputs/supplementary_tables/Sup_Table_3_ArChIPelago_performance.xlsx",
    "t4": "outputs/supplementary_tables/Sup_Table_4_operational_metrics.xlsx", "t5": "outputs/supplementary_tables/Sup_Table_5_cross_species_and_mouse_trained.xlsx",
    "hm": os.path.relpath(H.IN["hm_mm"], w2.R2), "gc": os.path.relpath(H.IN["gc_v3"], w2.R2), "cm": os.path.relpath(H.IN["counts_M"], w2.R2),
    "ch": os.path.relpath(H.IN["counts_H"], w2.R2), "b8": "results/datasets/split_leakage/b8w1g2_summary.csv",
    "b8top": "3_analysis/integration/top1000_check_h3.csv", "s4": "results/datasets/round1_fig_s4/Figure_S4_numbers.csv",
    "s4n": "4_visualisation/figures/output/Fig4/Figure_S4_stats.csv", "tok": "4_visualisation/figures/output/Fig4/text_tokens_h3b.csv",
    "negst": os.path.relpath(H.IN["neg_v3_stats"], w2.R2), "negck": os.path.relpath(H.IN["neg_v3_checks"], w2.R2),
    "bp": os.path.relpath(H.IN["build_positives"], w2.R2), "abl": "3_analysis/feature_importance/tables/B4_ablation_summary.csv",
    "t10": "outputs/supplementary_tables/Sup_Table_10_transfer_drivers_and_headroom.xlsx", "plan": "results/datasets/mouse_rebuild/mouse_rebuild_cpu_time.csv",
    "widths": "3_analysis/top_peaks/b3_peak_widths_summary.csv", "topn": "3_analysis/top_peaks/b3_topN_summary.csv",
    "rho": "3_analysis/top_peaks/b3_gain_vs_npos_spearman.csv", "b2": "3_analysis/low_fpr_metrics/B2_summary.csv",
    "b2c": "3_analysis/low_fpr_metrics/B2_per_chromosome_summary.csv", "drv": "3_analysis/transfer_drivers/B5_A_spearman_drivers.csv",
    "head": "3_analysis/transfer_drivers/B5_B_ceiling_spearman.csv", "clouds": "3_analysis/transfer_drivers/B5_B_cloud_members.txt",
    "cs": "3_analysis/package/analysis/cross_species/cross_species_table.csv", "tab1": "3_analysis/integration/table1_rows.csv",
    "fig1": "3_analysis/integration/figure1_values.csv", "sat": "3_analysis/package/analysis/saturation/saturation_summary_by_k.csv",
    "mask": "results/datasets/repeat_mask/repeat_mask_counts.csv", "subrank": "4_visualisation/figures/FigS7/tf_ranking_h3.csv",
    "subverd": "4_visualisation/figures/FigS7/subtype_verdicts_v11_h3.csv",
    "v11": "4_visualisation/figures/FigS7/hocomoco11_check/subtype_vs_hocomoco11.csv", "cv11": "4_visualisation/figures/FigS7/hocomoco11_check/cluster_vs_v11.csv", "subsrc": "4_visualisation/figures/FigS7/Figure_S7_source_data.csv",
    "v11n": "4_visualisation/figures/FigS7/hocomoco11_check/subtype_vs_named_v11.csv",
    "f4pair": "4_visualisation/figures/output/Fig4/Figure_4_S4_paired_RF_vs_Slim.csv",
    "gate": "3_analysis/low_fpr_metrics/B2_gate_h3_summary.csv", "porig": "3_analysis/integration/pwm_origin_h3.csv",
    "mvd": "3_analysis/integration/mono_vs_monodi_h3.csv", "gct": "3_analysis/integration/gc_per_tf_h3.csv", "gcs": "3_analysis/integration/gc_vs_gain_h3.csv",
    "ties": "3_analysis/integration/positive_count_ties_h3.csv", "abf": "3_analysis/feature_importance/tables/B4_ablation_all_fits.csv",
    "rt": "3_analysis/package/results_table.csv", "b2w": "3_analysis/low_fpr_metrics/B2_metrics_wide.csv", "s6": "3_analysis/integration/sup6_methods_values_h3.csv",
    "t6": "outputs/supplementary_tables/Sup_Table_6_runtime_memory.xlsx", "hpt": "3_analysis/transfer_drivers/B5_B_headroom_per_TF.csv",
    "conc": "3_analysis/feature_importance/tables/B4_importance_concentration.csv", "s1p": "4_visualisation/figures/output/S1/Figure_S1_plotdata.csv",
    "f3c": "4_visualisation/figures/Fig3_S3/Figure_3_C_plotdata.csv", "b5t": "3_analysis/transfer_drivers/B5_per_TF_table.csv",
    "f4src": "3_analysis/package/Figures/source_data/Figure_4_source_data.csv"}.items()}
H4_REF = "10feef6"    # frozen copy results/previous_values/10feef6 of the earlier NUMBERS_R2_35TF.csv (git commit 10feef6); its values fill the h4 column
H3A_REF = "ffcef4d"   # frozen copy results/previous_values/ffcef4d of the earlier NUMBERS_R2_35TF.csv (git commit ffcef4d); its values fill the h3a column
rows = []
def add(item, old, new, src, inv, note=""):
    rows.append(dict(item=item, **{"old text value": old, "new value": new, "source file": src, "inventory row": inv, "statement_note": note}))
r = lambda x, k=3: ("%." + str(k) + "f") % x
pc = lambda x, k=1: ("%." + str(k) + "f %%") % (100 * x)
e = lambda p: "%.0e" % p if p < 1e-3 else "%.2g" % p
rel = lambda k: os.path.relpath(SRC[k], R2)

t1 = pd.read_excel(SRC["t1"], sheet_name=None); t2 = pd.read_excel(SRC["t2"], sheet_name="per_TF").set_index("TF")
t2s = pd.read_excel(SRC["t2"], sheet_name="summary").set_index("Statistic")
t3 = pd.read_excel(SRC["t3"], sheet_name="summary"); t3p = pd.read_excel(SRC["t3"], sheet_name="performance")
t4 = pd.read_excel(SRC["t4"], sheet_name="summary"); t5 = pd.read_excel(SRC["t5"], sheet_name="per_TF").set_index("TF")
t5s = pd.read_excel(SRC["t5"], sheet_name="summary")
hm = pd.read_csv(SRC["hm"]).set_index("TF"); NTF = PD.check_tfs(hm.index, "hm_mm_results")
cm = pd.read_csv(SRC["cm"]).set_index("TF"); ch = pd.read_csv(SRC["ch"]).set_index("TF")
def _num(v):
    try:
        return float(v)
    except ValueError:
        return v
f1 = pd.read_csv(SRC["fig1"]).set_index("quantity").value.map(_num)

# ---- datasets, TFs
nH, nM = len(t1["human"]), len(t1["mouse"])
add("ChIP-Seq datasets used, human + mouse", "685", "%d (%d + %d)" % (nH + nM, nH, nM), rel("t1"), "M01, M02, M15, M16, L06, R01, F01")
add("ChIP-Seq datasets, human", "395", str(nH), rel("t1"), "M03, M16, L01, F01")
add("ChIP-Seq datasets, mouse", "290", str(nM), rel("t1"), "M03, M16, L01, F01")
add("orthologous TF pairs", "36", str(NTF), rel("t1"), "M01, M02, M06, M15, M16, F01, L02-L04, L08, L10, L11, R04", "P53 removed: no mouse P53 peak set passes the dataset rule")
add("minimum datasets per TF, mouse", "at least three per TF in mouse", "%d (%s)" % (f1["min datasets per TF, mouse"], f1["TFs at the mouse minimum"]), rel("t1"), "M03, F01",
    "statement false as written; true statement: at least %d datasets per TF in mouse (%s) and at least %d in human (%s)" % (
    f1["min datasets per TF, mouse"], f1["TFs at the mouse minimum"].replace(",", ", "), f1["min datasets per TF, human"], f1["TFs at the human minimum"].replace(",", ", ")))
add("minimum datasets per TF, human", "-", "%d (%s)" % (f1["min datasets per TF, human"], f1["TFs at the human minimum"]), rel("t1"), "M03, F01")
add("Sup. Table 1 rows (sets used)", "414 human / 290 mouse listed (395 human with a peak file)", "%d / %d, exactly the sets used" % (nH, nM), rel("t1"), "L01, R01, T01, T02")
# PWM origin: producer 3_analysis/integration/pwm_origin_h3.py over the four Zenodo PWM archives and Sup. Table 1
_po = pd.read_csv(SRC["porig"]); _pg = _po.groupby("species")[["n_PWMs", "n_from_listed_sets", "n_source_sets"]].sum()
_pwm_old = "802 human and 652 mouse datasets"   # old text value (quote of Methods 2.1)
add("PWM source datasets (Methods 2.1)", _pwm_old, "%d human and %d mouse datasets" % (_pg.loc["human", "n_source_sets"], _pg.loc["mouse", "n_source_sets"]),
    rel("porig"), "M06, D9 (PLAN_20261106_mouse_symmetric.md)",
    "statement changes: %s counts the PWMs of all TFs in the archives including the excluded TF; for the %d TFs: %d human and %d mouse" % (
        _pwm_old, _po.TF.nunique(), _pg.loc["human", "n_source_sets"], _pg.loc["mouse", "n_source_sets"]))
add("PWMs derived from one of the TF's listed positive sets (same species; D9)", "-",
    "; ".join("%s %s of %s (%.1f %%)" % (sp_, format(int(_pg.loc[sp_, "n_from_listed_sets"]), ","), format(int(_pg.loc[sp_, "n_PWMs"]), ","),
                                         100 * _pg.loc[sp_, "n_from_listed_sets"] / _pg.loc[sp_, "n_PWMs"]) for sp_ in ("human", "mouse")),
    rel("porig"), "D9 (PLAN_20261106_mouse_symmetric.md), R1-F3 (REVIEW_FINDINGS_20261107.md)",
    "Methods (D9): disclose that these PWMs come from peak sets that also give the positives; per TF %.1f-%.0f %% (human), %.1f-%.0f %% (mouse)" % tuple(
        100 * v for sp_ in ("human", "mouse") for v in (_po[_po.species == sp_].frac_from_listed_sets.min(), _po[_po.species == sp_].frac_from_listed_sets.max())))
# ---- PWMs
for sp in ("Human", "Mouse"):
    old = {"Human": ("1,495", "780", "2,275", "34.5"), "Mouse": ("1,192", "624", "1,816", "37")}[sp]
    for i, c in enumerate(("%s monoPWMs" % sp, "%s diPWMs" % sp, "%s PWMs (mono + di)" % sp)):
        add("Sup. Table 2 total, " + c, old[i], format(int(t2[c].sum()), ","), rel("t2"), "M15, M16, F01, R04")
    add("median PWMs per TF, " + sp.lower(), old[3], "%g" % t2["%s PWMs (mono + di)" % sp].median(), rel("t2"), "M06, M15, F01, R04")
add("minimum monoPWMs / diPWMs per TF, both species", ">= 3 monoPWMs and >= 3 diPWMs", "%d / %d" % (f1["min monoPWMs per TF (both species)"], f1["min diPWMs per TF (both species)"]), rel("t2"), "F01, R04")
# ---- pooled peaks
add("median pooled peak set per TF and species", "14,099", "removed; Methods describes the pooling without a count", "-", "M07, F01, R08", "count removed from Methods, Fig. 1, letter R1-L6 and Sup. Table 2 (authors' decision)")
add("median pooled peak set per TF, human", "9,641 (Sup. Table 2)", "removed; Methods describes the pooling without a count", "-", "M07, T03",
    "Sup. Table 2 column and summary rows removed")
add("median pooled peak set per TF, mouse", "18,654", "removed; Methods describes the pooling without a count", "-", "M07, T03", "Sup. Table 2 column and summary rows removed")
# ---- positives, cap, split
add("human training / test positives", "229,081 / 88,384", "%s / %s" % (format(int(f1["human training positives"]), ","), format(int(f1["human test positives"]), ",")), rel("t2"), "F01, T03", "human positives of the rebuild (D6): each peak set read once; at most 1,000 peaks per dataset (%d of %d TFs at or below 1,000 x datasets)" % (
    int(((t2[H.COL["HTR"]] + t2[H.COL["HTE"]]) <= 1000 * t2["Human ChIP-Seq experiments (Sup. Table 1)"]).sum()), NTF))
add("mouse training / test positives", "301,994 / 144,779", "%s / %s" % (format(int(f1["mouse training positives"]), ","), format(int(f1["mouse test positives"]), ",")), rel("t2"), "F01, T03")
add("median mouse test positives per TF", "2,771.5", "%g" % t2["Mouse test positives (chr1, 8, 19)"].median(), rel("t2"), "T03")
add("smallest mouse test sets", "-", ", ".join("%s %d" % (k, v) for k, v in t2["Mouse test positives (chr1, 8, 19)"].nsmallest(5).items()), rel("t2"), "M09")
add("training sets at the 10,000 cap", "12 human and 24 mouse", "%d human and %d mouse" % ((hm.n_train_pos_H == 10000).sum(), (hm.n_train_pos_M == 10000).sum()), rel("hm"), "M08, R12")
add("test:train ratio of positives", "approximately 30:70", "%.0f:%.0f human, %.0f:%.0f mouse" % (f1["test share of positives, human (%)"], 100 - f1["test share of positives, human (%)"],
    f1["test share of positives, mouse (%)"], 100 - f1["test share of positives, mouse (%)"]), rel("t2"), "M12, F01", "Methods: give the realised test:train shares per species (%.0f:%.0f human, %.0f:%.0f mouse) instead of 'about 30:70'" % (
    f1["test share of positives, human (%)"], 100 - f1["test share of positives, human (%)"], f1["test share of positives, mouse (%)"], 100 - f1["test share of positives, mouse (%)"]))
q = hm.n_test_neg_M / hm.n_test_pos_M
add("mouse test class ratio, median negatives per positive (range)", "1:44, range 1:14 to 1:100", "1:%.0f, range 1:%.0f to 1:%.0f" % (q.median(), q.min(), q.max()), rel("hm"), "M09")
# ---- masking, negatives, GC check
_mk = pd.read_csv(SRC["mask"], dtype=str).set_index("quantity")["value"]
_mH, _tH, _mM, _tM = [int(_mk[k].replace(",", "")) for k in ("removed_human", "total_human", "removed_mouse", "total_mouse")]
add("mouse repeat mask", "not applied (1-line placeholder)", "UCSC mm10 RepeatMasker, f = 0.7; %s of %s mouse top-1,000 windows removed (%.1f %%)" % (format(_mM, ","), format(_tM, ","), 100 * _mM / _tM),
    rel("mask"), "M04", "statement now true for mouse: name the hg38 and mm10 UCSC rmsk tracks")
add("human windows removed by the mask (B8 rerun)", "-", "%s of %s (%.1f %%)" % (format(_mH, ","), format(_tH, ","), 100 * _mH / _tH), rel("mask"), "M04")
add("negative pool, both species", "TFs of other structural classes (39 human, 41 mouse TFs pooled)", "the %d modelled TFs in both species (D5, D6), other TFClass classes only (D8)" % NTF,
    rel("negst"), "M10, R18", "statement changes: one pool rule for both species; exclusion is by TFClass class, not family (D8)")
# negatives v3: own-positive overlaps removed; realised negatives per positive
nc = pd.read_csv(SRC["negck"]); ns = pd.read_csv(SRC["negst"])
SPL = {("HUMAN", "train"): "human training", ("HUMAN", "control"): "human test", ("MOUSE", "control"): "mouse training", ("MOUSE", "train"): "mouse test"}
for (sp_, fs), lab in SPL.items():
    g = nc[(nc.species == sp_) & (nc.file_split == fs)]
    v2, v3 = g[g.version == "v2"], g[g.version == "v3"]
    # the v3 lists were built directly (no v2 lists): the removed count does not exist there, so the
    # row reports the v3 negatives and their own-positive overlaps (0 by construction) instead of "0 of 0 removed"
    _rem = ("%s of %s removed" % (format(int(v2.overlap_own_pos.sum()), ","), format(int(v2.n.sum()), ","))) if len(v2) else \
        ("%s overlaps with own positives (H5 builds the final lists directly)" % format(int(v3.overlap_own_pos.sum()), ","))
    add("negatives v3, %s: v2 negatives within 300 bp of the TF's own positives (removed); v3 negatives" % lab, "-",
        "%s; %s negatives" % (_rem, format(int(v3.n.sum()), ",")), rel("negck"), "M10, R1-F2, D7",
        "statement changes (D7): negatives overlapping the TF's own positives of the same split are removed before GC matching")
    r = ns[(ns.species == sp_) & (ns.split == fs)].neg_per_pos
    add("realised negatives per positive, %s (min / median / max)" % lab, "1:100 (design)", "%.1f / %.1f / %.1f" % (r.min(), r.median(), r.max()), rel("negst"), "M09, M10, R1-F7",
        "Methods: BiasAway is asked for 100 negatives per positive; the realised ratio is lower where the pool is small")
for sp_, lab in (("H", "human"), ("M", "mouse")):
    q = pd.concat([hm["n_train_neg_" + sp_] / hm["n_train_pos_" + sp_], hm["n_test_neg_" + sp_] / hm["n_test_pos_" + sp_]])
    add("realised negatives per positive used by the models, %s (range over TFs and splits)" % lab, "1:100", "1:%.1f to 1:%.1f" % (q.min(), q.max()), rel("hm"), "M09, R1-F7",
        "Methods: BiasAway is asked for 1:100; the realised ratio of the training and test sets ranges from 1:%.1f to 1:%.1f" % (q.min(), q.max()))
# training chromosome lists as the code builds them: the "11" "12" literal concatenates to "1112", so chr11 and chr12 are unused
_bp = open(SRC["bp"]).read()
if H.STAGE == "H5":   # the split is defined once in h5_splits.LOGICAL (every autosome in one set; sex chromosomes in neither)
    _ns = {}
    exec(compile(_bp, SRC["bp"], "exec"), _ns)
    for SP, lab, n_auto in (("HUMAN", "human", 22), ("MOUSE", "mouse", 19)):
        tr, test = _ns["LOGICAL"][SP]["train"], _ns["LOGICAL"][SP]["test"]
        _pp = _ns["partition_problems"](SP, tr, test)
        if _pp:
            raise SystemExit("CHECK FAILED split %s: %s" % (SP, _pp))
        _k = lambda c: int(c[3:])
        add("chromosomes used, %s (R1-F6)" % lab, "training: all autosomes except the test chromosomes", "training %s; test %s; unused %s" % (
            ", ".join(c[3:] for c in sorted(tr, key=_k)), ", ".join(c[3:] for c in sorted(test, key=_k)), "X, Y"), rel("bp"), "M12, R1-F6, D16",
            "Methods: name the training chromosomes; every autosome is in exactly one set, the sex chromosomes in neither")
_split_rows = () if H.STAGE == "H5" else (("HUMAN", "human", 22, "TEST_CHR"), ("MOUSE", "mouse", 19, "TRAIN_CHR"))
for SP, lab, n_auto, test_on_disk in _split_rows:
    blk = _bp[_bp.index('"%s": dict(' % SP):]
    lists = {k: ["chr" + x for x in ast.literal_eval(re.search(k + r'=\["chr" \+ x for x in (\[[^\]]+\])\]', blk).group(1))]
             for k in ("TRAIN_CHR", "TEST_CHR")}
    test = lists[test_on_disk]; train = lists["TEST_CHR" if test_on_disk == "TRAIN_CHR" else "TRAIN_CHR"]
    valid = ["chr%d" % i for i in range(1, n_auto + 1)]
    tr = [c for c in train if c in valid]; bad = [c for c in train if c not in valid]
    unused = [c for c in valid if c not in tr and c not in test] + ["chrX", "chrY"]
    add("chromosomes used, %s (R1-F6)" % lab, "training: all autosomes except the test chromosomes", "training %s; test %s; unused %s (code list entry %s)" % (
        ", ".join(c[3:] for c in tr), ", ".join(c[3:] for c in test), ", ".join(c[3:] for c in unused), ",".join(bad)), rel("bp"), "M12, R1-F6",
        "Methods: name the training chromosomes; chr11 and chr12 are not used (concatenated list entry %s)" % ",".join(bad))
gc = pd.read_csv(SRC["gc"])
for pool, _gp in (("human_sym_v3", H.GC_POOL["H"]), ("mouse_sym_v3", H.GC_POOL["M"])):   # item: logical pool name; rows: pool of the final run
    for chrom, g in gc[gc.pool == _gp].groupby("chromosomes"):   # negatives v3
        add("GC check %s, %s: median over TFs of (median GC positives - median GC negatives); median KS D; TFs with |diff| > 0.05" % (pool, chrom), "-",
            "%+.4f; %.3f; %d of %d" % (g.median_gc_diff.median(), g.ks_D.median(), (g.median_gc_diff.abs() > 0.05).sum(), len(g)), rel("gc"), "M10, D5")
    g = gc[gc.pool == _gp]
    add("GC check %s: largest |GC difference| (TF, split)" % pool, "-", "%+.4f (%s, %s)" % (g.loc[g.median_gc_diff.abs().idxmax(), "median_gc_diff"], g.loc[g.median_gc_diff.abs().idxmax(), "TF"],
        g.loc[g.median_gc_diff.abs().idxmax(), "chromosomes"]), rel("gc"), "M10, D5")
# ---- peak widths
w = pd.read_csv(SRC["widths"]).set_index("species")
add("median of dataset median peak widths, human / mouse (bp)", "575 / 551", "%g / %g" % (w.loc["H", "width_median_of_dataset_medians"], w.loc["M", "width_median_of_dataset_medians"]), rel("widths"), "M05, R02, R09")
add("median width of the top 1,000 peaks, human / mouse (bp)", "898 / 813", "%g / %g" % (w.loc["H", "width_median_top1000_median"], w.loc["M", "width_median_top1000_median"]), rel("widths"), "M05")
add("top 1,000 vs other peaks, both species (bp; datasets)", "858 vs 553 bp, 683 datasets", "%g vs %g bp, %d datasets" % (w.loc["all", "width_median_top1000_median"], w.loc["all", "width_median_rest_median"], w.loc["all", "n_datasets_ge1000_peaks"]),
    rel("widths"), "M32, R02, R09")
add("peak-width datasets", "685", str(int(w.loc["all", "n_datasets"])), rel("widths"), "L06, R01, R02")
# ---- top 1,000 check, split leakage
tp = pd.read_csv(SRC["b8top"])
for sp, lab in (("H", "human"), ("M", "mouse")):
    g = tp[tp.species == sp]
    add("%s positives (training and test) among the 1,000 most significant peaks of a listed set" % lab,
        "98.9 % (test)" if sp == "H" else "63 % found in a peak file, 99.8 % of these in the top 1,000",
        "%s (%d of %d)" % (pc(g.n_in_top1000.sum() / g.n_positives.sum()), g.n_in_top1000.sum(), g.n_positives.sum()), rel("b8top"), "L05, R02, R08, T10",
        "check C of the rebuild covers every positive of both splits (was: test positives of the published tables)")
b8 = pd.read_csv(SRC["b8"]).set_index(["version", "species"])
for sp, lab, old, oldiqr in (("H", "human", "63 %", "34-71 %"), ("M", "mouse", "57 %", "47-71 %")):
    s = b8.loc[("masked", sp)]
    add("random split: median share of %s test windows overlapping a training window (masked)" % lab, old, pc(s.median_overlap), rel("b8"), "M11, R01, R07, R10")
    add("random split: IQR, %s" % lab, oldiqr, "%.0f-%.0f %%" % (100 * s.q25_overlap, 100 * s.q75_overlap), rel("b8"), "R07")
    add("random split: maximum, %s" % lab, "up to 95 % for REST", "%s (%s)" % (pc(s.max_overlap), s.max_TF), rel("b8"), "R07")
    add("random split: overlaps from another dataset, pooled, %s" % lab, "more than 99.9 %", pc(s.share_other_ds_pooled, 2), rel("b8"), "R07")
    add("chromosome split overlaps, %s" % lab, "none, by construction", str(int(s.chromosome_split_overlaps)), rel("b8"), "R01, R07")
add("Sup. Table 7 random_split_overlap: repeat removal", "repeat removal not applied", "repeat mask applied in both species (pipeline rule)", rel("b8"), "L05", "legend must say the mask is applied")
# ---- human-trained on mouse test (Sup. Tables 3 and 5)
# the H>M rows give the mean of seeds 0 and 1, as the mouse-trained rows do (Sup. Table 5 columns and
# summary row '(mean of 2 seeds)'); the seed 0 values (Sup. Table 3 summary, Sup. Table 7) go into the statement note
HM_SEED = "mean of 2 seeds"
_t5r = lambda lab: t5s[t5s.iloc[:, 0] == lab].iloc[0]   # Sup. Table 5 summary row by its label
_HML = "Human-trained model on the mouse test set vs the best monoPWM (selected on human training data), "
hmr, hm0 = _t5r(_HML + HM_SEED), _t5r(_HML + "seed 0")
mt = t3[t3["Test set"].str.startswith("mouse")].iloc[0]   # seed 0
_s0 = "seed 0 (Sup. Table 3, Sup. Table 7): "
add("H>M RF median auROC / best monoPWM", "0.8745 / 0.8613", "%.4f / %.4f" % (hm.HM_rf_test_M_auroc_mean.median(), hm.HM_broc_test_M_auroc.median()),
    rel("hm") + "; " + rel("t5"), "T05, M28", "mean of 2 seeds (V006); " + _s0 + "%.4f / %.4f" % (mt["ArChIPelago auROC (median)"], mt["baseline auROC (median)"]))
add("H>M RF median auPRC / best monoPWM", "0.2481 / 0.1833", "%.4f / %.4f" % (hm.HM_rf_test_M_auprc_mean.median(), hm.HM_bprc_test_M_auprc.median()),
    rel("hm") + "; " + rel("t5"), "T05, M28", "mean of 2 seeds (V006); " + _s0 + "%.4f / %.4f" % (mt["ArChIPelago auPRC (median)"], mt["baseline auPRC (median)"]))
add("H>M median per-TF gain auROC / auPRC", "0.015 / 0.022", "%.4f / %.4f" % (hmr["median dauROC"], hmr["median dauPRC"]), rel("t5"), "M28, T05",
    "mean of 2 seeds (V006); " + _s0 + "%.4f / %.4f" % (hm0["median dauROC"], hm0["median dauPRC"]))
add("H>M TFs improved, auROC", "35 of 36", "%d of %d" % (hmr["TFs with dauROC > 0"], NTF), rel("t5"), "M28, T05",
    "mean of 2 seeds (V006); %s%d of %d; two-sided Wilcoxon P %.1e (mean) / %.1e (seed 0)" % (_s0, hm0["TFs with dauROC > 0"], NTF,
                                                                                          float(hmr["Wilcoxon P (auROC)"]), float(hm0["Wilcoxon P (auROC)"])))
add("H>M TFs improved, auPRC", "34 of 36", "%d of %d" % (hmr["TFs with dauPRC > 0"], NTF), rel("t5"), "M28, T05",
    "mean of 2 seeds (V006); %s%d of %d; two-sided Wilcoxon P %.1e (mean) / %.1e (seed 0)" % (_s0, hm0["TFs with dauPRC > 0"], NTF,
                                                                                          float(hmr["Wilcoxon P (auPRC)"]), float(hm0["Wilcoxon P (auPRC)"])))
dr = t5["H>M: dauROC (%s)" % HM_SEED]; dp = t5["H>M: dauPRC (%s)" % HM_SEED]; below = sorted(set(dr.index[dr < 0]) | set(dp.index[dp < 0]))
dr0, dp0 = t5["H>M: dauROC (seed 0)"], t5["H>M: dauPRC (seed 0)"]
mmr, mmp = t5["M>M: dauROC"], t5["M>M: dauPRC"]
add("H>M TFs below the baseline on >= 1 metric", "ERG, P53 and STA5A (by up to 0.02)", "%s (by up to %.1g)" % (" and ".join(below), max(-dr.min(), -dp.min())), rel("t5"), "M28, L09, L12, R16, T05",
    "true statement: only %s fell below the baseline on one metric (%s), by up to %.3f; both exceeded the baseline when trained on mouse data (%s)" % (
        " and ".join(below), "; ".join("%s dauROC %+.4f, dauPRC %+.4f" % (t, dr[t], dp[t]) for t in below), max(-dr.min(), -dp.min()),
        "; ".join("%s %+.4f / %+.4f" % (t, mmr[t], mmp[t]) for t in below)) + "; mean of 2 seeds (V006); seed 0: %s" % (
        "; ".join("%s %+.4f / %+.4f" % (t, dr0[t], dp0[t]) for t in below)) if all((mmr[t] > 0) and (mmp[t] > 0) for t in below) else "CHECK: not all recover")
mms = _t5r("Mouse-trained model on the mouse test set vs the best mouse monoPWM (selected on mouse training data)")
dd = _t5r("Mouse-trained minus human-trained model on the same mouse test set (positive = the mouse-trained model is better)")
add("M>M median gain auROC / auPRC; TFs improved", "0.0224 / 0.0357; 35 / 36 of 36", "%.4f / %.4f; %d / %d of %d" % (mms["median dauROC"], mms["median dauPRC"], mms["TFs with dauROC > 0"], mms["TFs with dauPRC > 0"], NTF), rel("t5"), "T07, R06")
_mmb = [t for t in str(mms["TFs below baseline on >= 1 metric"]).split(", ") if t != "none"]   # none below
add("M>M TFs below the baseline", "USF2", mms["TFs below baseline on >= 1 metric"], rel("t5"), "T07",
    "; ".join("mouse-trained %s: dauROC %+.4f, dauPRC %+.4f" % (t, mmr[t], mmp[t]) for t in _mmb) or "no TF below the baseline on either metric")
add("M>M minus H>M median auROC / auPRC; TFs > 0", "0.0086 / 0.0255; 31 / 32", "%.4f / %.4f; %d / %d" % (dd["median dauROC"], dd["median dauPRC"], dd["TFs with dauROC > 0"], dd["TFs with dauPRC > 0"]), rel("t5"), "T07")
add("M>M minus H>M TFs below 0", "ANDR, E2F4, MAX, REST, SRF, USF2", dd["TFs below baseline on >= 1 metric"], rel("t5"), "T07")
sim = t5["Similarity of the baseline human monoPWM to the nearest mouse monoPWM (Pearson r of aligned columns)"]
_simo = pd.read_csv(os.path.join(H.PKG, "analysis/cross_species/motif_similarity.csv")).set_index("TF").topH_vs_nearestM_pcc
add("motif similarity of the baseline monoPWMs (range, median)", "%.2f-%.2f, median %.2f" % (_simo.min(), _simo.max(), _simo.median()),
    "%.2f-%.2f, median %.2f (min %s)" % (sim.min(), sim.max(), sim.median(), sim.idxmin()), rel("t5") + "; old: PKG/analysis/cross_species/motif_similarity.csv", "M28, R06, R16",
    "minimum was %s (%.2f) in round 1, now %s (%.2f); the baseline monoPWM changed for %d TFs with the retrained data" % (_simo.idxmin(), _simo.min(), sim.idxmin(), sim.min(),
    int((pd.read_csv(SRC["cs"]).set_index("TF").topH_PWM != pd.read_csv(os.path.join(H.PKG, "analysis/cross_species/motif_similarity.csv")).set_index("TF").topH_PWM.reindex(pd.read_csv(SRC["cs"]).TF).values).sum())))
add("motif similarity of the TFs with reduced transfer", ", ".join("%s %.2f" % (t, _simo[t]) for t in ("ERG", "P53", "STA5A")), ", ".join("%s %.2f" % (t, sim[t]) for t in below), rel("t5"), "R16")
def _t1row(r, b):
    d = r.values - b.values
    return "%.3f %.3f %d/%d %.0e" % (np.median(b), np.median(r), (d > 0).sum(), len(d), wilcoxon(d).pvalue)
_pk = pd.read_csv(os.path.join(H.PKG, "results_table.csv"), sep="\t"); _pk = _pk[(_pk.Model == "RandomForestClassifier") & (_pk.PWM == "mono+di")]
_h1 = pd.read_csv(J("results", "round1_package", "analysis_mouse_swap", "hm_mm_results.csv"))
OLD_T1 = {"Human, human-trained, auROC": _t1row(_pk.roc_auc_test_H, _pk.roc_auc_test_H_PWM), "Human, human-trained, auPRC": _t1row(_pk.pr_auc_test_H, _pk.pr_auc_test_H_PWM)}
for x in pd.read_csv(SRC["tab1"]).itertuples():
    if x.row not in OLD_T1:   # the first submission printed the human-trained rows from seed 0 (only the new value changes)
        OLD_T1[x.row] = _t1row(_h1[x.rf_column.replace("HM_rf_test_M_auroc_mean", "HM_rf_test_M_auroc_s0").replace(
            "HM_rf_test_M_auprc_mean", "HM_rf_test_M_auprc_s0")], _h1[x.baseline_column])
for _, x in pd.read_csv(SRC["tab1"]).iterrows():
    add("Table 1, " + x.row, OLD_T1[x.row], "%.3f %.3f %s %s" % (x.baseline_median, x.archipelago_median, x.tfs_improved, x.wilcoxon_p), rel("tab1") + "; old: PKG/results_table.csv, round-1 hm_mm_results.csv", "M26, R06",
        ("RF: %s; below the baseline: %s" % (x.rf_column, x.below if isinstance(x.below, str) else "none")) if x.row.startswith("Mouse") else "")
# Table 1 rows of precision at 40 % recall and false positives at 50 % recall, with the exact two-sided
# Wilcoxon signed-rank P on the non-zero paired differences (scipy's default for n <= 50 with ties falls back to the normal
# approximation: 4.0e-07 / 4.8e-07). Precision: Sup. Table 7 per_TF (HH, RF seed 0 vs the best monoPWM by auPRC); false positives:
# Sup. Table 4 per_TF (human, per 10,000 regions).
def _exact(d):
    nz = d[d != 0]
    return wilcoxon(nz, method="exact").pvalue, len(nz)
_bw = pd.read_csv(SRC["b2w"]); _bw = _bw[_bw.scenario == "HH"]
_prf, _pmo = (_bw[_bw.model == m_].set_index("TF")["prec_at_tpr0.4"] for m_ in ("RF_s0", "mono_byprc"))
_pd = _prf - _pmo.loc[_prf.index]
_pp, _pn = _exact(_pd)
add("Table 1, Human, precision at 40 % recall (RF seed 0 vs best monoPWM by auPRC): medians, TFs improved, exact two-sided Wilcoxon P (V192)",
    "-", "%.3f %.3f %d/%d %s" % (_pmo.median(), _prf.median(), (_pd > 0).sum(), len(_pd), e(_pp)), rel("b2w"), "Table 1 (V192)",
    "exact test on the %d non-zero differences; tie %s, lower %s; normal approximation of Sup. Table 7 'summary': %.1e" % (
        _pn, ", ".join(_pd.index[_pd == 0]) or "none", ", ".join(_pd.index[_pd < 0]) or "none",
        pd.read_csv(SRC["b2"]).query("scenario == 'HH' and rf_seed == 0 and metric == 'prec_at_tpr0.4'").wilcoxon_p.iloc[0]))
_o4 = pd.read_excel(SRC["t4"], sheet_name="per_TF"); _o4 = _o4[_o4["Test set"].str.startswith("human")].set_index("TF")
_fm, _fr = _o4["FP per 10,000 regions at 50 % recall: best monoPWM"], _o4["FP per 10,000 regions at 50 % recall: ArChIPelago RF"]
_fd = _fr - _fm
_fp_, _fn = _exact(_fd)
# Sup. Table 4 'summary' holds the same test (the source the Table 1 caption names)
_s4 = t4[t4["Test set"].str.startswith("human") & (t4.Quantity == "FP reduction at 50 % recall (%)")].iloc[0]
if _s4["TFs in the test (non-zero differences)"] != _fn or float("%.2g" % _fp_) != _s4["Wilcoxon P (exact, two-sided; paired per TF, equal values left out)"]:
    raise SystemExit("CHECK FAILED Sup. Table 4 summary Wilcoxon P/n %s/%s vs recomputed %.2g/%d" % (
        _s4["Wilcoxon P (exact, two-sided; paired per TF, equal values left out)"], _s4["TFs in the test (non-zero differences)"], _fp_, _fn))
add("Table 1, Human, false positives per 10,000 regions at 50 % recall (RF seed 0 vs best monoPWM by auPRC): medians, TFs with fewer, exact two-sided Wilcoxon P (V193)",
    "-", "%.0f %.0f %d/%d %s" % (_fm.median(), _fr.median(), (_fd < 0).sum(), len(_fd), e(_fp_)), rel("t4"), "Table 1 (V193)",
    "exact test on the %d non-zero differences; tie %s, more %s; median reduction %.1f %% (Sup. Table 4 'summary')" % (
        _fn, ", ".join(_fd.index[_fd == 0]) or "none", ", ".join(_fd.index[_fd > 0]) or "none",
        t4[t4["Test set"].str.startswith("human") & (t4.Quantity == "FP reduction at 50 % recall (%)")]["Median over TFs"].iloc[0]))
# ---- mouse operational (Sup. Table 4)
for _, x in t4.iterrows():
    sp = "human" if x["Test set"].startswith("human") else "mouse"
    add("Sup. Table 4 %s: %s, median (IQR; TFs improved)" % (sp, x.Quantity), "-", "%g (%s; %d of %d)" % (x["Median over TFs"], x.IQR, x["TFs with improvement (of %d)" % NTF], NTF), rel("t4"), "M20, T06, H03",
        "median and IQR over %d TFs (Sup. Table 4 'TFs in the median'; V188)" % x["TFs in the median"] if x["TFs in the median"] != NTF else "")
m50 = t4[(t4["Test set"].str.startswith("mouse")) & (t4.Quantity.str.startswith("FP reduction at 50"))].iloc[0]
m80 = t4[(t4["Test set"].str.startswith("mouse")) & (t4.Quantity.str.startswith("FP reduction at 80"))].iloc[0]
add("mouse FP reduction at 50 % / 80 % recall", "also above 20 % (20.5 % and 20.7 %)", "%.1f %% / %.1f %%" % (m50["Median over TFs"], m80["Median over TFs"]), rel("t4"), "M20",
    ("statement holds" if min(m50["Median over TFs"], m80["Median over TFs"]) > 20 else "statement flips") +
    ": on the mouse test set the reductions were %.1f %% and %.1f %% (50 %% and 80 %% recall; %d and %d of %d TFs)" % (m50["Median over TFs"], m80["Median over TFs"],
    m50["TFs with improvement (of %d)" % NTF], m80["TFs with improvement (of %d)" % NTF], NTF))
# ---- low-FPR metrics, mouse (Sup. Table 7)
b2 = pd.read_csv(SRC["b2"]); b2 = b2[b2.rf_seed == 0]
for sc in ("HM", "MM"):
    for m in ("auroc", "auprc", "tpr_at_fpr0.01", "pauc_std_fpr0.01"):
        s = b2[(b2.scenario == sc) & (b2.metric == m)].iloc[0]
        add("Sup. Table 7 %s seed 0: %s RF / baseline median; TFs improved; Wilcoxon P" % (sc, m), "-", "%.4f / %.4f; %d of %d; %s" % (s.rf_median, s.base_median, s.n_improved, s.n_TF, e(s.wilcoxon_p)), rel("b2"), "T08")
b2c = pd.read_csv(SRC["b2c"])
for sc in ("HM", "MM"):
    s = b2c[(b2c.scenario == sc) & (b2c.metric == "auroc")]
    add("per mouse test chromosome %s: median dauROC (TFs improved) chr1, chr8, chr19" % sc, "-", "; ".join("%s %.4f (%d of %d)" % (c.chrom, c.diff_median, c.n_improved, c.n_TF) for c in s.itertuples()), rel("b2c"), "R07, R10")
# ---- top-N mouse (Sup. Table 8)
tn = pd.read_csv(SRC["topn"]); tn = tn[(tn.species == "M") & (tn.variant == "allneg") & (tn.subset == "avail")]
for c in tn.itertuples():
    add("Sup. Table 8 mouse top %s: TFs; RF vs PWM auROC median; d auROC" % c.N, "-", "%d; %.4f vs %.4f; %+.4f" % (c.n_tf, c.rf_auroc_med, c.pwm_auroc_med, c.d_auroc_med), rel("topn"), "T12")
# ---- gain vs positives
rho = pd.read_csv(SRC["rho"]).set_index(["test_set", "gain", "predictor_col"])
for k, old in ((("human test", "dauroc_H", "n_train_pos_H"), "rho = -0.33 (P = 0.053)"), (("human test", "dauprc_H", "n_train_pos_H"), "rho = +0.25 (P = 0.14)"),
               (("human test", "dauprc_H", "n_datasets_H"), "rho = +0.43 (P = 0.009)"), (("human test", "dauroc_H", "n_datasets_H"), "rho = -0.06 (P = 0.71)")):
    s = rho.loc[k]
    p_old = float(re.search(r"P = ([\d.]+)", old).group(1))
    add("gain vs predictor, %s %s vs %s" % k, old, "rho = %+.2f (P = %.2g)" % (s.rho, s.p), rel("rho"), "R12, T13",
        ("statement flips: %s vs %s %s nominally significant (P %.2g, was %s)" % (k[1], k[2], "is now" if s.p < 0.05 else "is no longer", s.p, p_old))
        if (s.p < 0.05) != (p_old < 0.05) else "")
# ---- drivers
d = pd.read_csv(SRC["drv"])
add("drivers of transfer: tests, n", "48 tests, n = 36", "%d tests, n = %d" % (len(d), d.n.iloc[0]), rel("drv"), "M28, R16, L15")
add("drivers: smallest BH-adjusted P", "0.15", "%.2f" % d.q_BH.min(), rel("drv"), "M28, R16, L15",
    "no correlation passes BH (smallest q %.2f)" % d.q_BH.min() if d.q_BH.min() >= 0.05 else "statement flips: %d tests pass BH at q < 0.05" % (d.q_BH < 0.05).sum())
for drv, old in (("IC_mean_per_pos_auROC", "-0.46"), ("base_auROC_M", "-0.45")):
    s = d[(d.driver == drv) & (d.outcome == "dROC_HM")].iloc[0]
    add("drivers: %s vs dauROC H>M (rho, nominal P)" % drv, old, "%+.2f (P = %.2g)" % (s.rho, s.p), rel("drv"), "M28, R16")
add("drivers: nominal P < 0.05 (daggers)", "-", "%d of %d" % ((d.p < 0.05).sum(), len(d)), rel("drv"), "L15")
hd = pd.read_csv(SRC["head"]).set_index("quantity")
clouds = open(SRC["clouds"]).read()
_cl = {k: next(l.split(":", 1)[1].strip() for l in clouds.splitlines() if l.startswith(k)) for k in ("top", "right")}
add("Fig. 3C groups (dauROC > 0.075; dauPRC > 0.15)", "4 TFs; 3 TFs (ANDR, MAFK and P53)", "; ".join(_cl[k] for k in ("top", "right")), rel("clouds"), "M21, L15, R13",
    "true statement: the TFs that stand out by their auPRC gain are %s" % _cl["right"])
add("headroom: dauROC vs baseline auROC (rho, P)", "-0.65 (2e-5)", "%+.2f (%s)" % (hd.loc["dauROC vs baseline auROC", "rho"], e(hd.loc["dauROC vs baseline auROC", "p"])), rel("head"), "M21, R13")
# ---- Fig. S4 mouse: Slim, LSlim, diChIPMunk and the RF2f family retrained on negatives v3; new values from the
#      Fig. S4 stats table of the figure producer, old from the first-submission numbers (36 TFs)
s4 = pd.read_csv(SRC["s4"]).set_index(["version", "model_key"])
s4n = pd.read_csv(SRC["s4n"]).set_index(["metric", "model"])
for mk, lab in (("RF2f_all5", "five-feature model"), ("RF_all", "RF on all PWMs"), ("slim1", "Slim m=1"), ("lslim5", "LSlim m=-5"), ("munk", "diChIPMunk"), ("RF2f", "RF2f")):
    o = s4.loc[("round1_36TF", mk)]
    nr, npr = s4n.loc[("auROC", o.model)], s4n.loc[("auPRC", o.model)]
    add("Fig. S4 mouse, %s: median d auROC / d auPRC; TFs > 0" % lab, "%+.3f / %+.3f; %d / %d of 36" % (o.median_d_auROC, o.median_d_auPRC, o.n_TF_d_auROC_gt0, o.n_TF_d_auPRC_gt0),
        "%+.3f / %+.3f; %d / %d of %d" % (nr.med, npr.med, nr.n_pos, npr.n_pos, nr.n), rel("s4n") + " (new); " + rel("s4") + " (old)", "M30, F08",
        "H3b: Fig. S4 models retrained on negatives v3; one-sided P %.1e / %.1e" % (nr.p, npr.p))
# ---- text tokens [PENDING_H3B:<name>] of the Abstract and Results 3.4 (fig4_tokens_h3b.py from the Fig. 4 / S4 stats)
for _, q in pd.read_csv(SRC["tok"]).iterrows():
    add("text token PENDING_H3B:%s" % q.token, q.round1_text, q.value, rel("tok"), "Results 3.4" if q.token != "Abstract_Slim_claim" else "Abstract",
        ("statement holds: " if q.statement_holds else "STATEMENT FLIPS: ") + q.note)
# ---- ablation (Methods): TFs per similarity cut and the size of the random control
ab = pd.read_csv(SRC["abl"])
_cut = ab.groupby("TF").sim_threshold.first()
_cuts = "; ".join("%d at cut %g (%s)" % ((_cut == c).sum(), c, ", ".join(sorted(_cut.index[_cut == c]))) for c in sorted(_cut.unique()))
_abl_old = "12 TFs (letter)"   # old text value (quote of the response letter)
add("ablation: TFs and similarity cut", _abl_old, "%d TFs: %s" % (_cut.size, _cuts), rel("abl"), "H11, M14, R2-F3",
    "Methods: the ablation ran on %d TFs (%s), not %s" % (_cut.size, _cuts, _abl_old.split(" (")[0]))
_rnd = ab[ab.cond == "drop_random"].set_index("TF"); _sec = ab[ab.cond == "drop_secondary"].set_index("TF")
_small = sorted(_rnd.index[_rnd.n_removed < _sec.n_removed])
add("ablation: random control smaller than the secondary removal", "-", "; ".join("%s: %d random PWMs removed vs %d secondary (%d primary PWMs)" % (
    t, _rnd.n_removed[t], _sec.n_removed[t], _rnd.n_prim[t]) for t in _small) or "none", rel("abl"), "H11, R2-F3",
    "Methods: the random control removes as many PWMs as the secondary removal, capped by the number of primary PWMs (%s)" % (
    ", ".join("%s %d of %d" % (t, _rnd.n_removed[t], _sec.n_removed[t]) for t in _small) or "none"))
# ---- Fig. S7 (motif subtypes): selection rule = TFs among the top 8 by auROC gain whose top PWMs show distinct subtypes
_s7r = pd.read_csv(SRC["subrank"]).set_index("TF"); _s7v = pd.read_csv(SRC["subverd"]).set_index("TF")
_s7o = O.csv("manuscript/B6_figures/FigS7/Figure_S7_source_data.csv").drop_duplicates("TF").set_index("TF")
_s7shown_old = ["E2F4", "RXRA", "TAL1", "TFE2"]    # old text value: the TFs of the first-submission legend (quote)
_ntop = int(_s7r.top8_dauROC.sum())
_top_new = list(_s7r.sort_values("rank_dauROC").index[_s7r.sort_values("rank_dauROC").top8_dauROC])
_top_old = list(_s7o[_s7o.rank_dauROC <= _ntop].sort_values("rank_dauROC").index)
_left = [t for t in _s7shown_old if t not in _top_new]; _gone = [t for t in _top_old if t not in _top_new]; _came = [t for t in _top_new if t not in _top_old]
add("Fig. S7: the TFs with the largest auROC gain (top %d, human test set)" % _ntop, ", ".join(_top_old) + " (36 TFs)", ", ".join(_top_new),
    rel("subrank") + "; old git %s:B6_figures/FigS7/Figure_S7_source_data.csv rank_dauROC" % O.TAG, "S7, R2-F2",
    "statement changes: %s left the top %d, %s entered" % (", ".join(_gone) or "none", _ntop, ", ".join(_came) or "none"))
add("Fig. S7: auROC gain rank of the TFs shown in round 1", ", ".join("%s %d" % (t, _s7o.rank_dauROC[t]) for t in _s7shown_old),
    ", ".join("%s %d" % (t, _s7r.rank_dauROC[t]) for t in _s7shown_old), rel("subrank"), "S7",
    "statement changes: %s no longer among the %d TFs with the largest auROC gain" % (", ".join("%s (rank %d)" % (t, _s7r.rank_dauROC[t]) for t in _left) or "none", _ntop))
# Fig. S7 on HOCOMOCO v11: TFs and roles from the v11 verdicts, tile labels from the v11 match
from decimal import Decimal, ROUND_HALF_UP   # noqa: E402  (r values are stored with 3 decimals; 0.975 -> 0.98, as in the text pass)
_f3 = lambda x: "%.3f" % x   # r is rebound above
_f2 = lambda x: str(Decimal(repr(round(float(x), 3))).quantize(Decimal("0.01"), rounding=ROUND_HALF_UP))
_dr = _s7v[_s7v.drawn].sort_values("figure_row")
_shown = list(_dr.index)
_byrole = lambda ro: ", ".join(_dr.index[_dr.role == ro])
add("Fig. S7: TFs shown (HOCOMOCO v11 subtype examples; own motif plus co-factor motif)",
    "%d of the %d (%s)" % (len(_s7shown_old), _ntop, ", ".join(_s7shown_old)),
    "%s (v11 subtypes); %s (own motif plus co-factor motif)" % (_byrole("v11 subtype example"), _byrole("own motif plus co-factor motif")),
    rel("subverd") + "; " + rel("subsrc"), "S7",
    "statement changes (legend): subtype labels come from the nearest HOCOMOCO v11 model, not from logo reading; %s dropped, "
    "%s added; panel letters change" % (", ".join(t for t in _s7shown_old if t not in _shown) or "none",
                                        ", ".join(t for t in _shown if t not in _s7shown_old) or "none"))
_v11 = pd.read_csv(SRC["v11"]).set_index(["TF", "panel"])
# Rows of the text pass (item names kept so that the text references stay valid), produced here from the v11 csv.
# Each row names its panels; the panels selected from the data (nearest own subtype / co-factor flag) must equal them.
_cof = lambda tf: set(p for (t, p), q in _v11.loc[[tf]].iterrows() if q.best_any_v11.split("_HUMAN")[0] != tf and q.r_best_any - q.r_best_own >= 0.2)
_sub = lambda tf, k: set(p for (t, p), q in _v11.loc[[tf]].iterrows() if q.best_own_subtype == k)
_rng = lambda v: "%s (%s)" % (_f3(v[0]), _f3(v[0])) if len(v) == 1 else "%s-%s (%s-%s)" % (_f2(min(v)), _f2(max(v)), _f3(min(v)), _f3(max(v)))
_rr = lambda tf, ps, c="r_best_own": [_v11.loc[(tf, p_), c] for p_ in sorted(ps)]
_rows_v11 = [  # item, TF, panels in the item, panels from the data, value, inventory, note
    ("Fig. S7 STAT1: GAS PWMs (mono 1, 2, 4, di 1-3), r with HOCOMOCO v11 STAT1 subtype 0 (V1-S7)", "STAT1", {1, 2, 4, 5, 6, 7}, _sub("STAT1", 0),
     lambda ps: _rng(_rr("STAT1", ps)), "V1-S7, Results 3.3, Fig. S7 legend, letter R1-B", ""),
    ("Fig. S7 STAT1: ISRE-like PWM (mono 3, K562), r with HOCOMOCO v11 STAT1 subtype 1 (V1-S7)", "STAT1", {3}, _sub("STAT1", 1),
     lambda ps: "%s (%s)" % (_f2(_rr("STAT1", ps)[0]), _f3(_rr("STAT1", ps)[0])), "V1-S7, Results 3.3, Fig. S7 legend", ""),
    ("Fig. S7 TAL1: E-box PWMs (mono 1-3; Jurkat E6-1, CCRF-CEM), r with HOCOMOCO v11 TAL1 subtype 1 (V1-S7)", "TAL1", {1, 2, 3}, _sub("TAL1", 1),
     lambda ps: _rng(_rr("TAL1", ps)), "V1-S7, Results 3.3, Fig. S7 legend", ""),
    ("Fig. S7 TAL1: GATA-E-box PWMs (mono 4, di 1-3; pro-erythroblasts, megakaryocytes), r with HOCOMOCO v11 TAL1 subtype 0 (V1-S7)", "TAL1",
     {4, 5, 6, 7}, _sub("TAL1", 0), lambda ps: _rng(_rr("TAL1", ps)), "V1-S7, Results 3.3, Fig. S7 legend", ""),
    ("Fig. S7 IRF4: own-motif PWMs (mono 1-4, di 1, 2), r with HOCOMOCO v11 IRF4 model (V1-S7)", "IRF4", {1, 2, 3, 4, 5, 6}, _sub("IRF4", 0) - _cof("IRF4"),
     lambda ps: _rng(_rr("IRF4", ps)), "V1-S7, Fig. S7 legend", ""),
    ("Fig. S7 IRF4: AP-1 PWM (di 3, GM12878), r with BATF_HUMAN.H11MO.0.A / with IRF4 model (V1-S7)", "IRF4", {7}, _cof("IRF4"),
     lambda ps: "%s / %s" % (_f3(_rr("IRF4", ps, "r_best_any")[0]), _f3(_rr("IRF4", ps)[0])), "V1-S7, Results 3.3, Fig. S7 legend, letter R1-B",
     "nearest model of another TF: co-factor motif"),
    ("Fig. S7 TFE2: ETS-like PWMs (mono 1, di 1; RPMI-8402), r with ETV7_HUMAN.H11MO.0.D / with TFE2 model (V1-S7)", "TFE2", {1, 5}, _cof("TFE2"),
     lambda ps: "%s / %s" % (" and ".join(map(_f3, _rr("TFE2", ps, "r_best_any"))), " and ".join(map(_f3, _rr("TFE2", ps)))),
     "V1-S7, Results 3.3, Fig. S7 legend, letter R1-B", "nearest model of another TF: co-factor motif"),
    ("Fig. S7 TFE2: E-box PWMs (mono 2-4, di 2, 3), r with HOCOMOCO v11 TFE2 model (V1-S7)", "TFE2", {2, 3, 4, 6, 7}, _sub("TFE2", 0) - _cof("TFE2"),
     lambda ps: _rng(_rr("TFE2", ps)), "V1-S7, Fig. S7 legend", ""),
]
for _it, _tf, _want, _got, _fv, _inv, _note in _rows_v11:
    if _want != _got:
        raise SystemExit("CHECK FAILED %s: panels from the v11 csv %s differ from the item %s" % (_it, sorted(_got), sorted(_want)))
    _bo = set(_v11.loc[(_tf, p_), "best_any_v11"] for p_ in _want) if "/ with" in _it else set()
    if _bo and not all(m_ in _it for m_ in _bo):
        raise SystemExit("CHECK FAILED %s: nearest any-TF model %s not the one named" % (_it, sorted(_bo)))
    add(_it, "-", _fv(_want), rel("v11"), _inv, _note)
# r of the TFE2 ETS-like PWMs to the best model of each reference TF the text names, with quality
# letters, and the co-factor gap r_ref - r_own (subtype_vs_named_v11.py)
_vn = pd.read_csv(SRC["v11n"]).set_index(["TF", "panel"])
_vt = [(1, "mono 1"), (5, "di 1")]
for _ref in ("ETV2", "ETV7"):
    _bm = set(_vn.loc[("TFE2", p_), "best_" + _ref] for p_, _ in _vt)
    _ql = set(_vn.loc[("TFE2", p_), "quality_" + _ref] for p_, _ in _vt)
    if len(_bm) != 1 or len(_ql) != 1:
        raise SystemExit("CHECK FAILED TFE2 mono 1 / di 1: different best %s models %s" % (_ref, sorted(_bm)))
    _bm, _ql = _bm.pop(), _ql.pop()
    add("Fig. S7 TFE2: ETS-like PWMs (mono 1, di 1; RPMI-8402), r with %s (quality %s)" % (_bm, _ql), "-",
        " and ".join(_f3(_vn.loc[("TFE2", p_), "r_" + _ref]) for p_, _ in _vt), rel("v11n"), "V1-S7, G1-5, G2-5, Fig. S7 legend, Results 3.3",
        "best model of the named reference TF; quality letter from the HOCOMOCO v11 annotation")
for p_, lab in _vt:
    add("Fig. S7 TFE2 %s: co-factor gap, r with ETV2 model minus r with TFE2 model / ETV7 model minus TFE2 model" % lab, "-",
        "%s / %s" % (_f3(_vn.loc[("TFE2", p_), "gap_ETV2"]), _f3(_vn.loc[("TFE2", p_), "gap_ETV7"])), rel("v11n"), "G2-4, G2-8, Fig. S7 legend",
        "the co-factor flag of Fig. S7 needs a gap of at least two tenths; with the ETV2 model mono one stays below it")
# Fig. S7 legend counts (top three diPWMs per TF): PWMs per TF, candidate TFs, drawn TFs and tiles
_ss = pd.read_csv(SRC["subsrc"])
_npc = _ss.groupby(["TF", "pwm_class"]).size().unstack()
if _npc.nunique().max() != 1:
    raise SystemExit("CHECK FAILED Fig. S7: TFs with different numbers of mono / di panels")
_drw = _ss[_ss.TF.isin(_shown)]
add("Fig. S7 legend counts: monoPWMs / diPWMs per TF; candidate TFs; TFs drawn; tiles drawn; co-factor tiles", "4 / 2; -; 4; 24; 3",
    "%d / %d; %d; %d; %d; %d" % (_npc["mono"].iloc[0], _npc["di"].iloc[0], _ss.TF.nunique(), len(_shown), len(_drw), int(_drw.cofactor.sum())),
    rel("subsrc"), "V1-S7, Fig. S7 legend", "authors' decision: the top three diPWMs of every TF, one rule for all TFs")
add("Fig. S7 subtype match: method (V1-S7)", "-", "pcc_max: max over offsets and strands of mean column Pearson r of z-scored PWMs, overlap >= 5; "
    "HOCOMOCO v11 HUMAN mono, full collection", "4_visualisation/figures/FigS7/hocomoco11_check/subtype_vs_hocomoco11.py", "V1-S7, Methods 2.10")
_cv0 = pd.read_csv(SRC["cv11"])
_cv = _cv0.assign(key=_cv0.TF.fillna(_cv0.summary)).set_index("key")
_cvt = _cv[_cv.TF.notna()]                        # the 11 TF rows (summary rows: TF empty, label in 'summary')
_cvi = _cv.index[_cv.index.str.startswith("median, informative")][0]
_cut9 = ", ".join(_cvt.index[_cvt.cut != 0.8])
add("Fig. S7 cluster subtype vs nearest HOCOMOCO v11 subtype, agreement (V1-S7)",
    "-", "median Rand index %s over the %d ablated TFs; median ARI %s over the %d informative TFs (both partitions with at least two "
         "groups: %s); ablation subtypes (cut 0.8; %s at 0.9), %d PWMs" % (
        _f2(_cv.rand_index["median, all TFs"]), len(_cvt), _f2(_cv.adjusted_rand[_cvi]), int(_cvt.informative.sum()),
        ", ".join(_cvt.index[_cvt.informative.astype(bool)]), _cut9, int(_cvt.n_PWMs.sum())),
    rel("cv11"), "V1-S7, Results 3.3, D12", "agreement on the ablation subtypes (Methods), median Rand over all TFs and median ARI over the informative TFs")
add("Fig. S7 cluster subtype vs nearest HOCOMOCO v11 subtype, per TF: Rand index / ARI (D12)", "-",
    ", ".join("%s %s / %s" % (t, _f2(q.rand_index), _f2(q.adjusted_rand)) for t, q in _cvt.iterrows()), rel("cv11"), "Sup. Table 9, D12",
    "per TF; ARI is zero by construction where one partition has a single group")
# ---- low-FPR step-0 gate (best-mono baseline)
_g = pd.read_csv(SRC["gate"]); _go = O.csv("manuscript/B2_metrics/B2_gate_summary.csv").set_index("check")
_gid = _g[_g.decision != "size only"]
add("B2 gate: checks passed (recomputed from the score npz vs hm_mm, b2_v2, results table)", "-",
    "%d of %d" % ((_gid.decision == "pass").sum(), len(_gid)), rel("gate"), "B2, R2-F1, R2-F2")
for _sc in ("auROC", "auPRC(integral)"):
    _row = _g[_g.check == "HH %s RF s0 vs results table RF (H3)" % _sc].iloc[0]
    add("B2 gate: %s RF seed 0 from the npz vs the RF of the results table, max |diff| per TF" % _sc,
        "%.2g (%s)" % (_go.loc[_row.round1_check, "max_abs_diff"], _row.round1_decision), "%.2g (%s)" % (_row.max_abs_diff, _row.decision),
        rel("gate") + "; old git %s:manuscript/B2_metrics/B2_gate_summary.csv" % O.TAG, "B2",
        "statement changes: the low-FPR metrics now use the forest of the results table (previous round: a retrained forest, gate failure accepted); "
        "drop the caveat 'per-TF RF values are from a re-trained forest' (%s)" % os.path.basename(SRC["gate"]))
# ---- verification rows
_mv = pd.read_csv(SRC["mvd"])
for _, _r in _mv.iterrows():
    add("RF on monoPWMs + diPWMs vs RF on monoPWMs only, %s test set, %s (V1-9)" % (_r.test_set, _r.metric), "-",
        "median %.4f vs %.4f (difference %+.4f); %d of %d TFs improved; two-sided Wilcoxon P %s" % (_r.median_monodi, _r.median_mono, _r.median_diff,
        _r.n_improved, _r.n_TF, e(_r.wilcoxon_p_two_sided)), rel("mvd"), "V1-9",
        "Results: the diPWMs add to the monoPWMs; TFs where mono+di is lower: %s" % _r.worse_TFs.replace(",", ", "))
_gc = pd.read_csv(SRC["gct"]); GC_T = 5   # rule of the Methods sentence: percentage points
for (_sp, _sl), _x in _gc.groupby(["species", "split"]):
    add("GC content, %s %s sets: per-TF median GC of positives minus negatives (percentage points), range; TFs with |difference| > %d points (V1-16)" % (_sp, _sl, GC_T),
        "-", "%+.1f to %+.1f; %d of %d TFs (%s)" % (_x.gc_diff_points.min(), _x.gc_diff_points.max(), (_x.gc_diff_points.abs() > GC_T).sum(), len(_x),
        ", ".join(_x.TF[_x.gc_diff_points.abs() > GC_T]) or "none"), rel("gct"), "V1-16, M10")
for _, _r in pd.read_csv(SRC["gcs"]).iterrows():
    add("GC control: Spearman of |GC difference| (test set) with the gain over the best monoPWM, %s, %s (V1-16)" % (_r.test_set, _r.metric), "-",
        "rho %+.2f, P %.2g (%d TFs)" % (_r.spearman_rho, _r.p, _r.n_TF), rel("gcs"), "V1-16",
        "control: the gain does not follow the residual GC difference" if _r.p >= 0.05 else "the gain follows the residual GC difference")
_ti = pd.read_csv(SRC["ties"]); _tn = _ti[~_ti.at_cap]
add("TFs that share a positive count below the cap (V1-17)", "-", "; ".join("%s %s: %s (%s)" % (r_.species, "training" if "training" in r_.column else "test",
    r_.TFs.replace(",", " and "), format(int(r_["count"]), ",")) for _, r_ in _tn.iterrows()) or "none", rel("ties"), "V1-17",
    "different pools after the mask (%s): equal counts by coincidence; the other shared counts are the cap" % "; ".join(_tn.pool_after_mask))
_rtt = pd.read_csv(SRC["rt"], sep="\t")
_rf3 = _rtt[(_rtt.Model == "RandomForestClassifier") & (_rtt.PWM == "mono+di")].set_index("TF_name")
_bel = _rf3[(_rf3.roc_auc_test_H < _rf3.roc_auc_test_H_PWM) | (_rf3.pr_auc_test_H < _rf3.pr_auc_test_H_PWM)]
add("Fig. 3C open diamond: RF (monoPWMs + diPWMs) below the best monoPWM on auROC or auPRC, human test set (V1-19)", "-",
    ", ".join("%s (dauROC %+.4f, dauPRC %+.4f)" % (t_, r_.roc_auc_test_H - r_.roc_auc_test_H_PWM, r_.pr_auc_test_H - r_.pr_auc_test_H_PWM) for t_, r_ in _bel.iterrows()) or "none",
    rel("rt"), "V1-19", "legend: open diamond = below the best monoPWM on at least one metric (the rule of the supplementary scatter plots)")
add("minimum PWMs per TF (V1-21)", ">= 5 monoPWMs", "; ".join("%s %d monoPWMs (%s), %d diPWMs (%s)" % (sp_.lower(), t2["%s monoPWMs" % sp_].min(),
    ", ".join(t2.index[t2["%s monoPWMs" % sp_] == t2["%s monoPWMs" % sp_].min()]), t2["%s diPWMs" % sp_].min(),
    ", ".join(t2.index[t2["%s diPWMs" % sp_] == t2["%s diPWMs" % sp_].min()])) for sp_ in ("Human", "Mouse")), rel("t2"), "V1-21, M06")
_ab = pd.read_csv(SRC["abf"])
_sd = _ab.groupby(["TF", "cond"]).seed.nunique()
add("ablation: seeds behind the 3 SD rule (random removals per TF; V1-25)", "-", "%d seeds per TF and condition (%s); %d TFs" % (
    _sd.min(), ", ".join(sorted(_ab.cond.unique())), _ab.TF.nunique()) if _sd.min() == _sd.max() else "seeds differ: %s" % _sd.to_dict(),
    rel("abf"), "V1-25", "Methods/legend: the SD of the random removals rests on %d fits per TF" % _sd.min())
# ---- Methods 2.11 / Sup. Table 6: benchmark rerun on the final training matrices; old = the released Sup. Table 6
for _, _r in pd.read_csv(SRC["s6"]).iterrows():
    add("Methods 2.11 / Sup. Table 6: %s (V1-3)" % _r.quantity, _r.round1_released, _r.h3, rel("s6") + "; " + rel("t6"), "V1-3, M2.11",
        "" if str(_r.round1_released) == str(_r.h3) else "statement changes: %s -> %s" % (_r.round1_released, _r.h3))
# ---- Fig. S5 mouse
sat = pd.read_csv(SRC["sat"])
for m in ("auroc_M", "auprc_M"):
    s = sat[(sat.design == "random") & (sat.metric == m)].set_index("k").median_delta
    add("Fig. S5 mouse random subsets, median delta %s at k = 16 / 32 / all (fraction of full)" % m, "-", "%.4f / %.4f / %.4f (%.2f / %.2f)" % (s["16"], s["32"], s["all"] if "all" in s.index else s.iloc[-1],
        s["16"] / s.iloc[-1], s["32"] / s.iloc[-1]), rel("sat"), "M32, F09")
# ---- rows for numbers that the text quotes (feature importance, drivers, response-letter analysis)
_g2 = lambda p: "%.2g" % p
for _q, _lab in (("dauPRC/(1-auPRC_base) vs baseline auPRC", "auPRC gain divided by the headroom (1 - baseline auPRC) vs baseline auPRC"),
                 ("baseline auPRC vs prevalence of positives (human test)", "baseline auPRC vs share of positives (human test)"),
                 ("dauROC/auROC_base (relative) vs baseline auROC", "relative auROC gain (dauROC / baseline auROC) vs baseline auROC")):
    add("headroom: %s (rho, P; N4-N6)" % _lab, "-", "%+.2f (P = %s)" % (hd.loc[_q, "rho"], _g2(hd.loc[_q, "p"])), rel("head"), "N4-N6, Results 3.2, letter R1-L14")
_hp = pd.read_csv(SRC["hpt"]).set_index("TF")
_lo, _hi = _hp.bROC[["IRF4", "SRF", "E2F4", "RXRA"]], _hp.bROC.nlargest(3).sort_index()   # the TFs the text names (the four of Fig. 3C, top auROC gain)
add("headroom: baseline auROC of IRF4, SRF, E2F4 and RXRA, human test (N7)", "-", "%s (range %.2f-%.2f)" % (", ".join("%s %.3f" % (t_, v_) for t_, v_ in _lo.items()),
    _lo.min(), _lo.max()), rel("hpt"), "N7, Results 3.2, letter R1-L14")
add("headroom: TFs with the highest baseline auROC, human test (N7)", "-", "%s (%.2f)" % (", ".join("%s %.4f" % (t_, v_) for t_, v_ in _hi.items()), _hi.median()),
    rel("hpt"), "N7, Results 3.2, letter R1-L14")
# the sentence names the three smallest auPRC gains, not JUND and COE1
_f3c = pd.read_csv(SRC["f3c"]).set_index("TF_name")
_sm = _f3c.dPRC.nsmallest(3)
add("Fig. 3C: the three smallest auPRC gains, human test (N9)", "-", ", ".join("%s %.4f" % (t_, v_) for t_, v_ in _sm.items()), rel("f3c"),
    "N9, Results 3.3", "replaces the JUND and COE1 row: the text names the three smallest auPRC gains")
_b5 = _f3c[(_f3c.dROC > 0.05) & (_f3c.dPRC > 0.05)]
add("Fig. 3C: TFs with both gains above 0.05, human test (dauROC / dauPRC)", "-", "%d: %s" % (len(_b5), ", ".join(
    "%s %.4f / %.4f" % (t_, q_.dROC, q_.dPRC) for t_, q_ in _b5.sort_index().iterrows())), rel("f3c"), "Results 3.3")
_tn = pd.read_csv(SRC["topn"]); _tn = _tn[(_tn.species == "H") & (_tn.subset == "avail")].set_index(["variant", "N"])
for (_v, _n, _m), _lab in ((("rorig", "500", "auroc"), "negatives at the original ratio, top 500, auROC"), (("allneg", "1000", "auprc"), "all negatives, top 1,000, auPRC"),
                           (("rorig", "1000", "auprc"), "negatives at the original ratio, top 1,000, auPRC"), (("rorig", "1000", "auroc"), "negatives at the original ratio, top 1,000, auROC")):
    _x = _tn.loc[(_v, _n)]
    add("top-N human: gain on the top N minus gain on all positives, %s (median; two-sided P; N8)" % _lab, "-",
        "%+.4f; P = %s" % (_x["topN_minus_allpos_%s_med" % _m], _g2(_x["p_topN_vs_allpos_%s" % _m])), rel("topn"), "N8, Results 3.2, letter R1-L6")
_sec = ab[ab.cond == "drop_secondary"].set_index("TF"); _rnd = ab[ab.cond == "drop_random"].set_index("TF")
add("ablation: d auROC after removing the secondary-subtype PWMs, IRF4, TFE2, TAL1, PRGR (mean over seeds; N10)", "-",
    ", ".join("%s %+.4f" % (t_, _sec.d_auroc_mean[t_]) for t_ in ("IRF4", "TFE2", "TAL1", "PRGR")), rel("abl"), "N10, Results 3.3, letter R1-B")
add("ablation: d auPRC of SOX2, secondary-subtype removal vs random removal (mean over seeds; N10)", "-",
    "%+.4f vs %+.4f" % (_sec.d_auprc_mean["SOX2"], _rnd.d_auprc_mean["SOX2"]), rel("abl"), "N10, Results 3.3, letter R1-B")
_lossm = lambda m: sorted(t_ for t_ in _sec.index if _sec.loc[t_, "d_%s_mean" % m] < _rnd.loc[t_, "d_%s_mean" % m] - 3 * _rnd.loc[t_, "d_%s_sd" % m])
add("ablation: TFs with a loss beyond 3 sd of the random removals, by auROC / by auPRC", "-",
    "auROC: %s; auPRC: %s" % (", ".join(_lossm("auroc")) or "none", ", ".join(_lossm("auprc")) or "none"), rel("abl"), "Results 3.3, letter R1-B")
_co = pd.read_csv(SRC["conc"])
add("importance concentration: PWMs that carry 50 % / 80 % of the Gini importance, median (range) over TFs (N11)", "-",
    "%g (%d-%d) / %g (%d-%d)" % (_co.n_for_50.median(), _co.n_for_50.min(), _co.n_for_50.max(), _co.n_for_80.median(), _co.n_for_80.min(), _co.n_for_80.max()),
    rel("conc"), "N11, Results 3.3, letter R1-B")
_dv = pd.read_csv(SRC["drv"]); _s = _dv[(_dv.driver == "base_auPRC_M") & (_dv.outcome == "dROC_HM")].iloc[0]
add("drivers: base_auPRC_M vs dauROC H>M (rho, nominal P; N12)", "-", "%+.2f (P = %.2g)" % (_s.rho, _s.p), rel("drv"), "N12, Results 3.3, letter R1-L15")
add("USF2 mouse-trained dauROC on the mouse test set: mean of two seeds / seed 0 (N12)", "-", "%+.4f / %+.4f" % (
    hm.MM_rf_test_M_auroc_mean["USF2"] - hm.MM_broc_test_M_auroc["USF2"], hm.MM_rf_test_M_auroc_s0["USF2"] - hm.MM_broc_test_M_auroc["USF2"]), rel("hm"), "N12, Results 3.3")
for _mo, _lab in (("Slim m=1", "Slim m=1"), ("diChIPMunk", "diChIPMunk"), ("ArChIPelago RF2f", "RF2f")):
    add("Fig. S4 mouse, %s: two-sided Wilcoxon P of d auROC / d auPRC (N13)" % _lab, "-", "%s / %s" % (_g2(s4n.loc[("auROC", _mo), "p_two_sided"]),
        _g2(s4n.loc[("auPRC", _mo), "p_two_sided"])), rel("s4n"), "N13, Results 3.4")
for _dr, _oc in (("base_auROC_H", "dROC_HM"), ("base_auPRC_M", "loss_prc")):
    _s = _dv[(_dv.driver == _dr) & (_dv.outcome == _oc)].iloc[0]
    add("drivers: %s vs %s (rho, nominal P, BH q)" % (_dr, _oc), "-", "%+.2f (P = %.2g, q = %.2f)" % (_s.rho, _s.p, _s.q_BH), rel("drv"), "Results 3.3")
_b5t = pd.read_csv(SRC["b5t"]).set_index("TF")
add("ERG: mouse-trained gain on the mouse test set, auROC / auPRC (mean of two seeds); human-trained dauPRC H>M (mean of two seeds)", "-",
    "%+.4f / %+.4f; %+.4f" % (_b5t.gain_MM_roc_2seed["ERG"], _b5t.gain_MM_prc_2seed["ERG"], _b5t.dPRC_HM["ERG"]), rel("b5t"), "Results 3.3")
_s1 = pd.read_csv(SRC["s1p"]); _s1 = _s1[_s1.PWM == "monoPWMs + diPWMs"]
add("Fig. S1 mouse test set, monoPWMs + diPWMs: TFs above the best monoPWM by auPRC, RF / XGBoost", "-", "%d / %d of %d" % tuple(
    [int((q_.auPRC > q_.base_auPRC).sum()) for q_ in (_s1[_s1.Model == "Random Forest (RF)"], _s1[_s1.Model == "XGBoost"])] + [NTF]),
    rel("s1p"), "Results 3.1")
_tnm = pd.read_csv(SRC["topn"]); _tnm = _tnm[(_tnm.species == "M") & (_tnm.variant == "allneg") & (_tnm.subset == "avail")].set_index("N")
add("Sup. Table 8 mouse top 500 / top 1,000: TFs improved by auROC of TFs available", "-", "; ".join("top %s %d of %d" % (
    n_, _tnm.loc[n_, "n_improved_auroc"], _tnm.loc[n_, "n_tf"]) for n_ in ("500", "1000")), rel("topn"), "Results 3.2, T12")
from scipy.stats import wilcoxon as _wx   # noqa: E402
_f4 = pd.read_csv(SRC["f4src"]); _f4 = _f4[_f4.test_set == "human"]
_rfa = _f4[_f4.model == "ArChIPelago RF on all PWMs"].set_index("TF")
_pp = []
for _m in ("auROC", "auPRC"):
    _cn = []
    for _sl in ("Slim m=0", "Slim m=1", "LSlim m=-5"):
        _dd = _rfa[_m] - _f4[_f4.model == _sl].set_index("TF").loc[_rfa.index, _m]
        _cn.append(int((_dd > 0).sum())); _pp.append(_wx(_dd).pvalue)
    add("Fig. 4 paired per TF: RF on all PWMs above Slim m=0, Slim m=1, LSlim m=-5, human test, %s (TFs)" % _m, "-", "%d, %d, %d of %d" % (
        tuple(_cn) + (len(_rfa),)), rel("f4src"), "Results 3.4")
add("Fig. 4 paired per TF: RF on all PWMs vs the three Slim models, largest two-sided Wilcoxon P (both metrics)", "-", "%.2g" % max(_pp),
    rel("f4src"), "Results 3.4", "text quotes the upper bound")
# paired counts of Sup. Table 3 'Fig4_S4_paired_counts' (fig4_paired_h4.py), all four comparators, both test sets
_fp = pd.read_csv(SRC["f4pair"])
for _ts, _g in _fp.groupby("test_set", sort=False):
    for _m in ("auROC", "auPRC"):
        add("Fig. 4 / S4 paired per TF (Sup. Table 3): RF on all PWMs above %s, %s test, %s (TFs)" % (", ".join(_g.comparator), _ts, _m), "-",
            "%s of %d" % (", ".join("%d" % v for v in _g["n_TF_RF_above_" + _m]), _g.n_TF.iloc[0]), rel("f4pair"), "Results 3.4, Fig. S4 legend",
            "paired per TF, strictly above; two-sided Wilcoxon P per comparator in the same sheet")
_t6 = pd.read_excel(SRC["t6"], sheet_name="RF_training").drop_duplicates("TF").set_index("TF")["Features (monoPWMs)"]
add("Sup. Table 6 benchmark: features (monoPWMs) of the RF matrices, STAT3 / HNF4A / CTCF (N17)", "-", " / ".join("%d" % _t6[t_] for t_ in ("STAT3", "HNF4A", "CTCF")),
    rel("t6"), "N17, Methods 2.11")
from scipy.stats import spearmanr   # noqa: E402
_sr = spearmanr(hm.HM_broc_train_auroc, hm.n_train_pos_H)
add("training auROC of the best monoPWM vs training positives of the TF, human (Spearman rho, P; letter R1-L13)", "-",
    "%+.2f (P = %.1g; %d TFs)" % (_sr[0], _sr[1], len(hm)), rel("hm"), "letter R1-L13",
    "the training auROC of the best monoPWM rises with the number of training positives; no scaler involved")
# ---- operational / process numbers
add("mouse rebuild CPU time (SARUS + negatives + positives)", "-", pd.read_csv(SRC["plan"], dtype=str).set_index("quantity")["value"]["mouse_rebuild_cpu_h"], rel("plan") + " (steps 3-5)", "operational")
add("Zenodo deposit", "10.5281/zenodo.23036475 (old mouse data)", "10.5281/zenodo.14927303 (concept DOI; version DOI of v2.0.0 added on upload)", "Zenodo record", "M33", "data of this version deposited as Zenodo version v2.0.0")
add("P53 sentence of Results 3.3 and the Vyas et al. 2017 reference", "present", "delete (P53 removed)", "-", "M23", "delete sentence and reference")
add("Sup. Table 3 performance rows", "900", str(len(t3p)), rel("t3"), "L03, T04")

# ---- text-pass rows: values the text took from source tables without a NUMBERS row. Each row is made
#      by one function from the current files and, for the h4 column, from the frozen copies at H4_REF.
import io as _io, subprocess as _sp   # noqa: E401,E402
from scipy.stats import wilcoxon as _wx   # noqa: E402
_XP = {"s4sd": ("3_analysis/package/Figures/source_data/Figure_S4_source_data.csv",) * 2, "drv": (SRC["drv"], "3_analysis/transfer_drivers/B5_A_spearman_drivers.csv"),
       "topn": (SRC["topn"], "3_analysis/top_peaks/b3_topN_summary.csv"), "b9pf": ("3_analysis/random_split/b9_random_split_per_fit.csv",) * 2,
       "f3d": ("4_visualisation/figures/Fig3_S3/Figure_3D_source_data.csv",) * 2, "v11": ("4_visualisation/figures/FigS7/hocomoco11_check/subtype_vs_hocomoco11.csv",) * 2,
       "hm": (SRC["hm"], "results/server_h4_test_sets/human/hm_mm_results.csv")}
def _xread(key, ref):
    cur, old = _XP[key]
    if ref is None:
        return pd.read_csv(cur if os.path.isabs(cur) else J(*cur.split("/")))
    rp = os.path.relpath(old, w2.R2) if os.path.isabs(old) else old
    return pd.read_csv(os.path.join(O.ROOT, ref, os.path.basename(rp)))
H4X = {}   # item -> value recomputed from the H4_REF sources (fills the h4 column of these rows)
def addx(item, fn, keys, inv, note=""):
    add(item, "-", fn(lambda k: _xread(k, None)), ", ".join(_XP[k][0] if not os.path.isabs(_XP[k][0]) else rel(k) for k in keys), inv, note)
    H4X[item] = fn(lambda k: _xread(k, H4_REF))
def _s4x(rd, mk, m):
    d = rd("s4sd"); x = d[d.model_key == mk]["d_" + m]
    return "%+.4f; %d of %d; P = %.2g" % (x.median(), (x > 0).sum(), len(x), _wx(x).pvalue)   # scipy default: two-sided, exact for n <= 50
for _mk, _lab in (("slim0", "Slim m=0"), ("slim1", "Slim m=1"), ("lslim5", "LSlim m=-5"), ("munk", "diChIPMunk"), ("RF2f", "RF2f")):
    for _m in ("auROC", "auPRC"):
        addx("Fig. S4 mouse test (human-trained), %s, %s: median delta; TFs above 0; two-sided Wilcoxon signed-rank P (computed directly)" % (_lab, _m),
             lambda rd, mk=_mk, m=_m: _s4x(rd, mk, m), ["s4sd"], "Results 3.4, Fig. S4 legend")
def _drvx(rd):
    d = rd("drv"); d = d[d.p < 0.05].sort_values("p")
    return "%d tests: " % len(d) + "; ".join("%s vs %s %+.2f (P = %.2g)" % (r_.driver, r_.outcome, r_.rho, r_.p) for _, r_ in d.iterrows())
addx("drivers: every test with nominal P < 0.05 (driver vs outcome, rho, P)", _drvx, ["drv"], "Results 3.3, Sup. Table 9")
for _v, _vl in (("allneg", "all negatives"), ("rorig", "negatives at the original ratio")):
    for _m in ("auroc", "auprc"):
        addx("top-N human: gain on the top 500 minus gain on as many random positives, %s, %s (median; two-sided P)" % (_vl, {"auroc": "auROC", "auprc": "auPRC"}[_m]),
             lambda rd, v=_v, m=_m: (lambda x: "%+.4f; P = %.2g" % (x["top_minus_random_%s_med" % m], x["p_top_vs_random_%s" % m]))(
                 rd("topn").astype({"N": str}).set_index(["species", "variant", "subset", "N"]).loc[("H", v, "avail", "500")]), ["topn"], "Results 3.2, letter R1-L6")
def _n5k(rd):
    h = rd("hm").set_index("TF").n_train_pos_H
    lo, hi = h[h < 5000].sort_values().iloc[-1:], h[h > 5000].sort_values().iloc[:1]
    return "%s %s / %s %s" % (lo.index[0], format(int(lo.iloc[0]), ","), hi.index[0], format(int(hi.iloc[0]), ","))
addx("human training positives nearest to 5,000 (below / above)", _n5k, ["hm"], "letter R1-L13")
def _b9x(rd):
    d = rd("b9pf"); win = lambda g: ((g.rf_auroc > g.pwm_auroc) & (g.rf_auprc > g.pwm_auprc))
    g = d.groupby(["TF", "split"])[["rf_auroc", "pwm_auroc", "rf_auprc", "pwm_auprc"]].mean()
    return "%d of %d (TF x split, mean of seeds); %d of %d (per seed)" % (win(g).sum(), len(g), win(d).sum(), len(d))
addx("random split retraining: fits where the RF beats the best monoPWM on auROC and auPRC", _b9x, ["b9pf"], "letter R1-F (all ten fits)")
def _f3dx(rd, tf):
    s = rd("f3d"); v = rd("v11")[["TF", "pwm_name", "best_own_v11", "r_best_own", "r_0.A"]]
    m = s[s.TF == tf].merge(v, on=["TF", "pwm_name"], how="left")
    if m.best_own_v11.isna().any():
        raise SystemExit("CHECK FAILED Fig. 3D %s: a PWM without its HOCOMOCO v11 match" % tf)
    return "; ".join("%s %d: %s %s, %s r %.3f (subtype 0 r %.3f)" % (r_.pwm_class, r_.rank_within_class, r_.feature, r_.dataset, r_.best_own_v11, r_.r_best_own, r_["r_0.A"])
                     for _, r_ in m.sort_values("panel").iterrows())
for _tf in ("RXRA", "SRF", "ANDR"):
    addx("Fig. 3D %s: PWMs shown (class rank: feature, peak set, nearest own HOCOMOCO v11 model, r; r to subtype 0)" % _tf, lambda rd, tf=_tf: _f3dx(rd, tf),
         ["f3d", "v11"], "Fig. 3D legend")
# ---- RF on monoPWMs + diPWMs vs the best single diPWM (selected on
#      the training set), human test set, from Sup. Table 3; Table 1 TPR cells rounded once from the unrounded low-FPR medians
def _vsdi(rd, m):
    p = rd("t3perf"); rf = p[(p.Model == "RandomForestClassifier") & (p["PWM type"] == "mono+di")].set_index("TF_name")
    sd = p[p.Model == "Single best di PWM"].drop_duplicates("TF_name").set_index("TF_name").loc[rf.index]
    if not np.allclose(sd[m + "_test_H"], rf[m + "_test_H_PWM_di"]):
        raise SystemExit("CHECK FAILED Sup. Table 3: best single diPWM differs between its own rows and the RF rows")
    x = rf[m + "_test_H"] - rf[m + "_test_H_PWM_di"]
    return "%d of %d; median %+.4f; P = %.2g" % ((x > 0).sum(), len(x), x.median(), _wx(x).pvalue)
_XP["t3perf"] = ("3_analysis/package/Sup_Tables/Sup_Table_3_ArChIPelago_performance.csv",) * 2   # = Sup. Table 3 sheet performance (AUC columns; tracked in git)
for _m, _ml in (("roc_auc", "auROC"), ("pr_auc", "auPRC")):
    addx("RF on monoPWMs + diPWMs vs the best single diPWM, human test set, %s: TFs improved; median delta; two-sided Wilcoxon P" % _ml,
         lambda rd, m=_m: _vsdi(rd, m), ["t3perf"], "Abstract, Results 3.2 (G4-9)")
def _t1tpr(rd):
    b = rd("b2s"); x = b[(b.scenario == "HH") & (b.rf_seed == 0) & (b.metric == "tpr_at_fpr0.01")].iloc[0]
    return "%.3f %.3f %d/%d %.0e" % (x.base_median, x.rf_median, x.n_improved, x.n_TF, x.wilcoxon_p)
_XP["b2s"] = ("3_analysis/low_fpr_metrics/B2_summary.csv",) * 2
addx("Table 1, Human, TPR at FPR 0.01 (RF seed 0): best monoPWM, RF median (3 decimals, rounded once), TFs improved, Wilcoxon P", _t1tpr, ["b2s"],
     "Table 1 (TPR row)")
w2r = pd.DataFrame(rows)
b12 = pd.read_csv(SRC["b12"]).rename(columns={"old value (36 TFs)": "old text value", "new value (35 TFs)": "new value", "source table": "source file"})
b12["statement_note"] = ""
for k in ("Sup. Table 2 total, Human pooled peak set", "Sup. Table 2 median per TF, Human pooled peak set"):
    m = b12.item.str.contains(k, regex=False)
    b12.loc[m, "new value"] = "removed; Methods describes the pooling without a count"
    b12.loc[m, "statement_note"] = "pooled peak count removed (authors' decision)"
pend = b12["new value"].astype(str).str.startswith("pending")
# Fig. 4 rows filled; a median whose sign differs from the first submission is noted
for i_ in b12.index[b12.item.str.startswith("Fig. 4") & b12.item.str.endswith("median delta")]:
    o_, n_ = float(b12.loc[i_, "old text value"]), float(b12.loc[i_, "new value"])
    if (o_ > 0) != (n_ > 0):
        b12.loc[i_, "statement_note"] = "sign of the median changes (H3b): %s -> %s" % (b12.loc[i_, "old text value"], b12.loc[i_, "new value"])
b12.loc[pend & b12.item.str.startswith(("Random", "Random-split")), "statement_note"] = "pending (H3b): B9 random split (Fig. S12) retrained on negatives v3"
b12.loc[pend & b12.item.str.startswith("Fig. 4"), "statement_note"] = "pending (H3b): Fig. 4 models retrained on negatives v3"
# random split: every number of Fig. S12 / Sup. Table 7 changes; the quoted shares change meaning (positives read once)
for i_ in b12.index[b12.item.str.startswith(("Random", "Random-split")) & ~pend]:
    b12.loc[i_, "statement_note"] = "statement changes (H3b, v3 negatives, positives read once): %s -> %s" % (b12.loc[i_, "old text value"], b12.loc[i_, "new value"])
if (b12.loc[pend, "statement_note"] == "").any():
    raise SystemExit("CHECK FAILED pending rows without a note: %s" % list(b12.item[pend & (b12.statement_note == "")]))
m = b12.item.str.startswith("Ablated TFs where the loss exceeds 3 sd")
for i_ in b12.index[m]:
    b12.loc[i_, "statement_note"] = "new: %s" % b12.loc[i_, "new value"]
# human rows: 'all TFs improved' that no longer holds (old n of N with n = N, new n < N)
for i_, r_ in b12[b12.item.str.contains("TFs improved") & (b12.statement_note == "")].iterrows():
    mo, mn = re.match(r"(\d+)/(\d+)$", str(r_["old text value"])), re.match(r"(\d+)/(\d+)$", str(r_["new value"]))
    if mo and mn and mo.group(1) == mo.group(2) and mn.group(1) != mn.group(2):
        b12.loc[i_, "statement_note"] = "statement changes: not every TF improved (%s of %s)" % (mn.group(1), mn.group(2))
b12["source file"] = "B12: " + b12["source file"].astype(str)
out = pd.concat([b12, w2r], ignore_index=True)[["item", "old text value", "new value", "source file", "inventory row", "statement_note"]]
# rows that no text quotes (the verifier's reverse check) are marked, not deleted
NOT_IN_TEXT = ["median mouse test positives per TF", "smallest mouse test sets", "mouse rebuild CPU time (SARUS + negatives + positives)",
               "Sup. Table 2 total, Human pooled peak set (peaks, all chromosomes)",
               "Sup. Table 2 median per TF, Human pooled peak set (peaks, all chromosomes)",
               "median pooled peak set per TF and species", "median pooled peak set per TF, human", "median pooled peak set per TF, mouse"]
_m = out.item.isin(NOT_IN_TEXT)
if _m.sum() != len(NOT_IN_TEXT):
    raise SystemExit("CHECK FAILED not-in-text items missing: %s" % sorted(set(NOT_IN_TEXT) - set(out.item)))
out.loc[_m, "statement_note"] = ("not in text; " + out.loc[_m, "statement_note"].fillna("").astype(str)).str.rstrip("; ")
# 'new value' = the final run; the value of the same item in the earlier table (frozen copy at H3A_REF) goes next to it to show
# which values moved. Items that occur twice (human and mouse rows) are matched by item and occurrence.
import io, subprocess   # noqa: E401,E402
_h3 = pd.read_csv(os.path.join(O.ROOT, H3A_REF, "NUMBERS_R2_35TF.csv"))
# items are matched without their chromosome list in parentheses, which follows the training split: labels of earlier runs name
# their own training chromosomes, labels of the final run name the final ones
_chrfree = lambda it: re.sub(r"\s*\((?:%; )?chr[^)]*\)", "", str(it))
_key = lambda x: list(zip(x.item.map(_chrfree), x.groupby(x.item.map(_chrfree)).cumcount()))
if out.item.map(_chrfree).duplicated().sum() != out.item.duplicated().sum():
    raise SystemExit("CHECK FAILED NUMBERS: two items differ only in their chromosome list")
_h3v = dict(zip(_key(_h3), _h3["new value"]))
# the value of the same item in the earlier table (frozen copy at H4_REF) next to the new value, as for the h3a column; for the final run
# 'h4 vs h3a' compares the h4 column with the h3a column, for the earlier run the new value with the h3a column
_h4 = pd.read_csv(os.path.join(O.ROOT, H4_REF, "NUMBERS_R2_35TF.csv"))
_h4v = dict(zip(_key(_h4), _h4["new value"]))
_cmp = lambda ref, k, v: "new row" if k not in ref else ("same" if str(ref[k]) == str(v) else "changed")
# items renamed when their H>M value moved from seed 0 to the two-seed mean keep the h4 value of the old label
H4_ALIAS = {"ERG: mouse-trained gain on the mouse test set, auROC / auPRC (mean of two seeds); human-trained dauPRC H>M (mean of two seeds)":
            "ERG: mouse-trained gain on the mouse test set, auROC / auPRC (mean of two seeds); human-trained dauPRC H>M (seed 0)"}
for _new, _old in H4_ALIAS.items():
    if (_chrfree(_old), 0) not in _h4v or _new not in set(out.item):
        raise SystemExit("CHECK FAILED H4_ALIAS: %r or %r not found" % (_old, _new))
    _h4v[(_chrfree(_new), 0)] = _h4v[(_chrfree(_old), 0)]
out.insert(3, "h4 value", [_h4v.get(k, "") for k in _key(out)])
for _i in out.index[out.item.isin(list(H4X))]:   # text-pass rows: h4 value recomputed from the H4_REF sources
    if out.loc[_i, "h4 value"] == "":
        out.loc[_i, "h4 value"] = H4X[out.loc[_i, "item"]]
        _h4v[(_chrfree(out.loc[_i, "item"]), 0)] = H4X[out.loc[_i, "item"]]
out.insert(4, "new vs h4", [_cmp(_h4v, k, v) for k, v in zip(_key(out), out["new value"])])
out.insert(5, "h3a value", [_h3v.get(k, "") for k in _key(out)])
_vs = out["h4 value"] if H.STAGE == "H5" else out["new value"]
out.insert(6, "h4 vs h3a", [_cmp(_h3v, k, v) for k, v in zip(_key(out), _vs)])
f = J("outputs/numbers/NUMBERS_R2_35TF.csv")
w2.backup(f); out.to_csv(f, index=False)
print("rows: B12 %d + W2 %d = %d" % (len(b12), len(w2r), len(out)))
w2.provenance(B13, "3_analysis/integration/make_numbers_w2.py", [f], list(SRC.values()))
