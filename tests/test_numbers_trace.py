"""Value tracing of outputs/numbers/NUMBERS_R2_35TF.csv: every 'new value' whose source is a table matches that table.

Rule per row:
  1. source paths are parsed from 'source file' (a '(_backup_<date>)' suffix names the old value's file, so it is removed; 'PKG/' is
     the released package results/round1_package/Manuscript_analysis, 'R1/' the first-submission folder
     results/round1_package/Resubmission; 'old: ...; new: ...' prefixes are ignored);
  2. rows recomputed here (DERIVED, the number is a median, count, range or test over table rows, not a cell) must equal the
     recomputation at the printed precision;
  3. every other number in 'new value' must occur in one of its source tables at the printed precision (half a unit of the
     last printed digit; percentages also match fractions x 100; numbers inside text cells count);
  4. rows without a table (NO_TABLE: flagged, removed or deleted items) and rows marked 'superseded' in statement_note
     (replaced by a later row, not quoted in the text) are not traced; both lists are asserted, so a new untraced row fails;
  5. rows whose statement_note starts with 'pending' wait for the Fig. 4 / Fig. S4 / random-split (Fig. S12)
     retraining on negatives v3 and have no new value yet; they are not traced, and only those three items may be pending.
"""
import os
import re

import numpy as np
import pandas as pd
import pytest

import t1_helpers as T

pytestmark = pytest.mark.r2
PKG = os.path.join(T.R2, "results/round1_package/Manuscript_analysis")
PATH = re.compile(r"([A-Za-z0-9_./()\-]+?\.(?:csv|xlsx|txt|md|tab))")
NUM = re.compile(r"(?<![\w.])[-+]?(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?(?:e[-+]?\d+)?(?![\w])")
PENDING_OK = ("pending (H3b): Fig. 4 ", "pending (H3b): Fig. S4 ", "pending (H3b): B9 random split")
NO_TABLE = {"median pooled peak set per TF and species", "median pooled peak set per TF, human",
            "median pooled peak set per TF, mouse", "P53 sentence of Results 3.3 and the Vyas et al. 2017 reference",
            "Fig. S7 subtype match: method (V1-S7)",   # a method description; its source is the producer script
            "Zenodo deposit"}                          # the DOI of the data deposit; no table holds it


def J(*p):
    return os.path.join(T.R2, *p)


def numbers_file():
    d = pd.read_csv(J("outputs/numbers/NUMBERS_R2_35TF.csv"))
    if T.mutate("numbers_value_changed"):
        i = d.index[d.item == "Human RF mono+di median auPRC, ArChIPelago RF"][0]
        d.loc[i, "new value"] = "0.2950"
    if T.mutate("numbers_source_missing"):
        d.loc[d.item == "Human RF mono+di median auPRC, ArChIPelago RF", "source file"] = "3_analysis/integration/inputs/no_such_table.csv"
    if T.mutate("numbers_pending_other"):
        d.loc[d.item == "Human RF mono+di median auPRC, ArChIPelago RF", "statement_note"] = "pending (H3b): something else"
    if T.mutate("numbers_derived_changed"):
        d.loc[d.item == "ChIP-Seq datasets, mouse", "new value"] = "275"
    return d


def source_files(src):
    out = []
    for m in PATH.findall(str(src)):
        m = m.replace("(_backup_20261106)", "")
        if m.startswith("PKG/"):
            out.append(os.path.join(PKG, m[4:]))
        elif m.startswith("R1/"):
            out.append(os.path.join(T.R2, "results", "round1_package", "Resubmission", m[3:]))
        else:
            out.append(J(m))
    return [f for f in out if os.path.exists(f)]


_cache = {}


def table_numbers(path):
    if path not in _cache:
        nums = []
        if path.endswith(".xlsx"):
            cells = [v for t in pd.read_excel(path, sheet_name=None, header=None).values() for v in t.values.ravel() if pd.notna(v)]
        elif path.endswith((".csv", ".tab")):
            sep = "\t" if "\t" in open(path).readline() else ","
            t = pd.read_csv(path, sep=sep, header=None, dtype=str, keep_default_na=False)
            cells = []
            for v in t.values.ravel():
                try:
                    cells.append(float(v))
                except ValueError:
                    cells.append(v)
        else:
            cells = [open(path).read()]
        for v in cells:
            if isinstance(v, (int, float, np.integer, np.floating)):
                nums.append(float(v))
            else:
                nums += [float(x.replace(",", "")) for x in NUM.findall(str(v))]
        _cache[path] = np.array(nums)
    return _cache[path]


def tol(tok):
    t = tok.lstrip("+-").replace(",", "")
    if "e" in t:
        mant, ex = t.split("e")
        dec = len(mant.split(".")[1]) if "." in mant else 0
        return 0.5 * 10 ** (-dec) * 10 ** int(ex)
    dec = len(t.split(".")[1]) if "." in t else 0
    return 0.5 * 10 ** (-dec) + 1e-12


def tokens(v):
    return [t for t in NUM.findall(str(v))]


def close(tok, x):
    return abs(float(tok.replace(",", "")) - float(x)) <= tol(tok)


# ------------------------------------------------------------------ recomputations of derived numbers
def _rf():
    t = pd.read_csv(J("3_analysis/integration", "inputs", "results_table.csv"), sep="\t")
    return t[(t.Model == "RandomForestClassifier") & (t.PWM == "mono+di")].set_index("TF_name")


def _wilcox(a):
    from scipy.stats import wilcoxon
    return wilcoxon(a).pvalue


def _op_h():   # the rebuilt Sup. Table 4 (per_TF, human rows)
    op = pd.read_excel(J("outputs/supplementary_tables", "Sup_Table_4_operational_metrics.xlsx"), sheet_name="per_TF")
    op = op[op["Test set"].str.startswith("human")]
    assert set(op.TF) == set(T.TF35), sorted(set(op.TF) ^ set(T.TF35))
    return op["FP reduction at 50 % recall, ArChIPelago vs best monoPWM (%)"]


def _b9(col):
    b = pd.read_csv(J("3_analysis", "random_split", "b9_summary.csv")).pivot(index="TF", columns="split", values=col)
    d = b["random"] - b["chromosome"]
    return [d.min(), d.max()]


def _b9d(col):
    b = pd.read_csv(J("3_analysis", "random_split", "b9_summary.csv")).pivot(index="TF", columns="split", values=col)
    return b["random"] - b["chromosome"]


def _t1():
    return pd.read_excel(J("outputs/supplementary_tables", "Sup_Table_1_ChIP-Seq_experiments.xlsx"), sheet_name=None)


def _top(sp):   # top-1,000 check of the training and test positives
    d = pd.read_csv(J("3_analysis/integration", "top1000_check_h3.csv"))
    d = d[d.species == sp]
    return [100 * d.n_in_top1000.sum() / d.n_positives.sum(), d.n_in_top1000.sum(), d.n_positives.sum()]


def _negv3(sp, fs):   # negatives v3: removed own-positive overlaps, v2 total, v3 total
    d = pd.read_csv(T.h3in().IN["neg_v3_checks"])   # not in the previous-run registry (no model; h3_inputs NOT_IN_H4)
    g = d[(d.species == sp) & (d.file_split == fs)]
    v2, v3 = g[g.version == "v2"], g[g.version == "v3"]
    if not len(v2):   # final run: v3 lists built directly; own-positive overlaps of v3, v3 total
        return [v3.overlap_own_pos.sum(), v3.n.sum()]
    return [v2.overlap_own_pos.sum(), v2.n.sum(), v3.n.sum()]


def _chroms(sp):
    """Chromosome lists of the split as build_positives_sym.py builds them (literal list of the code, parsed with ast):
    [training..., test..., unused autosomes...] as chromosome numbers."""
    import ast
    H = T.h3in()
    if H.STAGE == "H5":   # final run: one split definition, 1_data_preparation/h5_splits.py (LOGICAL, never swapped)
        ns = {}
        exec(open(H.IN["build_positives"]).read(), ns)
        tr, te = ns["LOGICAL"][sp]["train"], ns["LOGICAL"][sp]["test"]
        n = 22 if sp == "HUMAN" else 19
        return [int(c[3:]) for c in tr + te] + [i for i in range(1, n + 1) if "chr%d" % i not in tr + te]
    txt = open(J("B14_human_sym", "code", "build_positives_sym.py")).read()
    blk = txt[txt.index('"%s": dict(' % sp):]
    li = {k: ast.literal_eval(re.search(k + r'=\["chr" \+ x for x in (\[[^\]]+\])\]', blk).group(1)) for k in ("TRAIN_CHR", "TEST_CHR")}
    n = 22 if sp == "HUMAN" else 19
    test, train = (li["TEST_CHR"], li["TRAIN_CHR"]) if sp == "HUMAN" else (li["TRAIN_CHR"], li["TEST_CHR"])   # mouse files swapped
    valid = [str(i) for i in range(1, n + 1)]
    tr = [c for c in train if c in valid]
    return [int(c) for c in tr + list(test) + [c for c in valid if c not in tr and c not in test]]


def _conc():
    return pd.read_csv(J("3_analysis/feature_importance", "tables", "B4_importance_concentration.csv"))


def _mask_line(sp):
    """Top-1,000 windows removed by the repeat mask in the human or mouse set (repeat_mask_counts.csv) -> [removed, total, percent]."""
    v = pd.read_csv(J("results/datasets/repeat_mask", "repeat_mask_counts.csv"), dtype=str).set_index("quantity")["value"]
    a, b = v["removed_" + sp], v["total_" + sp]
    a, b = int(a.replace(",", "")), int(b.replace(",", ""))
    return [a, b, 100 * a / b]


def _drv():
    return pd.read_csv(J("3_analysis/transfer_drivers", "B5_A_spearman_drivers.csv"))


def _s1_mouse_prc():   # Fig. S1: monoPWMs + diPWMs model above the best monoPWM by auPRC, mouse test set, from the results table
    t = pd.read_csv(J("3_analysis/integration", "inputs", "results_table.csv"), sep="\t")
    n = lambda m: int(((lambda x: x.pr_auc_test_M > x.pr_auc_test_M_PWM)(t[(t.Model == m) & (t.PWM == "mono+di")])).sum())
    return [n("RandomForestClassifier"), n("XGBClassifier"), int(t.TF_name.nunique())]


def _f4_paired(metric=None):   # Fig. 4: RF on all PWMs vs Slim m=0, Slim m=1, LSlim m=-5 per TF, human test set (Fig. 4 source data)
    from scipy.stats import wilcoxon
    d = pd.read_csv(T.h3in().IN["fig4_src"]); d = d[d.test_set == "human"]
    rf = d[d.model == "ArChIPelago RF on all PWMs"].set_index("TF")
    out, ps = {}, []
    for m in ("auROC", "auPRC"):
        dd = [rf[m] - d[d.model == s_].set_index("TF").loc[rf.index, m] for s_ in ("Slim m=0", "Slim m=1", "LSlim m=-5")]
        out[m] = [int((x > 0).sum()) for x in dd] + [len(rf)]
        ps += [wilcoxon(x).pvalue for x in dd]
    return out[metric] if metric else [max(ps)]


DERIVED = {
    # counts and a test bound recomputed from the results table and the Fig. 4 source data
    "Fig. S1 mouse test set, monoPWMs + diPWMs: TFs above the best monoPWM by auPRC, RF / XGBoost": _s1_mouse_prc,
    "Fig. 4 paired per TF: RF on all PWMs above Slim m=0, Slim m=1, LSlim m=-5, human test, auROC (TFs)": lambda: _f4_paired("auROC"),
    "Fig. 4 paired per TF: RF on all PWMs above Slim m=0, Slim m=1, LSlim m=-5, human test, auPRC (TFs)": lambda: _f4_paired("auPRC"),
    "Fig. 4 paired per TF: RF on all PWMs vs the three Slim models, largest two-sided Wilcoxon P (both metrics)": _f4_paired,
    # rows that are differences or a test over hm_mm columns (h3_inputs hm_mm)
    "USF2 mouse-trained dauROC on the mouse test set: mean of two seeds / seed 0 (N12)": lambda: (lambda h: [
        h.MM_rf_test_M_auroc_mean["USF2"] - h.MM_broc_test_M_auroc["USF2"], h.MM_rf_test_M_auroc_s0["USF2"] - h.MM_broc_test_M_auroc["USF2"]])(
        pd.read_csv(T.h3in().IN["hm_mm"]).set_index("TF")),
    "training auROC of the best monoPWM vs training positives of the TF, human (Spearman rho, P; letter R1-L13)": lambda: (lambda h, r: [
        r[0], r[1], len(h)])(*(lambda h: (h, __import__("scipy.stats", fromlist=["spearmanr"]).spearmanr(h.HM_broc_train_auroc, h.n_train_pos_H)))(
        pd.read_csv(T.h3in().IN["hm_mm"]))),
    "Human Wilcoxon signed-rank P (two-sided), auROC gain": lambda: [_wilcox(_rf().roc_auc_test_H - _rf().roc_auc_test_H_PWM)],
    "Human Wilcoxon signed-rank P (two-sided), auPRC gain": lambda: [_wilcox(_rf().pr_auc_test_H - _rf().pr_auc_test_H_PWM)],
    "Human median per-TF gain in auPRC": lambda: [(_rf().pr_auc_test_H - _rf().pr_auc_test_H_PWM).median()],
    "Human median per-TF gain in auROC": lambda: [(_rf().roc_auc_test_H - _rf().roc_auc_test_H_PWM).median()],
    "Human FP reduction at 50 % recall, median (%)": lambda: [_op_h().median()],
    "Human FP reduction at 50 % recall, IQR low (%)": lambda: [_op_h().quantile(.25)],
    "Human FP reduction at 50 % recall, IQR high (%)": lambda: [_op_h().quantile(.75)],
    "Human test class ratio, median negatives per positive": lambda: [
        (lambda b: (b.n_neg / b.n_pos).median())(pd.read_csv(J("3_analysis", "low_fpr_metrics", "B2_metrics_wide.csv")).query("scenario == 'HH' and model == 'RF_s0'"))],
    "Top PWM is a diPWM (n of N)": lambda: [(_conc().top1_kind == "di").sum(), len(_conc())],
    "Feature PWMs (gini_all_PWMs rows)": lambda: [_conc().n_features.sum()],
    "TFs with a non-primary subtype (>= 2 subtypes at cut 0.8; Fig. S10D)": lambda: [
        (pd.read_csv(J("3_analysis/feature_importance", "tables", "B4_subtype_summary_per_TF.csv")).n_subtypes_80 >= 2).sum()],
    "Random minus chromosome split, RF auROC (range)": lambda: _b9("rf_auroc"),
    "Random minus chromosome split, best monoPWM auROC (range)": lambda: _b9("pwm_auroc"),
    "Random minus chromosome split, auROC gain (range)": lambda: _b9("gain_auroc"),
    "Random split: largest rise of the auPRC gain": lambda: [max(_b9("gain_auprc"))],
    "Random split raised the auPRC gain (n of N TFs)": lambda: [int((_b9d("gain_auprc") > 0).sum()), len(_b9d("gain_auprc"))],
    "Human datasets per TF, median (used sets)": lambda: [
        (lambda h: h[h.TF != "P53"].groupby("TF").size().median())(pd.read_csv(J("results/datasets/peak_ids", "human_used_peak_ids.csv")))],
    "ChIP-Seq datasets used, human + mouse": lambda: [len(_t1()["human"]) + len(_t1()["mouse"]), len(_t1()["human"]), len(_t1()["mouse"])],
    "ChIP-Seq datasets, human": lambda: [len(_t1()["human"])],
    "Human ChIP-Seq datasets used (Sup. Table 1)": lambda: [len(_t1()["human"])],
    "Datasets of both species (manuscript total)": lambda: [len(_t1()["human"]) + len(_t1()["mouse"])],
    "mouse repeat mask": lambda: [0.7] + _mask_line("mouse")[:2] + [1000] + _mask_line("mouse")[2:],
    "human windows removed by the mask (B8 rerun)": lambda: _mask_line("human"),
    "ChIP-Seq datasets, mouse": lambda: [len(_t1()["mouse"])],
    "orthologous TF pairs": lambda: [_t1()["human"].TF.nunique()],
    "minimum datasets per TF, mouse": lambda: [_t1()["mouse"].groupby("TF").size().min()],
    "minimum datasets per TF, human": lambda: [_t1()["human"].groupby("TF").size().min()],
    "Sup. Table 1 rows (sets used)": lambda: [len(_t1()["human"]), len(_t1()["mouse"])],
    "human positives (training and test) among the 1,000 most significant peaks of a listed set": lambda: _top("H"),
    "mouse positives (training and test) among the 1,000 most significant peaks of a listed set": lambda: _top("M"),
    "chromosomes used, human (R1-F6)": lambda: _chroms("HUMAN"),
    "chromosomes used, mouse (R1-F6)": lambda: _chroms("MOUSE"),
    "drivers of transfer: tests, n": lambda: [len(_drv()), _drv().n.max()],
    "drivers: nominal P < 0.05 (daggers)": lambda: [(_drv().p < 0.05).sum(), len(_drv())],
    "Sup. Table 3 performance rows": lambda: [len(pd.read_excel(J("outputs/supplementary_tables", "Sup_Table_3_ArChIPelago_performance.xlsx"), sheet_name="performance"))],
}
for _sp, _fs, _lab in (("HUMAN", "train", "human training"), ("HUMAN", "control", "human test"),
                       ("MOUSE", "control", "mouse training"), ("MOUSE", "train", "mouse test")):   # mouse files swapped
    DERIVED["negatives v3, %s: v2 negatives within 300 bp of the TF's own positives (removed); v3 negatives" % _lab] = (
        lambda sp=_sp, fs=_fs: _negv3(sp, fs))


def _gate_pass():
    g = pd.read_csv(J("3_analysis", "low_fpr_metrics", "B2_gate_h3_summary.csv"))
    g = g[g.decision != "size only"]
    return [(g.decision == "pass").sum(), len(g)]


DERIVED["B2 gate: checks passed (recomputed from the score npz vs hm_mm, b2_v2, results table)"] = _gate_pass


def _porig():   # PWM origin table (3_analysis/integration/pwm_origin_h3.csv, checked against the archives in test_pwm_origin.py)
    return pd.read_csv(J("3_analysis/integration", "pwm_origin_h3.csv")).groupby("species").sum(numeric_only=True)


def _gc_row(sp, sl, t=5):   # per-TF GC difference (points) recomputed from the v3 GC check, not from gc_per_tf_h3.csv
    g = pd.read_csv(T.h3in().IN["gc_v3"])   # not in the previous-run registry (no model; h3_inputs NOT_IN_H4)
    g = g[g.pool == T.h3in().GC_POOL[sp[0].upper()]]
    g = g[g.chromosomes.isin(["chr1/8/21", "chr1/8/19"]) == (sl == "test")]
    d = 100 * g.median_gc_diff
    return [d.min(), d.max(), int((d.abs() > t).sum()), len(d)]


for _sp in ("human", "mouse"):
    for _sl in ("test", "training"):
        DERIVED["GC content, %s %s sets: per-TF median GC of positives minus negatives (percentage points), range; TFs with |difference| > 5 points (V1-16)"
                % (_sp, _sl)] = (lambda sp=_sp, sl=_sl: _gc_row(sp, sl))


def _below():
    rf = _rf()
    b = rf[(rf.roc_auc_test_H < rf.roc_auc_test_H_PWM) | (rf.pr_auc_test_H < rf.pr_auc_test_H_PWM)]
    return [x for _, r in b.iterrows() for x in (r.roc_auc_test_H - r.roc_auc_test_H_PWM, r.pr_auc_test_H - r.pr_auc_test_H_PWM)]


DERIVED["Fig. 3C open diamond: RF (monoPWMs + diPWMs) below the best monoPWM on auROC or auPRC, human test set (V1-19)"] = _below
DERIVED["ablation: seeds behind the 3 SD rule (random removals per TF; V1-25)"] = lambda: [
    pd.read_csv(J("3_analysis/feature_importance", "tables", "B4_ablation_all_fits.csv")).query("cond == 'drop_random'").groupby("TF").seed.nunique().min()]
def _s7_counts():   # Fig. S7 legend counts recomputed from the source data and the verdicts
    s = pd.read_csv(J("4_visualisation/figures", "FigS7", "Figure_S7_source_data.csv"))
    v = pd.read_csv(J("4_visualisation/figures", "FigS7", "subtype_verdicts_v11_h3.csv"))
    d = s[s.TF.isin(v.TF[v.drawn])]
    return [int((s[s.TF == s.TF.iloc[0]].pwm_class == "mono").sum()), int((s[s.TF == s.TF.iloc[0]].pwm_class == "di").sum()),
            s.TF.nunique(), int(v.drawn.sum()), len(d), int(d.cofactor.sum())]


DERIVED["Fig. S7 legend counts: monoPWMs / diPWMs per TF; candidate TFs; TFs drawn; tiles drawn; co-factor tiles"] = _s7_counts
DERIVED["PWM source datasets (Methods 2.1)"] = lambda: [_porig().loc["human", "n_source_sets"], _porig().loc["mouse", "n_source_sets"]]
DERIVED["PWMs derived from one of the TF's listed positive sets (same species; D9)"] = lambda: [
    _porig().loc["human", "n_from_listed_sets"], _porig().loc["human", "n_PWMs"], 100 * _porig().loc["human", "n_from_listed_sets"] / _porig().loc["human", "n_PWMs"],
    _porig().loc["mouse", "n_from_listed_sets"], _porig().loc["mouse", "n_PWMs"], 100 * _porig().loc["mouse", "n_from_listed_sets"] / _porig().loc["mouse", "n_PWMs"]]
# rows of values the text took from source tables, recomputed here from the files
def _s4w(mk, m):   # Fig. S4 mouse test: median delta, TFs > 0, N, two-sided Wilcoxon signed-rank P computed on the deltas
    d = pd.read_csv(J("3_analysis", "package", "Figures", "source_data", "Figure_S4_source_data.csv"))
    x = d.loc[d.model_key == mk, "d_" + m].to_numpy()
    if T.mutate("s4_onesided"):   # red run: the one-sided P (the 'p' column of Figure_S4_stats.csv) quoted as two-sided
        from scipy.stats import wilcoxon
        return [np.median(x), int((x > 0).sum()), len(x), wilcoxon(x, alternative="greater").pvalue]
    return [np.median(x), int((x > 0).sum()), len(x), _wilcox(x)]


def _drv5():
    d = pd.read_csv(J("3_analysis/transfer_drivers", "B5_A_spearman_drivers.csv"))
    d = d[d.p < 0.05].sort_values("p")
    return [len(d)] + [v for _, r in d.iterrows() for v in (r.rho, r.p)]


def _drv5_names():
    d = pd.read_csv(J("3_analysis/transfer_drivers", "B5_A_spearman_drivers.csv"))
    return ["%s vs %s" % (r.driver, r.outcome) for _, r in d[d.p < 0.05].iterrows()]


def _top_rnd(v, m):
    t = pd.read_csv(J("3_analysis", "top_peaks", "b3_topN_summary.csv"), dtype={"N": str})
    x = t[(t.species == "H") & (t.variant == v) & (t.subset == "avail") & (t.N == "500")].iloc[0]
    return [x["top_minus_random_%s_med" % m], x["p_top_vs_random_%s" % m]]


def _near5k():
    h = pd.read_csv(T.h3in().IN["hm_mm"]).set_index("TF").n_train_pos_H
    lo, hi = h[h < 5000].idxmax(), h[h > 5000].idxmin()
    return (lo, hi), [h[lo], h[hi]]


def _b9wins():
    d = pd.read_csv(J("3_analysis", "random_split", "b9_random_split_per_fit.csv"))
    g = d.groupby(["TF", "split"]).mean(numeric_only=True)
    w = lambda x: int(((x.rf_auroc > x.pwm_auroc) & (x.rf_auprc > x.pwm_auprc)).sum())
    return [w(g), len(g), w(d), len(d)]


def _f3d(tf):
    s = pd.read_csv(J("4_visualisation/figures", "Fig3_S3", "Figure_3D_source_data.csv"))
    v = pd.read_csv(J("4_visualisation/figures", "FigS7", "hocomoco11_check", "subtype_vs_hocomoco11.csv")).set_index(["TF", "pwm_name"])
    out, names = [], []
    for _, r in s[s.TF == tf].sort_values("panel").iterrows():
        x = v.loc[(tf, r.pwm_name)]
        out += [r.rank_within_class, x.r_best_own, 0, x["r_0.A"]]
        names += [r.feature + " ", r.dataset, x.best_own_v11]
    return out, names


def _vs_di(m):   # RF mono+di vs the best single diPWM, human test; read from the published Sup. Table 3
    p = pd.read_excel(J("outputs/supplementary_tables", "Sup_Table_3_ArChIPelago_performance.xlsx"), sheet_name="performance")
    rf = p[(p.Model == "RandomForestClassifier") & (p["PWM type"] == "mono+di")].set_index("TF_name")
    ref = "Single best mono PWM" if T.mutate("vs_di_mono") else "Single best di PWM"   # red run: the best monoPWM as reference
    sd = p[p.Model == ref].drop_duplicates("TF_name").set_index("TF_name").loc[rf.index]
    x = (rf[m + "_test_H"] - sd[m + "_test_H"]).to_numpy()
    return [int((x > 0).sum()), len(x), np.median(x), _wilcox(x)]


def _t1_tpr():   # Table 1 TPR row: medians rounded once from B2_summary.csv
    b = pd.read_csv(J("3_analysis", "low_fpr_metrics", "B2_summary.csv"))
    x = b[(b.scenario == "HH") & (b.rf_seed == 0) & (b.metric == "tpr_at_fpr0.01")].iloc[0]
    return [round(x.base_median, 3), round(x.rf_median, 3), x.n_improved, x.n_TF, x.wilcoxon_p]


_S4X = "Fig. S4 mouse test (human-trained), %s, %s: median delta; TFs above 0; two-sided Wilcoxon signed-rank P (computed directly)"
for _mk, _lab in (("slim0", "Slim m=0"), ("slim1", "Slim m=1"), ("lslim5", "LSlim m=-5"), ("munk", "diChIPMunk"), ("RF2f", "RF2f")):
    for _m in ("auROC", "auPRC"):
        DERIVED[_S4X % (_lab, _m)] = lambda mk=_mk, m=_m: _s4w(mk, m)
DERIVED["drivers: every test with nominal P < 0.05 (driver vs outcome, rho, P)"] = _drv5
for _v, _vl in (("allneg", "all negatives"), ("rorig", "negatives at the original ratio")):
    for _m, _ml in (("auroc", "auROC"), ("auprc", "auPRC")):
        DERIVED["top-N human: gain on the top 500 minus gain on as many random positives, %s, %s (median; two-sided P)" % (_vl, _ml)] = (
            lambda v=_v, m=_m: _top_rnd(v, m))
DERIVED["human training positives nearest to 5,000 (below / above)"] = lambda: _near5k()[1]
DERIVED["random split retraining: fits where the RF beats the best monoPWM on auROC and auPRC"] = _b9wins
for _m, _ml in (("roc_auc", "auROC"), ("pr_auc", "auPRC")):
    DERIVED["RF on monoPWMs + diPWMs vs the best single diPWM, human test set, %s: TFs improved; median delta; two-sided Wilcoxon P" % _ml] = (
        lambda m=_m: _vs_di(m))
DERIVED["Table 1, Human, TPR at FPR 0.01 (RF seed 0): best monoPWM, RF median (3 decimals, rounded once), TFs improved, Wilcoxon P"] = _t1_tpr


def _hm_mean_medians(m, b):   # H>M medians of the RF (mean of seeds 0 and 1) and of the best monoPWM, from hm_mm of the stage
    h = pd.read_csv(T.h3in().IN["hm_mm"])
    return [h["HM_rf_test_M_%s_mean" % m].median(), h["HM_%s_test_M_%s" % (b, m)].median()]


def _exact_row(base, rf, better):   # medians, TFs better, n, exact two-sided Wilcoxon P on the non-zero differences
    from scipy.stats import wilcoxon
    d = rf - base.loc[rf.index]
    nz = d[d != 0]
    p = wilcoxon(nz, method="exact").pvalue
    if T.mutate("table1_exact_p"):   # red run: scipy's default (normal approximation with the zero difference)
        p = wilcoxon(d).pvalue
    return [base.median(), rf.median(), int(better(d).sum()), len(d), p]


def _t1_prec40():   # recomputed from Sup. Table 7 per_TF (the sheet the Table 1 row cites)
    w = pd.read_excel(J("outputs/supplementary_tables", "Sup_Table_7_low_FPR_metrics.xlsx"), sheet_name="per_TF")
    w = w[w.scenario == "HH"]
    rf, mo = (w[w.model == m_].set_index("TF")["prec_at_tpr0.4"] for m_ in ("RF_s0", "mono_byprc"))
    return _exact_row(mo, rf, lambda d: d > 0)


def _t1_fp50():   # recomputed from Sup. Table 4 per_TF, human rows
    o = pd.read_excel(J("outputs/supplementary_tables", "Sup_Table_4_operational_metrics.xlsx"), sheet_name="per_TF")
    o = o[o["Test set"].str.startswith("human")].set_index("TF")
    return _exact_row(o["FP per 10,000 regions at 50 % recall: best monoPWM"], o["FP per 10,000 regions at 50 % recall: ArChIPelago RF"],
                      lambda d: d < 0)


DERIVED["H>M RF median auROC / best monoPWM"] = lambda: _hm_mean_medians("auroc", "broc")
DERIVED["H>M RF median auPRC / best monoPWM"] = lambda: _hm_mean_medians("auprc", "bprc")
DERIVED["Table 1, Human, precision at 40 % recall (RF seed 0 vs best monoPWM by auPRC): medians, TFs improved, exact two-sided Wilcoxon P (V192)"] = _t1_prec40
DERIVED["Table 1, Human, false positives per 10,000 regions at 50 % recall (RF seed 0 vs best monoPWM by auPRC): medians, TFs with fewer, exact two-sided Wilcoxon P (V193)"] = _t1_fp50
_F3D = "Fig. 3D %s: PWMs shown (class rank: feature, peak set, nearest own HOCOMOCO v11 model, r; r to subtype 0)"
for _tf in ("RXRA", "SRF", "ANDR"):
    DERIVED[_F3D % _tf] = lambda tf=_tf: _f3d(tf)[0]
DERIVED_TEXT = {   # also check the TF names printed with a derived count
    "minimum datasets per TF, mouse": lambda: sorted(k for k, v in _t1()["mouse"].groupby("TF").size().items() if v == _t1()["mouse"].groupby("TF").size().min()),
    "minimum datasets per TF, human": lambda: sorted(k for k, v in _t1()["human"].groupby("TF").size().items() if v == _t1()["human"].groupby("TF").size().min()),
    "drivers: every test with nominal P < 0.05 (driver vs outcome, rho, P)": _drv5_names,
    "human training positives nearest to 5,000 (below / above)": lambda: list(_near5k()[0]),
}
for _tf in ("RXRA", "SRF", "ANDR"):
    DERIVED_TEXT[_F3D % _tf] = lambda tf=_tf: _f3d(tf)[1]


def _classify(d):
    traced, derived, notab, superseded = [], [], [], []
    for i, r in d.iterrows():
        note = str(r.statement_note)
        if note.startswith("pending"):
            continue
        if note.startswith("superseded"):
            superseded.append(r["item"])
        elif r["item"] in DERIVED:
            derived.append(i)
        elif r["item"] in NO_TABLE:
            notab.append(r["item"])
        else:
            traced.append(i)
    return traced, derived, notab, superseded


def test_pending_rows_are_only_the_h3b_items():
    d = numbers_file()
    note = d.statement_note.astype(str)
    pend = d[note.str.startswith("pending") | d["new value"].astype(str).str.startswith("pending")]
    bad = [(i, r["item"], str(r.statement_note)[:60]) for i, r in pend.iterrows() if not str(r.statement_note).startswith(PENDING_OK)]
    assert not bad, bad
    print("pending rows %d" % len(pend))


def test_untraced_rows_are_the_listed_ones():
    d = numbers_file()
    traced, derived, notab, superseded = _classify(d)
    assert sorted(set(notab)) == sorted(NO_TABLE)
    for i in traced:
        assert source_files(d.loc[i, "source file"]), "row %d (%s): no existing source table in '%s'" % (i, d.loc[i, "item"], d.loc[i, "source file"])
    assert len(superseded) <= 12, superseded
    print("traced %d, derived %d, no table %d, superseded %d" % (len(traced), len(derived), len(notab), len(superseded)))


def test_table_values_occur_in_their_source():
    d = numbers_file()
    traced, _, _, _ = _classify(d)
    bad = []
    for i in traced:
        r = d.loc[i]
        files = source_files(r["source file"])
        toks = tokens(r["new value"])
        miss = [t for t in toks if not any((np.abs(table_numbers(f) - float(t.replace(",", ""))) <= tol(t)).any()
                                           or (np.abs(100 * table_numbers(f) - float(t.replace(",", ""))) <= tol(t)).any() for f in files)]
        if miss:
            bad.append("row %d %s = %r: %s not in %s" % (i, r["item"], r["new value"], miss, [os.path.basename(f) for f in files]))
    assert not bad, "\n".join(bad)


def test_derived_values_match_recomputation():
    d = numbers_file()
    _, derived, _, _ = _classify(d)
    pend = set(d.item[d.statement_note.astype(str).str.startswith("pending")])   # the random-split rows can be pending
    missing = set(DERIVED) - set(d.item[derived]) - pend
    assert len(missing) <= 2, "DERIVED items no longer in the NUMBERS file: %s" % sorted(missing)
    bad = []
    for i in derived:
        r = d.loc[i]
        exp = DERIVED[r["item"]]()
        toks = tokens(r["new value"])[:len(exp)]
        if len(toks) != len(exp) or not all(close(t, x) for t, x in zip(toks, exp)):
            bad.append("row %d %s: file %r, recomputed %s" % (i, r["item"], r["new value"], exp))
        if r["item"] in DERIVED_TEXT:
            names = DERIVED_TEXT[r["item"]]()
            if not all(n in str(r["new value"]) for n in names):
                bad.append("row %d %s: %r lacks %s" % (i, r["item"], r["new value"], names))
    assert not bad, "\n".join(bad)


# ------------------------------------------------------------------ the *_PWM baseline is the train-selected best monoPWM
@pytest.mark.parametrize("rel", ["3_analysis/integration/inputs/results_table.csv", "3_analysis/package/results_table.csv"])
def test_pwm_baseline_is_best_mono(rel):
    t = pd.read_csv(J(rel), sep="\t")
    rf = t[(t.Model == "RandomForestClassifier") & (t.PWM == "mono+di")].set_index("TF_name").sort_index()
    mono = t[(t.Model == "Single best mono PWM") & (t.PWM == "mono")].drop_duplicates("TF_name").set_index("TF_name").loc[rf.index]
    if T.mutate("pwm_baseline_maxdi"):   # red run: the earlier evaluation convention, max(best mono, best di) for one TF
        rf = rf.copy()
        rf.iloc[0, rf.columns.get_loc("pr_auc_test_H_PWM")] += 0.01
    assert set(rf.index) == set(T.TF35)
    bad = [(tf, m) for m in ("roc_auc", "pr_auc") for tf in rf.index
           if rf.loc[tf, m + "_test_H_PWM"] != mono.loc[tf, m + "_test_H_PWM_mono"]]
    assert not bad, bad
