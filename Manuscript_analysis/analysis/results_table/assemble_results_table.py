"""
assemble_results_table.py -- build HUMAN_MOUSE_total_100k_recomputed.csv (identical schema / row semantics to
the notebook 4 table read by the R figure scripts) from the per-combination json rows written by
run_results_table.py, plus best_pwm_names.csv and a per-TF comparison with the notebook 4 table.  Then
validates against the notebook 4 table and prints the headline numbers.

usage: python assemble_results_table.py --rows-dir rows \
         --published ../inputs/HUMAN_MOUSE_total_100k_notebook4.csv --out-dir .
Never overwrites: refuses to run if an output file already exists.
"""
import os, sys, json, glob, argparse
import numpy as np
import pandas as pd

MODELS = ["RandomForestClassifier", "LogisticRegression", "XGBClassifier",
          "BaggingClassifier_XGBClassifier", "BaggingClassifier_LogisticRegression"]
MODE2PWM = {"mono": "mono", "di": "di", "mono_di": "mono+di"}
METRIC_BLOCKS = ["auc_test_M", "auc_test_H", "auc_train_H", "pr_auc_test_M", "pr_auc_test_H", "pr_auc_train_H"]


def load_rows(rows_dir):
    base, rows, col0 = {}, {}, {}
    for p in sorted(glob.glob(os.path.join(rows_dir, "*.json"))):
        with open(p) as f:
            j = json.load(f)
        if p.endswith("__baseline.json"):
            base[j["tf"]] = j
        elif p.endswith("__col0.json"):      # first mono / di PWM (published 'Single best di PWM' convention)
            col0[j["tf"]] = j
        else:
            rows[(j["tf"], j["model"], j["mode"])] = j
    return base, rows, col0


def base_vals(b, kind):
    """(roc: {train_H,test_H,test_M} of the by-auROC pick, pr: {...} of the by-auPRC pick) for kind mono|di."""
    r, p = b["%s_by_roc" % kind], b["%s_by_prc" % kind]
    return ({k: r[k][0] for k in ("train_H", "test_H", "test_M")},
            {k: p[k][1] for k in ("train_H", "test_H", "test_M")})


def fill_stats(row, roc, pr):
    """model/point values + the published mean/median/std columns (mean = median = value, std = 0)."""
    for blk, val in zip(METRIC_BLOCKS, [roc["test_M"], roc["test_H"], roc["train_H"], pr["test_M"], pr["test_H"], pr["train_H"]]):
        row[("roc_" + blk) if not blk.startswith("pr_") else blk] = val
        row["mean_" + blk] = val
        row["median_" + blk] = val
        row["std_" + blk] = 0.0


def make_rows(tf, b, rows, tf_index):
    """25 rows for one TF in the published block order:
    for each model: [Single best mono PWM, model mono, Single best di PWM, model di, model mono+di]."""
    bm_roc, bm_pr = base_vals(b, "mono")
    bd_roc, bd_pr = base_vals(b, "di")
    n_mono, n_di = b["n_mono"], b["n_di"]
    common = {"TF_name": tf, "Unnamed: 0": tf_index, "Names": tf, "Seq_count": b["seq_count_train_pos"]}
    out = []

    def baseline_row(kind):
        roc, pr = (bm_roc, bm_pr) if kind == "mono" else (bd_roc, bd_pr)
        row = dict(common, Count=n_mono if kind == "mono" else n_di,
                   Model="Single best %s PWM" % kind, PWM=kind)
        for k in ("train_H", "test_H", "test_M"):
            row["roc_auc_%s_PWM_%s" % (k, kind)] = roc[k]
            row["pr_auc_%s_PWM_%s" % (k, kind)] = pr[k]
            other = "di" if kind == "mono" else "mono"
            row["roc_auc_%s_PWM_%s" % (k, other)] = 0.0
            row["pr_auc_%s_PWM_%s" % (k, other)] = 0.0
        fill_stats(row, roc, pr)
        row["roc_auc_test_H_PWM"], row["roc_auc_test_M_PWM"] = roc["test_H"], roc["test_M"]
        row["pr_auc_test_H_PWM"], row["pr_auc_test_M_PWM"] = pr["test_H"], pr["test_M"]
        return row

    def model_row(model, mode):
        j = rows[(tf, model, mode)]
        pwm = MODE2PWM[mode]
        # published quirk: Count of the mono+di rows = number of mono PWMs (features_c counted mono files)
        row = dict(common, Count={"mono": n_mono, "di": n_di, "mono_di": n_mono}[mode], Model=model, PWM=pwm)
        if mode == "mono":
            roc, pr = bm_roc, bm_pr
            for k in ("train_H", "test_H", "test_M"):
                for kind in ("mono", "di"):          # published quirk: *_PWM_di columns repeat the mono values
                    row["roc_auc_%s_PWM_%s" % (k, kind)] = roc[k]
                    row["pr_auc_%s_PWM_%s" % (k, kind)] = pr[k]
        elif mode == "di":
            roc, pr = bd_roc, bd_pr
            for k in ("train_H", "test_H", "test_M"):
                for kind in ("mono", "di"):
                    row["roc_auc_%s_PWM_%s" % (k, kind)] = roc[k]
                    row["pr_auc_%s_PWM_%s" % (k, kind)] = pr[k]
        else:   # mono+di: test baselines = max(best mono, best di) per metric and test set; train = 0
            roc = {"train_H": 0.0, "test_H": max(bm_roc["test_H"], bd_roc["test_H"]), "test_M": max(bm_roc["test_M"], bd_roc["test_M"])}
            pr = {"train_H": 0.0, "test_H": max(bm_pr["test_H"], bd_pr["test_H"]), "test_M": max(bm_pr["test_M"], bd_pr["test_M"])}
            for k in ("train_H", "test_H", "test_M"):
                for kind in ("mono", "di"):
                    row["roc_auc_%s_PWM_%s" % (k, kind)] = roc[k]
                    row["pr_auc_%s_PWM_%s" % (k, kind)] = pr[k]
        m_roc = {k: j[k][0] for k in ("train_H", "test_H", "test_M")}
        m_pr = {k: j[k][1] for k in ("train_H", "test_H", "test_M")}
        fill_stats(row, m_roc, m_pr)
        row["roc_auc_test_H_PWM"], row["roc_auc_test_M_PWM"] = roc["test_H"], roc["test_M"]
        row["pr_auc_test_H_PWM"], row["pr_auc_test_M_PWM"] = pr["test_H"], pr["test_M"]
        return row

    for model in MODELS:
        out += [baseline_row("mono"), model_row(model, "mono"), baseline_row("di"),
                model_row(model, "di"), model_row(model, "mono_di")]
    return out


def headline(df, label):
    """Paper-style numbers from a table: RF mono+di vs single best mono PWM (test_H, test_M)."""
    rf = df[(df.Model == "RandomForestClassifier") & (df.PWM == "mono+di")].set_index("TF_name").sort_index()
    bm = df[df.Model == "Single best mono PWM"].groupby("TF_name").first().sort_index()   # 5 identical rows per TF (published: not identical for test_M)
    bm_all = df[df.Model == "Single best mono PWM"]
    print("\n== %s: RF (mono+di) vs single best mono PWM ('Single best mono PWM' rows; median of 36 TFs)" % label)
    res = {}
    for sp in ("H", "M"):
        for met, col in (("auROC", "roc_auc_test_%s" % sp), ("auPRC", "pr_auc_test_%s" % sp)):
            v_rf, v_b = rf[col].median(), bm_all[col].median()
            impr = int((rf[col].values > bm.loc[rf.index, col].values).sum())
            print("  test_%s %s: RF %.4f vs PWM %.4f (diff of medians %+.4f); TFs improved %d/36" % (sp, met, v_rf, v_b, v_rf - v_b, impr))
            res[(sp, met)] = (v_rf, v_b, impr)
    print("  (baseline taken from the *_PWM_mono columns of the RF mono+di rows [= max(best mono, best di) for test]):")
    for sp in ("H", "M"):
        for met, col, bcol in (("auROC", "roc_auc_test_%s" % sp, "roc_auc_test_%s_PWM_mono" % sp),
                               ("auPRC", "pr_auc_test_%s" % sp, "pr_auc_test_%s_PWM_mono" % sp)):
            print("    test_%s %s: RF %.4f vs PWM col %.4f; TFs improved %d/36" % (sp, met, rf[col].median(), rf[bcol].median(), int((rf[col] > rf[bcol]).sum())))
    worse = sorted(rf.index[(rf["roc_auc_test_M"] < bm.loc[rf.index, "roc_auc_test_M"].values) |
                            (rf["pr_auc_test_M"] < bm.loc[rf.index, "pr_auc_test_M"].values)])
    worse_col = sorted(rf.index[(rf["roc_auc_test_M"] < rf["roc_auc_test_M_PWM_mono"]) | (rf["pr_auc_test_M"] < rf["pr_auc_test_M_PWM_mono"])])
    print("  H->M RF worse than best mono PWM on >=1 metric (vs 'Single best mono PWM' rows): %d: %s" % (len(worse), " ".join(worse)))
    print("  H->M RF worse than *_PWM_mono column on >=1 metric: %d: %s" % (len(worse_col), " ".join(worse_col)))
    return rf, bm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows-dir", required=True)
    ap.add_argument("--published", required=True)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    out_csv = os.path.join(a.out_dir, "HUMAN_MOUSE_total_100k_recomputed.csv")
    out_names = os.path.join(a.out_dir, "best_pwm_names.csv")
    out_cmp = os.path.join(a.out_dir, "per_tf_notebook4_vs_recomputed.csv")
    for p in (out_csv, out_names, out_cmp):
        if os.path.exists(p):
            sys.exit("refusing to overwrite existing %s" % p)

    pub = pd.read_csv(a.published, sep="\t")
    cols = list(pub.columns)
    tf_index = pub.groupby("TF_name")["Unnamed: 0"].first().to_dict()
    tfs = sorted(pub["TF_name"].unique())
    base, rows, col0 = load_rows(a.rows_dir)
    missing = [(tf, m, mo) for tf in tfs for m in MODELS for mo in MODE2PWM if (tf, m, mo) not in rows]
    missing_b = [tf for tf in tfs if tf not in base]
    if missing or missing_b:
        sys.exit("missing rows: %d %s ... ; missing baselines: %s" % (len(missing), missing[:5], missing_b))

    recs = []
    for tf in tfs:
        recs += make_rows(tf, base[tf], rows, int(tf_index[tf]))
    df = pd.DataFrame(recs)[cols]
    assert list(df.columns) == cols, "column mismatch"
    assert len(df) == len(pub) == 900, (len(df), len(pub))
    key = lambda d: sorted(zip(d.TF_name, d.Model, d.PWM))
    assert key(df) == key(pub), "key set / multiplicity differs from the published table"
    assert (df.groupby("TF_name")["Seq_count"].first().sort_index().values ==
            pub.groupby("TF_name")["Seq_count"].first().sort_index().values).all(), "Seq_count differs"
    assert list(zip(df.TF_name, df.Model, df.PWM)) == list(zip(pub.TF_name, pub.Model, pub.PWM)), "row order differs"
    cnt_f = df.groupby(["TF_name", "PWM"])["Count"].first(); cnt_p = pub.groupby(["TF_name", "PWM"])["Count"].first()
    print("Count columns identical to published: %s" % (cnt_f == cnt_p.loc[cnt_f.index]).all())
    df.to_csv(out_csv, sep="\t", index=False)
    print("wrote %s (%d rows x %d cols)" % (out_csv, *df.shape))

    # side table with the selected PWM names
    nrec = []
    for tf in tfs:
        b = base[tf]
        r = {"TF_name": tf, "Seq_count": b["seq_count_train_pos"], "n_mono": b["n_mono"], "n_di": b["n_di"],
             "n_test_H_pos": b["n_test_H_pos"], "n_test_M_pos": b["n_test_M_pos"]}
        for kind in ("mono", "di"):
            for key_, mi in (("roc", 0), ("prc", 1)):
                x = b["%s_by_%s" % (kind, key_)]
                r["best_%s_by_%s" % (kind, key_)] = x["name"]
                for k in ("train_H", "test_H", "test_M"):
                    r["best_%s_by_%s_%s_%s" % (kind, key_, k, "auroc" if mi == 0 else "auprc")] = x[k][mi]
        if tf in col0:      # published convention for the di baseline: column 0 of the di matrix (first di PWM)
            for kind in ("mono", "di"):
                x = col0[tf][kind + "_0"]
                r["first_%s" % kind] = x["name"]
                for k, kk in (("train", "train_H"), ("test_H", "test_H"), ("test_M", "test_M")):
                    r["first_%s_%s_auroc" % (kind, kk)] = x[k][0]
                    r["first_%s_%s_auprc" % (kind, kk)] = x[k][1]
        nrec.append(r)
    pd.DataFrame(nrec).to_csv(out_names, sep="\t", index=False)
    print("wrote %s" % out_names)

    # headline numbers
    rf_p, bm_p = headline(pub, "PUBLISHED")
    rf_f, bm_f = headline(df, "FIXED")

    # per-TF comparison (RF mono+di minus single best mono PWM), published vs fixed
    cmp = pd.DataFrame({"TF_name": tfs})
    cmp = cmp.set_index("TF_name")
    for sp in ("H", "M"):
        for met, col in (("auROC", "roc_auc_test_%s" % sp), ("auPRC", "pr_auc_test_%s" % sp)):
            cmp["RF_%s_%s_pub" % (sp, met)] = rf_p.loc[tfs, col].values
            cmp["RF_%s_%s_fix" % (sp, met)] = rf_f.loc[tfs, col].values
            # published mono baseline for the mouse differs between the 5 model blocks (scrambled); use the RF block (first row)
            cmp["PWM_%s_%s_pub" % (sp, met)] = bm_p.loc[tfs, col].values
            cmp["PWM_%s_%s_fix" % (sp, met)] = bm_f.loc[tfs, col].values
            cmp["delta_%s_%s_pub" % (sp, met)] = cmp["RF_%s_%s_pub" % (sp, met)] - cmp["PWM_%s_%s_pub" % (sp, met)]
            cmp["delta_%s_%s_fix" % (sp, met)] = cmp["RF_%s_%s_fix" % (sp, met)] - cmp["PWM_%s_%s_fix" % (sp, met)]
    cmp.reset_index().to_csv(out_cmp, sep="\t", index=False, float_format="%.4f")
    print("wrote %s" % out_cmp)
    pd.set_option("display.width", 250)
    print("\nper-TF mouse deltas (RF mono+di minus best mono PWM), published vs fixed:")
    print(cmp[["delta_M_auROC_pub", "delta_M_auROC_fix", "delta_M_auPRC_pub", "delta_M_auPRC_fix"]].round(4).to_string())
    print("\nhuman: |RF_fix - RF_pub| median/max auROC %.4f/%.4f, auPRC %.4f/%.4f" % (
        (cmp.RF_H_auROC_fix - cmp.RF_H_auROC_pub).abs().median(), (cmp.RF_H_auROC_fix - cmp.RF_H_auROC_pub).abs().max(),
        (cmp.RF_H_auPRC_fix - cmp.RF_H_auPRC_pub).abs().median(), (cmp.RF_H_auPRC_fix - cmp.RF_H_auPRC_pub).abs().max()))
    print("human baseline identical to published (by TF): auROC max|diff| %.2e, auPRC max|diff| %.2e" % (
        (cmp.PWM_H_auROC_fix - cmp.PWM_H_auROC_pub).abs().max(), (cmp.PWM_H_auPRC_fix - cmp.PWM_H_auPRC_pub).abs().max()))


if __name__ == "__main__":
    main()
