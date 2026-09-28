"""
run_mouse_transfer.py -- H->M (human-trained) and M->M (mouse-trained) RF vs correct single-mono-PWM baselines
for the 36 ArChIPelago TFs (mouse test set = chr1/8/19).  Writes ONLY into --rows-dir, --scores-dir and --out-csv.

usage: python run_mouse_transfer.py --tfs SRF P53 ... --workers 3 --rf-jobs 4 --out-csv <path> --rows-dir <dir> --scores-dir <dir>
Existing per-TF row json / npz are never overwritten: a TF whose row json exists is skipped.
"""
import os, sys, json, time, argparse, traceback
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import numpy as np
import pandas as pd
from joblib import parallel_backend

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
import pipeline_data as ad

PRC = "integral"          # auPRC estimator of the manuscript (PRROC)
SEEDS = (0, 1)


def col_metrics(X, y, cols):
    """(auROC, auPRC-integral) of every column j in cols."""
    return np.array([ad.auroc_auprc(y, X[:, j], prc=PRC) for j in cols])


def mono_baselines(d, tag):
    """Best mono PWM by TRAIN auROC and (independently) by TRAIN auPRC, evaluated on test_H/test_M;
    plus the oracle (best mono PWM selected on each test set itself).  Returns (row dict, chosen cols)."""
    cols = np.flatnonzero(d["is_mono"])
    tr = col_metrics(d["X_train"], d["y_train"], cols)
    te = {k: col_metrics(d["X_" + k], d["y_" + k], cols) for k in ("test_H", "test_M")}
    row, chosen = {}, {}
    for key, m in (("roc", 0), ("prc", 1)):
        i = int(np.argmax(tr[:, m]))
        j = int(cols[i])
        chosen[key] = j
        row["%s_b%s_col" % (tag, key)] = d["feature_names"][j]
        row["%s_b%s_train_auroc" % (tag, key)] = float(tr[i, 0])
        row["%s_b%s_train_auprc" % (tag, key)] = float(tr[i, 1])
        for k in ("test_H", "test_M"):
            row["%s_b%s_%s_auroc" % (tag, key, k)] = float(te[k][i, 0])
            row["%s_b%s_%s_auprc" % (tag, key, k)] = float(te[k][i, 1])
    for k in ("test_H", "test_M"):
        i_roc, i_prc = int(np.argmax(te[k][:, 0])), int(np.argmax(te[k][:, 1]))
        row["%s_oracle_%s_auroc" % (tag, k)] = float(te[k][i_roc, 0])
        row["%s_oracle_%s_auroc_col" % (tag, k)] = d["feature_names"][int(cols[i_roc])]
        row["%s_oracle_%s_auprc" % (tag, k)] = float(te[k][i_prc, 1])
        row["%s_oracle_%s_auprc_col" % (tag, k)] = d["feature_names"][int(cols[i_prc])]
    return row, chosen


def raw_col(d, key, j):
    """Un-standardised log-odds of column j of matrix `key` (per-matrix scaler inverted)."""
    sc = d["scalers"][key]
    return (d["X_" + key][:, j] * sc.scale_[j] + sc.mean_[j]).astype(np.float32)


def fit_eval(d, tag, n_jobs):
    """Fit RF with each seed; return (row dict, {seed: {matrix: scores}})."""
    row, scores = {}, {}
    for seed in SEEDS:
        t0 = time.time()
        # loky refuses to parallelise inside a daemonic multiprocessing worker (falls back to n_jobs=1);
        # the threading backend gives the intended n_jobs threads per worker (trees release the GIL)
        with parallel_backend("threading", n_jobs=n_jobs):
            rf = ad.fit_paper_rf(d["X_train"], d["y_train"], n_jobs=n_jobs, random_state=seed)
        row["%s_rf_fit_sec_s%d" % (tag, seed)] = round(time.time() - t0, 1)
        scores[seed] = {}
        for k in ("train", "test_H", "test_M"):
            with parallel_backend("threading", n_jobs=n_jobs):
                s = rf.predict_proba(d["X_" + k])[:, 1].astype(np.float32)
            scores[seed][k] = s
            a, p = ad.auroc_auprc(d["y_" + k], s, prc=PRC)
            row["%s_rf_%s_auroc_s%d" % (tag, k, seed)] = a
            row["%s_rf_%s_auprc_s%d" % (tag, k, seed)] = p
    for k in ("train", "test_H", "test_M"):
        for m in ("auroc", "auprc"):
            row["%s_rf_%s_%s_mean" % (tag, k, m)] = float(np.mean(
                [row["%s_rf_%s_%s_s%d" % (tag, k, m, s)] for s in SEEDS]))
    return row, scores


def all_columns():
    """Fixed csv column order (rows without M->M leave the MM_* fields empty)."""
    cols = ["TF", "P_mono_H", "P_di_H", "n_train_pos_H", "n_train_neg_H", "n_test_pos_H", "n_test_neg_H",
            "n_test_pos_M", "n_test_neg_M", "MM_available", "P_mono_M", "P_di_M", "n_train_pos_M", "n_train_neg_M"]
    for tag in ("HM", "MM"):
        for k in ("train", "test_H", "test_M"):
            for m in ("auroc", "auprc"):
                cols += ["%s_rf_%s_%s_mean" % (tag, k, m)] + ["%s_rf_%s_%s_s%d" % (tag, k, m, s) for s in SEEDS]
        for key in ("roc", "prc"):
            cols += ["%s_b%s_col" % (tag, key), "%s_b%s_train_auroc" % (tag, key), "%s_b%s_train_auprc" % (tag, key)]
            for k in ("test_H", "test_M"):
                cols += ["%s_b%s_%s_auroc" % (tag, key, k), "%s_b%s_%s_auprc" % (tag, key, k)]
        for k in ("test_H", "test_M"):
            cols += ["%s_oracle_%s_auroc" % (tag, k), "%s_oracle_%s_auroc_col" % (tag, k),
                     "%s_oracle_%s_auprc" % (tag, k), "%s_oracle_%s_auprc_col" % (tag, k)]
        cols += ["%s_load_sec" % tag] + ["%s_rf_fit_sec_s%d" % (tag, s) for s in SEEDS]
    cols += ["total_sec"]
    return cols


def mm_available(tf):
    base = os.path.join(ad.OUT, tf)
    return (os.path.isdir(os.path.join(base, "MOUSE_seq_MOUSE_pwm_mono"))
            and os.path.isfile(os.path.join(base, "train_neg_id_out_MOUSE_result_M_train_full.csv")))


def run_tf(args):
    tf, n_jobs, rows_dir, scores_dir = args
    t_start = time.time()
    row = {"TF": tf}
    npz = {}
    # ---------------- A. human-trained, human PWMs (H->H and H->M)
    t0 = time.time()
    d = ad.load_tf(tf, "HUMAN", "HUMAN", pwm_set="mono_di", scale="per_matrix")
    row["HM_load_sec"] = round(time.time() - t0, 1)
    row["P_mono_H"] = int(d["is_mono"].sum()); row["P_di_H"] = int((~d["is_mono"]).sum())
    row["n_train_pos_H"] = int(d["y_train"].sum()); row["n_train_neg_H"] = int((d["y_train"] == 0).sum())
    row["n_test_pos_H"] = int(d["y_test_H"].sum()); row["n_test_neg_H"] = int((d["y_test_H"] == 0).sum())
    row["n_test_pos_M"] = int(d["y_test_M"].sum()); row["n_test_neg_M"] = int((d["y_test_M"] == 0).sum())
    brow, chosen = mono_baselines(d, "HM")
    row.update(brow)
    rrow, sc = fit_eval(d, "HM", n_jobs)
    row.update(rrow)
    npz["y_test_H"] = d["y_test_H"]; npz["y_test_M"] = d["y_test_M"]
    npz["names_test_H"] = d["meta_test_H"]["name"].values.astype(str)
    npz["names_test_M"] = d["meta_test_M"]["name"].values.astype(str)
    for seed in SEEDS:
        npz["rf_H_s%d" % seed] = sc[seed]["test_H"]        # human model on human test
        npz["rf_M_s%d" % seed] = sc[seed]["test_M"]        # human model on mouse test
    for key in ("roc", "prc"):
        j = chosen[key]
        npz["pwmH_by%s_test_H" % key] = raw_col(d, "test_H", j)
        npz["pwmH_by%s_test_M" % key] = raw_col(d, "test_M", j)
    npz["pwmH_byroc_col"] = d["feature_names"][chosen["roc"]]
    npz["pwmH_byprc_col"] = d["feature_names"][chosen["prc"]]
    yM_H = d["y_test_M"].copy(); namesM_H = npz["names_test_M"].copy()
    del d, sc
    # ---------------- B. mouse-trained, mouse PWMs (M->M; test_H of this loader = human seqs x mouse PWMs = M->H)
    row["MM_available"] = mm_available(tf)
    if row["MM_available"]:
        t0 = time.time()
        dm = ad.load_tf(tf, "MOUSE", "MOUSE", pwm_set="mono_di", scale="per_matrix")
        row["MM_load_sec"] = round(time.time() - t0, 1)
        row["P_mono_M"] = int(dm["is_mono"].sum()); row["P_di_M"] = int((~dm["is_mono"]).sum())
        row["n_train_pos_M"] = int(dm["y_train"].sum()); row["n_train_neg_M"] = int((dm["y_train"] == 0).sum())
        if not (np.array_equal(dm["y_test_M"], yM_H) and np.array_equal(dm["meta_test_M"]["name"].values.astype(str), namesM_H)):
            raise RuntimeError("%s: mouse test rows differ between HUMAN-PWM and MOUSE-PWM loaders" % tf)
        brow, chosen = mono_baselines(dm, "MM")
        row.update(brow)
        rrow, sc = fit_eval(dm, "MM", n_jobs)
        row.update(rrow)
        for seed in SEEDS:
            npz["rfMM_M_s%d" % seed] = sc[seed]["test_M"]    # mouse model on mouse test
            npz["rfMM_H_s%d" % seed] = sc[seed]["test_H"]    # mouse model on human test (M->H)
        for key in ("roc", "prc"):
            j = chosen[key]
            npz["pwmM_by%s_test_M" % key] = raw_col(dm, "test_M", j)
            npz["pwmM_by%s_test_H" % key] = raw_col(dm, "test_H", j)
        npz["pwmM_byroc_col"] = dm["feature_names"][chosen["roc"]]
        npz["pwmM_byprc_col"] = dm["feature_names"][chosen["prc"]]
        del dm, sc
    row["total_sec"] = round(time.time() - t_start, 1)
    # ---------------- C. write (never overwrite)
    npz_path = os.path.join(scores_dir, tf + ".npz")
    json_path = os.path.join(rows_dir, tf + ".json")
    for p in (npz_path, json_path):
        if os.path.exists(p):
            raise RuntimeError("refusing to overwrite existing %s" % p)
    np.savez_compressed(npz_path, **npz)
    with open(json_path, "w") as fh:
        json.dump(row, fh, indent=1)
    return row


def worker(args):
    tf = args[0]
    try:
        row = run_tf(args)
        return tf, row, None
    except Exception:
        return tf, None, traceback.format_exc()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tfs", nargs="+", required=True)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--rf-jobs", type=int, default=4)
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--rows-dir", required=True)
    ap.add_argument("--scores-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.rows_dir, exist_ok=True); os.makedirs(a.scores_dir, exist_ok=True)
    todo = [tf for tf in a.tfs if not os.path.exists(os.path.join(a.rows_dir, tf + ".json"))]
    skipped = [tf for tf in a.tfs if tf not in todo]
    print("start %s  workers=%d rf_jobs=%d  todo=%d skipped(done)=%s" % (
        time.strftime("%F %T"), a.workers, a.rf_jobs, len(todo), skipped), flush=True)
    # big TFs first so the tail is short
    size = {tf: len(ad.feature_list(tf, "HUMAN", "mono_di")) for tf in todo}
    todo.sort(key=lambda t: -size[t])
    print("order:", [(t, size[t]) for t in todo], flush=True)
    jobs = [(tf, a.rf_jobs, a.rows_dir, a.scores_dir) for tf in todo]
    failures = {}
    COLS = all_columns()
    header_written = os.path.exists(a.out_csv)
    from multiprocessing import get_context
    with get_context("fork").Pool(a.workers, maxtasksperchild=1) as pool:
        for tf, row, err in pool.imap_unordered(worker, jobs):
            if err is not None:
                failures[tf] = err
                print("FAILED %s\n%s" % (tf, err), flush=True)
                continue
            extra = set(row) - set(COLS)
            if extra:
                raise RuntimeError("row keys not in fixed column list: %s" % sorted(extra))
            df = pd.DataFrame([row]).reindex(columns=COLS)
            df.to_csv(a.out_csv, mode="a", header=not header_written, index=False)
            header_written = True
            print("%s done %s  %.0fs  HM test_H %.4f/%.4f test_M %.4f/%.4f | bHroc test_M %.4f bHprc test_M %.4f | MM test_M %s" % (
                time.strftime("%T"), tf, row["total_sec"],
                row["HM_rf_test_H_auroc_mean"], row["HM_rf_test_H_auprc_mean"],
                row["HM_rf_test_M_auroc_mean"], row["HM_rf_test_M_auprc_mean"],
                row["HM_broc_test_M_auroc"], row["HM_bprc_test_M_auprc"],
                ("%.4f/%.4f" % (row["MM_rf_test_M_auroc_mean"], row["MM_rf_test_M_auprc_mean"])) if row["MM_available"] else "NA"),
                flush=True)
    print("finished %s  failures=%s" % (time.strftime("%F %T"), list(failures)), flush=True)
    if failures:
        sys.exit(1)


if __name__ == "__main__":
    main()
