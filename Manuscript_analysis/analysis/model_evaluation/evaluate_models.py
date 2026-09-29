"""
evaluate_models.py - compute the rows of the evaluation table (schema of the results table) from the
per-TF feature files, with aligned human/mouse feature columns and the single-best-PWM baseline selected on
the human training set by PWM identity; mouse test set = chr1/8/19.

For every TF (dirs in $ARCHI_RELEASE_DIR/outputdir), pwm_set in {mono, di, mono_di} and the 5 paper models,
fits on human train (random_state=0) and records auROC / auPRC(PRROC integral) on train_H, test_H, test_M.
Baselines (best mono PWM, best di PWM; by train auROC and independently by train auPRC) once per TF.

Writes ONLY into --rows-dir (one json per TF/model/mode + one <TF>__baseline.json).  Existing row files
are never overwritten: combinations whose json exists are skipped, so a second instance can be started on
the remaining TFs (use --tfs ... or --reverse) and reruns resume.

usage: python evaluate_models.py --workers N --threads-per-worker T [--tfs SRF P53 ...] [--reverse]
       [--models ...] [--modes mono di mono_di] [--rows-dir rows]
Thread use = workers x threads-per-worker (RF/XGB n_jobs, Bagging n_jobs with single-threaded members,
OMP/BLAS caps).  Always launch with nice -n 10.
"""
import os, sys, json, time, argparse, traceback


def _threads_from_argv():
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument("--threads-per-worker", type=int, default=4)
    a, _ = ap.parse_known_args()
    return a.threads_per_worker


T_DEFAULT = _threads_from_argv()
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = str(T_DEFAULT)   # must precede numpy import

import numpy as np
from sklearn.ensemble import RandomForestClassifier, BaggingClassifier
from sklearn.linear_model import LogisticRegression
from xgboost import XGBClassifier

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
import pipeline_data as ad

PRC = "integral"     # auPRC estimator of the manuscript (PRROC pr.curve()$auc.integral)
SEED = 0
MODELS = ["RandomForestClassifier", "LogisticRegression", "XGBClassifier",
          "BaggingClassifier_XGBClassifier", "BaggingClassifier_LogisticRegression"]
MODES = ["mono", "di", "mono_di"]
# exact hyper-parameters of the release notebook (2_ArChIPelago_and_Slim_training.ipynb, model_building)
XGB_P = {'colsample_bytree': 0.8, 'gamma': 0.3, 'max_depth': 3, 'min_child_weight': 1,
         'n_estimators': 100, 'reg_alpha': 0.01, 'subsample': 0.8}
LR_P = {'C': 0.1, 'penalty': 'l2', 'solver': 'liblinear'}
RF_P = {'max_depth': 6, 'max_samples': 0.8, 'n_estimators': 100}
BAG_N = 10           # notebook: BaggingClassifier(...) without n_estimators -> sklearn default 10


def make_model(name, T):
    if name == "RandomForestClassifier":
        return RandomForestClassifier(n_jobs=T, random_state=SEED, **RF_P)
    if name == "XGBClassifier":
        return XGBClassifier(n_jobs=T, random_state=SEED, **XGB_P)
    if name == "LogisticRegression":
        return LogisticRegression(random_state=SEED, **LR_P)          # liblinear: single-threaded
    if name == "BaggingClassifier_XGBClassifier":
        return BaggingClassifier(estimator=XGBClassifier(n_jobs=1, random_state=SEED, **XGB_P),
                                 n_estimators=BAG_N, n_jobs=T, random_state=SEED)
    if name == "BaggingClassifier_LogisticRegression":
        return BaggingClassifier(estimator=LogisticRegression(random_state=SEED, **LR_P),
                                 n_estimators=BAG_N, n_jobs=T, random_state=SEED)
    raise ValueError(name)


def write_json(path, obj):
    """Create `path` only if it does not exist (never overwrite)."""
    if os.path.exists(path):
        return False
    tmp = "%s.tmp.%d" % (path, os.getpid())
    with open(tmp, "w") as f:
        json.dump(obj, f, indent=1, sort_keys=True)
    if os.path.exists(path):        # another instance finished the same combination meanwhile
        return False
    os.replace(tmp, path)           # renames the file created 1 line above; target did not exist
    return True


def col_metrics(X, y, cols):
    return np.array([ad.auroc_auprc(y, X[:, j], prc=PRC) for j in cols])


def baselines(d):
    """Best single PWM chosen on human TRAIN (by auROC, and independently by auPRC), separately among
    mono and among di PWMs, evaluated by NAME (= same column index; columns are aligned) on test_H/test_M."""
    out = {"tf": d["tf"], "n_mono": int(d["is_mono"].sum()), "n_di": int((~d["is_mono"]).sum()),
           "seq_count_train_pos": int(d["y_train"].sum()), "n_train": int(len(d["y_train"])),
           "n_test_H_pos": int(d["y_test_H"].sum()), "n_test_H": int(len(d["y_test_H"])),
           "n_test_M_pos": int(d["y_test_M"].sum()), "n_test_M": int(len(d["y_test_M"])),
           "feature_names": list(d["feature_names"]), "prc": PRC}
    for kind, mask in (("mono", d["is_mono"]), ("di", ~d["is_mono"])):
        cols = np.flatnonzero(mask)
        tr = col_metrics(d["X_train"], d["y_train"], cols)
        te = {k: col_metrics(d["X_" + k], d["y_" + k], cols) for k in ("test_H", "test_M")}
        for key, m in (("roc", 0), ("prc", 1)):
            i = int(np.argmax(tr[:, m]))       # first max wins, as in the notebook loop
            out["%s_by_%s" % (kind, key)] = {
                "name": d["feature_names"][int(cols[i])], "col": int(cols[i]),
                "train_H": [float(tr[i, 0]), float(tr[i, 1])],
                "test_H": [float(te["test_H"][i, 0]), float(te["test_H"][i, 1])],
                "test_M": [float(te["test_M"][i, 0]), float(te["test_M"][i, 1])]}
    return out


def run_tf(tf, models, modes, rows_dir, T):
    def log(msg):
        print("[%s] %s %s" % (time.strftime("%H:%M:%S"), tf, msg), flush=True)
    todo = [(mo, m) for mo in modes for m in models
            if not os.path.exists(os.path.join(rows_dir, "%s__%s__%s.json" % (tf, m, mo)))]
    bpath = os.path.join(rows_dir, "%s__baseline.json" % tf)
    if not todo and os.path.exists(bpath):
        log("all rows present, skipping")
        return
    t0 = time.time()
    d = ad.load_tf(tf, "HUMAN", "HUMAN", pwm_set="mono_di", scale="per_matrix")
    log("loaded %d features (%d mono, %d di), train %d (%d pos), test_H %d, test_M %d in %.0fs" % (
        len(d["feature_names"]), d["is_mono"].sum(), (~d["is_mono"]).sum(), len(d["y_train"]),
        d["y_train"].sum(), len(d["y_test_H"]), len(d["y_test_M"]), time.time() - t0))
    if not os.path.exists(bpath):
        t1 = time.time()
        b = baselines(d)
        write_json(bpath, b)
        log("baselines in %.0fs: mono by_roc %s test_H %.4f test_M %.4f; by_prc %s test_H %.4f test_M %.4f" % (
            time.time() - t1, b["mono_by_roc"]["name"], b["mono_by_roc"]["test_H"][0], b["mono_by_roc"]["test_M"][0],
            b["mono_by_prc"]["name"], b["mono_by_prc"]["test_H"][1], b["mono_by_prc"]["test_M"][1]))
    for mode in modes:
        mask = {"mono": d["is_mono"], "di": ~d["is_mono"], "mono_di": np.ones(len(d["is_mono"]), bool)}[mode]
        cols = np.flatnonzero(mask)
        # per-matrix StandardScaler is column-wise, so slicing the scaled mono_di matrices == loading pwm_set=mode
        X = {k: np.ascontiguousarray(d["X_" + k][:, cols]) for k in ("train", "test_H", "test_M")}
        for model_name in models:
            path = os.path.join(rows_dir, "%s__%s__%s.json" % (tf, model_name, mode))
            if os.path.exists(path):
                continue
            t1 = time.time()
            try:
                m = make_model(model_name, T)
                m.fit(X["train"], d["y_train"])
                fit_s = time.time() - t1
                res = {}
                for k in ("train", "test_H", "test_M"):
                    s = m.predict_proba(X[k])[:, 1]
                    res[k] = list(ad.auroc_auprc(d["y_" + k], s, prc=PRC))
                    res[k + "_auprc_trapz"] = ad.auroc_auprc(d["y_" + k], s, prc="trapz")[1]
                row = {"tf": tf, "model": model_name, "mode": mode, "n_features": int(len(cols)),
                       "feature_names": [d["feature_names"][j] for j in cols], "seed": SEED, "prc": PRC,
                       "scale": "per_matrix", "fit_seconds": fit_s, "total_seconds": time.time() - t1,
                       "threads": T, "train_H": res["train"], "test_H": res["test_H"], "test_M": res["test_M"],
                       "auprc_trapz": {k: res[k + "_auprc_trapz"] for k in ("train", "test_H", "test_M")}}
                write_json(path, row)
                log("%s %s: train %.4f/%.4f test_H %.4f/%.4f test_M %.4f/%.4f (fit %.0fs)" % (
                    mode, model_name, res["train"][0], res["train"][1], res["test_H"][0], res["test_H"][1],
                    res["test_M"][0], res["test_M"][1], fit_s))
            except Exception:
                log("FAILED %s %s\n%s" % (mode, model_name, traceback.format_exc()))
    log("done in %.0fs" % (time.time() - t0))


def _worker(args):
    tf, models, modes, rows_dir, T = args
    try:
        run_tf(tf, models, modes, rows_dir, T)
    except Exception:
        print("[%s] %s FAILED\n%s" % (time.strftime("%H:%M:%S"), tf, traceback.format_exc()), flush=True)
    return tf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tfs", nargs="*", default=None)
    ap.add_argument("--order", choices=["cost", "name"], default="cost",
                    help="cost: largest TFs (most PWM features) first, so the long jobs start early")
    ap.add_argument("--reverse", action="store_true", help="reverse the order (for a second instance: it then meets the first in the middle)")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--threads-per-worker", type=int, default=4)
    ap.add_argument("--models", nargs="*", default=MODELS)
    ap.add_argument("--modes", nargs="*", default=MODES)
    ap.add_argument("--rows-dir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "rows"))
    a = ap.parse_args()
    os.makedirs(a.rows_dir, exist_ok=True)
    tfs = a.tfs or sorted(x for x in os.listdir(ad.OUT) if os.path.isdir(os.path.join(ad.OUT, x)) and not x.startswith("."))
    if a.order == "cost":
        tfs = sorted(tfs, key=lambda tf: -len(ad.feature_list(tf, "HUMAN", "mono_di")))
    if a.reverse:
        tfs = tfs[::-1]
    for m in a.models:
        make_model(m, 1)          # validate names early
    print("evaluate_models: %d TFs, workers=%d threads/worker=%d models=%s modes=%s rows=%s" % (
        len(tfs), a.workers, a.threads_per_worker, a.models, a.modes, a.rows_dir), flush=True)
    jobs = [(tf, a.models, a.modes, a.rows_dir, a.threads_per_worker) for tf in tfs]
    if a.workers == 1:
        for j in jobs:
            _worker(j)
    else:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(a.workers, maxtasksperchild=1) as pool:
            for tf in pool.imap_unordered(_worker, jobs):
                pass
    print("evaluate_models: finished", flush=True)


if __name__ == "__main__":
    main()
