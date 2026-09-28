#!/usr/bin/env python
"""PWM-subsampling saturation curves (Fig. S5): RF performance vs number of PWMs k.

Design A ("random"): for each k in {1,2,4,...,128} (k < P) draw R random subsets of k columns
(mono+di), fit the paper's RF (max_depth=6, max_samples=0.8, n_estimators=100, random_state=r);
k = P is fitted with 2 seeds.  Design B ("topk"): the k columns with the highest TRAIN auROC
(deterministic, one fit per k < P).  Metrics: auROC and auPRC (PRROC integral) on test_H and
test_M, the human-train-selected best single mono-PWM baseline, and the best single PWM inside
the subset.  Uses the loader ../common/pipeline_data.py.

Outputs (this folder):
  per_tf/<TF>.csv            rows of this TF (written when the TF finishes; a TF with an existing
                             file is skipped -> restartable, nothing is ever overwritten)
  per_tf/<TF>_subsets.json   column indices / names of every drawn subset and the top-k ranking
  per_tf/<TF>_columns.csv    per-column single-PWM auROC/auPRC on train/test_H/test_M
  saturation_results.csv     tidy long table, appended as TFs finish
"""
import os
import sys
import json
import time
import zlib
import argparse
import traceback
import multiprocessing as mp

import numpy as np
import pandas as pd
from joblib import parallel_backend

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
import pipeline_data as ad  # noqa: E402

OUTDIR = os.path.dirname(os.path.abspath(__file__))
PER_TF = os.path.join(OUTDIR, "per_tf")
MASTER = os.path.join(OUTDIR, "saturation_results.csv")
CSV36 = os.path.join(OUTDIR, "..", "inputs", "notebook4_results_table.csv")
K_GRID = [1, 2, 4, 8, 16, 32, 64, 128]
COLS = ["tf", "design", "k", "P", "rep", "seed", "auroc_H", "auprc_H", "auroc_M", "auprc_M",
        "base_auroc_H", "base_auprc_H", "base_auroc_M", "base_auprc_M",
        "sub_auroc_H", "sub_auprc_H", "sub_auroc_M", "sub_auprc_M",
        "fit_sec"]


def tf_list():
    csv = pd.read_csv(CSV36, sep="\t")
    rf = csv[(csv.Model == "RandomForestClassifier") & (csv.PWM == "mono+di")]
    tfs = sorted(rf.TF_name.unique())
    if len(tfs) != 36:
        raise RuntimeError("expected 36 TFs, got %d" % len(tfs))
    return tfs


def subset_seed(tf, k, r):
    return zlib.crc32(("%s|%d|%d" % (tf, k, r)).encode()) % (2 ** 31)


def column_table(d):
    """Single-PWM auROC / auPRC(integral) of every column on train, test_H, test_M."""
    rows = []
    for j, name in enumerate(d["feature_names"]):
        row = dict(col=j, feature=name, is_mono=bool(d["is_mono"][j]))
        for split in ("train", "test_H", "test_M"):
            a, p = ad.auroc_auprc(d["y_" + split], d["X_" + split][:, j], prc="integral")
            row["auroc_" + split] = a
            row["auprc_" + split] = p
        rows.append(row)
    return pd.DataFrame(rows)


def fit_eval(d, cols, random_state, n_jobs):
    t0 = time.time()
    out = {}
    # inside a daemonic Pool worker joblib's default loky backend silently falls back to n_jobs=1;
    # sklearn forests parallelise over trees with threads anyway, so force the threading backend
    with parallel_backend("threading", n_jobs=n_jobs):
        rf = ad.fit_paper_rf(d["X_train"][:, cols], d["y_train"], n_jobs=n_jobs, random_state=random_state)
        for split in ("test_H", "test_M"):
            s = rf.predict_proba(d["X_" + split][:, cols])[:, 1]
            a, p = ad.auroc_auprc(d["y_" + split], s, prc="integral")
            out["auroc_" + split[-1]] = a
            out["auprc_" + split[-1]] = p
    out["fit_sec"] = time.time() - t0
    return out


def subset_best(ct, cols):
    """Best single PWM inside the subset: auROC columns via best train auROC, auPRC columns via
    best train auPRC (same rule as the paper's baseline)."""
    sub = ct.iloc[list(cols)]
    jr = sub.auroc_train.idxmax()
    jp = sub.auprc_train.idxmax()
    return dict(sub_auroc_H=ct.at[jr, "auroc_test_H"], sub_auroc_M=ct.at[jr, "auroc_test_M"],
                sub_auprc_H=ct.at[jp, "auprc_test_H"], sub_auprc_M=ct.at[jp, "auprc_test_M"])


def run_tf(args):
    tf, R_small, R_large, n_jobs = args
    out_csv = os.path.join(PER_TF, tf + ".csv")
    if os.path.exists(out_csv):
        return tf, "skipped (exists)"
    t_start = time.time()
    try:
        d = ad.load_tf(tf, species_train="HUMAN", pwm_species="HUMAN", pwm_set="mono_di", scale="per_matrix")
        P = d["X_train"].shape[1]
        R = R_large if P > 100 else R_small
        t_load = time.time() - t_start

        ct = column_table(d)
        ct_path = os.path.join(PER_TF, tf + "_columns.csv")
        if not os.path.exists(ct_path):
            ct.to_csv(ct_path, index=False)

        # paper baseline: best single MONO PWM selected on human train (column choice as in the
        # release code), evaluated with the PRROC integral auPRC estimator
        tests = {"test_H": (d["X_test_H"], d["y_test_H"]), "test_M": (d["X_test_M"], d["y_test_M"])}
        b = ad.best_pwm_baseline(d["X_train"], d["y_train"], tests, mask=d["is_mono"])
        jr, jp = b["by_roc"]["col"], b["by_prc"]["col"]
        base = dict(base_auroc_H=ct.at[jr, "auroc_test_H"], base_auroc_M=ct.at[jr, "auroc_test_M"],
                    base_auprc_H=ct.at[jp, "auprc_test_H"], base_auprc_M=ct.at[jp, "auprc_test_M"])

        rows, subsets = [], {"tf": tf, "P": P, "feature_names": d["feature_names"],
                             "baseline_cols": {"by_roc": int(jr), "by_prc": int(jp)},
                             "random": {}, "topk": {}}
        ks = [k for k in K_GRID if k < P]

        def add(design, k, rep, seed, cols, rs):
            res = fit_eval(d, cols, rs, n_jobs)
            row = dict(tf=tf, design=design, k=k, P=P, rep=rep, seed=seed, **res, **base, **subset_best(ct, cols))
            rows.append(row)
            return row

        # design A: random subsets
        for k in ks:
            for r in range(R):
                seed = subset_seed(tf, k, r)
                cols = np.sort(np.random.RandomState(seed).choice(P, size=k, replace=False))
                subsets["random"]["%d_%d" % (k, r)] = [int(c) for c in cols]
                add("random", k, r, seed, cols, r)
        # k = P, two RF seeds
        allc = np.arange(P)
        for r in range(2):
            add("random", P, r, -1, allc, r)
        # design B: top-k by train auROC
        order = ct.sort_values("auroc_train", ascending=False, kind="mergesort").col.to_numpy()
        subsets["topk_ranking"] = [int(c) for c in order]
        for k in ks:
            cols = np.sort(order[:k])
            subsets["topk"]["%d" % k] = [int(c) for c in cols]
            add("topk", k, 0, -1, cols, 0)

        df = pd.DataFrame(rows)[COLS]
        df.to_csv(out_csv, index=False)
        subsets["load_sec"] = t_load
        subsets["total_sec"] = time.time() - t_start
        js_path = os.path.join(PER_TF, tf + "_subsets.json")
        if not os.path.exists(js_path):
            with open(js_path, "w") as fh:
                json.dump(subsets, fh)
        return tf, "done P=%d R=%d fits=%d load=%.0fs total=%.0fs" % (P, R, len(rows), t_load, time.time() - t_start)
    except Exception:
        return tf, "FAILED\n" + traceback.format_exc()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tfs", nargs="*", default=None)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--rf-jobs", type=int, default=4)
    ap.add_argument("--R", type=int, default=5)
    ap.add_argument("--R-large", type=int, default=5, help="reps for TFs with P > 100")
    a = ap.parse_args()
    os.makedirs(PER_TF, exist_ok=True)
    tfs = a.tfs or tf_list()
    # biggest first for load balancing
    Ps = {tf: len(ad.feature_list(tf, "HUMAN", "mono_di")) for tf in tfs}
    tfs = sorted(tfs, key=lambda t: -Ps[t])
    print("[%s] %d TFs, workers=%d rf_jobs=%d R=%d R_large=%d" % (time.ctime(), len(tfs), a.workers, a.rf_jobs, a.R, a.R_large),
          {t: Ps[t] for t in tfs}, flush=True)
    todo = [(tf, a.R, a.R_large, a.rf_jobs) for tf in tfs]
    ctx = mp.get_context("spawn")
    with ctx.Pool(a.workers) as pool:
        for tf, msg in pool.imap_unordered(run_tf, todo):
            print("[%s] %s: %s" % (time.ctime(), tf, msg), flush=True)
            if msg.startswith("done"):
                df = pd.read_csv(os.path.join(PER_TF, tf + ".csv"))
                df.to_csv(MASTER, mode="a", index=False, header=not os.path.exists(MASTER))
    print("[%s] all finished" % time.ctime(), flush=True)


if __name__ == "__main__":
    main()
