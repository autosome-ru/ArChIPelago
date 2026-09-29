"""
fit_RF2f_models.py - the models of Fig. 4 / S4 that contain the two single-PWM features (RF2f family).

Loader = ../common/pipeline_data.py; mouse test set = chr1/8/19 (logical ("MOUSE", "control")).
  * Slim test_M scores = scans of the mouse test-set fasta with the Slim models of the TF
    ($ARCHI_SLIM_SCANS_DIR/slim_scans/<TF>/m<m>_<MOUSE_SET>/model1_predictions.txt, written by scan_slim.sh;
    checked by row count and by sequence-name order against the tab file); ANDR m=1 from the model in
    slim_ANDR_m1/ (its human scores from slim_scans/ANDR/m1_HUMAN_*);
  * diChIPMunk test_M scores = SARUS scans $ARCHI_SLIM_SCANS_DIR/munk/<TF>_full_train_M_1_ChIPMunk_no_repeats_0.tab
    (scan_dichipmunk.sh);
  * every TF must have Slim m=0, m=1, m=-5 and diChIPMunk on all three sets, otherwise the script fails.

The two single-PWM features of the RF2f models are the best monoPWM and the best diPWM of the TF, each
selected by training auROC among its class. The Slim and diChIPMunk single-model scores are read as saved.

What is fitted per TF (paper RF: max_depth 6, max_samples 0.8, 100 trees, random_state 0,
each matrix standardised on itself as in the release code):
  RF2f            best monoPWM + best diPWM (both selected by TRAIN auROC among their class)
  RF2f_munk       + diChIPMunk (de novo diPWM of the training positives, notebook 2 scan files)
  RF2f_slim1      + Slim m=1
  RF2f_lslim5     + LSlim m=-5
  RF2f_all5       + Slim m=1 + LSlim m=-5 + diChIPMunk   (the five-feature model of the text)
Also written: single-model metrics of Slim m=0, m=1, m=-5, m=-7 and diChIPMunk on test_H/test_M,
and the metrics of the two selected PWMs.  Slim scores are read from the saved Slim predictions on the global sequence sets
(<TF>_SlimModel_<m>/*_predictions/model1_predictions.txt, column 4 = maxscore, one row per line of
out_tab_<SP>_10000_<split>.tab; the newest complete file per set is used), diChIPMunk scores from
HUMAN_seq_HUMAN_pwm_mono/full_{train_H,control_H}_1_ChIPMunk_no_repeats_0.tab and
MOUSE_seq_HUMAN_pwm_mono/full_control_M_1_ChIPMunk_no_repeats_0.tab (global order).

usage: python fit_RF2f_models.py --tfs A B ... --rows-dir rows     (never overwrites)
"""
import os, sys, json, time, glob, argparse
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
import pipeline_data as ad

PRC = "integral"
N_GLOBAL = {("HUMAN", "train"): 241362, ("HUMAN", "control"): 91074, ("MOUSE", "control"): 170571}
SWAP = os.environ.get("ARCHI_SLIM_SCANS_DIR", os.path.dirname(os.path.abspath(__file__)))
MOUSE_SET = "MOUSE_10000_" + ad._file_split("MOUSE", "control")   # on-disk name of the chr1/8/19 mouse table / fasta
SLIM = {"slim0": 0, "slim1": 1, "lslim5": 5}


def global_rows(tf):
    out = {}
    for key, (sp, split) in {"train": ("HUMAN", "train"), "test_H": ("HUMAN", "control"), "test_M": ("MOUSE", "control")}.items():
        gt = ad.global_table(sp, split)
        if gt.shape[0] != N_GLOBAL[(sp, split)]:
            raise ValueError("global table %s %s has %d rows" % (sp, split, gt.shape[0]))
        rows, y = ad.select_rows(gt, tf, sp, split)
        out[key] = (rows, y)
    return out


def _mouse_names():
    """sequence names of out_tab_MOUSE_10000_train.tab (column 2 before the '$', as global_table strips it), global order."""
    gt = ad.global_table("MOUSE", "control")
    if gt.shape[0] != N_GLOBAL[("MOUSE", "control")]:
        raise ValueError("mouse test table has %d rows" % gt.shape[0])
    return gt["name"].astype(str).to_numpy()


def slim_scores(tf, m, mouse_names):
    """global-order maxscore vectors {train, test_H, test_M} of Slim model order m."""
    d = os.path.join(ad.OUT, tf, "%s_SlimModel_%d" % (tf, m))
    files = sorted(glob.glob(os.path.join(d, "*_predictions", "model1_predictions.txt")), key=os.path.getmtime)
    # fallback: predictions regenerated from the saved Slim model (ANDR m=0 has none in the pipeline output)
    files += sorted(glob.glob(os.path.join(SWAP, "slim_predictions_%s" % tf, "m%d_*" % m, "model1_predictions.txt")),
                    key=os.path.getmtime)
    # human scans of the ANDR m=1 model of slim_ANDR_m1/
    files += sorted(glob.glob(os.path.join(SWAP, "slim_scans", tf, "m%d_HUMAN_*" % m, "model1_predictions.txt")), key=os.path.getmtime)
    by_len = {}
    for f in files:                                # newest wins
        with open(f) as fh:
            first = fh.readline().rstrip("\n").split("\t")[-1]
        n = sum(1 for _ in open(f))
        by_len[n] = (f, first)
    out = {}
    for key, (sp, split) in {"train": ("HUMAN", "train"), "test_H": ("HUMAN", "control")}.items():
        n = N_GLOBAL[(sp, split)]
        if n not in by_len:
            raise FileNotFoundError("%s Slim m=%d: no prediction file with %d rows" % (tf, m, n))
        f, first = by_len[n]
        tag = "_%s_" % sp
        if tag not in first or (split == "train") != ("_macs_train_" in first):
            raise ValueError("%s Slim m=%d: %s does not look like %s %s (%s)" % (tf, m, f, sp, split, first))
        v = pd.read_csv(f, sep="\t", header=None, usecols=[3])[3].to_numpy(dtype=np.float64)
        out[key] = (v, f)
    # mouse test set: the scan of the chr1/8/19 fasta, checked by row count and sequence-name order
    f = os.path.join(SWAP, "slim_scans", tf, "m%d_%s" % (m, MOUSE_SET), "model1_predictions.txt")
    if not os.path.isfile(f):
        raise FileNotFoundError("%s Slim m=%d: %s missing" % (tf, m, f))
    t = pd.read_csv(f, sep="\t", header=None, usecols=[3, 6])
    if t.shape[0] != N_GLOBAL[("MOUSE", "control")]:
        raise ValueError("%s Slim m=%d: %s has %d rows, expected %d" % (tf, m, f, t.shape[0], N_GLOBAL[("MOUSE", "control")]))
    if not np.array_equal(t[6].astype(str).to_numpy(), mouse_names):
        raise ValueError("%s Slim m=%d: sequence names of %s are not in the order of out_tab_%s.tab" % (tf, m, f, MOUSE_SET))
    out["test_M"] = (t[3].to_numpy(dtype=np.float64), f)
    return out


def munk_scores(tf):
    base = os.path.join(ad.OUT, tf)
    paths = {"train": os.path.join(base, "HUMAN_seq_HUMAN_pwm_mono", "full_train_H_1_ChIPMunk_no_repeats_0.tab"),
             "test_H": os.path.join(base, "HUMAN_seq_HUMAN_pwm_mono", "full_control_H_1_ChIPMunk_no_repeats_0.tab"),
             "test_M": os.path.join(SWAP, "munk", "%s_full_train_M_1_ChIPMunk_no_repeats_0.tab" % tf)}
    out = {}
    for k, p in paths.items():
        v = pd.read_csv(p, sep="\t", header=None, usecols=[0])[0].to_numpy(dtype=np.float64)
        sp, split = {"train": ("HUMAN", "train"), "test_H": ("HUMAN", "control"), "test_M": ("MOUSE", "control")}[k]
        if len(v) != N_GLOBAL[(sp, split)]:
            raise ValueError("%s diChIPMunk %s: %d rows, expected %d" % (tf, k, len(v), N_GLOBAL[(sp, split)]))
        out[k] = v
    return out


def metrics(y, s):
    return ad.auroc_auprc(y, s, prc=PRC)


def fit_eval(Xs, ys, seed=0):
    """Xs: dict split -> raw feature matrix; per-matrix standardisation as in the release code."""
    Z = {k: StandardScaler().fit_transform(X) for k, X in Xs.items()}
    m = ad.fit_paper_rf(Z["train"], ys["train"], n_jobs=1, random_state=seed)
    out = {}
    for k in ("train", "test_H", "test_M"):
        out[k] = metrics(ys[k], m.predict_proba(Z[k])[:, 1])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tfs", nargs="+", required=True)
    ap.add_argument("--rows-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.rows_dir, exist_ok=True)
    for tf in a.tfs:
        out_path = os.path.join(a.rows_dir, tf + ".json")
        if os.path.exists(out_path):
            print("skip", tf, flush=True); continue
        t0 = time.time()
        d = ad.load_tf(tf, "HUMAN", "HUMAN", pwm_set="mono_di", scale=None)     # raw log-odds
        rows = global_rows(tf)
        mouse_names = _mouse_names()
        for k in ("train", "test_H", "test_M"):
            if not np.array_equal(rows[k][1], d["y_" + k]):
                raise ValueError("%s: row selection differs from load_tf for %s" % (tf, k))
        ys = {k: d["y_" + k] for k in ("train", "test_H", "test_M")}
        # best mono / di by TRAIN auROC (and, for the record, by train auPRC)
        tr = np.array([metrics(ys["train"], d["X_train"][:, j]) for j in range(d["X_train"].shape[1])])
        mono = np.flatnonzero(d["is_mono"]); di = np.flatnonzero(~d["is_mono"])
        j_mono = int(mono[np.argmax(tr[mono, 0])]); j_di = int(di[np.argmax(tr[di, 0])])
        row = {"TF": tf, "n_mono": int(len(mono)), "n_di": int(len(di)),
               "best_mono_col": d["feature_names"][j_mono], "best_di_col": d["feature_names"][j_di],
               "best_mono_train_auroc": float(tr[j_mono, 0]), "best_di_train_auroc": float(tr[j_di, 0])}
        for k in ("test_H", "test_M"):
            row["best_mono_%s" % k] = metrics(ys[k], d["X_" + k][:, j_mono])
            row["best_di_%s" % k] = metrics(ys[k], d["X_" + k][:, j_di])
        # Slim and diChIPMunk feature vectors on the selected rows
        feats = {}
        for name, m in SLIM.items():
            sc = slim_scores(tf, m, mouse_names)          # raises if any of the three sets is missing
            feats[name] = {k: sc[k][0][rows[k][0]] for k in rows}
            row["slim_file_%s" % name] = {k: os.path.relpath(sc[k][1], ad.OUT) for k in rows}
        mk = munk_scores(tf)
        feats["munk"] = {k: mk[k][rows[k][0]] for k in rows}
        for name in list(feats):
            for k in ("test_H", "test_M"):
                row["single_%s_%s" % (name, k)] = metrics(ys[k], feats[name][k])
        # RF models
        pwm2 = {k: np.column_stack([d["X_" + k][:, j_mono], d["X_" + k][:, j_di]]) for k in rows}
        models = {"RF2f": [], "RF2f_munk": ["munk"], "RF2f_slim1": ["slim1"], "RF2f_lslim5": ["lslim5"],
                  "RF2f_all5": ["slim1", "lslim5", "munk"]}
        for mname, extra in models.items():
            if any(e not in feats for e in extra):
                raise RuntimeError("%s: feature(s) %s missing for %s" % (tf, extra, mname))
            Xs = {k: np.column_stack([pwm2[k]] + [feats[e][k][:, None] for e in extra]) for k in rows}
            row[mname] = fit_eval(Xs, ys)
            row[mname]["n_features"] = 2 + len(extra)
        row["sec"] = round(time.time() - t0, 1)
        with open(out_path, "w") as fh:
            json.dump(row, fh, indent=1)
        print("%s done %s %.0fs  mono %s di %s | RF2f test_H %.4f (best mono %.4f, best di %.4f) | slim1 %.4f munk %.4f" % (
            time.strftime("%T"), tf, row["sec"], row["best_mono_col"], row["best_di_col"],
            row["RF2f"]["test_H"][0], row["best_mono_test_H"][0], row["best_di_test_H"][0],
            row.get("single_slim1_test_H", [float("nan")])[0], row["single_munk_test_H"][0]), flush=True)


if __name__ == "__main__":
    main()
