"""train_models_v2d.py - the models of Models.tar.gz, fitted exactly as evaluate_models_h5.py fits them:
load_tf(SP, SP, mono_di, scale="train") with training = all autosomes outside the test chromosomes, columns sliced by
mode, RandomForestClassifier(max_depth 6, max_samples 0.8, n_estimators 100, random_state 0).
- HUMAN mono_di: the model saved by run_hm_mm_h5 (H5_alltrain/human_sym/models) is only checked; json written.
- HUMAN mono, di; MOUSE mono, di, mono_di: not saved by the evaluation; trained here, saved (n_jobs=1), reloaded and checked.
Every model is scored again and compared with the evaluation rows (rows_HM / rows_MM, RandomForestClassifier).
Writes <out>/Models/<TF>/ArChIPelago_<TF>[_MOUSE]_<mode>.{sav,json} (new files only).
Usage: train_models_v2d.py <out_dir> <threads_per_worker> <workers>"""
import os, sys
OUT, T, W = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[v] = str(T)
import json, time
from multiprocessing import Pool
import numpy as np, pandas as pd, joblib
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")); import archi_env as AE
sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "1_data_preparation"))
import archi_data_v5 as ad
assert ad.NEG_VERSION == "v3"
RB = AE.work()
RF_P = {"max_depth": 6, "max_samples": 0.8, "n_estimators": 100}
ROWS = {"HUMAN": RB + "/H5_alltrain/human_sym/eval_models/rows_HM", "MOUSE": RB + "/H5_alltrain/mouse_sym/eval_models/rows_MM"}
H5_MODEL = RB + "/H5_alltrain/human_sym/models/{tf}_HUMAN_RandomForestClassifier_mono_di_seed0_human_sym_negv3.sav"
CHR = {"HUMAN": ("chr2-7, 9-20, 22", "chr1, 8, 21"), "MOUSE": ("chr2-7, 9-18", "chr1, 8, 19")}
MODES = {"HUMAN": ["mono", "di", "mono_di"], "MOUSE": ["mono", "di", "mono_di"]}

def name(tf, sp, mode):
    return "ArChIPelago_%s_%s" % (tf, mode) if sp == "HUMAN" else "ArChIPelago_%s_MOUSE_%s" % (tf, mode)

def run(tf):
    t0 = time.time(); res = []
    for sp in ("HUMAN", "MOUSE"):
        d = ad.load_tf(tf, sp, sp, pwm_set="mono_di", scale="train")
        sc = d["scalers"]["train"]
        for mode in MODES[sp]:
            mask = {"mono": d["is_mono"], "di": ~d["is_mono"], "mono_di": np.ones(len(d["is_mono"]), bool)}[mode]
            cols = np.flatnonzero(mask)
            X = {k: np.ascontiguousarray(d["X_" + k][:, cols]) for k in ("train", "test_H", "test_M")}
            ref = json.load(open("%s/%s__RandomForestClassifier__%s.json" % (ROWS[sp], tf, mode)))
            assert ref["feature_names"] == [d["feature_names"][j] for j in cols] and ref["seed"] == 0 and ref["scale"] == "train", (tf, sp, mode)
            stem = os.path.join(OUT, "Models", tf, name(tf, sp, mode)); os.makedirs(os.path.dirname(stem), exist_ok=True)
            if sp == "HUMAN" and mode == "mono_di":
                m = joblib.load(H5_MODEL.format(tf=tf)); saved = "H5 run_hm_mm_h5 model (linked into the archive)"
            else:
                m = RandomForestClassifier(n_jobs=T, random_state=0, **RF_P).fit(X["train"], d["y_train"])
                m.n_jobs = 1
                if os.path.exists(stem + ".sav"): raise FileExistsError(stem)
                joblib.dump(m, stem + ".sav"); m = joblib.load(stem + ".sav"); saved = stem + ".sav"
            out = {"TF": tf, "species_train": sp, "mode": mode, "k": len(cols), "saved": saved}
            for k, rk in (("train", "train_H"), ("test_H", "test_H"), ("test_M", "test_M")):
                s = m.predict_proba(X[k])[:, 1]
                a = roc_auc_score(d["y_" + k], s)
                out["auroc_" + k] = a; out["auroc_" + k + "_h5"] = ref[rk][0]; out["absdiff_" + k] = abs(a - ref[rk][0])
                out["auprc_" + k] = ad.auprc_integral(d["y_" + k], s)
            own, other = ("test_H", "test_M") if sp == "HUMAN" else ("test_M", "test_H")
            sp_other = "MOUSE" if sp == "HUMAN" else "HUMAN"
            spec = {"tf": tf, "pwm_set": mode, "training_species": sp.lower(), "pwm_species": sp,
                    "model": "RandomForestClassifier(max_depth=6, max_samples=0.8, n_estimators=100, random_state=0), scikit-learn 1.3",
                    "features": ref["feature_names"],
                    "pwm_files": "mono_<k> = %s/<TF>/<k>.pwm, di_<k> = %s/<TF>/<k>.dpwm" % (
                        ("PWMs_mono_HUMAN", "PWMs_di_HUMAN") if sp == "HUMAN" else ("PWMs_MOUSE_mono", "PWMs_MOUSE_di")),
                    "scaler_mean": [float(x) for x in sc.mean_[cols]], "scaler_scale": [float(x) for x in sc.scale_[cols]],
                    "training_set": {"species": sp.lower(), "chromosomes": CHR[sp][0], "negatives": "v3",
                                     "n_positives": int(d["y_train"].sum()), "n_negatives": int((d["y_train"] == 0).sum())},
                    "test_set_same_species": {"species": sp.lower(), "chromosomes": CHR[sp][1],
                        "n_positives": int(d["y_" + own].sum()), "n_negatives": int((d["y_" + own] == 0).sum()),
                        "auROC": round(out["auroc_" + own], 4), "auPRC": round(out["auprc_" + own], 4)},
                    "test_set_other_species": {"species": sp_other.lower(), "chromosomes": CHR[sp_other][1],
                        "n_positives": int(d["y_" + other].sum()), "n_negatives": int((d["y_" + other] == 0).sum()),
                        "auROC": round(out["auroc_" + other], 4), "auPRC": round(out["auprc_" + other], 4)},
                    "scaling": "every matrix (training and test) is standardised with scaler_mean / scaler_scale, the "
                               "StandardScaler fitted on the training matrix; this is the scaling of the manuscript's evaluation",
                    "sarus": "--skipn --show-non-matching --output-scoring-mode score besthit", "sequence_length": 300}
            with open(stem + ".json", "x") as fh: json.dump(spec, fh, indent=1)
            res.append(out)
    print("%s done in %.0fs: %s" % (tf, time.time() - t0, " ".join("%s/%s %.1e" % (r["species_train"][0], r["mode"], max(r["absdiff_train"], r["absdiff_test_H"], r["absdiff_test_M"])) for r in res)), flush=True)
    return res

if __name__ == "__main__":
    tfs = sorted(os.listdir(RB + "/human_sym_h5/TF-ML/outputdir")); tfs = [t for t in tfs if os.path.isdir(RB + "/human_sym_h5/TF-ML/outputdir/" + t)]
    assert len(tfs) == 35 and "P53" not in tfs, tfs
    os.makedirs(OUT)
    with Pool(W) as p: rows = [r for rr in p.imap_unordered(run, tfs) for r in rr]
    R = pd.DataFrame(rows).sort_values(["species_train", "mode", "TF"])
    R.to_csv(os.path.join(OUT, "models_v2d_checks.csv"), index=False)
    for (sp, mode), g in R.groupby(["species_train", "mode"]):
        print(sp, mode, "n", len(g), "max |auROC - H5 eval| train %.2e test_H %.2e test_M %.2e" % (g.absdiff_train.max(), g.absdiff_test_H.max(), g.absdiff_test_M.max()))
    print("ALL <= 1e-6:", bool((R[["absdiff_train", "absdiff_test_H", "absdiff_test_M"]].max(axis=1) <= 1e-6).all()))
    print("MODELS_V2D_DONE")
