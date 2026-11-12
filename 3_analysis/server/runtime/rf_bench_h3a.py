"""RF runtime/memory benchmark for ArChIPelago.
Usage: rf_bench_h3a.py <TF> <step: load|train|predict> <n_jobs>
Data: $ARCHI_RELEASE/outputdir/<TF>/result_H_train_full.csv (the pipeline's human training
matrix; columns: index, 0-3 = id/name/threshold/TF, then k PWM best-hit scores, 'ind' = label).
Model: StandardScaler + RandomForestClassifier(max_depth=6, max_samples=0.8, n_estimators=100),
exactly as in the notebook (random_state=0 added for reproducibility).
"""
import sys, os, time, resource, json
import numpy as np, pandas as pd, joblib
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier

TF, step, n_jobs = sys.argv[1], sys.argv[2], int(sys.argv[3])
R = os.environ["BENCH_R"]  # output root set by the driver
csv = os.path.expanduser(os.environ["BENCH_CSV"].format(TF=TF))  # matrix set by the driver
out = {"TF": TF, "step": step, "n_jobs": n_jobs}
def rss_gb(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6  # kB -> GB
def rss_children_gb(): return resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1e6

t0 = time.perf_counter()
df = pd.read_csv(csv, index_col=0)
out["load_s"] = round(time.perf_counter() - t0, 2)
X, y = df.iloc[:, 4:-1], df["ind"].values
assert list(df.columns[-1:]) == ["ind"] and X.shape[1] > 0
out.update(n_rows=int(X.shape[0]), k_features=int(X.shape[1]), n_pos=int(y.sum()), n_neg=int((y == 0).sum()),
           csv_mb=round(os.path.getsize(csv) / 1e6, 1), rss_after_load_gb=round(rss_gb(), 3))
X = X.values.astype(np.float64)

if step == "train":
    t0 = time.perf_counter()
    scaler = StandardScaler(); Xs = scaler.fit_transform(X)
    out["scaler_s"] = round(time.perf_counter() - t0, 2)
    model = RandomForestClassifier(max_depth=6, max_samples=0.8, n_estimators=100, random_state=0, n_jobs=n_jobs)
    t0 = time.perf_counter(); model.fit(Xs, y); out["fit_s"] = round(time.perf_counter() - t0, 2)
    os.makedirs(f"{R}/models", exist_ok=True)
    mp = f"{R}/models/model_{TF}_HUMAN_RandomForestClassifier_mono_njobs{n_jobs}.sav"
    joblib.dump(model, mp); joblib.dump(scaler, mp.replace("model_", "scaler_"))
    out["model_bytes"] = os.path.getsize(mp); out["model_path"] = mp
    out["train_auc_sanity"] = round(float(__import__("sklearn.metrics", fromlist=["roc_auc_score"]).roc_auc_score(y, model.predict_proba(Xs)[:, 1])), 4)
elif step == "predict":
    mp = f"{R}/models/model_{TF}_HUMAN_RandomForestClassifier_mono_njobs1.sav"
    model = joblib.load(mp); scaler = joblib.load(mp.replace("model_", "scaler_"))
    model.n_jobs = n_jobs
    rng = np.random.default_rng(0)
    N = 1_000_000
    X1m = X[rng.integers(0, X.shape[0], N)]            # 1e6 x k raw best-hit scores (resampled rows)
    out["rss_after_1M_matrix_gb"] = round(rss_gb(), 3)
    t0 = time.perf_counter(); Xs = scaler.transform(X1m); out["scaler_transform_1M_s"] = round(time.perf_counter() - t0, 2)
    t0 = time.perf_counter(); p = model.predict_proba(Xs)[:, 1]; out["predict_proba_1M_s"] = round(time.perf_counter() - t0, 2)
    out["n_pred_rows"] = int(len(p)); out["mean_p"] = round(float(p.mean()), 4)
out["max_rss_gb_self"] = round(rss_gb(), 3); out["max_rss_gb_children"] = round(rss_children_gb(), 3)
print(json.dumps(out))
with open(f"{R}/rf_timings.jsonl", "a") as f: f.write(json.dumps(out) + "\n")
