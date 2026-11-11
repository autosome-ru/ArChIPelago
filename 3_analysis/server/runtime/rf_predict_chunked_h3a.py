"""Chunked inference benchmark: 1,000,000 rows scored in chunks of 100,000 (scaler.transform +
RandomForest.predict_proba), as a deployment/scanning tool would do. Usage: rf_predict_chunked_h3a.py <TF> <n_jobs>"""
import sys, os, time, resource, json
import numpy as np, pandas as pd, joblib
TF, n_jobs = sys.argv[1], int(sys.argv[2])
R = os.environ["BENCH_R"]  # output root set by the driver
csv = os.path.expanduser(os.environ["BENCH_CSV"].format(TF=TF))  # matrix set by the driver
rss = lambda: round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6, 3)
df = pd.read_csv(csv, index_col=0); X = df.iloc[:, 4:-1].values.astype(np.float64); del df
mp = f"{R}/models/model_{TF}_HUMAN_RandomForestClassifier_mono_njobs1.sav"
t0 = time.perf_counter(); model = joblib.load(mp); scaler = joblib.load(mp.replace("model_", "scaler_")); load_model_s = time.perf_counter() - t0
model.n_jobs = n_jobs
rng = np.random.default_rng(0); N, CH = 1_000_000, 100_000
idx = rng.integers(0, X.shape[0], N)
rss_before = rss()
t_tr = t_pr = 0.0; out = np.empty(N)
t_all = time.perf_counter()
for s in range(0, N, CH):
    xb = X[idx[s:s + CH]]
    t0 = time.perf_counter(); xs = scaler.transform(xb); t_tr += time.perf_counter() - t0
    t0 = time.perf_counter(); out[s:s + CH] = model.predict_proba(xs)[:, 1]; t_pr += time.perf_counter() - t0
t_all = time.perf_counter() - t_all
res = dict(TF=TF, n_jobs=n_jobs, k_features=int(X.shape[1]), n_rows=N, chunk=CH, model_load_s=round(load_model_s, 3),
           scaler_transform_s=round(t_tr, 2), predict_proba_s=round(t_pr, 2), total_s=round(t_all, 2),
           rows_per_s=int(N / t_all), rss_before_gb=rss_before, max_rss_gb=rss(), mean_p=round(float(out.mean()), 4))
print(json.dumps(res))
with open(f"{R}/rf_predict_chunked.jsonl", "a") as f: f.write(json.dumps(res) + "\n")
