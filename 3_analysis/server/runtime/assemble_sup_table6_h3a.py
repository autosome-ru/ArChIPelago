"""Assemble Sup. Table 6 (RF training, RF inference, PWM scanning single-PWM rows) from the runtime rerun, same columns and
rounding as Manuscript_analysis/scripts/make_sup_tables_3_to_6.py, plus old-vs-new comparison tables.
Usage: assemble_sup_table6_h3a.py <runtime dir> <published runtime dir>"""
import os, sys, json
import numpy as np, pandas as pd
H, R1 = sys.argv[1], sys.argv[2]
TFS = ["STAT3", "HNF4A", "CTCF"]
def jl(p): return [json.loads(l) for l in open(p)]
def proc(p): return pd.read_csv(p, sep="\t")
out = os.path.join(H, "tables"); os.makedirs(out, exist_ok=False)

def train_table(J, P, label):
    rows = []
    for tf in TFS:
        for nj in (1, 4):
            j = [r for r in J if r["TF"] == tf and r["step"] == "train" and r["n_jobs"] == nj]
            assert len(j) == 1, (label, tf, nj, len(j)); j = j[0]
            p = P[(P.TF == tf) & (P.step == "train") & (P.n_jobs == nj)]; assert len(p) == 1; p = p.iloc[0]
            rows.append([tf, j["k_features"], "{:,} ({:,} + {:,})".format(j["n_rows"], j["n_pos"], j["n_neg"]) if nj == 1 else "{:,}".format(j["n_rows"]),
                         nj, round(j["fit_s"], 1), j["scaler_s"], round(p.wall_s, 1), round(p.user_s, 1), round(p.max_rss_kb / 1e6, 2),
                         round(j["model_bytes"] / 1e6, 2), j["train_auc_sanity"]])
    return pd.DataFrame(rows, columns=["TF", "Features (monoPWMs)", "Training rows (positives + negatives)", "Threads (n_jobs)",
        "RandomForest fit (s, in-process)", "StandardScaler (s)", "Whole process wall-clock incl. CSV reading (s)",
        "Whole process user CPU (s)", "Peak memory (GB)", "Model size on disk (MB)", "Training auROC (sanity)"])

trH = train_table(jl(os.path.join(H, "rf_timings.jsonl")), proc(os.path.join(H, "rf_process_timings.tsv")), "H3a")
trO = train_table(jl(os.path.join(R1, "rf_timings_v2.jsonl"))[:15], proc(os.path.join(R1, "rf_process_timings_v2.tsv")), "round1")
trC = train_table(jl(os.path.join(H, "calib_round1", "rf_timings.jsonl")), proc(os.path.join(H, "calib_round1", "rf_process_timings.tsv")), "calib")
trH.to_csv(os.path.join(out, "Sup_Table_6_runtime_RF_training_H3a.csv"), index=False)

# inference: chunked (rf_predict_chunked*) + single call (rf_bench predict)
def infer_table(CJ, CP, J):
    rows = []
    for tf in TFS:
        for nj in (1, 4):
            c = [r for r in CJ if r["TF"] == tf and r["n_jobs"] == nj]; assert len(c) == 1, (tf, nj); c = c[0]
            p = CP[(CP.TF == tf) & (CP.n_jobs == nj)].iloc[0]
            rows.append(["chunked, 10 x 100,000 rows", tf, c["k_features"], nj, round(c["scaler_transform_s"], 2), round(c["predict_proba_s"], 1),
                         round(c["total_s"], 1), round(p.wall_s, 1), round(p.max_rss_kb / 1e6, 2), "%dk" % round(1e6 / c["total_s"] / 1e3), c["mean_p"]])
    for tf in TFS:
        for nj in (1, 4):
            j = [r for r in J if r["TF"] == tf and r["step"] == "predict" and r["n_jobs"] == nj]
            if not j: continue
            j = j[0]
            rows.append(["single call on the full 1e6 x k matrix", tf, j["k_features"], nj, j["scaler_transform_1M_s"], j["predict_proba_1M_s"], np.nan,
                         np.nan, round(j["max_rss_gb_self"], 2), "%dk" % round(1e6 / j["predict_proba_1M_s"] / 1e3), j["mean_p"]])
    return pd.DataFrame(rows, columns=["Variant (1,000,000 sequences)", "TF", "Features (monoPWMs)", "Threads (n_jobs)", "scaler.transform (s)",
        "predict_proba (s)", "Total in-process (s)", "Whole process wall-clock (s)", "Peak memory (GB)", "Sequences per second", "mean p"])
JH = jl(os.path.join(H, "rf_timings.jsonl"))
inH = infer_table(jl(os.path.join(H, "rf_predict_chunked.jsonl")), proc(os.path.join(H, "rf_predict_chunked_process.tsv")), JH)
PH = proc(os.path.join(H, "rf_process_timings.tsv"))
for i, r in inH.iterrows():   # whole-process wall of the single-call runs
    if r.iloc[0].startswith("single"):
        inH.loc[i, "Whole process wall-clock (s)"] = round(PH[(PH.TF == r.TF) & (PH.step == "predict") & (PH.n_jobs == r["Threads (n_jobs)"])].wall_s.iloc[0], 1)
inH.to_csv(os.path.join(out, "Sup_Table_6_runtime_RF_inference_H3a.csv"), index=False)
inO = infer_table(jl(os.path.join(R1, "rf_predict_chunked.jsonl")), proc(os.path.join(R1, "rf_predict_chunked_process.tsv")), jl(os.path.join(R1, "rf_timings_v2.jsonl"))[:15])

# scanning: single-PWM rows (medians of 3), same FASTA
def scan_med(p):
    s = proc(p); s = s[s.task == "scan_1pwm"]
    if "jvm_flags" in s: s = s[s.jvm_flags == "AP1"]
    g = s.groupby(["TF", "mode"]).agg(wall=("wall_s", "median"), user=("user_s", "median"), rss=("max_rss_kb", "median"), n=("wall_s", "size")).reset_index()
    return g
sH, sO = scan_med(os.path.join(H, "sarus_timings.tsv")), scan_med(os.path.join(R1, "sarus_timings.tsv"))
sc = sH.merge(sO, on=["TF", "mode"], suffixes=("_H3a", "_round1"))
sc.insert(0, "Task", sc["mode"].map({"mono": "1 PWM", "di": "1 diPWM"}))
sc.to_csv(os.path.join(out, "Sup_Table_6_runtime_PWM_scanning_single_H3a_vs_round1.csv"), index=False)

# old vs new comparison
cmp = trO[["TF", "Threads (n_jobs)", "Training rows (positives + negatives)", "RandomForest fit (s, in-process)", "Whole process wall-clock incl. CSV reading (s)", "Peak memory (GB)"]].copy()
cmp.columns = ["TF", "n_jobs", "rows_round1", "fit_s_round1_Sep", "wall_s_round1_Sep", "rss_gb_round1_Sep"]
cmp["fit_s_round1_today"] = trC["RandomForest fit (s, in-process)"].values
cmp["wall_s_round1_today"] = trC["Whole process wall-clock incl. CSV reading (s)"].values
cmp["rows_H3a"] = trH["Training rows (positives + negatives)"].values
cmp["fit_s_H3a"] = trH["RandomForest fit (s, in-process)"].values
cmp["wall_s_H3a"] = trH["Whole process wall-clock incl. CSV reading (s)"].values
cmp["rss_gb_H3a"] = trH["Peak memory (GB)"].values
cmp["fit_ratio_H3a_vs_round1_today"] = (cmp.fit_s_H3a / cmp.fit_s_round1_today).round(2)
cmp.to_csv(os.path.join(out, "compare_training_old_vs_new.csv"), index=False)
ic = inO[["Variant (1,000,000 sequences)", "TF", "Threads (n_jobs)", "predict_proba (s)", "Total in-process (s)", "Peak memory (GB)", "mean p"]].merge(
     inH[["Variant (1,000,000 sequences)", "TF", "Threads (n_jobs)", "predict_proba (s)", "Total in-process (s)", "Peak memory (GB)", "mean p"]],
     on=["Variant (1,000,000 sequences)", "TF", "Threads (n_jobs)"], suffixes=("_round1", "_H3a"))
ic.to_csv(os.path.join(out, "compare_inference_old_vs_new.csv"), index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
for n, t in (("TRAIN", cmp), ("INFER", ic), ("SCAN", sc)): print("---", n); print(t.to_string(index=False))
