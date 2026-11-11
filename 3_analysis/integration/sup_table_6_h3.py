"""Sup. Table 6 (runtime and memory) on the final training data.
Sources
  h3_inputs runtime = the server mirror human/runtime/ (benchmark rerun on the human training matrices; README.md there):
    tables/Sup_Table_6_runtime_RF_training_H3a.csv, tables/Sup_Table_6_runtime_RF_inference_H3a.csv,
    tables/Sup_Table_6_runtime_PWM_scanning_single_H3a_vs_round1.csv (medians of 3; outputs byte-identical to round 1),
    rf_timings.jsonl (matrix load times), load_during_run.log (1-min load), machine.txt
  the released Sup. Table 6 (results/round1_package/Manuscript_analysis, read-only): JVM start-up, the pipeline-flags row and the k-PWM loops, measured for the first submission and
  kept (scan time does not depend on the training data); marked in the column "Measured"
  hm_mm_results.csv: training rows and positives of the 3 TFs must equal the benchmark matrices
Every number in the notes is computed from these tables. Also writes 3_analysis/integration/sup6_methods_values_h3.csv: the values that
Methods 2.11 quotes, first submission (released table) and final run, for the NUMBERS file.
The companion-scanner formula of RF_inference quotes the per-Mb rate of the single scans (mean) and of the loops (lowest),
both from sheet PWM_scanning, instead of the first-submission 0.15 s/Mb.
usage: import sup_table_6_h3 as S6; S6.build(out_xlsx)   (called by 3_analysis/supplementary_tables/make_sup_tables_1_3_R2.py 6)
"""
import json
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w2  # noqa: E402
import h3_inputs as H  # noqa: E402

R2, B13 = w2.R2, w2.B13
RT = H.IN["runtime"]   # not rerun for the previous training run (timing only; h3_inputs NOT_IN_H4)
IN = {"train": os.path.join(RT, "tables", "Sup_Table_6_runtime_RF_training_H3a.csv"),
      "inf": os.path.join(RT, "tables", "Sup_Table_6_runtime_RF_inference_H3a.csv"),
      "scan": os.path.join(RT, "tables", "Sup_Table_6_runtime_PWM_scanning_single_H3a_vs_round1.csv"),
      "jsonl": os.path.join(RT, "rf_timings.jsonl"), "load": os.path.join(RT, "load_during_run.log"), "machine": os.path.join(RT, "machine.txt"),
      "pkg6": os.path.join(H.PKG, "Sup_Tables", "Sup_Table_6_runtime_memory.xlsx"), "hm": os.path.join(R2, H.IN["hm_mm"])}
NEW, OLD = "10 Nov 2026", "10 Sep 2026"   # the runtime folder (h3_inputs "runtime") holds the rerun of 10 Nov
MB = 15.0
TFS = ["STAT3", "HNF4A", "CTCF"]
rng = lambda v, f="%.2f": (f + "-" + f) % (min(v), max(v)) if f % min(v) != f % max(v) else f % min(v)


def build(out):
    old = pd.read_excel(IN["pkg6"], sheet_name=None)
    tr, inf, sc = pd.read_csv(IN["train"]), pd.read_csv(IN["inf"]), pd.read_csv(IN["scan"])
    js = [json.loads(l) for l in open(IN["jsonl"]) if l.strip()]
    load = {d["TF"]: d for d in js if d["step"] == "load"}
    trj = {(d["TF"], d["n_jobs"]): d for d in js if d["step"] == "train"}
    hm = pd.read_csv(IN["hm"]).set_index("TF")
    for tf in TFS:   # the benchmark matrices are the training sets
        if (load[tf]["n_pos"], load[tf]["n_neg"]) != (hm.loc[tf, "n_train_pos_H"], hm.loc[tf, "n_train_neg_H"]):
            raise SystemExit("CHECK FAILED %s: benchmark matrix %d + %d rows, hm_mm %d + %d" % (
                tf, load[tf]["n_pos"], load[tf]["n_neg"], hm.loc[tf, "n_train_pos_H"], hm.loc[tf, "n_train_neg_H"]))
    if not (sc.n_H3a == 3).all():
        raise SystemExit("CHECK FAILED single-PWM scans: not 3 repeats each")
    lv = np.array([float(l.split()[1]) for l in open(IN["load"]) if re.match(r"^\d\d:\d\d:\d\d ", l)])
    mach = open(IN["machine"]).read()
    nthr = int(re.search(r"CPU\(s\):\s+(\d+)", mach).group(1))

    # ---- PWM_scanning: single-PWM rows from the rerun, the rest of the released sheet as measured on 10 Sep
    ps_old = old["PWM_scanning"]
    task = ps_old.columns[0]
    keep = ps_old[ps_old[task].notna() & ps_old.TF.notna()].copy()
    single = keep[task].isin(["1 PWM", "1 diPWM"])
    new_rows = []
    sc = sc.assign(_o=sc.TF.map({t: k for k, t in enumerate(TFS)}), _m=(sc["mode"] == "di").astype(int)).sort_values(["_m", "_o"])
    for _, r in sc.iterrows():
        new_rows.append({task: r.Task, "TF": r.TF, "PWM type": r["mode"], "Number of PWMs": 1, "Scanned sequence (Mb)": MB,
                         "Wall-clock (s)": r.wall_H3a, "User CPU (s)": r.user_H3a, "Peak memory (GB)": round(r.rss_H3a / 1e6, 2), "Measured": NEW})
    kept = keep[~single].assign(Measured=OLD)
    ps = pd.concat([kept.iloc[:1], pd.DataFrame(new_rows), kept.iloc[1:]], ignore_index=True)
    ps["Number of PWMs"] = ps["Number of PWMs"].astype(int)
    s1 = ps[ps.Measured == NEW]
    rate = lambda d: d["Wall-clock (s)"] / (d["Number of PWMs"] * d["Scanned sequence (Mb)"])
    mono1, di1 = rate(s1[s1["PWM type"] == "mono"]), rate(s1[s1["PWM type"] == "di"])
    loops = ps[ps[task].str.startswith("all ")]
    monoL, diL = rate(loops[loops["PWM type"] == "mono"]), rate(loops[loops["PWM type"] == "di"])
    cpu1 = s1["User CPU (s)"] / s1["Wall-clock (s)"]
    jvm = float(ps.loc[ps[task].str.startswith("JVM"), "Wall-clock (s)"].iloc[0])
    loopmin = loops.assign(m=loops["Wall-clock (s)"] / 60)
    notes_ps = [
        "Derived quantities",
        "Single PWM over 15 Mb (%s, medians of 3 runs; SARUS outputs byte-identical to the %s runs): %s s wall-clock per monoPWM per Mb "
        "(mean %.2f) and %s s per diPWM per Mb (mean %.2f), JVM start-up included." % (NEW, OLD, rng(mono1, "%.3f"), mono1.mean(), rng(di1, "%.3f"), di1.mean()),
        "Sequential loops over all PWMs of a TF (%s, one run each; scan time does not depend on the training data): %s s per monoPWM per Mb, "
        "%s s per diPWM per Mb; JVM start-up %.2f s per invocation." % (OLD, rng(monoL, "%.3f"), rng(diL, "%.3f"), jvm),
        "CPU time is higher than wall-clock because the JVM keeps JIT/GC helper threads with -XX:ActiveProcessorCount=1: user CPU / wall-clock "
        "= %s for the single scans of %s." % (rng(cpu1, "%.1f"), NEW),
        "Scan time is linear in the number of PWMs k; the full feature step of the loops took %s." % ", ".join(
            "%s (%s, %d %sPWMs)" % ("%.0f s" % (60 * r.m) if r.m < 1 else "%.1f min" % r.m, r.TF, r["Number of PWMs"], r["PWM type"]) for _, r in loopmin.iterrows()),
        "Memory is flat (%s GB per mono scan, %s GB per di scan) because PWMs are scanned one at a time; the scans are independent, so wall "
        "time divides by the number of cores used." % (rng(ps.loc[ps["PWM type"] == "mono", "Peak memory (GB)"][1:]), rng(ps.loc[ps["PWM type"] == "di", "Peak memory (GB)"])),
        "Extrapolation (not measured): one monoPWM over the human genome (3.1 Gb) ~ %.0f min wall-clock at the mean single-scan rate of %s." % (
            3100 * mono1.mean() / 60, NEW),
    ]
    ps = pd.concat([ps, pd.DataFrame([{task: None}] * 2), pd.DataFrame({task: notes_ps})], ignore_index=True)

    # ---- RF_training
    cols_tr = list(old["RF_training"].columns)
    t = tr[cols_tr].copy()
    sp = {tf: trj[(tf, 1)]["fit_s"] / trj[(tf, 4)]["fit_s"] for tf in TFS}
    notes_tr = [
        "Notes",
        "RandomForestClassifier(max_depth = 6, max_samples = 0.8, n_estimators = 100, random_state = 0) on StandardScaler-transformed "
        "features, as in the paper; training matrices of the final models (%s). Training auROC (sanity check): %s." % (
            NEW, ", ".join("%s %.2f" % (tf, trj[(tf, 1)]["train_auc_sanity"]) for tf in TFS)),
        "The pickled model is %s MB for every TF. Peak memory is dominated by the feature matrix, not by the forest. Speed-up of the fit "
        "with 4 threads: %s x. Each configuration was run once; the matrices hold the monoPWM features only (the mono+di models have about "
        "twice as many features)." % (rng(t["Model size on disk (MB)"]), rng(list(sp.values()), "%.1f")),
        "Reading the feature matrix takes %s." % ", ".join("%.2f s (%s, %.0f MB)" % (load[tf]["load_s"], tf, load[tf]["csv_mb"]) for tf in TFS),
    ]
    t = pd.concat([t, pd.DataFrame([{cols_tr[0]: None}] * 2), pd.DataFrame({cols_tr[0]: notes_tr})], ignore_index=True)

    # ---- RF_inference
    cols_in = list(old["RF_inference"].columns)
    OLD_RATE = "0.15 s/Mb x Mb"
    scanner_note = old["RF_inference"].iloc[:, 0].dropna().iloc[-1]
    if scanner_note.count(OLD_RATE) != 1:
        raise SystemExit("CHECK FAILED Sup. Table 6 companion scanner note: %r not found once in %r" % (OLD_RATE, scanner_note))
    scanner_note = scanner_note.replace(OLD_RATE, "r x Mb, with r = %.2f-%.2f s/Mb (single scans and sequential loops, sheet PWM_scanning)" % (
        mono1.mean(), monoL.min()))
    i = inf[cols_in].copy()
    ch, sg = inf[inf[cols_in[0]].str.startswith("chunked")], inf[inf[cols_in[0]].str.startswith("single")]
    p1, p4 = ch[ch["Threads (n_jobs)"] == 1], ch[ch["Threads (n_jobs)"] == 4]
    s1i = sg[sg["Threads (n_jobs)"] == 1]
    notes_in = [
        "Notes",
        "predict_proba takes %s s per 10^6 sequences single-threaded in 100,000-row chunks (%s s on 4 threads; %s s as one call on the "
        "full matrix); with scaler.transform included the whole single-threaded chunked loop takes %s s. 'Sequences per second' = 10^6 / total in-process "
        "time for the chunked rows and 10^6 / predict_proba for the single-call rows." % (rng(p1["predict_proba (s)"], "%.1f"),
        rng(p4["predict_proba (s)"], "%.1f"), rng(s1i["predict_proba (s)"], "%.1f"), rng(p1["Total in-process (s)"], "%.1f")),
        "Feeding one monolithic 10^6 x k matrix raises peak memory to %s GB (CTCF); the chunked numbers are the representative deployment "
        "figures (identical results)." % rng(sg.loc[sg.TF == "CTCF", "Peak memory (GB)"]),
        scanner_note,   # companion scanner note: the per-Mb rate of the single scans and of the loops
    ]
    i = pd.concat([i, pd.DataFrame([{cols_in[0]: None}] * 2), pd.DataFrame({cols_in[0]: notes_in})], ignore_index=True)

    # ---- setup
    su = old["setup"].set_index("Item")
    su.loc["CPU", "Value"] = ("Intel Xeon E5-4607 v2 @ 2.60 GHz (4 sockets x 6 cores, 2 threads per core = %d threads; shared server). "
                              "1-minute load during the %s runs %.2f-%.2f (median %.2f) of %d; during the %s PWM loops 3-8." % (
                                  nthr, NEW, lv.min(), lv.max(), np.median(lv), nthr, OLD))
    su.loc["OS", "Value"] = "Ubuntu 24.04.5 LTS, kernel 6.8.0 (%s runs); the %s rows ran on Ubuntu 20.04.6 LTS, kernel 5.4.0" % (NEW, OLD)
    su.loc["Java / SARUS", "Value"] = su.loc["Java / SARUS", "Value"].replace(
        "OpenJDK 1.8.0_232;", "OpenJDK 21.0.12 (%s runs; %s rows: OpenJDK 1.8.0_232);" % (NEW, OLD))
    su.loc["Training matrices", "Value"] = ("the human training matrices of the final models (human positives with each peak set read "
        "once, negatives v3, chromosomes %s): %s; features = best-hit log-odds of the monoPWMs" % (H.SPLIT["H_train"][3:], ", ".join(
            "%s %s rows (%s positives)" % (tf, format(load[tf]["n_rows"], ","), format(load[tf]["n_pos"], ",")) for tf in TFS)))
    rf_cpu = [trj[(tf, 1)] for tf in TFS]
    su.loc["Thread policy", "Value"] = (
        "SARUS with the JVM limited to one active processor (JIT/GC helper threads still run: user CPU = %s x wall-clock); Random Forest "
        "with n_jobs = 1 and n_jobs = 4 (OMP_NUM_THREADS = 4); all jobs at nice 10. Repeats: single-PWM scans 3 (medians reported), PWM "
        "loops and RF runs 1. Feature matrices of the RF benchmarks contain the monoPWMs only (%s features); the mono+di models of the "
        "paper have about twice as many features and scale about linearly." % (rng(cpu1, "%.1f"), " / ".join(str(int(load[tf]["k_features"])) for tf in TFS)))
    if "Python" not in su.index or "scikit-learn 1.3.0" not in su.loc["Python", "Value"]:
        raise SystemExit("CHECK FAILED setup: Python line changed")
    su = su.reset_index()
    w2.backup(out)
    with pd.ExcelWriter(out) as xw:
        for name, d in (("setup", su), ("PWM_scanning", ps), ("RF_training", t), ("RF_inference", i)):
            d.to_excel(xw, sheet_name=name, index=False)

    # ---- values of Methods 2.11 (round 1 = released table; h3 = this table)
    otr, oin, ops = old["RF_training"], old["RF_inference"], old["PWM_scanning"]
    otr1 = otr[otr["Threads (n_jobs)"] == 1]
    o1 = ops[ops[task].isin(["1 PWM", "1 diPWM"])]
    oc = oin[oin[cols_in[0]].astype(str).str.startswith("chunked") & (oin["Threads (n_jobs)"] == 1)]
    trn1 = tr[tr["Threads (n_jobs)"] == 1]
    rows_n = [load[tf]["n_rows"] for tf in TFS]
    rows_o = [int(str(x).split()[0].replace(",", "")) for x in otr1["Training rows (positives + negatives)"]]
    V = [
        ("SARUS wall-clock per monoPWM per Mb (mean of the single-PWM scans of STAT3, HNF4A, CTCF)", "%.2f" % rate(o1[o1["PWM type"] == "mono"]).mean(), "%.2f" % mono1.mean()),
        ("SARUS wall-clock per diPWM per Mb (mean of the single-PWM scans)", "%.2f" % rate(o1[o1["PWM type"] == "di"]).mean(), "%.2f" % di1.mean()),
        ("RF training matrices, rows (STAT3 / HNF4A / CTCF)", " / ".join(format(x, ",") for x in rows_o), " / ".join(format(x, ",") for x in rows_n)),
        ("RF training matrices, range of rows", "%s-%s" % (format(min(rows_o), ","), format(max(rows_o), ",")), "%s-%s" % (format(min(rows_n), ","), format(max(rows_n), ","))),
        ("RF fit on one core, s (STAT3 / HNF4A / CTCF)", " / ".join("%.1f" % x for x in otr1["RandomForest fit (s, in-process)"]),
         " / ".join("%.1f" % x for x in trn1["RandomForest fit (s, in-process)"])),
        ("RF training, maximum peak memory incl. the feature table (GB)", "%.2f" % otr["Peak memory (GB)"].max(), "%.2f" % tr["Peak memory (GB)"].max()),
        ("RF inference, predict_proba for 10^6 sequences on one core, chunked (s, range)", rng(oc["predict_proba (s)"].astype(float), "%.1f"), rng(p1["predict_proba (s)"], "%.1f")),
        ("RF model size on disk (MB, range)", rng(otr["Model size on disk (MB)"].dropna()), rng(tr["Model size on disk (MB)"])),
        ("benchmark system (OS; Java)", "Ubuntu 20.04.6, kernel 5.4.0; OpenJDK 1.8.0_232",
         "Ubuntu 24.04.5, kernel 6.8.0; OpenJDK 21.0.12 (PWM loops: %s on the old system)" % OLD),
        ("benchmark CPU", "Intel Xeon E5-4607 v2, 2.6 GHz, %d threads" % nthr, "Intel Xeon E5-4607 v2, 2.6 GHz, %d threads" % nthr),
    ]
    mv = pd.DataFrame(V, columns=["quantity", "round1_released", "h3"])
    dst = os.path.join(B13, "sup6_methods_values_h3.csv")
    w2.backup(dst)
    mv.to_csv(dst, index=False)
    print(mv.to_string(index=False))
    w2.provenance(os.path.dirname(out), "3_analysis/integration/sup_table_6_h3.py (Sup. Table 6)", [out, dst], list(IN.values()),
                  note="V1-3: RF rows and single-PWM scans rerun on the H3a training matrices (server, %s); k-PWM loops of %s kept" % (NEW, OLD))
    return out
