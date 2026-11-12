"""b3_topN_h5.py - top-N re-evaluation of RF (human-trained, seed 0) vs best monoPWM on the strongest test peaks.
Usage: python b3_topN_h5.py OUTDIR (environment variable ARCHI_WORK)
Inputs (read-only): the score npz of run_hm_mm_h5.py, TF-ML/out_tab_{HUMAN_10000_control,MOUSE_10000_train}.tab
(human test = file 'control', mouse test = file 'train'; names are inverted on disk, see archi_data_v5.py).
Selection: top-N positives by MACS signal (header 'signal:' = -10log10 P of the source peak), ties broken by file order.
Negatives: variant allneg = all; r100 = min(100*n_pos_kept, all) random, seed 0; rorig = n_pos_kept*(neg/pos of the full test set), seed 0.
Control 'random': N random positives (mean of NDRAW draws, seeds 0..NDRAW-1) with the same negatives rule.
Baselines: pwmH_byroc (for auROC) and pwmH_byprc (for auPRC), selected on the human training set.
"""
import os, re, sys, json, glob
import numpy as np, pandas as pd

OUT = sys.argv[1]
sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "..")); import archi_env as AE
sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "..", "1_data_preparation")); import prroc_metrics
REL = AE.work("human_sym_h5", "TF-ML")  # human data build
REL_M = AE.work("mouse_sym_h5", "TF-ML")  # mouse data build
SC = AE.work("H5_alltrain", "human_sym", "scores")
NS = [500, 1000, 2000]
NDRAW = 10
SEED = 0


def read_tab(path):
    d = {}
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            nm, hdr = p[1].split("$", 1)
            m = re.search(r"chrom: (\w+); center: ([\d.]+); peak: \d+; signal: ([\d.]+)", hdr)
            assert m, line[:200]
            assert nm not in d, "duplicate name " + nm
            d[nm] = (m.group(1), float(m.group(2)), float(m.group(3)))
    return d


TAB = {"H": read_tab(f"{REL}/out_tab_HUMAN_10000_control.tab"),
       "M": read_tab(f"{REL_M}/out_tab_MOUSE_10000_train.tab")}
print("out_tab rows", {k: len(v) for k, v in TAB.items()}, flush=True)


def pick_neg(neg_idx, n, seed):
    if n >= len(neg_idx):
        return neg_idx
    return np.random.RandomState(seed).choice(neg_idx, n, replace=False)


def metrics(y, rf, pr_roc, pr_prc, idx):
    yy = y[idx]
    return (prroc_metrics.auroc(yy, rf[idx]), prroc_metrics.auroc(yy, pr_roc[idx]),
            prroc_metrics.auprc(yy, rf[idx]), prroc_metrics.auprc(yy, pr_prc[idx]))


def neg_n(variant, nk, npos, nneg):
    if variant == "allneg":
        return nneg
    if variant == "r100":
        return min(100 * nk, nneg)
    return int(round(nk * nneg / npos))  # rorig


rows, sigrows = [], []
log = open(f"{OUT}/b3_topN_selection_check.txt", "w")
for fn in sorted(glob.glob(SC + "/*.npz")):
    tf = os.path.basename(fn)[:-4]
    z = np.load(fn, allow_pickle=True)
    for sp in ("H", "M"):
        y = z[f"y_test_{sp}"].astype(int); names = z[f"names_test_{sp}"]
        rf = z[f"rf_{sp}_s0"].astype(float)
        pr_roc = z[f"pwmH_byroc_test_{sp}"].astype(float); pr_prc = z[f"pwmH_byprc_test_{sp}"].astype(float)
        pos = np.flatnonzero(y == 1); neg = np.flatnonzero(y == 0)
        sig = np.array([TAB[sp][n][2] for n in names[pos]])  # KeyError = join failure
        assert all(n.startswith(tf + "_") for n in names[pos]), tf
        sigrows.append(dict(TF=tf, species=sp, n_pos=len(pos), n_pos_in_table=len(sig), sig_min=sig.min(), sig_median=np.median(sig), sig_max=sig.max()))
        order = np.argsort(-sig, kind="mergesort")  # strongest first, stable
        for N in NS + ["all"]:
            nk = len(pos) if N == "all" else N
            if N != "all" and nk > len(pos):
                continue
            keep = pos[order[:nk]]
            if N != "all":
                dropped = sig[order[nk:]]
                kept = sig[order[:nk]]
                if len(dropped):
                    assert kept.min() >= dropped.max(), (tf, sp, N)
                if tf in ("ANDR", "ZEB1") and N == 500 or (tf == "AP2A" and sp == "H" and N == 500):
                    log.write(f"{tf} {sp} N={N}: first kept signal {kept[0]:.2f}, last kept {kept[-1]:.2f}, first dropped {dropped[0] if len(dropped) else float('nan'):.2f}; kept min {kept.min():.2f} >= dropped max {dropped.max() if len(dropped) else float('nan'):.2f}\n")
            for variant in ("allneg", "r100", "rorig"):
                nn = neg_n(variant, nk, len(pos), len(neg))
                ng = pick_neg(neg, nn, SEED)
                idx = np.r_[keep, ng]
                m = metrics(y, rf, pr_roc, pr_prc, idx)
                rows.append(dict(TF=tf, species=sp, selection="top", variant=variant, N=N, n_pos=nk, n_neg=len(ng), neg_per_pos=len(ng) / nk,
                                 min_signal_kept=float(sig[order[:nk]].min()),
                                 rf_auroc=m[0], pwm_auroc=m[1], rf_auprc=m[2], pwm_auprc=m[3]))
                if N == "all":
                    continue
                ms = []
                for dr in range(NDRAW):
                    kp = np.random.RandomState(1000 + dr).choice(pos, nk, replace=False)
                    ng2 = pick_neg(neg, nn, SEED + dr)
                    ms.append(metrics(y, rf, pr_roc, pr_prc, np.r_[kp, ng2]))
                ms = np.mean(ms, axis=0)
                rows.append(dict(TF=tf, species=sp, selection="random", variant=variant, N=N, n_pos=nk, n_neg=len(ng), neg_per_pos=len(ng) / nk,
                                 min_signal_kept=np.nan, rf_auroc=ms[0], pwm_auroc=ms[1], rf_auprc=ms[2], pwm_auprc=ms[3]))
    print(tf, "done", flush=True)
log.close()
df = pd.DataFrame(rows)
df["d_auroc"] = df.rf_auroc - df.pwm_auroc
df["d_auprc"] = df.rf_auprc - df.pwm_auprc
df.to_csv(f"{OUT}/b3_topN_results.csv", index=False)
pd.DataFrame(sigrows).to_csv(f"{OUT}/b3_topN_signal_summary.csv", index=False)
print("written", len(df))
