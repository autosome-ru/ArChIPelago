"""b9_random_split_h5.py - does a random split of the scored sequences inflate the test performance compared with the
chromosome hold-out?
Human data, mono+di features, paper RF (max_depth 6, max_samples 0.8, 100 trees), training scaler applied to the test rows.
Pool = the human training set (train chromosomes) + human test set (chr1/8/21), positives and negatives.
  chromosome: the published assignment, RF seeds 0-2.
  random:     the same pool re-assigned at random with the same numbers of training/test positives and negatives
              (seed s shuffles the rows and seeds the RF), s = 0-2.
Per fit: RF and best single monoPWM (chosen on the training set by auROC / by auPRC, ad.best_pwm_baseline) on the test set.
Random split only: test positives split into those with a training positive within 300 bp on the same chromosome
('overlap') and the rest; RF and PWM auROC/auPRC on each subset against all test negatives; exact duplicates (same centre).
Centres: out_tab headers ('chrom: ..; center: ..').
Refuses to overwrite its output.
usage (server): nice python b9_random_split_h5.py OUTDIR TF [TF ...]   env: N_JOBS (default 8)
"""
import os, sys, re, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.expanduser("~/archi_work/human_sym_h5/code"))
import archi_data_v5 as ad
assert ad.NEG_VERSION == "v3"
from sklearn.preprocessing import StandardScaler

OUT, TFS = sys.argv[1], sys.argv[2:]
if os.path.exists(os.path.join(OUT, "b9_random_split_per_fit.csv")):
    raise SystemExit("refusing to overwrite " + os.path.join(OUT, "b9_random_split_per_fit.csv"))
N_JOBS = int(os.environ.get("N_JOBS", 8))
SEEDS = (0, 1, 2)
WIN = 300


def centres(split):
    path = os.path.join(ad.HUMAN_REL, "out_tab_HUMAN_10000_%s.tab" % split)
    d = {}
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2 or "$" not in p[1]:
                continue
            nm, hdr = p[1].split("$", 1)
            m = re.search(r"chrom: (chr\w+); center: ([\d.]+);", hdr)
            d[nm] = (m.group(1), int(float(m.group(2))))
    return d


CEN = {**centres("train"), **centres("control")}
sc = lambda X: StandardScaler().fit_transform(X).astype(np.float32)


def overlap_mask(tr_meta, te_meta):
    """Test rows with a training row within WIN bp on the same chromosome; and exact same centre."""
    tr = pd.DataFrame([CEN[n] for n in tr_meta], columns=["c", "p"])
    te = pd.DataFrame([CEN[n] for n in te_meta], columns=["c", "p"])
    ov = np.zeros(len(te), bool); ex = np.zeros(len(te), bool)
    for c, g in te.groupby("c"):
        tp = np.sort(tr.p[tr.c == c].values)
        if len(tp) == 0:
            continue
        lo = np.searchsorted(tp, g.p.values - WIN); hi = np.searchsorted(tp, g.p.values + WIN, side="right")
        ov[g.index.values] = hi > lo
        e = np.searchsorted(tp, g.p.values); e2 = np.searchsorted(tp, g.p.values, side="right")
        ex[g.index.values] = e2 > e
    return ov, ex


rows = []
for tf in TFS:
    t0 = time.time()
    feats = ad.feature_list(tf, "HUMAN", "mono_di")
    is_mono = np.array([f.startswith("mono_") for f in feats])
    Xa, ya, ma = ad.load_matrix(tf, "HUMAN", "HUMAN", "train", feats)
    Xb, yb, mb = ad.load_matrix(tf, "HUMAN", "HUMAN", "control", feats)
    X = np.vstack([Xa, Xb]); y = np.concatenate([ya, yb]); names = np.concatenate([ma.name.values, mb.name.values])
    is_tr_chr = np.r_[np.ones(len(ya), bool), np.zeros(len(yb), bool)]
    n_tr_pos, n_tr_neg = int(ya.sum()), int((ya == 0).sum())
    pos, neg = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    for split in ("chromosome", "random"):
        for s in SEEDS:
            if split == "chromosome":
                tr = np.flatnonzero(is_tr_chr); te = np.flatnonzero(~is_tr_chr)
            else:
                r = np.random.RandomState(s)
                pp, nn = r.permutation(pos), r.permutation(neg)
                tr = np.r_[pp[:n_tr_pos], nn[:n_tr_neg]]; te = np.r_[pp[n_tr_pos:], nn[n_tr_neg:]]
            scl = StandardScaler().fit(X[tr])  # training scaler, applied to the test rows
            Xtr, Xte = scl.transform(X[tr]).astype(np.float32), scl.transform(X[te]).astype(np.float32); ytr, yte = y[tr], y[te]
            m = ad.fit_paper_rf(Xtr, ytr, n_jobs=N_JOBS, random_state=s)
            srf = m.predict_proba(Xte)[:, 1]
            b = ad.best_pwm_baseline(Xtr, ytr, {"te": (Xte, yte)}, mask=is_mono)
            jr, jp = b["by_roc"]["col"], b["by_prc"]["col"]
            rf_roc, rf_prc = ad.auroc_auprc(yte, srf)
            row = dict(TF=tf, split=split, seed=s, n_train_pos=int(ytr.sum()), n_train_neg=int((ytr == 0).sum()),
                       n_test_pos=int(yte.sum()), n_test_neg=int((yte == 0).sum()),
                       rf_auroc=rf_roc, rf_auprc=rf_prc,
                       pwm_auroc=b["by_roc"]["tests"]["te"][0], pwm_auprc=b["by_prc"]["tests"]["te"][1])
            tp_names = names[te][yte == 1]
            ov, ex = overlap_mask(names[tr][ytr == 1], tp_names)
            row.update(frac_test_pos_overlap=float(ov.mean()), frac_test_pos_exact=float(ex.mean()))
            negi = np.flatnonzero(yte == 0); posi = np.flatnonzero(yte == 1)
            for lab, msk in (("ov", ov), ("noov", ~ov)):
                if msk.sum() >= 50:
                    idx = np.r_[posi[msk], negi]; yy = yte[idx]
                    row["rf_auroc_" + lab], row["rf_auprc_" + lab] = ad.auroc_auprc(yy, srf[idx])
                    row["pwm_auroc_" + lab] = ad.auroc_auprc(yy, Xte[idx, jr])[0]
                    row["pwm_auprc_" + lab] = ad.auroc_auprc(yy, Xte[idx, jp])[1]
                    row["n_pos_" + lab] = int(msk.sum())
            rows.append(row)
            print(tf, split, s, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items() if k not in ("TF", "split", "seed")}, flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "b9_random_split_per_fit.csv"), index=False)
    print(tf, "done %.0f s" % (time.time() - t0), flush=True)
print("ALL DONE", flush=True)
