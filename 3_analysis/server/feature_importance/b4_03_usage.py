"""b4_03_usage.py - feature importance, step 3: subtype usage per TF (human test set chr1/8/21, saved mono+di RF).
(a) Activation: per feature, P_j(x) = fraction of the test NEGATIVES with raw score >= x (rank-normalised against negatives).
    Strong hit = P_j <= 0.01 (top 1 % of negatives). Best PWM of a sequence = feature with the smallest P_j (ties broken by the
    larger z-score against the negatives). Reported for positives and, as null, for negatives:
      - fraction of sequences whose best PWM is in a non-primary subtype (primary = subtype 1 = largest summed Gini importance);
        compare with the fraction of features that are non-primary;
      - groups by strong hits: primary only / secondary only / both / none, with the RF score (median, and fraction above the
        99th percentile of the negatives' RF score = recall at 1 % FPR).
    'secondary' is defined twice: any non-primary subtype, and subtype 2 (= top secondary subtype) alone.
(b) Grouped permutation: all columns of a subtype are permuted with ONE shared row permutation (5 repeats, seed 0); drop in auROC /
    auPRC of the saved RF. Done for the subtypes carrying >= 10 % of the importance or the top 2 subtypes.
Thresholds for subtypes: 0.8 (main), 0.7 and 0.9 for the best-PWM fraction.
Outputs: usage/<TF>_best_pwm.csv, usage/<TF>_groups.csv, usage/<TF>_grouped_perm.csv
"""
import argparse, joblib
from b4_common import *

NREP, SEED, STRONG = 5, 0, 0.01


def pvals(X, y):
    Xn = X[y == 0]
    n = len(Xn)
    P = np.empty(X.shape, dtype=np.float32); Z = np.empty(X.shape, dtype=np.float32)
    for j in range(X.shape[1]):
        s = np.sort(Xn[:, j])
        P[:, j] = (n - np.searchsorted(s, X[:, j], side="left")) / n    # fraction of negatives >= x
        sd = Xn[:, j].std() or 1.0
        Z[:, j] = (X[:, j] - Xn[:, j].mean()) / sd
    return P, Z


def run(tf, n_jobs):
    C = np.load(os.path.join(CACHE, "%s_test.npz" % tf), allow_pickle=True)
    X, y, rf, feats = C["X"], C["y"], C["rf"], list(C["feats"])
    S = pd.read_csv(os.path.join(OUT, "subtypes", "%s_subtypes.csv" % tf)).set_index("feature").loc[feats]
    P, Z = pvals(X, y)
    key = P - 1e-9 * np.clip(Z, -1000, 1000)
    best = key.argmin(axis=1)
    strong = P <= STRONG
    rf_thr = np.quantile(rf[y == 0], 0.99)
    os.makedirs(os.path.join(OUT, "usage"), exist_ok=True)
    rows, grows = [], []
    for thr in (70, 80, 90):
        st = S["subtype_%02d" % thr].values
        nonprim_feat = float((st != 1).mean())
        for cls, m in (("pos", y == 1), ("neg", y == 0)):
            bsub = st[best[m]]
            ms = m & strong.any(axis=1)
            bs = st[best[ms]]
            rows.append(dict(TF=tf, thr=thr / 100, cls=cls, n=int(m.sum()), n_with_strong_hit=int(ms.sum()),
                             frac_best_nonprimary_strong=float((bs != 1).mean()) if ms.any() else np.nan,
                             frac_best_nonprimary=float((bsub != 1).mean()),
                             frac_best_subtype2=float((bsub == 2).mean()), frac_features_nonprimary=nonprim_feat,
                             n_subtypes=int(len(set(st)))))
        if thr != 80:
            continue
        for defn, sec in (("any_nonprimary", st != 1), ("subtype2", st == 2)):
            if sec.sum() == 0:
                continue
            ph = strong[:, st == 1].any(axis=1); sh = strong[:, sec].any(axis=1)
            grp = np.where(ph & ~sh, "primary_only", np.where(~ph & sh, "secondary_only", np.where(ph & sh, "both", "none")))
            for cls, m in (("pos", y == 1), ("neg", y == 0)):
                for g in ("primary_only", "secondary_only", "both", "none"):
                    mm = m & (grp == g)
                    grows.append(dict(TF=tf, defn=defn, cls=cls, group=g, n=int(mm.sum()), frac=float(mm.sum() / m.sum()),
                                      median_rf=float(np.median(rf[mm])) if mm.any() else np.nan,
                                      frac_rf_above_neg99=float((rf[mm] > rf_thr).mean()) if mm.any() else np.nan,
                                      n_sec_features=int(sec.sum()), n_prim_features=int((st == 1).sum())))
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "usage", "%s_best_pwm.csv" % tf), index=False)
    pd.DataFrame(grows).to_csv(os.path.join(OUT, "usage", "%s_groups.csv" % tf), index=False)
    # grouped permutation
    model = joblib.load(model_path(tf))  # human RF
    model.n_jobs = n_jobs
    Xs = scale_test(X, tf, feats)
    auroc, auprc = ad.auroc_auprc(y, model.predict_proba(Xs)[:, 1])
    st = S["subtype_80"].values
    share = pd.Series(S["gini"].values).groupby(st).sum()
    sel = sorted(set(share.index[share >= 0.10]) | {1, 2})
    prow = []
    for s in sel:
        cols = np.flatnonzero(st == s)
        if len(cols) == 0:
            continue
        rng = np.random.RandomState(SEED)
        for r in range(NREP):
            Xp = Xs.copy(); perm = rng.permutation(len(Xp)); Xp[:, cols] = Xs[perm][:, cols]
            ro, pr = ad.auroc_auprc(y, model.predict_proba(Xp)[:, 1])
            prow.append(dict(TF=tf, subtype=s, n_features=len(cols), gini_share=float(share[s]), repeat=r,
                             d_auroc=auroc - ro, d_auprc=auprc - pr))
    pd.DataFrame(prow).to_csv(os.path.join(OUT, "usage", "%s_grouped_perm.csv" % tf), index=False)
    print(tf, "done", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tfs", nargs="+", default=TFS)
    ap.add_argument("--n-jobs", type=int, default=2)
    a = ap.parse_args()
    for tf in a.tfs:
        run(tf, a.n_jobs)
