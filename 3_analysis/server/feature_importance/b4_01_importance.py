"""b4_01_importance.py - feature importance, step 1: Gini importance, test-set metrics, permutation importance (top 20 Gini features), PWM origin check.
Per TF writes (in $OUT): gini/<TF>_gini.csv, perm/<TF>_perm.csv, cache/<TF>_test.npz, and appends metrics via
metrics/<TF>_metrics.csv. Model: finalized_model_<TF>_HUMAN_RandomForestClassifier_mono_di_full_MODEL_all_features.sav
Permutation: top-20 Gini features, 5 repeats, RandomState(seed 0), drop in auROC and auPRC (PRROC-style integral, ad.auprc_integral)
on the full human test set (chr1/8/21; negatives would be subsampled to 100,000 with seed 0 only if more - never triggered).
usage: python b4_01_importance.py --tfs SRF P53 ... [--n-jobs 2]
"""
import argparse, time, re, json
import joblib
from sklearn.metrics import roc_auc_score
from b4_common import *

TOPN, NREP, SEED = 20, 5, 0


def run(tf, n_jobs, retrain=False):
    t0 = time.time()
    path = model_path(tf)  # human RF, mono+di, seed 0
    if not retrain:
        model = joblib.load(path)
        feats = feats_in_model_order(tf, model)
        model.n_jobs = n_jobs
        source = "h3a_human_sym_negv3_seed0"
    else:
        # JUND: the .sav has 50 features but the file code_H_base_1_RandomForestClassifier_mono_di.sh was overwritten with a
        # 25-feature run, so the column order of the .sav is lost. Seed-0 retrain, same hyperparameters, release scaling.
        from sklearn.ensemble import RandomForestClassifier
        feats = ad.feature_list(tf, "HUMAN", "mono_di")
        Xtr, ytr, _ = ad.load_matrix(tf, "HUMAN", "HUMAN", "train", feats)
        model = RandomForestClassifier(n_jobs=n_jobs, random_state=0, **ad.RF_PARAMS).fit(scale_test(Xtr, tf, feats), ytr)
        source = "retrained_seed0"
        os.makedirs(os.path.join(OUT, "models"), exist_ok=True)
        joblib.dump(model, os.path.join(OUT, "models", "%s_retrained_seed0.sav" % tf))
    gini = model.feature_importances_
    # --- PWM origin check
    rows = []
    for f, g in zip(feats, gini):
        header, _ = read_pwm_file(tf, f)
        fn = feature_file_name(tf, f)
        m = re.match(r"(\w+)_HUMAN~C([MD])~(\w+)_HUMAN\.(PEAKS\d+)\.macs\.([MS])$", header)
        assert m, header
        fm = re.match(r"\d+_(\w+)_HUMAN~C([MD])~(\w+)_HUMAN_(PEAKS\d+)_macs_([MS])_feature_compare_table", fn)
        assert fm, fn
        rows.append(dict(TF=tf, feature=f, kind=f.split("_")[0], idx=int(f.split("_")[1]), pwm_header=header,
                         pwm_tf=m.group(1), dataset=m.group(4), variant=m.group(5), feature_file=fn,
                         header_tf_is_target=(m.group(1) == tf and m.group(3) == tf),
                         file_matches_header=(fm.groups() == (m.group(1), m.group(2), m.group(3), m.group(4), m.group(5))),
                         gini=g))
    G = pd.DataFrame(rows)
    G["rank"] = G.gini.rank(ascending=False, method="first").astype(int)
    G = G.sort_values("rank").reset_index(drop=True)
    G["cum_gini"] = G.gini.cumsum()
    os.makedirs(os.path.join(OUT, "gini"), exist_ok=True)
    G.to_csv(os.path.join(OUT, "gini", "%s_gini.csv" % tf), index=False)
    # --- test matrix
    X, y, meta = load_test(tf, feats)
    if (y == 0).sum() > NEG_CAP:
        rng = np.random.RandomState(SEED)
        neg = np.flatnonzero(y == 0); keep = np.sort(np.r_[np.flatnonzero(y == 1), rng.choice(neg, NEG_CAP, replace=False)])
        X, y, meta = X[keep], y[keep], meta.iloc[keep].reset_index(drop=True)
    Xs = scale_test(X, tf, feats)
    p = model.predict_proba(Xs)[:, 1]
    auroc, auprc = ad.auroc_auprc(y, p)
    os.makedirs(CACHE, exist_ok=True)
    np.savez_compressed(os.path.join(CACHE, "%s_test.npz" % tf), X=X, y=y, rf=p.astype(np.float32),
                        names=meta["name"].values.astype(str), feats=np.array(feats))
    # --- permutation importance on the top-20 Gini features
    fidx = {f: i for i, f in enumerate(feats)}
    prow = []
    for f in G.feature[:TOPN]:
        j = fidx[f]
        rng = np.random.RandomState(SEED + j)  # seed 0 offset by column so repeats differ between features
        for r in range(NREP):
            Xp = Xs.copy(); Xp[:, j] = Xp[rng.permutation(len(Xp)), j]
            pp = model.predict_proba(Xp)[:, 1]
            ro, pr = ad.auroc_auprc(y, pp)
            prow.append(dict(TF=tf, feature=f, repeat=r, d_auroc=auroc - ro, d_auprc=auprc - pr))
    P = pd.DataFrame(prow)
    os.makedirs(os.path.join(OUT, "perm"), exist_ok=True)
    P.to_csv(os.path.join(OUT, "perm", "%s_perm.csv" % tf), index=False)
    os.makedirs(os.path.join(OUT, "metrics"), exist_ok=True)
    pd.DataFrame([dict(TF=tf, n_features=len(feats), n_test=len(y), n_pos=int(y.sum()), n_neg=int((y == 0).sum()),
                       model_source=source, sav_md5=md5(path), auroc_sav=auroc, auprc_sav=auprc, seconds=time.time() - t0)]
                 ).to_csv(os.path.join(OUT, "metrics", "%s_metrics.csv" % tf), index=False)
    print(tf, len(feats), len(y), "auROC %.4f auPRC %.4f" % (auroc, auprc), "%.0fs" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tfs", nargs="+", default=TFS)
    ap.add_argument("--n-jobs", type=int, default=2)
    ap.add_argument("--retrain", nargs="*", default=[], help="TFs whose .sav column order is unrecoverable (JUND)")
    a = ap.parse_args()
    for tf in a.tfs:
        run(tf, a.n_jobs, tf in a.retrain)
