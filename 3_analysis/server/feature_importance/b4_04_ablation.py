"""b4_04_ablation.py - feature importance, step 4: retrain-without-subtype ablation.
Per TF: training set = human train chromosomes (positives <= 10,000, negatives <= 1,000,000, loader archi_data_v5.load_matrix, no extra
subsampling), test set = human chr1/8/21 (cache/<TF>_test.npz). Both matrices are standardised with the training scaler.
RF = RandomForestClassifier(max_depth=6, max_samples=0.8, n_estimators=100, random_state=seed) (hyperparameters of the saved model).
Conditions (subtypes at similarity --thr, default 0.8; 0.9 for TFs with a single subtype at 0.8, subtype 1 = primary, subtype 2 = top secondary, mono + di features alike):
  full            all features, seeds 0,1,2
  drop_secondary  all features of subtype 2 removed, seeds 0,1,2
  drop_random     the same NUMBER of randomly chosen primary-subtype features removed (capped at n_primary - 1), draw d = 0,1,2
                  (feature draw RandomState(d), RF seed d)
  drop_primary    all primary-subtype features removed, seeds 0,1,2 (positive control: the model must degrade when a used subtype goes)
Output: ablation/<TF>_ablation.csv (one row per fit: test auROC, auPRC-integral, n features kept)
usage: python b4_04_ablation.py --tfs E2F4 ... [--n-jobs 4]
"""
import argparse, time
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from b4_common import *

SEEDS = (0, 1, 2)


def fit_eval(Xtr, ytr, Xte, yte, cols, seed, n_jobs):
    m = RandomForestClassifier(n_jobs=n_jobs, random_state=seed, **ad.RF_PARAMS).fit(Xtr[:, cols], ytr)
    return ad.auroc_auprc(yte, m.predict_proba(Xte[:, cols])[:, 1])


def run(tf, n_jobs, thr):
    t0 = time.time()
    C = np.load(os.path.join(CACHE, "%s_test.npz" % tf), allow_pickle=True)
    feats = list(C["feats"]); yte = C["y"]
    Xtr, ytr, _ = ad.load_matrix(tf, "HUMAN", "HUMAN", "train", feats)
    sc_tr = StandardScaler().fit(Xtr)  # training scaler, applied to the test matrix
    Xte = sc_tr.transform(C["X"]).astype(np.float32)
    Xtr = sc_tr.transform(Xtr).astype(np.float32)
    S = pd.read_csv(os.path.join(OUT, "subtypes", "%s_subtypes.csv" % tf)).set_index("feature").loc[feats]
    st = S["subtype_%02d" % thr].values
    prim = np.flatnonzero(st == 1); sec = np.flatnonzero(st == 2)
    if len(sec) == 0:
        print(tf, "no secondary subtype at %.1f - skipped" % (thr / 100) + "", flush=True); return
    allc = np.arange(len(feats))
    rows = []
    def add(cond, seed, cols, **kw):
        ro, pr = fit_eval(Xtr, ytr, Xte, yte, cols, seed, n_jobs)
        rows.append(dict(TF=tf, sim_threshold=thr / 100, cond=cond, seed=seed, n_features=len(cols), n_removed=len(feats) - len(cols), auroc=ro, auprc=pr,
                         n_train_pos=int(ytr.sum()), n_train_neg=int((ytr == 0).sum()), n_prim=len(prim), n_sec=len(sec), **kw))
    for s in SEEDS:
        add("full", s, allc)
        add("drop_secondary", s, np.setdiff1d(allc, sec))
        add("drop_primary", s, np.setdiff1d(allc, prim))
    k = min(len(sec), len(prim) - 1)
    for d in SEEDS:
        rem = np.random.RandomState(d).choice(prim, k, replace=False) if k > 0 else np.array([], int)
        add("drop_random", d, np.setdiff1d(allc, rem), n_random_requested=len(sec))
    os.makedirs(os.path.join(OUT, "ablation"), exist_ok=True)
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "ablation", "%s_ablation.csv" % tf), index=False)
    print(tf, "done %.0fs" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tfs", nargs="+", required=True)
    ap.add_argument("--n-jobs", type=int, default=4)
    ap.add_argument("--thr", type=int, default=80, help="similarity cut x100 defining the subtypes (80 = 0.8)")
    a = ap.parse_args()
    for tf in a.tfs:
        run(tf, a.n_jobs, a.thr)
