"""b4_02_subtypes.py - feature importance, step 2: PWM similarity and subtypes per TF.
Similarity of two PWMs = the Fig. S6C definition (cross_species/motif_similarity.py): maximum over both strands and
offsets (overlap >= 5 columns) of the mean per-column Pearson correlation of the row-z-scored log-odds matrices.
monoPWM: file log-odds. diPWM (16 columns per row): P(xy) = exp(w)/sum, reshaped [row, first nt, second nt]; mono marginal =
first-nt marginal of rows 1..L, plus second-nt marginal of row L (L+1 positions); log(p/0.25) is z-scored like a monoPWM.
Subtypes: average-linkage hierarchical clustering on 1 - similarity of ALL features of the model (mono + di marginal), cut at similarity
0.7 / 0.8 / 0.9. Subtype 1 = the cluster with the largest summed Gini importance ("primary"), numbered by decreasing importance.
Outputs: subtypes/<TF>_subtypes.csv (feature, gini, subtype at 0.7/0.8/0.9), cache/<TF>_sim.npy, subtypes/<TF>_sim_summary.csv
usage: python b4_02_subtypes.py --tfs ... [--workers 3]
"""
import argparse
from multiprocessing import Pool
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from b4_common import *

THRS = (0.7, 0.8, 0.9)


def zrows(m):
    m = m - m.mean(axis=1, keepdims=True)
    sd = m.std(axis=1, keepdims=True)
    sd[sd == 0] = 1.0
    return m / sd


def pcc_max(a, b, min_overlap=5):
    best = -1.0
    for bb in (b, b[::-1, ::-1]):
        R = a @ bb.T / 4.0
        la, lb = R.shape
        for k in range(-(lb - min_overlap), la - min_overlap + 1):
            d = np.diagonal(R, offset=-k)
            if len(d) >= min_overlap:
                best = max(best, d.mean())
    return best


def mono_matrix(kind, w):
    if kind == "mono":
        return w
    L = len(w)
    P = np.exp(w - w.max(axis=1, keepdims=True)); P = (P / P.sum(axis=1, keepdims=True)).reshape(L, 4, 4)
    first = P.sum(axis=2)                       # L x 4
    last = P[-1].sum(axis=0)[None, :]           # second nt of row L
    p = np.vstack([first, last])
    return np.log(p / 0.25 + 1e-6)


_M = None
def _row(i):
    n = len(_M)
    return [(i, j, pcc_max(_M[i], _M[j])) for j in range(i + 1, n)]


def run(tf, workers):
    global _M
    G = pd.read_csv(os.path.join(OUT, "gini", "%s_gini.csv" % tf))
    mats = []
    for f in G.feature:
        _, w = read_pwm_file(tf, f)
        mats.append(zrows(mono_matrix(f.split("_")[0], w)))
    _M = mats
    n = len(mats)
    S = np.eye(n)
    with Pool(workers) as pool:
        for part in pool.imap_unordered(_row, range(n - 1), chunksize=4):
            for i, j, v in part:
                S[i, j] = S[j, i] = v
    np.save(os.path.join(CACHE, "%s_sim.npy" % tf), S)
    D = np.clip(1 - S, 0, None); np.fill_diagonal(D, 0)
    Z = linkage(squareform(D, checks=False), method="average") if n > 1 else None
    out = G[["TF", "feature", "kind", "dataset", "variant", "gini", "rank"]].copy()
    for thr in THRS:
        lab = fcluster(Z, t=1 - thr, criterion="distance") if n > 1 else np.array([1])
        tot = pd.Series(G.gini.values).groupby(lab).sum().sort_values(ascending=False)
        remap = {c: i + 1 for i, c in enumerate(tot.index)}
        out["subtype_%02d" % int(thr * 100)] = [remap[c] for c in lab]
    os.makedirs(os.path.join(OUT, "subtypes"), exist_ok=True)
    out.to_csv(os.path.join(OUT, "subtypes", "%s_subtypes.csv" % tf), index=False)
    ism = (G.kind == "mono").values
    iu = np.triu_indices(n, 1)
    mm = S[np.ix_(ism, ism)]; mu = np.triu_indices(ism.sum(), 1)
    top = int(np.argmax(G.gini.values))
    summ = dict(TF=tf, n_features=n, n_mono=int(ism.sum()), n_di=int((~ism).sum()),
                median_pairwise_all=float(np.median(S[iu])), median_pairwise_mono=float(np.median(mm[mu])) if ism.sum() > 1 else np.nan,
                min_pairwise_all=float(S[iu].min()), median_sim_to_top_feature=float(np.median(np.delete(S[top], top))),
                top_feature=G.feature[top])
    for thr in THRS:
        c = "subtype_%02d" % int(thr * 100)
        sh = out.groupby(c).gini.sum()
        summ["n_subtypes_%02d" % int(thr * 100)] = int(out[c].nunique())
        summ["n_subtypes_ge10pct_%02d" % int(thr * 100)] = int((sh >= 0.10).sum())
        summ["primary_share_%02d" % int(thr * 100)] = float(sh.max())
        summ["n_features_primary_%02d" % int(thr * 100)] = int((out[c] == 1).sum())
    pd.DataFrame([summ]).to_csv(os.path.join(OUT, "subtypes", "%s_sim_summary.csv" % tf), index=False)
    print(tf, n, {k: v for k, v in summ.items() if k.startswith("n_subtypes")}, flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tfs", nargs="+", default=TFS)
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    os.makedirs(CACHE, exist_ok=True)
    for tf in a.tfs:
        run(tf, a.workers)
