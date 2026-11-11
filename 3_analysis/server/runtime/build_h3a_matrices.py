"""Write the human training matrices (negatives v3, monoPWM features) for the runtime benchmark,
in the same layout as ~/archi_release/TF-ML/outputdir/<TF>/result_H_train_full.csv read by rf_bench_h3a.py:
index, columns 0-3 = position / name / threshold (NaN, not used) / TF label, then k monoPWM best-hit scores, 'ind'.
Usage: build_h3a_matrices.py <out_dir> TF [TF ...]"""
import os, sys, json, hashlib
import numpy as np, pandas as pd
sys.path.insert(0, os.path.expanduser("~/archi_work/human_sym_h5/code"))
import archi_data_v5 as ad
assert ad.NEG_VERSION == "v3"
out_dir = sys.argv[1]; os.makedirs(out_dir, exist_ok=True)
prov = {}
for tf in sys.argv[2:]:
    feats = ad.feature_list(tf, "HUMAN", "mono")
    X, y, meta = ad.load_matrix(tf, "HUMAN", "HUMAN", "train", feats)
    # scores carry 2 decimals in the score tables; write them with 2 decimals like the published matrices
    dev = float(np.abs(np.round(X.astype(np.float64), 2) - X.astype(np.float64)).max())
    assert dev < 1e-4, (tf, dev)
    df = pd.DataFrame({"0": np.arange(len(y)), "1": meta["name"].values, "2": np.nan, "3": meta["tf"].values})
    F = pd.DataFrame(np.round(X.astype(np.float64), 2), columns=[str(4 + i) for i in range(X.shape[1])])
    df = pd.concat([df, F], axis=1); df["ind"] = y.astype(int)
    p = os.path.join(out_dir, "%s_H3a_human_train_mono.csv" % tf)
    if os.path.exists(p): raise FileExistsError(p)
    df.to_csv(p, float_format="%.2f")
    prov[tf] = dict(path=p, n_rows=int(len(y)), n_pos=int(y.sum()), n_neg=int((y == 0).sum()), k_mono=len(feats),
                    md5=hashlib.md5(open(p, "rb").read()).hexdigest(), max_round_dev=dev,
                    chroms=sorted(meta["chrom"].dropna().unique().tolist()))
    print(json.dumps({k: v for k, v in prov[tf].items() if k != "chroms"}), "chroms:", ",".join(prov[tf]["chroms"]))
json.dump(prov, open(os.path.join(out_dir, "matrices_provenance.json"), "x"), indent=1)
