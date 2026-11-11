"""b2_compute.py - low-FPR metrics, RF (ArChIPelago) vs best single monoPWM, from saved per-sequence test scores.
Inputs: <TF>.npz in --scores (written by run_hm_mm_h5.py).  Writes only into --out.
Scenarios: HH = human-trained RF on human test (chr1,8,21); HM = human-trained RF on mouse test (chr1,8,19);
MM = mouse-trained RF on mouse test.  Baseline = best monoPWM chosen on the training set by auROC (model
'mono_byroc', used for ROC-type metrics) or by auPRC ('mono_byprc', PR-type metrics).
usage: python b2_compute.py --scores DIR --tables DIR --out DIR"""
import os, sys, glob, argparse
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import b2_lib as L
ap = argparse.ArgumentParser()
ap.add_argument("--scores", required=True)  # no defaults
ap.add_argument("--tables", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()

def chrom_map(fname):
    t = pd.read_csv(os.path.join(a.tables, fname), sep="\t", header=None, dtype=str, usecols=[1]).dropna()
    s = t[1]
    nm = s.str.split("$").str[0]
    ch = s.str.extract(r"chrom: (chr\w+);")[0]
    d = dict(zip(nm, ch))
    return d
CH = {"H": chrom_map("out_tab_HUMAN_10000_control.tab"),     # human test = file 'control'
      "M": chrom_map("out_tab_MOUSE_10000_train.tab")}  # mouse test = file 'train' (names are inverted on disk, see archi_data_v5)
TEST_CH = {"H": ["chr1", "chr8", "chr21"], "M": ["chr1", "chr8", "chr19"]}

rows, crows = [], []
for f in sorted(glob.glob(a.scores + "/*.npz")):
    tf = os.path.basename(f)[:-4]
    z = np.load(f, allow_pickle=True)
    keys = set(z.files)
    scen = {"HH": ("y_test_H", "names_test_H", "H", "rf_H_s%d", "pwmH_by%s_test_H"),
            "HM": ("y_test_M", "names_test_M", "M", "rf_M_s%d", "pwmH_by%s_test_M")}
    if "rfMM_M_s0" in keys:
        scen["MM"] = ("y_test_M", "names_test_M", "M", "rfMM_M_s%d", "pwmM_by%s_test_M")
    for sc, (ky, kn, sp, krf, kb) in scen.items():
        y = z[ky]; names = z[kn]
        chs = np.array([CH[sp].get(n.split("$")[0], "NA") for n in names])
        bad = int((chs == "NA").sum())
        if bad or set(np.unique(chs)) - set(TEST_CH[sp]):
            raise RuntimeError("%s %s: %d unmapped; chroms %s" % (tf, sc, bad, np.unique(chs)))
        models = {"RF_s0": z[krf % 0], "RF_s1": z[krf % 1], "mono_byroc": z[kb % "roc"], "mono_byprc": z[kb % "prc"]}
        for mn, s in models.items():
            m = L.all_metrics(y, s)
            m.update(TF=tf, scenario=sc, model=mn)
            rows.append(m)
            for c in TEST_CH[sp]:
                k = chs == c
                mc = L.all_metrics(y[k], s[k])
                crows.append({"TF": tf, "scenario": sc, "model": mn, "chrom": c, "n_pos": mc["n_pos"], "n_neg": mc["n_neg"],
                              "auroc": mc["auroc"], "auprc": mc["auprc"], "tpr_at_fpr0.01": mc["tpr_at_fpr0.01"],
                              "pauc_std_fpr0.01": mc["pauc_std_fpr0.01"]})
    print(tf, "done", flush=True)
wide = pd.DataFrame(rows)
idc = ["TF", "scenario", "model"]
wide = wide[idc + [c for c in wide.columns if c not in idc]]
wide.to_csv(os.path.join(a.out, "B2_metrics_wide.csv"), index=False)
long = wide.melt(id_vars=idc, var_name="metric", value_name="value")
long.to_csv(os.path.join(a.out, "B2_metrics_per_TF_long.csv"), index=False)
pd.DataFrame(crows).to_csv(os.path.join(a.out, "B2_per_chromosome.csv"), index=False)
print("wide", wide.shape, "chrom", len(crows))
