"""
gc_check_v3_h5.py - GC content of the 300 bp positives vs their BiasAway negatives (v3), per TF and split.
Pools: the human and mouse data builds of this folder (neg-id file suffix "_v3").
Rows exactly as the loader selects them (archi_data_v5.select_rows): positives = rows of <TF>_<SP> (.sample(10000,
random_state=0) if more), negatives = rows whose name is in the neg-id file (.sample(1e6, random_state=0) if more).
GC = (G+C)/(A+C+G+T) over the 300 bp sequence (case-insensitive). KS = scipy.stats.ks_2samp(pos, neg).
File split names as on disk; chromosomes: HUMAN train = chr2-7, 9-20, 22, control = chr1/8/21;
MOUSE train = chr1/8/19 (the mouse test set), control = chr2-7, 9-18.
usage: python gc_check_v3_h5.py --tfs A B ... --out OUT.csv   (refuses to overwrite)
"""
import os, sys, argparse
import numpy as np, pandas as pd
from scipy.stats import ks_2samp

H = os.path.expanduser("~")
POOLS = {"human_sym_h5_v3": ("HUMAN", f"{H}/archi_work/human_sym_h5/TF-ML", "_v3"),  # v3 files only
         "mouse_sym_h5_v3": ("MOUSE", f"{H}/archi_work/mouse_sym_h5/TF-ML", "_v3")}
CHROMS = {("HUMAN", "train"): "chr2-7,9-20,22", ("HUMAN", "control"): "chr1/8/21",  # split labels
          ("MOUSE", "train"): "chr1/8/19", ("MOUSE", "control"): "chr2-7,9-18"}
POS_CLASS, NEG_CLASS = 10000, 1000000


def gc_of(seq):
    s = seq.upper()
    acgt = sum(s.count(c) for c in "ACGT")
    return (s.count("G") + s.count("C")) / acgt if acgt else np.nan


assert gc_of("GCGCAT") == 4 / 6 and gc_of("atatat") == 0.0          # two-line example of the GC function


def read_gc(fa):
    out, name = {}, None
    with open(fa) as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].strip()
            else:
                out[name] = gc_of(line.strip())
    return out


ap = argparse.ArgumentParser()
ap.add_argument("--tfs", nargs="+", required=True); ap.add_argument("--out", required=True)
a = ap.parse_args()
if os.path.exists(a.out):
    raise SystemExit("refusing to overwrite " + a.out)
rows = []
for pool, (sp, root, sfx) in POOLS.items():
    for split in ("train", "control"):
        gt = pd.read_csv(f"{root}/out_tab_{sp}_10000_{split}.tab", sep="\t", header=None, dtype={1: str, 3: str})
        gt = gt.dropna(subset=[1])
        gt["name"] = gt[1].str.split("$").str[0]
        gc = read_gc(f"{root}/all_mfa_file_{sp}_10000_{split}_no_NF.fasta")
        if set(gt["name"]) - set(gc):
            raise SystemExit(f"{pool} {split}: table names missing from fasta")
        gt["gc"] = gt["name"].map(gc)
        mode, flag = ("train", f"result_{sp[0]}_train_full") if split == "train" else ("test", f"result_{sp[0]}_control_full")
        for tf in a.tfs:
            pos = gt[gt[3] == f"{tf}_{sp}"]
            if pos.shape[0] > POS_CLASS:
                pos = pos.sample(n=POS_CLASS, replace=False, random_state=0)
            ids = set(pd.read_csv(f"{root}/outputdir/{tf}/{mode}_neg_id_out_{sp}_{flag}{sfx}.csv", header=None)[0])
            neg = gt[gt["name"].isin(ids)]
            if neg.shape[0] > NEG_CLASS:
                neg = neg.sample(n=NEG_CLASS, replace=False, random_state=0)
            ks = ks_2samp(pos["gc"].values, neg["gc"].values)
            rows.append(dict(pool=pool, species=sp, TF=tf, file_split=split, chromosomes=CHROMS[(sp, split)],
                             n_pos=pos.shape[0], n_neg=neg.shape[0], median_gc_pos=round(float(pos["gc"].median()), 4),
                             median_gc_neg=round(float(neg["gc"].median()), 4),
                             median_gc_diff=round(float(pos["gc"].median() - neg["gc"].median()), 4),
                             ks_D=round(float(ks.statistic), 4), ks_p=float(ks.pvalue)))
            print(rows[-1], flush=True)
df = pd.DataFrame(rows)
df.to_csv(a.out, index=False)
s = df.groupby(["pool", "file_split"]).agg(n_TF=("TF", "size"), median_ks_D=("ks_D", "median"), max_ks_D=("ks_D", "max"),
                                           median_abs_gc_diff=("median_gc_diff", lambda x: float(np.median(np.abs(x)))))
print(s.to_string())
print("GC_DONE")
