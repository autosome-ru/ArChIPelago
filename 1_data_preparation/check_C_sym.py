"""
check_C_sym.py - independent check of the rebuilt positives (does not import build_positives_v5.py), --species HUMAN|MOUSE.
(1) every positive window (bed rows of <TFML>/{train,test}/<TF>_{train,test}_set_<SPECIES>.bed) is the +-500 bp window of
    a peak that ranks in the top 1,000 (tags >= 10; rank = 1 + number of peaks with a strictly larger -10log10 P, so
    ties at the cut-off count as inside) of at least one peak set of the TF in --peaks; key = chrom, summit, P (2 dp);
(2) union coverage of each window by the repeat track (bedtools coverage) is < 0.7.
usage: python check_C_sym.py --tfml DIR --peaks CSV --mask BED --tmp DIR --species HUMAN|MOUSE [--tfs A B]
(prints one line per TF and a pooled line)
"""
import os, sys, argparse, subprocess, tempfile
import numpy as np, pandas as pd

H = os.path.expanduser("~")
ap = argparse.ArgumentParser()
ap.add_argument("--tfml", required=True); ap.add_argument("--peaks", required=True); ap.add_argument("--mask", required=True)
ap.add_argument("--tfs", nargs="*"); ap.add_argument("--tmp", required=True)
ap.add_argument("--species", choices=["HUMAN", "MOUSE"], required=True)
a = ap.parse_args()
pk = pd.read_csv(a.peaks)
tfs = a.tfs or sorted(pk.TF.unique())
tot = dict(n=0, not_top=0, rep=0)
for tf in tfs:
    keys = set()
    for pid in sorted(set(pk.peak_id[pk.TF == tf])):
        d = pd.read_csv(f"{H}/macs/{pid}.interval", sep="\t")
        d = d[d["tags"] >= 10]
        p = d["-10*log10(pvalue)"].values
        sp = np.sort(p)
        rank = 1 + (len(sp) - np.searchsorted(sp, p, side="right"))
        d = d[rank <= 1000]
        keys |= set(zip("chr" + d["#CHROM"].astype(str), (d["START"] + d["summit"]).astype(int), np.round(d["-10*log10(pvalue)"].values, 2)))
    beds = []
    for s, t in (("train", "train"), ("test", "test")):
        b = pd.read_csv(f"{a.tfml}/{s}/{tf}_{t}_set_{a.species}.bed", sep="\t", index_col=0)
        beds.append(b)
    b = pd.concat(beds, ignore_index=True)
    pv = b["name"].str.split(":").str[0].astype(float).round(2)
    k = list(zip(b["chrom"], b["score"].astype(int), pv))
    not_top = sum(1 for x in k if x not in keys)
    win = b[["chrom", "start", "end"]]
    fd, fn = tempfile.mkstemp(dir=a.tmp, suffix=".bed"); os.close(fd)
    win.to_csv(fn, sep="\t", header=False, index=False)
    cov = subprocess.run(f"bedtools coverage -a {fn} -b {a.mask} | awk '$NF>=0.7' | wc -l", shell=True,
                         capture_output=True, text=True, check=True).stdout.strip()
    os.unlink(fn)
    rep = int(cov)
    print(f"{tf}\tpositives={len(k)}\tnot_in_top1000_of_listed_sets={not_top}\trepeat_cov>=0.7={rep}", flush=True)
    tot["n"] += len(k); tot["not_top"] += not_top; tot["rep"] += rep
print(f"POOLED\tpositives={tot['n']}\tnot_in_top1000_of_listed_sets={tot['not_top']}\trepeat_cov>=0.7={tot['rep']}", flush=True)
