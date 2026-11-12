"""
make_negatives_v3.py - GC-matched negatives of one TF on one split (negatives v3), --species HUMAN|MOUSE.

Follows group_selection_GC of the release notebook (notebook_code/RELEASE_SLIM_MONO_..._sav_model_generation.py) for the "_full" flags that the loader reads, plus one rule before GC matching: every candidate negative whose
300 bp window overlaps (>= 1 bp, same chromosome) any 300 bp window of the TF's own positives of the same split (all rows
<TF>_<SP> of the split table) is dropped. 300 bp window = [center - 150, center + 150) with center = abs. summit from
the out_tab header (the 300 bp sequence is seq[350:650] of the +-500 bp window), so two windows overlap iff
|center1 - center2| < 300.
  positives = rows of the pooled out_tab with column 3 == <TF>_<SP> (.sample(n=10000, random_state=0) if more);
  pool      = rows whose family code (column 2) differs from the TF's family (TFClass code j.split("{")[1][:3] of
              core_annotation_HUMAN_mono.csv), or every other TF's row when the TF has no family code;
  seqtk subseq from the pooled _no_NF fasta, then biasaway g -n 100 -e 1 (-r/-f/-b/-p as in the notebook), then
  awk 'sub(/^>/, "")' -> <mode>_neg_id_out_<SP>_<flag>_v3.csv
HUMAN: names <TF>_HUMAN, flags result_H_{train,control}_full. MOUSE: names <TF>_MOUSE, flags result_M_{train,control}_full.
Split "train" -> mode train; split "control" -> mode test. New files only: working dir <outdir>/<TF>/neg_v3_work_<split>/,
final list <mode>_neg_id_out_<SP>_<flag>_v3.csv (hard link in <outdir>/<TF>/), side file neg_v3_stats_<split>.tsv.
The neg-id list holds the names only (seqtk matches names only) and the BiasAway stderr log is per flag.
usage: python make_negatives_v3.py --tfml TFMLDIR --outdir OUTPUTDIR --tf SRF --split train --species HUMAN [--biasaway PATH]
"""
import os, sys, argparse, subprocess
from collections import defaultdict
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")); import archi_env as AE
FAMILY = AE.family("core_annotation_HUMAN_mono.csv")
SEQTK = AE.var("ARCHI_SEQTK", "seqtk executable")
POS_CLASS, NEG_CLASS = 10000, 1000000

ap = argparse.ArgumentParser()
ap.add_argument("--tfml", required=True); ap.add_argument("--outdir", required=True)
ap.add_argument("--tf", required=True); ap.add_argument("--split", choices=["train", "control"], required=True)
ap.add_argument("--biasaway", default=None, help="default: $ARCHI_BIASAWAY")
ap.add_argument("--species", choices=["HUMAN", "MOUSE"], required=True)
a = ap.parse_args()
a.biasaway = a.biasaway or AE.var("ARCHI_BIASAWAY", "biasaway executable")
SP = a.species
mode, flag = ("train", f"result_{SP[0]}_train_full") if a.split == "train" else ("test", f"result_{SP[0]}_control_full")
TF_name = a.tf + "_" + SP
a.tfml = os.path.abspath(a.tfml)
tfdir = os.path.abspath(os.path.join(a.outdir, a.tf))
d = os.path.join(tfdir, f"neg_v3_work_{a.split}")   # working files; the final id list is hard-linked into tfdir
os.makedirs(d, exist_ok=True)
final = os.path.join(tfdir, f"{mode}_neg_id_out_{SP}_{flag}_v3.csv")
stats = os.path.join(tfdir, f"neg_v3_stats_{a.split}.tsv")
for p in (final, stats):
    if os.path.exists(p):
        raise SystemExit("refusing to overwrite " + p)
outs = [f"{mode}_pos_id_{SP}.csv", f"{mode}_neg_id_{SP}.csv", f"{mode}_pos_id_{SP}.fasta", f"{mode}_neg_id_{SP}.fasta",
        f"out_{SP}_{flag}.fasta", f"{mode}_neg_id_out_{SP}_{flag}.csv", f"BiasAway_results_{SP}_{flag}"]
for o in outs:
    if os.path.exists(os.path.join(d, o)):
        raise SystemExit("refusing to overwrite " + os.path.join(d, o))

fam = pd.read_csv(FAMILY, sep='\t')
TF_family_dict = defaultdict(None)
for i, j in zip(fam["Model"], fam["TF family"]):
    if type(j) == str:
        TF_family_dict[i.split('_')[0]] = j.split('{')[1][:3]

fasta = os.path.join(a.tfml, f"all_mfa_file_{SP}_10000_{a.split}_no_NF.fasta")
df = pd.read_csv(os.path.join(a.tfml, f"out_tab_{SP}_10000_{a.split}.tab"), header=None, sep='\t')
hdr = df[1].str.split("$").str[1]
df["chrom"] = hdr.str.extract(r"chrom: (\w+);")[0]
df["center"] = hdr.str.extract(r"center: ([\d.]+);")[0].astype(float)
if df["chrom"].isna().any() or df["center"].isna().any():
    raise SystemExit("unparsed out_tab headers")
df[1] = [x[0] for x in df[1].str.split("$")]
df_group = df.dropna(subset=[1]).copy()
columns_name = df_group.columns[3]
train_pos = df_group[[x == TF_name for x in df_group[columns_name]]]
df1 = train_pos.sample(n=POS_CLASS, replace=False, random_state=0) if train_pos.shape[0] > POS_CLASS else train_pos
TF_mask = TF_family_dict.get(a.tf, None)
if TF_mask is not None:
    train_neg = df_group[df_group[2] != float(TF_mask)]
else:
    train_neg = df_group[[x != TF_name for x in df_group[columns_name]]]
# v3 rule: drop candidates whose 300 bp window overlaps one of the TF's own positive 300 bp windows (same split)
n_cand = train_neg.shape[0]
own = {c: np.sort(g["center"].values) for c, g in train_pos.groupby("chrom")}
def overlaps(chrom, cen):
    s = own.get(chrom)
    if s is None:
        return False
    i = np.searchsorted(s, cen)
    return (i < len(s) and s[i] - cen < 300) or (i > 0 and cen - s[i - 1] < 300)
ov = np.array([overlaps(c, x) for c, x in zip(train_neg["chrom"].values, train_neg["center"].values)], dtype=bool)
train_neg = train_neg[~ov]
assert overlaps("chrZ", 0) is False and (len(own) == 0 or overlaps(next(iter(own)), own[next(iter(own))][0] + 299))
os.chdir(d)
df1[1].to_csv(f"{mode}_pos_id_{SP}.csv", sep="\n", index=False)
train_neg[1].to_csv(f"{mode}_neg_id_{SP}.csv", sep="\n", index=False)
cmds = [f"{SEQTK} subseq {fasta} {mode}_pos_id_{SP}.csv > {mode}_pos_id_{SP}.fasta",
        f"{SEQTK} subseq {fasta} {mode}_neg_id_{SP}.csv > {mode}_neg_id_{SP}.fasta",
        f"{a.biasaway} g -r BiasAway_results_{SP}_{flag} -f {mode}_pos_id_{SP}.fasta -b {mode}_neg_id_{SP}.fasta "
        f"-n {NEG_CLASS // POS_CLASS} -p BiasAway_results_{SP}_{mode} -e 1 > out_{SP}_{flag}.fasta 2> BiasAway_log_{SP}_{flag}.txt",
        "awk 'sub(/^>/, \"\")' < out_" + SP + "_" + flag + f".fasta > {mode}_neg_id_out_{SP}_{flag}.csv"]
for c in cmds:
    subprocess.run(c, shell=True, check=True, executable="/bin/bash")
os.link(f"{mode}_neg_id_out_{SP}_{flag}.csv", final)
n_neg = sum(1 for _ in open(final))
with open(stats, "x") as fh:
    fh.write("TF\tsplit\tn_pos\tclass\tcandidates_after_class\tremoved_overlap_own\tcandidates_v3\tneg_out\tneg_per_pos\n")
    fh.write(f"{a.tf}\t{a.split}\t{df1.shape[0]}\t{TF_mask}\t{n_cand}\t{int(ov.sum())}\t{train_neg.shape[0]}\t{n_neg}\t{n_neg / df1.shape[0]:.3f}\n")
print(f"{a.tf}\t{a.split}\tpos={df1.shape[0]}\tclass={TF_mask}\tcandidates={n_cand}\tremoved_overlap={int(ov.sum())}\tpool={train_neg.shape[0]}\tneg_out={n_neg}", flush=True)
