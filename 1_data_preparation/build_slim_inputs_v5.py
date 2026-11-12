"""build_slim_inputs_v5.py - Slim training inputs of one TF, built standalone from the RF loader archi_data_v5.
  train_H_1_SLIM.fasta : 1,000 bp SLIM sequences of the selected training positives (result_H_train ind == 1),
                         header = the '$' part of the SLIM fasta name without '>' ('chrom: ..; center: ..; peak: 501; signal: ..')
  train_H_0_SLIM.fasta : 1,000 bp SLIM sequences of the selected training negatives (ind == 0), header '> peak: 0; signal: 0'
  Selection = loader select_rows (group_selection_GC replica: positives capped at 10,000, rs 0; negatives = neg-id file).
  Biopython SeqIO fasta layout: 60 bases per line.
--export also writes the selection (row positions into out_tab_HUMAN_10000_train.tab and y) as <OUTDIR>/selection.npz, so
that the same FASTAs can be rebuilt from the pooled SLIM fasta without the loader (2_training/slim_chipmunk/build_slim_inputs_from_selection.py).
usage: python build_slim_inputs_v5.py LOADER_DIR TF OUTDIR [--export]   (LOADER_DIR holds archi_data_v5.py; refuses to overwrite)"""
import os, sys
loader, tf, out = sys.argv[1:4]
EXPORT = "--export" in sys.argv
sys.path.insert(0, os.path.expanduser(loader))
import numpy as np
import archi_data_v5 as ad
assert ad.NEG_VERSION == "v3"

os.makedirs(out, exist_ok=True)
f1, f0 = os.path.join(out, "train_H_1_SLIM.fasta"), os.path.join(out, "train_H_0_SLIM.fasta")
for f in (f1, f0):
    if os.path.exists(f):
        raise SystemExit("refusing to overwrite " + f)
gt = ad.global_table("HUMAN", "train")
rows, y = ad.select_rows(gt, tf, "HUMAN", "train")
if EXPORT:
    p = os.path.join(out, "selection.npz")
    if os.path.exists(p):
        raise SystemExit("refusing to overwrite " + p)
    np.savez_compressed(p, rows=rows.astype(np.int64), y=y.astype(np.int8), n_table=np.int64(gt.shape[0]))
want = gt["name"].to_numpy()[rows]
need = set(want)
slim = os.path.join(ad._rel("HUMAN") if hasattr(ad, "_rel") else ad.REL, "all_mfa_file_HUMAN_SLIM_10000_train_no_NF_no_N.fasta")
seq, hdr = {}, {}
with open(slim) as fh:
    name = None
    for line in fh:
        if line.startswith(">"):
            n, h = line[1:].rstrip("\n").split("$", 1)
            name = n if n in need else None
            if name is not None:
                hdr[name] = h[1:]
                seq[name] = []
        elif name is not None:
            seq[name].append(line.strip())
missing = need - set(seq)
if missing:
    raise SystemExit("%d selected names absent from %s" % (len(missing), slim))


def write(path, names, header):
    with open(path, "w") as fo:
        for n in names:
            s = "".join(seq[n])
            fo.write(">" + header(n) + "\n" + "".join(s[i:i + 60] + "\n" for i in range(0, len(s), 60)))


write(f1, want[y == 1], lambda n: hdr[n])
write(f0, want[y == 0], lambda n: "> peak: 0; signal: 0")
print(tf, "positives", int((y == 1).sum()), "background", int((y == 0).sum()), "->", out, flush=True)
