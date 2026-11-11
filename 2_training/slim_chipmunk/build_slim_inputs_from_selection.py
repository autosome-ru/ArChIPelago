"""build_slim_inputs_from_selection.py (server) - the Slim inputs of build_slim_inputs_v5.py rebuilt from a shipped selection, without the loader.
Inputs (shipped from server, md5 checked): OUT_TAB = out_tab_HUMAN_10000_train.tab, SLIM =
all_mfa_file_HUMAN_SLIM_10000_train_no_NF_no_N.fasta, SEL = <TF>/selection.npz (row positions into OUT_TAB and y, written by
build_slim_inputs_v5.py --export on server). Same selection order, headers and 60-base lines as on server.
usage: python build_slim_inputs_from_selection.py OUT_TAB SLIM SEL OUTDIR   (refuses to overwrite)"""
import os, sys
import numpy as np
tab, slim, sel, out = sys.argv[1:5]
os.makedirs(out, exist_ok=True)
f1, f0 = os.path.join(out, "train_H_1_SLIM.fasta"), os.path.join(out, "train_H_0_SLIM.fasta")
for f in (f1, f0):
    if os.path.exists(f):
        raise SystemExit("refusing to overwrite " + f)
z = np.load(sel)
rows, y = z["rows"], z["y"]
names = []
with open(tab) as fh:
    for line in fh:
        names.append(line.split("\t")[1].split("$")[0])
if len(names) != int(z["n_table"]):
    raise SystemExit("out_tab has %d rows, selection was made on %d" % (len(names), int(z["n_table"])))
want = np.array(names, dtype=object)[rows]
need = set(want)
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
def write(path, ns, header):
    with open(path, "w") as fo:
        for n in ns:
            s = "".join(seq[n])
            fo.write(">" + header(n) + "\n" + "".join(s[i:i + 60] + "\n" for i in range(0, len(s), 60)))
write(f1, want[y == 1], lambda n: hdr[n])
write(f0, want[y == 0], lambda n: "> peak: 0; signal: 0")
print(os.path.basename(os.path.dirname(sel)), "positives", int((y == 1).sum()), "background", int((y == 0).sum()), "->", out, flush=True)
