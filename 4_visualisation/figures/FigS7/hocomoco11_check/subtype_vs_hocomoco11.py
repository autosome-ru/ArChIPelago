"""Fig. S7 subtype labels from HOCOMOCO v11 (subtypes come from the PWM collection, not from logos).
Each Fig. S7 matrix (matrices/<TF>_<kind>_<k>.txt, position probabilities; di collapsed to mono by the S7 builder)
is matched to every HOCOMOCO v11 HUMAN mono model (full collection, pcm -> probabilities) with pcc_max of the
S7 builder (max over offsets and strands of mean column Pearson r of z-scored rows, overlap >= 5).
Output: subtype_vs_hocomoco11.csv (one row per Fig. S7 matrix): best own-TF v11 model (its subtype index = rank field),
r, and the best model over all 769 human v11 models (shows co-factor motifs)."""
import os, sys, hashlib, datetime, subprocess
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
S7 = os.path.dirname(HERE)
R2 = os.path.dirname(os.path.dirname(os.path.dirname(S7)))
sys.path.insert(0, os.path.join(R2, "3_analysis", "package", "analysis", "motif_subtypes"))
from build_Figure_S7_source_data import pcc_max, zrows  # noqa: E402

PCM = os.path.join(R2, "results/datasets/hocomoco11/HOCOMOCOv11_full_pcms_HUMAN_mono.txt")
SRC = os.path.join(S7, "Figure_S7_source_data.csv")
MIN_OVERLAP = 5


def read_pcms(path):
    mats, name, rows = {}, None, []
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name:
                mats[name] = np.array(rows)
            name, rows = line[1:].split()[0], []
        else:
            rows.append([float(v) for v in line.split()])
    mats[name] = np.array(rows)
    out = {}
    for k, m in mats.items():
        if m.shape[1] != 4:
            raise SystemExit("pcm %s has %d columns" % (k, m.shape[1]))
        out[k] = m / m.sum(axis=1, keepdims=True)
    return out


def read_matrix(path):
    t = pd.read_csv(path, sep="\t", index_col=0)
    if list(t.index) != ["A", "C", "G", "T"]:
        raise SystemExit("unexpected rows in %s" % path)
    return t.values.T


def md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()


def main():
    v11 = read_pcms(PCM)
    z11 = {k: zrows(m) for k, m in v11.items()}
    src = pd.read_csv(SRC)
    tfs = src.TF.unique().tolist()   # all S7 candidate TFs
    rows = []
    for _, r in src[src.TF.isin(tfs)].iterrows():
        za = zrows(read_matrix(os.path.join(S7, r.matrix_file)))
        own = {k: pcc_max(za, z, MIN_OVERLAP) for k, z in z11.items() if k.startswith(r.TF + "_HUMAN.")}
        if not own:
            raise SystemExit("no HOCOMOCO v11 model for %s" % r.TF)
        allm = {k: pcc_max(za, z, MIN_OVERLAP)[0] for k, z in z11.items()}
        bo = max(own, key=lambda k: own[k][0])
        ba = max(allm, key=allm.get)
        rows.append(dict(TF=r.TF, panel=r.panel, pwm_class=r.pwm_class, pwm_name=r.pwm_name, cell_line=r.cell_line,
                         consensus=r.consensus, n_v11_subtypes=len(own),
                         **{"r_" + k.split(".")[-2] + "." + k.split(".")[-1]: round(v[0], 3) for k, v in sorted(own.items())},
                         best_own_v11=bo, best_own_subtype=int(bo.split(".")[-2]), r_best_own=round(own[bo][0], 3),
                         best_any_v11=ba, r_best_any=round(allm[ba], 3)))
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(HERE, "subtype_vs_hocomoco11.csv"), index=False)
    git = subprocess.run(["git", "-C", R2, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    with open(os.path.join(HERE, "_provenance.yaml"), "w") as fh:
        fh.write("script: subtype_vs_hocomoco11.py\ndate: %s\ngit: %s\nmin_overlap: %d\ninputs:\n" % (
            datetime.datetime.now().isoformat(timespec="seconds"), git, MIN_OVERLAP))
        for p in (PCM, SRC):
            fh.write("  - {path: %s, md5: %s}\n" % (os.path.relpath(p, R2), md5(p)))
        fh.write("source: server ~/hocomoco11/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_pcms_HUMAN_mono.txt\n")
    print(out.drop(columns=["pwm_name"]).to_string())


if __name__ == "__main__":
    main()
