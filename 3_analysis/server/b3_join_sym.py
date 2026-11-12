"""b3_join_sym.py - peak widths of the source peak sets and the (chrom, center, signal) join of the scored test positives
to their source peaks, --species H|M.
H: the 361 human peak sets and out_tab_HUMAN_10000_control.tab (chr1/8/21), names_test_H of the score npz.
M: the 276 mouse peak sets and out_tab_MOUSE_10000_train.tab (chr1/8/19), names of the score npz.
Datasets = <sp>_sym_peak_ids.csv (TF, peak_id); interval columns are read by NAME (7- and 8-column files);
--signal-offset adds a constant to the table signal (control run: any offset > 0.0051 must drop the match rate to ~0).
usage: python b3_join_sym.py OUTDIR --species H|M [--signal-offset 0] [--suffix _mouse]   (refuses to overwrite)"""
import os, re, sys, argparse
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("out"); ap.add_argument("--signal-offset", type=float, default=0.0); ap.add_argument("--suffix", default="_mouse")
ap.add_argument("--species", choices=["H", "M"], required=True)
a = ap.parse_args()
sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "..")); import archi_env as AE
MACS = AE.macs()
MS = AE.work("mouse_sym")
SP = a.species
ROOT = {"H": AE.work("human_sym"), "M": MS}[SP]
PEAKS = {"H": f"{ROOT}/inputs/human_sym_peak_ids.csv", "M": f"{MS}/inputs/mouse_sym_peak_ids.csv"}[SP]
TABF = {"H": f"{ROOT}/TF-ML/out_tab_HUMAN_10000_control.tab", "M": f"{MS}/TF-ML/out_tab_MOUSE_10000_train.tab"}[SP]
SCD = AE.work("human_sym", "H3a", "scores")
HALF, TOPK = 150, 1000
outs = {k: f"{a.out}/b3_{k}{a.suffix}.csv" for k in ("peak_widths_per_dataset", "join_test_per_tf", "join_test_per_sequence")}
for p in outs.values():
    if os.path.exists(p):
        raise SystemExit("refusing to overwrite " + p)
t1 = pd.read_csv(PEAKS)
ds_rows, store = [], {}
for r in t1.itertuples(index=False):
    d = pd.read_csv(f"{MACS}/{r.peak_id}.interval", sep="\t")
    d = d.rename(columns={"#CHROM": "chrom", "START": "start", "END": "end", "-10*log10(pvalue)": "p"})
    w = (d.end - d.start).values
    s = d.start.values + d.summit.values
    ov = np.minimum(d.end.values, s + HALF) - np.maximum(d.start.values, s - HALF)
    cover = np.clip(ov, 0, None) / w
    o = np.argsort(-d.p.values, kind="mergesort")
    wt, wr = (w[o[:TOPK]], w[o[TOPK:]]) if len(d) > TOPK else (None, None)
    ds_rows.append(dict(TF=r.TF, species=SP, peak_set=r.peak_id, n_peaks=len(d), n_peaks_table=np.nan,
                        width_median=np.median(w), width_q25=np.percentile(w, 25), width_q75=np.percentile(w, 75),
                        frac_width_le300=float((w <= 2 * HALF).mean()), frac_window_covers_ge50=float((cover >= 0.5).mean()),
                        median_cover=float(np.median(cover)),
                        width_median_top1000=np.median(wt) if wt is not None else np.nan,
                        width_median_rest=np.median(wr) if wr is not None and len(wr) else np.nan,
                        cover_ge50_top1000=float((cover[o[:TOPK]] >= 0.5).mean()) if wt is not None else np.nan,
                        cover_ge50_rest=float((cover[o[TOPK:]] >= 0.5).mean()) if wt is not None and len(wr) else np.nan))
    store[(r.TF, SP, r.peak_id)] = (d.chrom.astype(str).values, s, np.round(d.p.values, 2), w)
pd.DataFrame(ds_rows).to_csv(outs["peak_widths_per_dataset"], index=False)
print("datasets", len(ds_rows), flush=True)


def read_tab(path):
    dd = {}
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            nm, hdr = p[1].split("$", 1)
            m = re.search(r"chrom: chr(\w+); center: ([\d.]+); peak: \d+; signal: ([\d.]+)", hdr)
            dd[nm] = (m.group(1), int(float(m.group(2))), round(float(m.group(3)), 2))
    return dd


TAB = read_tab(TABF)
jrows, perseq = [], []
for tf in sorted({k[0] for k in store}):
    z = np.load(f"{SCD}/{tf}.npz", allow_pickle=True)  # own npz written by run_hm_mm_h5.py
    names = z["names_test_" + SP][z["y_test_" + SP] == 1]
    idx = {}
    for (t, s_, pid), (ch, sm, p, w) in store.items():
        if t == tf:
            for c, m, pp, ww in zip(ch, sm, p, w):
                idx.setdefault((c, int(m)), []).append((pp, ww, pid))
    n0 = n1 = nm = 0
    for nmz in names:
        c, cen, sg = TAB[nmz]
        sg = sg + a.signal_offset
        cand = [(ww, pid) for pp, ww, pid in idx.get((c, cen), []) if abs(pp - sg) < 0.0051]
        if not cand:
            n0 += 1; perseq.append((tf, SP, nmz, c, cen, sg, 0, np.nan)); continue
        ws = {ww for ww, _ in cand}
        if len(ws) == 1:
            n1 += 1; perseq.append((tf, SP, nmz, c, cen, sg, len(cand), ws.pop()))
        else:
            nm += 1; perseq.append((tf, SP, nmz, c, cen, sg, len(cand), np.nan))
    jrows.append(dict(TF=tf, species=SP, n_test_pos=len(names), n_matched_unique_width=n1, n_matched_ambiguous_width=nm,
                      n_unmatched=n0, match_rate=(n1 + nm) / len(names)))
    print(tf, jrows[-1], flush=True)
j = pd.DataFrame(jrows)
j.to_csv(outs["join_test_per_tf"], index=False)
pd.DataFrame(perseq, columns=["TF", "species", "name", "chrom", "center", "signal", "n_candidates", "width"]).to_csv(
    outs["join_test_per_sequence"], index=False)
print("POOLED", "n_test_pos", int(j.n_test_pos.sum()), "matched", int((j.n_matched_unique_width + j.n_matched_ambiguous_width).sum()),
      "ambiguous_width", int(j.n_matched_ambiguous_width.sum()), "unmatched", int(j.n_unmatched.sum()), flush=True)
