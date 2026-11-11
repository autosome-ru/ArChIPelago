"""check_neg_v3_h5.py - checks of the negative-ID files (v3).
Per species x TF x split: n, overlapping own positives (300 bp windows, >= 1 bp, same chromosome, i.e.
|center - own center| < 300), negatives that are the TF itself, negatives in the TF's class (TFClass code of
core_annotation_HUMAN_mono.csv, column 2 of the out_tab), missing from the split table, duplicates.
Also aggregates neg_v3_stats_<split>.tsv. Writes <out>/neg_v3_checks.csv and neg_v3_stats_all.csv (refuse to overwrite).
usage: python check_neg_v3_h5.py OUTDIR"""
import os, sys
import numpy as np, pandas as pd
H = os.path.expanduser("~"); R = f"{H}/archi_work"
OUT = sys.argv[1]
for f in ("neg_v3_checks.csv", "neg_v3_stats_all.csv"):
    if os.path.exists(f"{OUT}/{f}"): raise SystemExit("refusing to overwrite " + f)
fam = pd.read_csv(f"{H}/TF-ML_last_v_22/core_annotation_HUMAN_mono.csv", sep="\t")
F = {i.split("_")[0]: float(j.split("{")[1][:3]) for i, j in zip(fam["Model"], fam["TF family"]) if type(j) == str}
rows, st = [], []
for SP, root in (("HUMAN", f"{R}/human_sym_h5/TF-ML"), ("MOUSE", f"{R}/mouse_sym_h5/TF-ML")):
    for split in ("train", "control"):
        mode, flag = ("train", f"result_{SP[0]}_train_full") if split == "train" else ("test", f"result_{SP[0]}_control_full")
        t = pd.read_csv(f"{root}/out_tab_{SP}_10000_{split}.tab", sep="\t", header=None)
        h = t[1].str.split("$").str[1]
        t["id"] = t[1].str.split("$").str[0]; t["chrom"] = h.str.extract(r"chrom: (\w+);")[0]
        t["center"] = h.str.extract(r"center: ([\d.]+);")[0].astype(float)
        info = t.set_index("id")
        for tf in sorted(os.listdir(f"{root}/outputdir")):
            own = {c: np.sort(g.center.values) for c, g in t[t[3] == f"{tf}_{SP}"].groupby("chrom")}
            for ver, fn in (("v3", f"{mode}_neg_id_out_{SP}_{flag}_v3.csv"),):  # v3 files only
                ids = [x.strip() for x in open(f"{root}/outputdir/{tf}/{fn}")]
                ok = [x for x in ids if x in info.index]
                sub = info.loc[ok]
                ovl = 0
                for c, g in sub.groupby("chrom"):
                    s = own.get(c)
                    if s is None: continue
                    x = g.center.values; i = np.searchsorted(s, x)
                    right = np.where(i < len(s), s[np.minimum(i, len(s) - 1)] - x, np.inf)
                    left = np.where(i > 0, x - s[np.maximum(i - 1, 0)], np.inf)
                    ovl += int(((right < 300) | (left < 300)).sum())
                rows.append(dict(species=SP, TF=tf, file_split=split, version=ver, n=len(ids), missing=len(ids) - len(ok),
                                 duplicates=len(ids) - len(set(ids)), overlap_own_pos=ovl, self=int((sub[3] == f"{tf}_{SP}").sum()),
                                 same_class=int((sub[2] == F[tf]).sum()) if tf in F else -1))
            st.append(pd.read_csv(f"{root}/outputdir/{tf}/neg_v3_stats_{split}.tsv", sep="\t").assign(species=SP))
        print(SP, split, "done", flush=True)
df = pd.DataFrame(rows); df.to_csv(f"{OUT}/neg_v3_checks.csv", index=False)
S = pd.concat(st, ignore_index=True); S.to_csv(f"{OUT}/neg_v3_stats_all.csv", index=False)
pd.set_option("display.width", 200)
print(df.groupby(["species", "file_split", "version"])[["n", "missing", "duplicates", "overlap_own_pos", "self", "same_class"]].sum().to_string())
print(S.groupby(["species", "split"]).agg(removed=("removed_overlap_own", "sum"), removed_median=("removed_overlap_own", "median"),
      neg=("neg_out", "sum"), ratio_min=("neg_per_pos", "min"), ratio_median=("neg_per_pos", "median"), ratio_max=("neg_per_pos", "max")).to_string())
