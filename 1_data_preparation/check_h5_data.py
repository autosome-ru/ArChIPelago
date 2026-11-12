"""check_h5_data.py - data checks after the build of the symmetric human and mouse tables (three checks and a split check):
 split  chromosomes present in the tables: logical train only training chromosomes, test only test chromosomes,
        no sex chromosome; every autosome with positives in exactly one table; training-only chromosomes listed.
 D      per TF and logical split: positives (out_tab rows) and v3 negatives of this build vs the earlier build
        (counts_h5_vs_h4.csv).
 super  per TF (training split): every training window of the earlier build is a training window of this build when it
        is not capped (pre-cap pool <= 10,000); key = (chrom, centre) from the out_tab header.
 B      SARUS: for 3 TFs x 4 feature folders x first 2 PWMs, scores of the training table on the windows shared with
        the earlier build (same key) are equal; red control: the same comparison with the scores shifted by one row.
usage: python check_h5_data.py OUTDIR   (refuses to overwrite); prints CHECK_H5_DATA_OK / CHECK_H5_DATA_FAILED"""
import os, sys, glob
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h5_splits as S

out = sys.argv[1]
sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")); import archi_env as AE; R = AE.work()
ROOT = {("H4", "HUMAN"): f"{R}/human_sym/TF-ML", ("H5", "HUMAN"): f"{R}/human_sym_h5/TF-ML",
        ("H4", "MOUSE"): f"{R}/mouse_sym/TF-ML", ("H5", "MOUSE"): f"{R}/mouse_sym_h5/TF-ML"}
FILE = {("HUMAN", "train"): "train", ("HUMAN", "test"): "control", ("MOUSE", "train"): "control", ("MOUSE", "test"): "train"}
for f in ("counts_h5_vs_h4.csv", "split_chromosomes.csv", "sarus_shared_rows.csv"):
    if os.path.exists(os.path.join(out, f)):
        raise SystemExit("refusing to overwrite " + f)
os.makedirs(out, exist_ok=True)
fails = []


def tab(v, sp, split):
    t = pd.read_csv(f"{ROOT[(v, sp)]}/out_tab_{sp}_10000_{FILE[(sp, split)]}.tab", sep="\t", header=None, dtype=str)
    h = t[1].str.split("$").str[1]
    t["chrom"] = h.str.extract(r"chrom: (\w+);")[0]
    t["center"] = h.str.extract(r"center: ([\d.]+);")[0]
    t["name"] = t[1].str.split("$").str[0]
    return t


def neg_n(v, sp, split, tf):
    fs = FILE[(sp, split)]
    mode = "train" if fs == "train" else "test"
    p = f"{ROOT[(v, sp)]}/outputdir/{tf}/{mode}_neg_id_out_{sp}_result_{sp[0]}_{fs}_full_v3.csv"
    return sum(1 for _ in open(p))


rows, chrom_rows, sar_rows = [], [], []
for sp in ("HUMAN", "MOUSE"):
    T = {(v, s): tab(v, sp, s) for v in ("H4", "H5") for s in ("train", "test")}
    # split
    seen = {}
    for s in ("train", "test"):
        cs = sorted(T[("H5", s)].chrom.unique(), key=lambda c: (len(c), c))
        seen[s] = set(cs)
        extra = set(cs) - set(S.LOGICAL[sp][s])
        if extra:
            fails.append("%s %s holds chromosomes outside its set: %s" % (sp, s, sorted(extra)))
        for c in cs:
            chrom_rows.append(dict(species=sp, logical_split=s, chrom=c, rows_h5=int((T[("H5", s)].chrom == c).sum()),
                                   rows_h4=int((T[("H4", s)].chrom == c).sum())))
    if set(S.SEX) & (seen["train"] | seen["test"]):
        fails.append("%s: sex chromosome in a table" % sp)
    missing = sorted(set(S.LOGICAL[sp]["train"]) - seen["train"])
    if missing:
        fails.append("%s: training chromosomes without positives %s" % (sp, missing))
    new_only = sorted(seen["train"] - set(T[("H4", "train")].chrom), key=lambda c: int(c[3:]))
    print("%s training chromosomes in H5: %d; not in H4: %s; test chromosomes: %s" % (sp, len(seen["train"]), new_only, sorted(seen["test"])))
    pc = pd.read_csv(f"{ROOT[('H5', sp)]}/per_tf_positive_counts.csv").set_index("TF")
    precap_col = "control_before_cap" if sp == "MOUSE" else "train_chr1_8_19_before_cap"   # logical training split
    for tf in sorted(pc.index):
        r = dict(species=sp, TF=tf)
        for s in ("train", "test"):
            for v in ("H4", "H5"):
                r["pos_%s_%s" % (s, v)] = int((T[(v, s)][3] == "%s_%s" % (tf, sp)).sum())
                r["neg_%s_%s" % (s, v)] = neg_n(v, sp, s, tf)
        r["train_precap_H5"] = int(pc.loc[tf, precap_col])
        k4 = set(zip(*[T[("H4", "train")].loc[T[("H4", "train")][3] == "%s_%s" % (tf, sp), c] for c in ("chrom", "center")]))
        k5 = set(zip(*[T[("H5", "train")].loc[T[("H5", "train")][3] == "%s_%s" % (tf, sp), c] for c in ("chrom", "center")]))
        r["h4_train_windows_in_h5"] = len(k4 & k5)
        r["h5_train_windows_new_chroms"] = sum(1 for c, _ in k5 if c in new_only)
        if r["train_precap_H5"] <= 10000 and not k4 <= k5:
            fails.append("%s %s: %d H4 training windows missing from the uncapped H5 training set" % (sp, tf, len(k4 - k5)))
        if r["pos_test_H4"] != r["pos_test_H5"] or r["neg_test_H4"] != r["neg_test_H5"]:
            fails.append("%s %s: test counts differ" % (sp, tf))
        rows.append(r)
    # B: SARUS on shared training windows
    t4, t5 = T[("H4", "train")], T[("H5", "train")]
    key4 = pd.Series(np.arange(len(t4)), index=pd.MultiIndex.from_arrays([t4[3], t4.chrom, t4.center]))
    key5 = pd.Series(np.arange(len(t5)), index=pd.MultiIndex.from_arrays([t5[3], t5.chrom, t5.center]))
    key4 = key4[~key4.index.duplicated()]; key5 = key5[~key5.index.duplicated()]
    common = key4.index.intersection(key5.index)
    i4, i5 = key4.loc[common].values, key5.loc[common].values
    fs = FILE[(sp, "train")]
    for tf in ("SRF", "CTCF", "STA5A"):
        for folder in ("%s_seq_HUMAN_pwm_mono" % sp, "%s_seq_HUMAN_pwm_di" % sp, "%s_seq_MOUSE_pwm_mono" % sp, "%s_seq_MOUSE_pwm_di" % sp):
            fl = sorted(glob.glob(f"{ROOT[('H5', sp)]}/outputdir/{tf}/{folder}/*_feature_compare_table_{fs}_cut.tab"),
                        key=lambda p: int(os.path.basename(p).split("_")[0]))[:2]
            for f5 in fl:
                f4 = f5.replace(ROOT[("H5", sp)], ROOT[("H4", sp)])
                a = np.loadtxt(f4, dtype=str); b = np.loadtxt(f5, dtype=str)
                if len(b) != len(t5) or len(a) != len(t4):
                    fails.append("SARUS rows %s" % f5)
                eq = float((a[i4] == b[i5]).mean())
                red = float((a[i4] == np.roll(b, 1)[i5]).mean())
                sar_rows.append(dict(species=sp, TF=tf, file=os.path.relpath(f5, ROOT[("H5", sp)]), shared_rows=len(common),
                                     equal_fraction=eq, red_shifted_equal_fraction=red))
                if eq != 1.0 or red > 0.9:
                    fails.append("SARUS shared rows %s: equal %.4f, red %.4f" % (f5, eq, red))
pd.DataFrame(rows).to_csv(os.path.join(out, "counts_h5_vs_h4.csv"), index=False)
pd.DataFrame(chrom_rows).to_csv(os.path.join(out, "split_chromosomes.csv"), index=False)
sr = pd.DataFrame(sar_rows); sr.to_csv(os.path.join(out, "sarus_shared_rows.csv"), index=False)
d = pd.DataFrame(rows)
for sp, g in d.groupby("species"):
    print("%s positives train H4 %d -> H5 %d; test H4 %d -> H5 %d; negatives train H4 %d -> H5 %d; test %d -> %d" % (
        sp, g.pos_train_H4.sum(), g.pos_train_H5.sum(), g.pos_test_H4.sum(), g.pos_test_H5.sum(),
        g.neg_train_H4.sum(), g.neg_train_H5.sum(), g.neg_test_H4.sum(), g.neg_test_H5.sum()))
print("SARUS shared-row files: %d, all equal: %s, red control max equal fraction %.4f" % (len(sr), bool((sr.equal_fraction == 1).all()), sr.red_shifted_equal_fraction.max()))
print("\n".join(fails) if fails else "no failures")
print("CHECK_H5_DATA_FAILED" if fails else "CHECK_H5_DATA_OK")
