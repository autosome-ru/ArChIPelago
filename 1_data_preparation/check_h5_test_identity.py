"""check_h5_test_identity.py - every TEST-set file of the current data build must be byte-identical to the earlier
symmetric build (human_sym / mouse_sym). Test side on disk: HUMAN files "test"/"control" (chr1, 8, 21); MOUSE files
"train" (chr1, 8, 19).
stage positives: per-TF test beds and mfa, all_bed_df test, pooled fasta (300 bp, SLIM, _no_NF, _no_NF_no_N), headerStr csv,
                 out_tab of the test split.
stage negatives: <TF> v3 negative-id list of the test split.
stage sarus:     every feature table of the test split (<SP>_seq_{HUMAN,MOUSE}_pwm_{mono,di}/*_<split>_cut.tab).
Red control in every stage: the TRAINING-side counterpart of the first file is compared with the earlier build and must
differ (the training sets changed), so a comparison that cannot see a difference is caught.
Writes OUTDIR/test_identity_<stage>.csv (path, md5_h4, md5_h5, identical) and prints the counts; exit 1 if any test file
differs or is missing, or if the red control does not differ.
usage: python check_h5_test_identity.py STAGE OUTDIR"""
import os, sys, glob, hashlib
import pandas as pd

stage, outdir = sys.argv[1:3]
H = os.path.expanduser("~"); R = f"{H}/archi_work"
ROOTS = {"HUMAN": (f"{R}/human_sym/TF-ML", f"{R}/human_sym_h5/TF-ML"), "MOUSE": (f"{R}/mouse_sym/TF-ML", f"{R}/mouse_sym_h5/TF-ML")}
TEST = {"HUMAN": ("test", "control"), "MOUSE": ("train", "train")}  # (per-TF folder/bed split, pooled split) of the test set
TRAIN = {"HUMAN": ("train", "train"), "MOUSE": ("test", "control")}
out = os.path.join(outdir, "test_identity_%s.csv" % stage)
if os.path.exists(out):
    raise SystemExit("refusing to overwrite " + out)


def md5(p):
    if not os.path.isfile(p):
        return "MISSING"
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def rel_files(sp, side):
    old, _ = ROOTS[sp]
    d, s = side[sp]
    if stage == "positives":
        fs = [f"{d}/{os.path.basename(f)}" for f in sorted(glob.glob(f"{old}/{d}/*_{d}_set_{sp}.bed"))]
        fs += [f"{d}/{os.path.basename(f)}" for f in sorted(glob.glob(f"{old}/{d}/*_{sp}.PEAKS000001.macs.{d}.mfa"))]
        fs += [f"all_bed_df_{d}_set_{sp}.bed", f"out_tab_{sp}_10000_{s}.tab"]
        for tag in ("", "SLIM_"):
            fs += [f"all_mfa_file_{sp}_{tag}10000_{s}{x}.fasta" for x in ("", "_no_NF", "_no_NF_no_N")]
            fs += [f"headerStr_names_{sp}_{tag}10000_{s}.csv"]
        return fs
    if stage == "negatives":
        mode = "train" if s == "train" else "test"
        return [os.path.relpath(f, old) for f in sorted(glob.glob(f"{old}/outputdir/*/{mode}_neg_id_out_{sp}_result_{sp[0]}_{s}_full_v3.csv"))]
    if stage == "sarus":
        return [os.path.relpath(f, old) for f in sorted(glob.glob(f"{old}/outputdir/*/{sp}_seq_*_pwm_*/*_feature_compare_table_{s}_cut.tab"))]
    raise SystemExit("unknown stage " + stage)


rows, bad, red_ok = [], 0, True
for sp in ("HUMAN", "MOUSE"):
    old, new = ROOTS[sp]
    fs = rel_files(sp, TEST)
    if not fs:
        raise SystemExit("no H4 files found for %s %s" % (sp, stage))
    for f in fs:
        a, b = md5(os.path.join(old, f)), md5(os.path.join(new, f))
        rows.append(dict(species=sp, path=f, md5_h4=a, md5_h5=b, identical=(a == b and a != "MISSING")))
    ctrl = rel_files(sp, TRAIN)[0]
    a, b = md5(os.path.join(old, ctrl)), md5(os.path.join(new, ctrl))
    red = a != b
    red_ok &= red
    print("%s %s: red control (training-side file %s) differs from H4: %s (%s vs %s)" % (sp, stage, ctrl, red, a[:8], b[:8]))
df = pd.DataFrame(rows)
os.makedirs(outdir, exist_ok=True)
df.to_csv(out, index=False)
for sp, g in df.groupby("species"):
    print("%s %s: %d test files, %d identical, %d differ or missing" % (sp, stage, len(g), int(g.identical.sum()), int((~g.identical).sum())))
    print("  e.g.", g.path.iloc[0], g.md5_h4.iloc[0], g.md5_h5.iloc[0])
if not df.identical.all() or not red_ok:
    print("TEST_IDENTITY_FAILED", stage)
    sys.exit(1)
print("TEST_IDENTITY_OK", stage)
