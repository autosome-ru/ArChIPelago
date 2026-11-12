"""
counts_D_mouse_h5.py - per-TF counts table of the mouse data build, with the published mouse numbers beside it.
Reads (all read-only): mouse_sym_h5/TF-ML/{per_tf_positive_counts.csv, out_tab_MOUSE_10000_{train,control}.tab,
n_pwms_per_folder.tsv, outputdir/<TF>/*_neg_id_out_MOUSE_*_full.csv} and the same out_tab / neg-id files of the published TF-ML folder ($ARCHI_RELEASE).
Writes mouse_sym_h5/TF-ML/counts_per_tf.csv (refuses to overwrite).
File naming: on disk, mouse files named "train" hold the test chromosomes 1, 8 and 19; files named "control" hold the
training chromosomes (2-7, 9-18).
usage: python counts_D_mouse_h5.py > counts_D.log
"""
import os
import pandas as pd

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")); import archi_env as AE
NEW = AE.work("mouse_sym_h5", "TF-ML")
OLD = AE.release()
out = f"{NEW}/counts_per_tf.csv"
if os.path.exists(out):
    raise SystemExit("refusing to overwrite " + out)
c = pd.read_csv(f"{NEW}/per_tf_positive_counts.csv")
npw = pd.read_csv(f"{NEW}/n_pwms_per_folder.tsv", sep="\t").pivot(index="TF", columns="folder", values="n_pwm")


def tab_counts(root, split):
    t = pd.read_csv(f"{root}/out_tab_MOUSE_10000_{split}.tab", sep="\t", header=None, usecols=[3])[3]
    return t.str.replace("_MOUSE", "", regex=False).value_counts()


def nlines(p):
    return sum(1 for _ in open(p)) if os.path.exists(p) else None


rows = []
tn, cn = tab_counts(NEW, "train"), tab_counts(NEW, "control")
to, co = tab_counts(OLD, "train"), tab_counts(OLD, "control")
for _, r in c.iterrows():
    tf = r.TF
    d = dict(TF=tf, peak_sets=r.peak_sets, pool_rows_before_mask=r.pool_rows, pool_unique_windows=r.pool_unique_windows,
             pool_after_mask=r.after_mask, removed_by_mask=r.removed_by_mask,
             share_removed=round(r.removed_by_mask / r.pool_rows, 4),
             pos_chr1_8_19_file_train=int(tn.get(tf, 0)), pos_other_file_control=int(cn.get(tf, 0)),
             neg_chr1_8_19=nlines(f"{NEW}/outputdir/{tf}/train_neg_id_out_MOUSE_result_M_train_full_v3.csv"),
             neg_other=nlines(f"{NEW}/outputdir/{tf}/test_neg_id_out_MOUSE_result_M_control_full_v3.csv"),
             orig_pos_chr1_8_19=int(to.get(tf, 0)), orig_pos_other=int(co.get(tf, 0)),
             orig_neg_chr1_8_19=nlines(f"{OLD}/outputdir/{tf}/train_neg_id_out_MOUSE_result_M_train_full.csv"),
             orig_neg_other=nlines(f"{OLD}/outputdir/{tf}/test_neg_id_out_MOUSE_result_M_control_full.csv"))
    for f in ("MOUSE_seq_HUMAN_pwm_mono", "MOUSE_seq_HUMAN_pwm_di", "MOUSE_seq_MOUSE_pwm_mono", "MOUSE_seq_MOUSE_pwm_di"):
        d["n_" + f.replace("MOUSE_seq_", "")] = int(npw.loc[tf, f])
    rows.append(d)
df = pd.DataFrame(rows)
tot = {k: (df[k].sum() if df[k].dtype != object else "") for k in df.columns}
tot["TF"] = "POOLED"
tot["share_removed"] = round(df.removed_by_mask.sum() / df.pool_rows_before_mask.sum(), 4)
df = pd.concat([df, pd.DataFrame([tot])], ignore_index=True)
df.to_csv(out, index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(df.to_string(index=False))
