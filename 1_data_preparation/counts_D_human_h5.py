"""counts_D_human_h5.py - per-TF counts of the human data build beside the published ones.
Human naming: file "train" = training set (chr2-7, 9-20, 22), file "control"/"test" = test set (chr1, 8, 21).
Positives are counted from out_tab_HUMAN_10000_{train,control}.tab (rows, and unique windows = unique 'chrom/center/
signal' header after '$'), negatives from outputdir/<TF>/{train_neg_id_out_HUMAN_result_H_train_full,
test_neg_id_out_HUMAN_result_H_control_full}.csv. Writes human_sym_h5/TF-ML/counts_per_tf.csv (refuses to overwrite).
usage: python counts_D_human_h5.py > counts_D.log"""
import os
import pandas as pd
H = os.path.expanduser("~")
NEW = f"{H}/archi_work/human_sym_h5/TF-ML"; OLD = f"{H}/Release/TF-ML"
out = f"{NEW}/counts_per_tf.csv"
if os.path.exists(out): raise SystemExit("refusing to overwrite " + out)
c = pd.read_csv(f"{NEW}/per_tf_positive_counts.csv")
npw = pd.read_csv(f"{NEW}/n_pwms_per_folder.tsv", sep="\t").pivot(index="TF", columns="folder", values="n_pwm")


def tab(root, split):
    t = pd.read_csv(f"{root}/out_tab_HUMAN_10000_{split}.tab", sep="\t", header=None, usecols=[1, 3])
    t["TF"] = t[3].str.replace("_HUMAN", "", regex=False); t["key"] = t[1].str.split("$").str[1]
    return t.groupby("TF")[3].size(), t.groupby("TF")["key"].nunique()


def nlines(p):
    return sum(1 for _ in open(p)) if os.path.exists(p) else None


(tn, tnu), (cn, cnu) = tab(NEW, "train"), tab(NEW, "control")
(to, tou), (co, cou) = tab(OLD, "train"), tab(OLD, "control")
rows = []
for _, r in c.iterrows():
    tf = r.TF
    d = dict(TF=tf, peak_sets=r.peak_sets, peak_set_reads=r.peak_set_reads, pool_rows_before_mask=r.pool_rows,
             pool_unique_windows=r.pool_unique_windows, pool_after_mask=r.after_mask, removed_by_mask=r.removed_by_mask,
             share_removed=round(r.removed_by_mask / r.pool_rows, 4),
             pos_train=int(tn.get(tf, 0)), pos_train_unique=int(tnu.get(tf, 0)),
             pos_test=int(cn.get(tf, 0)), pos_test_unique=int(cnu.get(tf, 0)),
             neg_train=nlines(f"{NEW}/outputdir/{tf}/train_neg_id_out_HUMAN_result_H_train_full_v3.csv"),
             neg_test=nlines(f"{NEW}/outputdir/{tf}/test_neg_id_out_HUMAN_result_H_control_full_v3.csv"),
             orig_pos_train=int(to.get(tf, 0)), orig_pos_train_unique=int(tou.get(tf, 0)),
             orig_pos_test=int(co.get(tf, 0)), orig_pos_test_unique=int(cou.get(tf, 0)),
             orig_neg_train=nlines(f"{OLD}/outputdir/{tf}/train_neg_id_out_HUMAN_result_H_train_full.csv"),
             orig_neg_test=nlines(f"{OLD}/outputdir/{tf}/test_neg_id_out_HUMAN_result_H_control_full.csv"))
    for f in ("HUMAN_seq_HUMAN_pwm_mono", "HUMAN_seq_HUMAN_pwm_di", "HUMAN_seq_MOUSE_pwm_mono", "HUMAN_seq_MOUSE_pwm_di"):
        d["n_" + f.replace("HUMAN_seq_", "")] = int(npw.loc[tf, f])
    rows.append(d)
df = pd.DataFrame(rows)
tot = {k: (df[k].sum() if df[k].dtype != object else "") for k in df.columns}
tot["TF"] = "POOLED"; tot["share_removed"] = round(df.removed_by_mask.sum() / df.pool_rows_before_mask.sum(), 4)
df = pd.concat([df, pd.DataFrame([tot])], ignore_index=True)
df.to_csv(out, index=False)
pd.set_option("display.width", 300); pd.set_option("display.max_columns", 40)
print(df.to_string(index=False))
