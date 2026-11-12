"""Every number of Fig. 1 recomputed from the final Sup. Tables and the final-data inputs, written to
3_analysis/integration/figure1_values.csv. Fig. 1 is drawn by hand by the authors and is not in the repository; this file holds
the numbers it shows (and feeds make_numbers_w2.py).
Sources: Sup. Tables 1-4 (outputs/supplementary_tables), hm_mm (the negatives and positives each model was trained and tested
on, for the realised negative:positive ratios).
usage: ARCHI_STAGE=H5 python3 3_analysis/integration/figure1_values.py
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w2  # noqa: E402
import h3_inputs as H  # noqa: E402

PD, R2 = w2.PD, w2.R2
D = os.path.join(R2, "outputs/supplementary_tables")
SRC = {"t1": os.path.join(D, "Sup_Table_1_ChIP-Seq_experiments.xlsx"), "t2": os.path.join(D, "Sup_Table_2_PWMs_and_datasets.xlsx"),
       "t3": os.path.join(D, "Sup_Table_3_ArChIPelago_performance.xlsx"), "t4": os.path.join(D, "Sup_Table_4_operational_metrics.xlsx"),
       "hm": H.IN["hm_mm"], "s6": os.path.join(w2.B13, "sup6_methods_values_h3.csv")}   # Fig. 1A cost = Methods 2.11
t1 = pd.read_excel(SRC["t1"], sheet_name=None)
t2 = pd.read_excel(SRC["t2"], sheet_name="per_TF").set_index("TF")
t3 = pd.read_excel(SRC["t3"], sheet_name="summary")
t4 = pd.read_excel(SRC["t4"], sheet_name="summary")
hm = pd.read_csv(SRC["hm"]).set_index("TF")
s6 = pd.read_csv(SRC["s6"]).set_index("quantity").h3
s6m = lambda prefix: next(v for k, v in s6.items() if k.startswith(prefix))
NTF = PD.check_tfs(t2.index, "Sup. Table 2")
PD.check_tfs(hm.index, "hm_mm_results")
nH, nM = len(t1["human"]), len(t1["mouse"])
minH, minM = t1["human"].groupby("TF").size(), t1["mouse"].groupby("TF").size()
trH, teH = t2[H.COL["HTR"]].sum(), t2[H.COL["HTE"]].sum()
trM, teM = t2[H.COL["MTR"]].sum(), t2[H.COL["MTE"]].sum()
if (trH, teH, trM, teM) != (hm.n_train_pos_H.sum(), hm.n_test_pos_H.sum(), hm.n_train_pos_M.sum(), hm.n_test_pos_M.sum()):
    raise SystemExit("CHECK FAILED positives: Sup. Table 2 vs hm_mm_results")
h = t3[t3["Test set"].str.startswith("human")].iloc[0]
op = t4[(t4["Test set"].str.startswith("human")) & (t4.Quantity == "FP reduction at 50 % recall (%)")].iloc[0]
# realised negatives per positive of every training and test set (BiasAway asks for up to 100 per positive)
rH = pd.concat([hm.n_train_neg_H / hm.n_train_pos_H, hm.n_test_neg_H / hm.n_test_pos_H])
rM = pd.concat([hm.n_train_neg_M / hm.n_train_pos_M, hm.n_test_neg_M / hm.n_test_pos_M])
V = [
    ("human datasets", nH), ("mouse datasets", nM), ("datasets total", nH + nM), ("TF pairs", NTF),
    ("min datasets per TF, human", int(minH.min())), ("TFs at the human minimum", ",".join(minH.index[minH == minH.min()])),
    ("min datasets per TF, mouse", int(minM.min())), ("TFs at the mouse minimum", ",".join(minM.index[minM == minM.min()])),
    ("human training positives", int(trH)), ("human test positives", int(teH)), ("mouse training positives", int(trM)), ("mouse test positives", int(teM)),
    ("test share of positives, human (%)", 100 * teH / (teH + trH)), ("test share of positives, mouse (%)", 100 * teM / (teM + trM)),
    ("negatives per positive, human, min", rH.min()), ("negatives per positive, human, median", rH.median()), ("negatives per positive, human, max", rH.max()),
    ("negatives per positive, mouse, min", rM.min()), ("negatives per positive, mouse, median", rM.median()), ("negatives per positive, mouse, max", rM.max()),
    ("human monoPWMs", int(t2["Human monoPWMs"].sum())), ("human diPWMs", int(t2["Human diPWMs"].sum())), ("human PWMs", int(t2["Human PWMs (mono + di)"].sum())),
    ("mouse monoPWMs", int(t2["Mouse monoPWMs"].sum())), ("mouse diPWMs", int(t2["Mouse diPWMs"].sum())), ("mouse PWMs", int(t2["Mouse PWMs (mono + di)"].sum())),
    ("median PWMs per TF, human", float(t2["Human PWMs (mono + di)"].median())), ("median PWMs per TF, mouse", float(t2["Mouse PWMs (mono + di)"].median())),
    ("min monoPWMs per TF (both species)", int(min(t2["Human monoPWMs"].min(), t2["Mouse monoPWMs"].min()))),
    ("min diPWMs per TF (both species)", int(min(t2["Human diPWMs"].min(), t2["Mouse diPWMs"].min()))),
    ("human best monoPWM median auROC", h["baseline auROC (median)"]), ("human RF median auROC", h["ArChIPelago auROC (median)"]),
    ("human best monoPWM median auPRC", h["baseline auPRC (median)"]), ("human RF median auPRC", h["ArChIPelago auPRC (median)"]),
    ("human TFs improved auROC", int(h["TFs with dauROC > 0"])), ("human TFs improved auPRC", int(h["TFs with dauPRC > 0"])),
    ("human FP reduction at 50 % recall (%)", float(op["Median over TFs"])),   # header without the TF count (34 TFs here)
    ("scan s per monoPWM per Mb", s6m("SARUS wall-clock per monoPWM per Mb")),
    ("inference s per 10^6 sequences, one core (range)", s6m("RF inference, predict_proba for 10^6 sequences on one core, chunked")),
]
out = pd.DataFrame(V, columns=["quantity", "value"])
dst = os.path.join(w2.B13, "figure1_values.csv")
w2.backup(dst)
out.to_csv(dst, index=False)
w2.provenance(w2.B13, "3_analysis/integration/figure1_values.py", [dst], list(SRC.values()))
print(out.to_string(index=False))
