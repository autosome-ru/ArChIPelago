"""Inputs of the top-peaks summary figure (Fig. S9, Sup. Table 8) in 3_analysis/top_peaks from the server mirrors. Every file
holds the 35 TFs; the peak widths hold exactly the used peak sets (361 + 276).
Inputs (h3_inputs.py): b3_topN (both species), b3_widths_H/M, b3_join_H/M and the per-sequence joins next to them.
usage: python3 3_analysis/integration/make_b3_h3.py
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w2  # noqa: E402
import h3_inputs as H  # noqa: E402

PD, R2, B13 = w2.PD, w2.R2, w2.B13
O = os.path.join(B13, "../top_peaks")
H.check()
ins, outs = [], []


def put(df, name, srcs):
    dst = os.path.join(O, name)
    w2.backup(dst)
    df.to_csv(dst, index=False)
    ins.extend(srcs); outs.append(dst)
    print("%s: %d rows" % (name, len(df)))


t = pd.read_csv(H.IN["b3_topN"])
for sp in ("H", "M"):
    PD.check_tfs(t[t.species == sp].TF, "b3_topN_results %s" % sp)
put(t, "b3_topN_results.csv", [H.IN["b3_topN"]])

w = pd.concat([pd.read_csv(H.IN["b3_widths_H"]), pd.read_csv(H.IN["b3_widths_M"])], ignore_index=True)
used = pd.concat([PD.human_ids().assign(species="H"), PD.mouse_ids().assign(species="M")]).rename(columns={"peak_id": "peak_set"})
key = ["peak_set", "TF", "species"]
if sorted(map(tuple, w[key].values)) != sorted(map(tuple, used[key].values)):
    raise SystemExit("CHECK FAILED widths: the H3a width tables do not list exactly the %d used peak sets" % len(used))
put(w, "b3_peak_widths_per_dataset.csv", [H.IN["b3_widths_H"], H.IN["b3_widths_M"]])

jt = pd.concat([pd.read_csv(H.IN["b3_join_H"]), pd.read_csv(H.IN["b3_join_M"])], ignore_index=True)
for sp in ("H", "M"):
    PD.check_tfs(jt[jt.species == sp].TF, "join per TF %s" % sp)
put(jt, "b3_join_test_per_tf.csv", [H.IN["b3_join_H"], H.IN["b3_join_M"]])

seqH, seqM = H.IN["b3_join_H"].replace("per_tf", "per_sequence"), H.IN["b3_join_M"].replace("per_tf", "per_sequence")
js = pd.concat([pd.read_csv(seqH), pd.read_csv(seqM)], ignore_index=True)
pos = {"H": int(pd.read_csv(H.IN["pos_H"]).control_file_pos.sum()), "M": int(pd.read_csv(H.IN["pos_M"]).train_file_pos.sum())}
for sp in ("H", "M"):
    if (js.species == sp).sum() != pos[sp]:
        raise SystemExit("CHECK FAILED per-sequence join %s: %d rows, %d test positives" % (sp, (js.species == sp).sum(), pos[sp]))
put(js, "b3_join_test_per_sequence.csv", [seqH, seqM, H.IN["pos_H"], H.IN["pos_M"]])
w2.provenance(O, "3_analysis/integration/make_b3_h3.py", outs, sorted(set(ins)), note="H3: B3 inputs from the H3a mirrors (both species)")
