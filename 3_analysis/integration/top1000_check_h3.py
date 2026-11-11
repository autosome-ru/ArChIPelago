"""Sup. Table 7 sheet top1000_check from the check C logs of the two data builds (every positive, training and test, lies
among the 1,000 most significant peaks of a listed peak set). The first-submission table counted the published (repeated-read)
human positives (88,144 test / 316,563 all) instead of the positives of the rebuild.
Inputs: the check_C.log of the human and the mouse data build (results/server/human_data/logs, mouse_data/logs); each TF's count
is asserted equal to the positives of the rebuild (pos_H, pos_M: training + test files), and the POOLED line to their sum.
Output: 3_analysis/integration/top1000_check_h3.csv (TF, species, n_positives, n_in_top1000, frac_in_top1000).
usage: ARCHI_STAGE=H5 python3 3_analysis/integration/top1000_check_h3.py
"""
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w2  # noqa: E402
import h3_inputs as H  # noqa: E402

PD, R2, B13 = w2.PD, w2.R2, w2.B13
LOGS = {"H": H.IN["checkC_H"], "M": H.IN["checkC_M"]}   # check C logs of the data build of the stage (h3_inputs)
POS = {"H": H.IN["pos_H"], "M": H.IN["pos_M"]}
rows = []
for sp, log in LOGS.items():
    txt = open(log).read()
    rec = {m.group(1): (int(m.group(2)), int(m.group(3))) for m in re.finditer(r"^(\S+)\tpositives=(\d+)\tnot_in_top1000_of_listed_sets=(\d+)", txt, re.M)}
    if "POOLED" not in rec:
        raise SystemExit("CHECK FAILED %s: no POOLED line (incomplete log)" % log)
    pooled = rec.pop("POOLED")
    PD.check_tfs(rec, "check C %s" % sp)
    pos = pd.read_csv(POS[sp]).set_index("TF")
    n = (pos.train_file_pos + pos.control_file_pos).reindex(sorted(rec))
    bad = [t for t in rec if rec[t][0] != n[t]]
    if bad or pooled[0] != n.sum():
        raise SystemExit("CHECK FAILED check C %s: positives differ from %s for %s (pooled %d vs %d)" % (sp, POS[sp], bad, pooled[0], n.sum()))
    for t in sorted(rec):
        rows.append(dict(TF=t, species=sp, n_positives=rec[t][0], n_in_top1000=rec[t][0] - rec[t][1], frac_in_top1000=(rec[t][0] - rec[t][1]) / rec[t][0]))
out = pd.DataFrame(rows)
dst = os.path.join(B13, "top1000_check_h3.csv")
w2.backup(dst)
out.to_csv(dst, index=False)
print(out.groupby("species")[["n_positives", "n_in_top1000"]].sum())
w2.provenance(B13, "3_analysis/integration/top1000_check_h3.py", [dst], list(LOGS.values()) + list(POS.values()))
