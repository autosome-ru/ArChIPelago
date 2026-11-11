"""Per-TF low-FPR metric tables (Sup. Table 7, Fig. S8) in 3_analysis/low_fpr_metrics from the server mirror.

The partial auROC (pauc_*) of the first b2 run was biased; the corrected run is b2_v2 (h3_inputs b2_dir). Without the b2_v2
mirror every pauc_* column is left out of the copies, so no table or number carries the biased values; with b2_v2 present
the full tables are copied from there.
Files: B2_metrics_wide.csv, B2_per_chromosome.csv (HH, HM, MM; RF seeds 0/1, best monoPWM by auROC / auPRC selected on training).
usage: python3 3_analysis/integration/make_b2_h3.py
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w2  # noqa: E402
import h3_inputs as H  # noqa: E402

PD, B13 = w2.PD, w2.B13
O = os.path.join(B13, "../low_fpr_metrics")
V2 = H.IN["b2_dir"]
PENDING = not os.path.isdir(V2)
SRC = {"B2_metrics_wide.csv": H.IN["b2_wide"], "B2_per_chromosome.csv": H.IN["b2_chrom"]}
if not PENDING:
    SRC = {k: os.path.join(V2, k) for k in SRC}
ins, outs = [], []
for name, src in SRC.items():
    d = pd.read_csv(src)
    for sc in ("HH", "HM", "MM"):
        PD.check_tfs(d[d.scenario == sc].TF, "%s %s" % (name, sc))
    if PENDING:
        d = d.drop(columns=[c for c in d.columns if c.startswith("pauc_")])
    dst = os.path.join(O, name)
    w2.backup(dst)
    d.to_csv(dst, index=False)
    ins.append(src); outs.append(dst)
    print("%s: %d rows from %s%s" % (name, len(d), os.path.relpath(src, w2.R2), " (pauc_* pending: b2_v2)" if PENDING else ""))
stale = os.path.join(O, "B2_metrics_per_TF_long.csv")   # old copy of the long table; not read by any producer
if os.path.exists(stale):
    w2.backup(stale); os.remove(stale)
w2.provenance(O, "3_analysis/integration/make_b2_h3.py", outs, ins, params="{pauc_pending: %s}" % PENDING,
              note="H3: B2 tables from the H3a mirror" + ("; pauc_* columns pending (R2-F1, b2_v2)" if PENDING else " (b2_v2)"))
