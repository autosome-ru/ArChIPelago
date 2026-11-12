"""b9_summary_h2b.py - b9_summary.csv from a b9_random_split_per_fit.csv: per TF and split, mean over seeds 0-2;
gain = RF minus best monoPWM. No figure (rendered with archi_style).
usage: python b9_summary_h2b.py PER_FIT_CSV OUT_CSV   (refuses to overwrite)"""
import os, sys
import pandas as pd
src, out = sys.argv[1:3]
if os.path.exists(out):
    raise SystemExit("exists: " + out)
d = pd.read_csv(src)
assert "P53" not in set(d.TF)
d["gain_auroc"] = d.rf_auroc - d.pwm_auroc; d["gain_auprc"] = d.rf_auprc - d.pwm_auprc
m = d.groupby(["TF", "split"])[["rf_auroc", "pwm_auroc", "gain_auroc", "rf_auprc", "pwm_auprc", "gain_auprc",
                                "frac_test_pos_overlap", "frac_test_pos_exact"]].mean().reset_index()
m.to_csv(out, index=False)
print(m.round(4).to_string())
