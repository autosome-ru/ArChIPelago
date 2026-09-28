#!/usr/bin/env python3
"""Assemble the Figure 4 / S4 table from the per-TF rows (analysis/fig4/rows/<TF>.json) and
the results table (Manuscript_analysis/HUMAN_MOUSE_total_100k.csv).

Mouse test set = chromosomes 1, 8 and 19. Every TF carries all five RF2f-family
models (n = 36 in every row); a missing model or a missing single-model value stops the script.

Rows of the figure (delta = model metric - the best single monoPWM, train-selected; the *_PWM columns of the results table):
  Slim m=0, Slim m=1, LSlim m=-5, diChIPMunk            single models: human test set = the notebook 2 Slim
                                                          tables (analysis/inputs/HUMAN_MOUSE_SLIM_*.txt), mouse test set = the scans of
                                                          the chr1/8/19 set (single_<model>_test_M of the rows)
  RF2f, RF2f+diChIPMunk, RF2f+Slim m=1, RF2f+LSlim m=-5,
  RF2f+Slim m=1+LSlim m=-5+diChIPMunk                     fitted on the best monoPWM + best diPWM (random_state 0)
  RF on all PWMs                                          the ArChIPelago model of Fig. 2/3 (results table)
Outputs: analysis/fig4/fig4_table.csv (wide), Figures/source_data/Figure_4_source_data.csv
and Figure_S4_source_data.csv (long), Sup_Tables/fig4_numbers.json (headline numbers).
"""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts"))
from archi_paths import INPUTS, BASIS, RES  # noqa: E402

ROWS = os.path.join(BASIS, "fig4", "rows")
SLIM_COLS_H = {"slim0": 11, "slim1": 12, "lslim5": 13, "munk": 14}   # 1-based columns of the notebook 2 Slim table
SINGLE = ("slim0", "slim1", "lslim5", "munk")
RF2F = ("RF2f", "RF2f_munk", "RF2f_slim1", "RF2f_lslim5", "RF2f_all5")
MODELS = list(SINGLE) + list(RF2F) + ["RF_all"]
LABEL = {"slim0": "Slim m=0", "slim1": "Slim m=1", "lslim5": "LSlim m=-5", "munk": "diChIPMunk",
         "RF2f": "ArChIPelago RF2f", "RF2f_munk": "ArChIPelago RF2f + diChIPMunk", "RF2f_slim1": "ArChIPelago RF2f + Slim m=1",
         "RF2f_lslim5": "ArChIPelago RF2f + LSlim m=-5",
         "RF2f_all5": "ArChIPelago RF2f + Slim m=1, LSlim m=-5 + diChIPMunk", "RF_all": "ArChIPelago RF on all PWMs"}


def main():
    rows = {os.path.basename(f)[:-5]: json.load(open(f)) for f in sorted(glob.glob(os.path.join(ROWS, "*.json")))}
    T = pd.read_csv(os.path.join(RES, "HUMAN_MOUSE_total_100k.csv"), sep="\t")
    rf = T[(T.Model == "RandomForestClassifier") & (T.PWM == "mono+di")].set_index("TF_name")
    tfs = sorted(rf.index)
    assert len(tfs) == 36
    missing = [tf for tf in tfs if tf not in rows]
    if missing:
        raise SystemExit("fig4 rows missing in %s for: %s" % (ROWS, ", ".join(missing)))
    for tf in tfs:
        r = rows[tf]
        bad = [m for m in RF2F if not r.get(m)] + [f"single_{m}_test_M" for m in SINGLE if not r.get(f"single_{m}_test_M")]
        if bad:
            raise SystemExit("%s.json lacks %s (every TF must carry all five RF2f-family models and the four single models on the mouse test set)" % (tf, ", ".join(bad)))
    slim_h = {}
    for met, fn in (("roc", "HUMAN_MOUSE_SLIM_roc_mono_di_RandomForestClassifier.txt"), ("pr", "HUMAN_MOUSE_SLIM_pr_mono_di_RandomForestClassifier.txt")):
        s = pd.read_csv(os.path.join(INPUTS, fn), sep="\t", header=None)
        s[0] = s[0].str[:-6]
        slim_h[met] = s.set_index(0)

    wide, long = [], []
    checks = []
    for tf in tfs:
        r = rows[tf]
        for sp in ("H", "M"):
            base_roc, base_prc = rf.loc[tf, f"roc_auc_test_{sp}_PWM"], rf.loc[tf, f"pr_auc_test_{sp}_PWM"]
            b_mono, b_di = r[f"best_mono_test_{sp}"], r[f"best_di_test_{sp}"]
            checks.append(abs(b_mono[0] - base_roc))          # the rows' best monoPWM must be the table's baseline
            rec = dict(TF=tf, test_set={"H": "human", "M": "mouse"}[sp], base_auroc=base_roc, base_auprc=base_prc,
                       best_mono_col=r["best_mono_col"], best_di_col=r["best_di_col"],
                       best_mono_auroc=b_mono[0], best_di_auroc=b_di[0], best_mono_auprc=b_mono[1], best_di_auprc=b_di[1])
            vals = {}
            for mname in SINGLE:
                if sp == "H":                     # human test set: the notebook 2 Slim tables
                    c = SLIM_COLS_H[mname]
                    vals[mname] = (slim_h["roc"].loc[tf, c - 1], slim_h["pr"].loc[tf, c - 1])
                else:                             # mouse test set chr1/8/19: scans of the rows
                    vals[mname] = tuple(r[f"single_{mname}_test_M"])
            for mname in RF2F:
                vals[mname] = tuple(r[mname][f"test_{sp}"])
            vals["RF_all"] = (rf.loc[tf, f"roc_auc_test_{sp}"], rf.loc[tf, f"pr_auc_test_{sp}"])
            for mname in MODELS:
                rec[f"{mname}_auroc"], rec[f"{mname}_auprc"] = vals[mname]
                rec[f"{mname}_d_auroc"] = vals[mname][0] - base_roc
                rec[f"{mname}_d_auprc"] = vals[mname][1] - base_prc
                long.append(dict(TF=tf, test_set=rec["test_set"], model=LABEL[mname], model_key=mname,
                                 auROC=vals[mname][0], auPRC=vals[mname][1], base_auROC=base_roc, base_auPRC=base_prc,
                                 d_auROC=vals[mname][0] - base_roc, d_auPRC=vals[mname][1] - base_prc))
            wide.append(rec)
    assert max(checks) < 1e-6, "best-PWM values of the rows disagree with the results table (max diff %g)" % max(checks)
    wide = pd.DataFrame(wide)
    long = pd.DataFrame(long)
    assert long[["auROC", "auPRC"]].notna().all().all() and len(long) == 36 * 2 * len(MODELS)
    os.makedirs(os.path.join(RES, "Figures", "source_data"), exist_ok=True)
    os.makedirs(os.path.join(RES, "Sup_Tables"), exist_ok=True)
    wide.to_csv(os.path.join(BASIS, "fig4", "fig4_table.csv"), index=False)
    long[long.test_set == "human"].to_csv(os.path.join(RES, "Figures", "source_data", "Figure_4_source_data.csv"), index=False)
    long[long.test_set == "mouse"].to_csv(os.path.join(RES, "Figures", "source_data", "Figure_S4_source_data.csv"), index=False)

    num = {}
    for mname in MODELS:
        num[mname] = {}
        for sp, lab in (("human", "H"), ("mouse", "M")):
            d = long[(long.model_key == mname) & (long.test_set == sp)]
            assert len(d) == 36
            num[mname][lab] = dict(n=int(len(d)), median_d_auroc=round(float(d.d_auROC.median()), 4), median_d_auprc=round(float(d.d_auPRC.median()), 4),
                                   n_pos_auroc=int((d.d_auROC > 0).sum()), n_pos_auprc=int((d.d_auPRC > 0).sum()),
                                   p_auroc=float(wilcoxon(d.d_auROC, alternative="greater").pvalue),
                                   p_auprc=float(wilcoxon(d.d_auPRC, alternative="greater").pvalue),
                                   median_auroc=round(float(d.auROC.median()), 4), median_auprc=round(float(d.auPRC.median()), 4))
    num["n_features"] = {m: rows["SRF"][m]["n_features"] for m in RF2F}
    with open(os.path.join(RES, "Sup_Tables", "fig4_numbers.json"), "w") as fh:
        json.dump(num, fh, indent=1)
    for sp, title in (("H", "human test set"), ("M", "mouse test set (chr1/8/19)")):
        print("\n%s: median delta vs best monoPWM / TFs above 0 / one-sided Wilcoxon P" % title)
        for mname in MODELS:
            v = num[mname][sp]
            print("  %-52s auROC %+.4f %2d/%d p=%.1e | auPRC %+.4f %2d/%d p=%.1e" % (LABEL[mname], v["median_d_auroc"], v["n_pos_auroc"], v["n"], v["p_auroc"],
                                                                                  v["median_d_auprc"], v["n_pos_auprc"], v["n"], v["p_auprc"]))


if __name__ == "__main__":
    main()
