"""s4_numbers.py - headline numbers of a Figure_S4_source_data CSV, as numbers_in_text_Figure_4.json of the package assembler
(per model: median delta vs best monoPWM, TFs with delta > 0, one-sided Wilcoxon P, median auROC/auPRC).
usage: python s4_numbers.py OUT_CSV LABEL=CSV[:EXCLUDE_TF,...] ...   (refuses to overwrite OUT_CSV)"""
import os, sys
import pandas as pd
from scipy.stats import wilcoxon
out = sys.argv[1]
if os.path.exists(out):
    raise SystemExit("exists: " + out)
rec = []
for arg in sys.argv[2:]:
    lab, rest = arg.split("=", 1)
    path, _, excl = rest.partition(":")
    d = pd.read_csv(path)
    d = d[~d.TF.isin([x for x in excl.split(",") if x])]
    for key, g in d.groupby("model_key", sort=False):
        rec.append(dict(version=lab, model_key=key, model=g.model.iloc[0], n_TF=len(g),
                        median_d_auROC=g.d_auROC.median(), n_TF_d_auROC_gt0=int((g.d_auROC > 0).sum()),
                        p_auROC=wilcoxon(g.d_auROC, alternative="greater").pvalue,
                        median_d_auPRC=g.d_auPRC.median(), n_TF_d_auPRC_gt0=int((g.d_auPRC > 0).sum()),
                        p_auPRC=wilcoxon(g.d_auPRC, alternative="greater").pvalue,
                        median_auROC=g.auROC.median(), median_auPRC=g.auPRC.median(),
                        median_base_auROC=g.base_auROC.median(), median_base_auPRC=g.base_auPRC.median()))
r = pd.DataFrame(rec)
r.to_csv(out, index=False)
pd.set_option("display.width", 250)
print(r.drop(columns=["model"]).round(4).to_string(index=False))
