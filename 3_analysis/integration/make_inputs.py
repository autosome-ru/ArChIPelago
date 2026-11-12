"""Inputs of the human figure scripts in 4_visualisation/main_figures (ARCHI_PKG=3_analysis/integration/inputs).

The human tables are the final ones of the package 3_analysis/package; this script copies them byte for byte into
3_analysis/integration/inputs/ (folder layout of the package) after the 35-TF check. Files:
  results_table.csv (tab)                          <- 3_analysis/package/results_table.csv
  analysis/inputs/training_positives_per_TF.csv    <- 3_analysis/package/analysis/inputs/training_positives_per_TF.csv
Figures/source_data/Figure_4_source_data.csv is not copied.
usage: python3 3_analysis/integration/make_inputs.py
"""
import os
import shutil

import pandas as pd

import p53drop as PD

SRC = os.path.join(PD.R2, "3_analysis", "package")
OUT = os.path.join(PD.B12, "inputs")
FILES = [("results_table.csv", "\t", "TF_name"), ("analysis/inputs/training_positives_per_TF.csv", "\t", "Names")]

outs, ins = [], []
for rel, sep, col in FILES:
    src, dst = os.path.join(SRC, rel), os.path.join(OUT, rel)
    PD.check_tfs(pd.read_csv(src, sep=sep)[col], rel)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    PD.backup(dst)
    shutil.copy2(src, dst)
    print("copied %s" % rel)
    outs.append(dst); ins.append(src)
PD.write_provenance(OUT, "3_analysis/integration/make_inputs.py", outs, ins, note="H3: byte copies of the integration package tables; ARCHI_PKG for B6")
