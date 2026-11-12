"""Fig. 3D source data on the final data: for RXRA, SRF and ANDR the two monoPWMs and the two diPWMs
with the highest Gini importance in the human-trained RF (monoPWMs + diPWMs, seed 0), as the legend says. Same inputs and
logo conversion as Fig. S7 (4_visualisation/supplementary_figures/Figure_S7/build_Figure_S7_source_data_h3.py, imported: importances(), the package
to_prob / consensus / pcc_max); logos of a TF are oriented on the strand of its top monoPWM.
Outputs (this folder): Figure_3D_source_data.csv, matrices_3D/<TF>_<feature>.txt, _provenance_h3.yaml.
usage: python3 4_visualisation/main_figures/Figure_3/build_Figure_3D_source_data_h3.py
"""
import glob
import importlib.util
import os
import re
import shutil

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
R2 = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
_s = importlib.util.spec_from_file_location("s7h3", os.path.join(R2, "4_visualisation", "supplementary_figures", "Figure_S7", "build_Figure_S7_source_data_h3.py"))
S7 = importlib.util.module_from_spec(_s)
_s.loader.exec_module(S7)
P, w2 = S7.P, S7.w2
TFS, N_MONO, N_DI = ["RXRA", "SRF", "ANDR"], 2, 2
out = os.path.join(HERE, "Figure_3D_source_data.csv")
mdir = os.path.join(HERE, "matrices_3D")
w2.backup(out)
if os.path.isdir(mdir):
    if not os.path.exists(mdir + "_backup_20261107"):
        shutil.copytree(mdir, mdir + "_backup_20261107")
    shutil.rmtree(mdir)
os.makedirs(mdir)
rk = S7.ranking()   # gains from the package results table (asserted)
rows = []
for tf in TFS:
    g = S7.importances(tf)
    bad = sum(P.read_pwm(P.pwm_path(tf, f))[0] != h for f, h in zip(g.feature, g.pwm_header))
    if bad:
        raise SystemExit("CHECK FAILED %s: %d PWM headers differ from the gini table" % (tf, bad))
    sel = pd.concat([g[g.kind == "mono"].head(N_MONO), g[g.kind == "di"].head(N_DI)])
    probs = [(P.read_pwm(P.pwm_path(tf, f))) for f in sel.feature]
    probs = [(h, w, P.to_prob(f.split("_")[0], w)) for (h, w), f in zip(probs, sel.feature)]
    ref = P.zrows(probs[0][2])
    row = rk[rk.TF == tf].iloc[0]
    for j, ((_, r), (header, w, p)) in enumerate(zip(sel.iterrows(), probs), start=1):
        r_max, strand, _, _ = P.pcc_max(ref, P.zrows(p))
        if strand == "-":
            p = p[::-1, ::-1]
        mfile = os.path.join("matrices_3D", "%s_%s.txt" % (tf, r.feature))
        pd.DataFrame(p.T, index=list(P.ACGT), columns=["pos%d" % (i + 1) for i in range(len(p))]).to_csv(os.path.join(HERE, mfile), sep="\t")
        rows.append(dict(TF=tf, panel=j, pwm_class=r.kind, rank_overall=int(r["rank"]), feature=r.feature, pwm_name=header,
                         dataset=re.search(r"PEAKS\d+", header).group(0), importance=float(r.importance), n_features=len(g),
                         motif_length=len(p), strand_for_display=strand, r_max_vs_panel1=round(r_max, 3), consensus=P.consensus(p),
                         delta_auROC=float(row.delta_auROC), delta_auPRC=float(row.delta_auPRC), matrix_file=mfile))
src = pd.DataFrame(rows)
src["rank_within_class"] = src.groupby(["TF", "pwm_class"]).cumcount() + 1
src.to_csv(out, index=False)
print(src[["TF", "panel", "feature", "dataset", "importance", "strand_for_display", "r_max_vs_panel1", "consensus"]].to_string(index=False))
w2.provenance(HERE, "4_visualisation/main_figures/Figure_3/build_Figure_3D_source_data_h3.py", [out],
              [S7.IN["rt"], S7.IN["t2"]] + [os.path.join(S7.IN["gini"], "%s_gini.csv" % t) for t in TFS],
              params="{TFs: [%s], n_mono: %d, n_di: %d}" % (", ".join(TFS), N_MONO, N_DI))
