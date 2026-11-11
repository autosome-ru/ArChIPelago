"""Copy every rebuilt figure to outputs/figures (pdf, png; Fig. 1 also svg), backing up each target first. Fig. 3 is written
into outputs/figures by 4_visualisation/figures/Fig3_S3/assemble_figure3.py (run_h3.py stage s7). Fig. 4 and S4 come from
4_visualisation/figures/output/Fig4. Fig. S7 depends on the training data (gain ranking, Gini importances of the models) and
is rebuilt by run_h3.py stage s7.
usage: ARCHI_STAGE=H5 python3 3_analysis/integration/copy_figures.py"""
import os, shutil, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w2
R2 = w2.R2
MAP = {"Figure_1": "4_visualisation/figure1/Figure_1_R2", "Figure_2": "4_visualisation/figures/Fig2_S1/Figure_2_restyled",
       "Figure_S1": "4_visualisation/figures/output/S1/Figure_S1_restyled", "Figure_S2": "4_visualisation/figures/output/S2/Figure_S2_restyled",
       "Figure_S3": "4_visualisation/figures/output/S3/Figure_S3_ABC_restyled", "Figure_S5": "4_visualisation/figures/output/S5/Figure_S5_restyled",
       "Figure_S6": "4_visualisation/figures/output/S6/Figure_S6_restyled", "Figure_S8": "3_analysis/low_fpr_metrics/B2_low_FPR_metrics",
       "Figure_S9": "3_analysis/top_peaks/B3_peaks_figure", "Figure_S10": "3_analysis/feature_importance/B4_importance_figure",
       "Figure_S11": "3_analysis/transfer_drivers/B5_figure", "Figure_S7": "4_visualisation/figures/FigS7/Figure_S7_motif_subtypes",
       "Figure_S12": "3_analysis/random_split/B9_random_split_figure",
       "Figure_4": "4_visualisation/figures/output/Fig4/Figure_4_restyled", "Figure_S4": "4_visualisation/figures/output/Fig4/Figure_S4_restyled"}
ins, outs = [], []
for tgt, src in MAP.items():
    for ext in ("pdf", "png") + (("svg",) if tgt == "Figure_1" else ()):
        s, t = os.path.join(R2, src + "." + ext), os.path.join(R2, "outputs/figures", tgt + "." + ext)
        w2.backup(t); shutil.copy2(s, t); ins.append(s); outs.append(t)
        assert w2.md5(s) == w2.md5(t)
outs += [os.path.join(R2, "outputs/figures", "Figure_3." + e) for e in ("pdf", "png")]
ins += [os.path.join(R2, "4_visualisation/figures", "Fig3_S3", "Figure_3_ABC_restyled.png")]
print("copied", len(outs) - 2)
w2.provenance(os.path.join(R2, "outputs/figures"), "3_analysis/integration/copy_figures.py", outs, ins,
              note="H3: Fig. 1, 2, 3 (assembled), S1-S3, S5, S6, S7, S8-S11, S12 rebuilt; Fig. 4, S4 pending (H3b)")
