#!/bin/bash
# Rebuild the results table, Sup. Tables 1-6, the figure source data and every figure panel of the manuscript
# from the per-TF analysis outputs in ../analysis/ (no cluster steps). Run from any directory:
#   bash Manuscript_analysis/scripts/run_local.sh
# Needs python3 (pandas, numpy, scipy, matplotlib, openpyxl) and Rscript (ggplot2, dplyr, tidyr, patchwork,
# cowplot, ggrepel, ggbeeswarm, ggseqlogo); the motif-similarity and Fig. S7 steps also need the Zenodo PWMs
# ($ARCHI_ZENODO_DIR, default <repository>/14927304 with PWMs_mono_HUMAN/, PWMs_di_HUMAN/, PWMs_mono_MOUSE.tar.gz).
# crossspecies_analysis.py refuses to overwrite its table: delete analysis/crossspecies/crossspecies_table.csv first
# to rebuild it.
set -euo pipefail
PY=${PYTHON:-python3}
RS=${RSCRIPT:-Rscript}
S=$(cd "$(dirname "$0")" && pwd)
A=$S/../analysis

$PY $S/make_results_table.py                       # -> ../HUMAN_MOUSE_total_100k.csv
$PY $A/saturation/analyze_saturation.py            # -> analysis/saturation/ summaries
$PY $A/crossspecies/motif_similarity.py            # -> analysis/crossspecies/motif_similarity.csv
[ -e $A/crossspecies/crossspecies_table.csv ] || $PY $A/crossspecies/crossspecies_analysis.py
$PY $S/make_sup_tables_1_2.py                      # -> ../Sup_Tables/Sup_Table_1*, Sup_Table_2*
$PY $S/make_sup_tables.py                          # -> ../Sup_Tables/Sup_Table_3..6*, headline_numbers.json
$PY $S/make_figures_S2_S5_S6.py                    # -> Figure_S2.pdf, source data of S2, S5, S6
$PY $S/make_figure_source_data.py                  # -> source data of Fig. 2, S1, 3, S3
$PY $A/fig4_refit/assemble_fig4.py                 # -> source data of Fig. 4, S4, fig4_refit_numbers.json
$PY $A/subtypes/build_source_data.py               # -> analysis/subtypes/ (Fig. S7 source data, matrices)
$PY -c "import pandas as pd; d = pd.read_csv('$A/subtypes/Figure_S7_source_data.csv'); \
d[d.TF.isin(['E2F4', 'RXRA', 'TAL1', 'TFE2'])].to_csv('$S/../Figures/source_data/Figure_S7_source_data.csv', index=False)"

for r in Figure_2_H_H Figure_S1_H_M Figure_3_H_H Figure_S3_H_M Figure_4_and_S4 Figure_S5_saturation Figure_S6_cross_species; do
  $RS $S/$r.R
done
$RS $A/subtypes/Figure_S7_motif_subtypes.R
