#!/bin/bash
# Rebuild the results table, Sup. Tables 1-6, the figure source data and every figure panel of the manuscript
# from the per-TF analysis outputs in ../analysis/ (no cluster steps). Run from any directory:
#   bash Manuscript_analysis/scripts/rebuild_tables_and_figures.sh
# Needs python3 (pandas, numpy, scipy, matplotlib, openpyxl) and Rscript (ggplot2, dplyr, tidyr, patchwork,
# cowplot, ggrepel, ggbeeswarm, ggseqlogo); the motif-similarity and Fig. S7 steps also need the Zenodo PWMs
# ($ARCHI_ZENODO_DIR, default <repository>/14927304 with PWMs_mono_HUMAN/, PWMs_di_HUMAN/, PWMs_mono_MOUSE.tar.gz).
# cross_species_table.py refuses to overwrite its table: delete analysis/cross_species/cross_species_table.csv first
# to rebuild it.
set -euo pipefail
PY=${PYTHON:-python3}
RS=${RSCRIPT:-Rscript}
S=$(cd "$(dirname "$0")" && pwd)
A=$S/../analysis

$PY $S/make_results_table.py                            # -> ../results_table.csv
$PY $A/saturation/analyze_saturation.py                 # -> analysis/saturation/ summaries
$PY $A/cross_species/motif_similarity.py                # -> analysis/cross_species/motif_similarity.csv
[ -e $A/cross_species/cross_species_table.csv ] || $PY $A/cross_species/cross_species_table.py
$PY $S/make_sup_tables_1_2.py                           # -> ../Sup_Tables/Sup_Table_1*, Sup_Table_2*
$PY $S/make_sup_tables_3_to_6.py                        # -> ../Sup_Tables/Sup_Table_3..6*, ../numbers_in_text.json
$PY $A/slim_dichipmunk/assemble_Figure_4_S4.py          # -> source data of Fig. 4, S4, ../numbers_in_text_Figure_4.json
$PY $A/motif_subtypes/build_Figure_S7_source_data.py    # -> analysis/motif_subtypes/ (Fig. S7 source data, matrices)
$PY $S/make_figure_source_data.py                       # -> source data of Fig. 2, S1, 3, S3, S5, S6, S7
$PY $S/Figure_S2_RF_vs_best_monoPWM.py                  # -> Fig. S2 panel and source data

for r in Figure_2_human_test Figure_S1_mouse_test Figure_3ABC_human_test Figure_S3ABC_mouse_test Figure_4_and_S4 \
         Figure_S5_saturation Figure_S6_cross_species Figure_S7_motif_subtypes; do
  $RS $S/$r.R
done
