# Figure 3

A-C: Random Forest against the best single monoPWM on the human test set (A auROC, B auPRC, C gain in auROC against gain in auPRC per TF, bubble size the number of training positives, fill colour the number of PWMs); D: logos of the two monoPWMs and two diPWMs with the highest Gini importance for RXRA, SRF and ANDR.

Producer: `Figure_3_and_S3.R` (A-C; also writes Fig. S3 to `../../supplementary_figures/Figure_S3/`), `build_Figure_3D_source_data_h3.py` and `Figure_3D_logos.R` (D), `assemble_figure3.py` (vector composite `Figure_3.pdf`).

Inputs (paths from the repository root): `3_analysis/integration/inputs/results_table.csv`, Gini importances in `results/server/human/B4_importance/gini/`, PWM files from `results/datasets/zenodo_v1.0.0/`, `outputs/supplementary_tables/Sup_Table_2_PWMs_and_datasets.xlsx`.

Rebuilt by stages `figs` (A-C) and `s7` (D, composite) of `3_analysis/integration/run_h3.py`.

Files: `Figure_3.pdf`, `Figure_3.png`, `Figure_3D_logos.R`, `Figure_3D_logos.pdf`, `Figure_3D_logos.png`, `Figure_3D_source_data.csv`, `Figure_3_ABC_restyled.pdf`, `Figure_3_ABC_restyled.png`, `Figure_3_AB_plotdata.csv`, `Figure_3_C_plotdata.csv`, `Figure_3_and_S3.R`, `_provenance_h3.yaml`, `assemble_figure3.py`, `build_Figure_3D_source_data_h3.py`, `matrices_3D/`.
