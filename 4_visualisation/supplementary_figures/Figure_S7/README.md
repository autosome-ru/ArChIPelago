# Figure S7

For the TFs with the largest gains over the best single monoPWM (human test set), logos of the PWMs with the highest Random Forest feature importance, with their HOCOMOCO v11 matches.

Producer: `build_Figure_S7_source_data_h3.py`, then `Figure_S7.R`; `hocomoco11_check/` holds the HOCOMOCO v11 comparisons (`subtype_vs_hocomoco11.py`, `subtype_vs_named_v11.py`, `cluster_vs_v11.py`).

Inputs (paths from the repository root): `3_analysis/package/results_table.csv`, Gini importances in `results/server/human/B4_importance/gini/`, PWM files from `results/datasets/zenodo_v1.0.0/`, `results/datasets/hocomoco11/`, `outputs/supplementary_tables/Sup_Table_2_PWMs_and_datasets.xlsx`.

Rebuilt by stage `s7` (`cluster_vs_v11.py`: stage `b234`) of `3_analysis/integration/run_h3.py`.

Files: `Figure_S7.R`, `Figure_S7.pdf`, `Figure_S7.png`, `Figure_S7_source_data.csv`, `_provenance_h3.yaml`, `build_Figure_S7_source_data_h3.py`, `hocomoco11_check/`, `matrices/`, `subtype_verdicts_v11_h3.csv`, `tf_ranking_h3.csv`.
