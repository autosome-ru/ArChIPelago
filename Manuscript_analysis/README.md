# Manuscript_analysis

Code and outputs behind every number, supplementary table and figure panel of the manuscript
(Kravchenko et al., *Classic machine learning on top of multiple position weight matrices improves genomic
prediction of transcription factor binding sites*).

## Basis of every number

| | human | mouse |
|---|---|---|
| test set | chromosomes 1, 8, 21 | chromosomes 1, 8, 19 |
| training set | chromosomes 2-7, 9, 10, 13-20 | chromosomes 2-7, 9, 10, 13-18 |

* Chromosomes 11, 12 and the sex chromosomes are used in neither set.
* The ArChIPelago models are trained on the human training set (five algorithms: Random Forest, Logistic
  Regression, XGBoost, Bagging with XGBoost, Bagging with Logistic Regression; three PWM sets: monoPWMs,
  diPWMs, monoPWMs + diPWMs) and evaluated on the human and on the mouse test set.
* Baseline: the best single monoPWM of each TF, selected on the human training set, independently by
  auROC and by auPRC. The mouse-trained control models (Fig. S6, Sup. Table 5) use the best mouse
  monoPWM selected on the mouse training set.
* auROC and auPRC as computed by the PRROC R package (auPRC = `pr.curve()$auc.integral`).

## Layout

```
Manuscript_analysis/
├── results_table.csv                  results of all models and baselines (36 TFs x 5 algorithms x 3 PWM sets)
├── numbers_in_text.json               every number quoted in the text
├── numbers_in_text_Figure_4.json      the numbers of Fig. 4 / S4 quoted in the text
├── Sup_Tables/                        Supplementary Tables 1-6 (xlsx, and one csv per sheet)
├── Figures/
│   ├── panels/                        figure panels as drawn by the scripts
│   └── source_data/                   one csv per figure with exactly the plotted values
├── scripts/                           builders of the tables, source data and panels (run on a laptop)
└── analysis/                          per-TF analyses (run on a compute server) and their outputs
```

`results_table.csv` is tab-separated, 900 rows x 47 columns: one row per TF x algorithm x PWM set plus the
single-PWM rows; `roc_auc_*` / `pr_auc_*` on the human training set (`train_H`) and the human (`test_H`) and
mouse (`test_M`) test sets; `*_PWM` = best single monoPWM (the baseline), `*_PWM_mono` / `*_PWM_di` = best
single monoPWM / diPWM.

The figures of the manuscript were laid out in Adobe Illustrator from the panels. Fig. 1 (scheme) and the logo
panel of Fig. 3D were drawn in Illustrator.

## Figures and tables

| item | drawn / written by (`scripts/`) | panel (`Figures/panels/`) | data |
|---|---|---|---|
| Fig. 2 | `Figure_2_human_test.R` | `Figure_2_human_test.pdf`, `Figure_2_and_S1_legend.pdf` | `results_table.csv` |
| Fig. S1 | `Figure_S1_mouse_test.R` | `Figure_S1_mouse_test.pdf` | `results_table.csv` |
| Fig. S2 | `Figure_S2_RF_vs_best_monoPWM.py` | `Figure_S2_RF_vs_best_monoPWM.pdf` | `results_table.csv` |
| Fig. 3A-C | `Figure_3ABC_human_test.R` | `Figure_3ABC_human_test.pdf` | `results_table.csv` |
| Fig. S3 | `Figure_S3ABC_mouse_test.R` | `Figure_S3ABC_mouse_test.pdf` | `results_table.csv` |
| Fig. 4, S4 | `Figure_4_and_S4.R` | `Figure_4_human_test.pdf`, `Figure_S4_mouse_test.pdf` | `analysis/slim_dichipmunk/` |
| Fig. S5 | `Figure_S5_saturation.R` | `Figure_S5_saturation.pdf` | `analysis/saturation/` |
| Fig. S6 | `Figure_S6_cross_species.R` | `Figure_S6_cross_species.pdf` | Sup. Table 5 |
| Fig. S7 | `Figure_S7_motif_subtypes.R` | `Figure_S7_motif_subtypes.pdf` | `analysis/motif_subtypes/` |
| Sup. Tables 1, 2 | `make_sup_tables_1_2.py` | | `analysis/inputs/`, `analysis/datasets/`, `analysis/mouse_transfer/`, `analysis/slim_dichipmunk/rows/` |
| Sup. Tables 3-6, `numbers_in_text.json` | `make_sup_tables_3_to_6.py` | | `results_table.csv`, `analysis/operational/`, `analysis/mouse_transfer/`, `analysis/cross_species/`, `analysis/saturation/`, `analysis/runtime/` |

The other scripts of `scripts/`: `make_results_table.py` writes `results_table.csv`; `make_figure_source_data.py`
writes the source data of Fig. 2, S1, 3, S3, S5, S6 and S7 (Fig. 4 / S4 by
`analysis/slim_dichipmunk/assemble_Figure_4_S4.py`, Fig. S2 by its drawing script); `archi_paths.py` holds the
shared locations; `rebuild_tables_and_figures.sh` runs everything in order.

## `analysis/`

| folder | what | producer | output |
|---|---|---|---|
| `inputs/` | tables written by the pipeline notebooks and the GTRD metadata (see below) | pipeline | |
| `common/` | `pipeline_data.py`: loader of the per-TF feature matrices of the pipeline output | | |
| `model_evaluation/` | all ArChIPelago models and single-PWM baselines evaluated on both test sets | `evaluate_models.py`, `assemble_evaluation_table.py` | `evaluation_table.csv`, `best_single_PWMs.csv` |
| `mouse_transfer/` | human-trained Random Forest on the mouse test set and mouse-trained control models (seeds 0 and 1) | `run_mouse_transfer.py` | `mouse_transfer_results.csv` (+ per-TF score vectors, not tracked) |
| `operational/` | false positives and precision at fixed recall, true sites in the top of the ranking (Sup. Table 4) | `operational_metrics.py` | `operational_metrics.csv` |
| `saturation/` | Random Forest on random and top-k subsets of k = 1, 2, 4, ..., 128 PWMs (Fig. S5) | `run_saturation.py`, `analyze_saturation.py` | `saturation_results.csv`, `per_tf/`, `saturation_summary_by_k.csv`, `saturation_per_tf.csv` |
| `slim_dichipmunk/` | Fig. 4 / S4: Random Forest on the best monoPWM + best diPWM (RF2f), alone and with Slim / diChIPMunk features, and the single Slim and diChIPMunk models | `scan_slim.sh`, `scan_dichipmunk.sh`, `train_slim_ANDR_m1.sh`, `fit_RF2f_models.py`, `assemble_Figure_4_S4.py` | `rows/<TF>.json`, `RF2f_Slim_diChIPMunk_table.csv`, `slim_ANDR_m1/SlimDimont_1.xml` |
| `cross_species/` | per-TF cross-species table; similarity of the baseline human monoPWM to the mouse monoPWMs (Fig. S6, Sup. Table 5) | `motif_similarity.py`, `cross_species_table.py` | `motif_similarity.csv`, `cross_species_table.csv` |
| `datasets/` | GC content of every positive and negative set, pooled peak sets per TF (Sup. Table 2) | `gc_content.py`, `pooled_peaks_human.py`, `pooled_peaks_mouse.py` | `gc_content.csv`, `pooled_peaks_human.csv`, `pooled_peaks_mouse.csv` |
| `runtime/` | raw timings of SPRY-SARUS scanning and Random Forest training (Sup. Table 6) | | `sarus_timings.tsv`, `rf_process_timings.tsv` |
| `motif_subtypes/` | Random Forest feature importances and probability matrices of the top-ranked PWMs (Fig. S7) | `build_Figure_S7_source_data.py` | `Figure_S7_source_data.csv`, `tf_ranking.csv`, `importances/`, `matrices/` |

`analysis/inputs/`:

| file | content |
|---|---|
| `notebook4_results_table.csv` | the results table written by notebook 4 (`HUMAN_MOUSE_total_100k.csv`); source of the human values of the ArChIPelago models |
| `notebook2_Slim_diChIPMunk_RF2f_auROC.txt`, `..._auPRC.txt` | the Slim / diChIPMunk result tables of notebook 2 (human values of the single Slim and diChIPMunk models of Fig. 4) |
| `GTRD_experiments.xlsx` | the ChIP-Seq experiments of every TF (Sup. Table 1) |
| `GTRD_metadata.csv` | GTRD metadata used by notebook 0 |
| `TF_table.xlsx` | TF table of the pipeline (numbers of experiments and PWMs) |
| `training_positives_per_TF.csv` | number of human training positives per TF |
| `single_PWM_features.csv` | feature index -> PWM name, with the single-PWM auROC / auPRC |

`results_table.csv` takes the human values of the ArChIPelago models and the human best-monoPWM baselines from
`notebook4_results_table.csv`, and the mouse columns and the best-diPWM baselines from
`model_evaluation/evaluation_table.csv` (`make_results_table.py`).

## Rebuilding the tables and figures

```bash
bash Manuscript_analysis/scripts/rebuild_tables_and_figures.sh
```

rebuilds the results table, Sup. Tables 1-6, the numbers of the text, the source data and all panels from the
outputs in `analysis/` (a few minutes). Requirements: Python 3.9+ with pandas, numpy, scipy, matplotlib and
openpyxl; R 4.4 with ggplot2, dplyr, tidyr, patchwork, cowplot, ggrepel, ggbeeswarm and ggseqlogo; the PWM
files of the Zenodo archive for `motif_similarity.py` and `build_Figure_S7_source_data.py` (`ARCHI_ZENODO_DIR`,
default `<repository>/14927304`, with `PWMs_mono_HUMAN/`, `PWMs_di_HUMAN/` and `PWMs_mono_MOUSE.tar.gz`).
`cross_species_table.py` does not overwrite its table; delete `analysis/cross_species/cross_species_table.csv`
to rebuild it. With pandas 2.3, numpy 1.26 and R 4.4.3 the run reproduces every csv and json file of this
folder byte for byte and every panel except the random horizontal jitter of the dots in Fig. 2 and S1.

## Running the per-TF analyses

These steps read the pipeline output directory written by notebooks 0-2 (`ARCHI_RELEASE_DIR`, default
`~/Release/TF-ML`: the global sequence tables `out_tab_<HUMAN|MOUSE>_10000_<train|control>.tab`, the fasta
files `all_mfa_file_*`, and `outputdir/<TF>/` with the per-PWM score files, negative ids, Slim models and the
diChIPMunk diPWMs). "train" and "control" name the training and the test chromosomes. In the pipeline
directory of the manuscript run the two mouse sets carry each other's names; for that directory set
`ARCHI_MOUSE_FILES_SWAPPED=1`. The loader checks the chromosomes of every table it reads and stops if a split
does not contain the chromosomes listed above.

Python 3.8, scikit-learn 1.3.0, xgboost, joblib; run with `nice -n 10`. From `Manuscript_analysis/analysis/`:

```bash
TFS="ANDR AP2A CEBPB COE1 CTCF E2F4 ERG ESR1 FLI1 GATA1 GATA2 GATA3 GCR HNF4A IRF1 IRF4 JUND MAFK MAX MYC P53 PPARG PRGR REST RUNX1 RXRA SOX2 SPI1 SRF STA5A STAT1 STAT3 TAL1 TF65 TFE2 USF2"

# all models and baselines on both test sets (5 algorithms x 3 PWM sets)
python model_evaluation/evaluate_models.py --workers 4 --threads-per-worker 4 --rows-dir model_evaluation/rows
python model_evaluation/assemble_evaluation_table.py --rows-dir model_evaluation/rows \
       --notebook4-table inputs/notebook4_results_table.csv --out-dir model_evaluation

# human-trained models on the mouse test set, mouse-trained control; score vectors for the operational metrics
python mouse_transfer/run_mouse_transfer.py --tfs $TFS --workers 2 --rf-jobs 4 \
       --out-csv mouse_transfer/mouse_transfer_results.csv --rows-dir mouse_transfer/rows --scores-dir mouse_transfer/scores
python operational/operational_metrics.py mouse_transfer/scores operational/operational_metrics.csv

# saturation curves (Fig. S5)
python saturation/run_saturation.py --workers 2 --rf-jobs 4

# Fig. 4 / S4: Slim (m = 0, 1, -5) and diChIPMunk scans of the mouse test set, then the RF2f models
for tf in $TFS; do
  for m in 0 1 5; do slim_dichipmunk/scan_slim.sh $tf $m MOUSE_10000_control; done
  slim_dichipmunk/scan_dichipmunk.sh $tf
done
python slim_dichipmunk/fit_RF2f_models.py --tfs $TFS --rows-dir slim_dichipmunk/rows

# datasets (Sup. Table 2): GC content of every set, pooled peak sets
python datasets/gc_content.py
python datasets/pooled_peaks_human.py
python datasets/pooled_peaks_mouse.py

# motif similarity and the cross-species table (Fig. S6, Sup. Table 5)
python cross_species/motif_similarity.py
python cross_species/cross_species_table.py
```

`train_slim_ANDR_m1.sh` trains the Slim m=1 model of ANDR (`slim_dichipmunk/slim_ANDR_m1/SlimDimont_1.xml`),
which is then scanned with `scan_slim.sh ANDR 1 <set> slim_dichipmunk/slim_ANDR_m1/SlimDimont_1.xml` for the three
sequence sets. The scripts never overwrite existing outputs; a TF whose output exists is skipped, so runs can be
resumed. `pooled_peaks_*.py` read the MACS peak files and the HOCOMOCO v11 curation tables from the home
directory (`~/macs`, `~/hocomoco11`, `~/Cleaning_data`).
