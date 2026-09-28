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

## Contents

| path | content |
|---|---|
| `HUMAN_MOUSE_total_100k.csv` | the results table, tab-separated, 900 rows x 47 columns: 36 TFs x 5 algorithms x 3 PWM sets plus the single-PWM rows; `roc_auc_*` / `pr_auc_*` on the human training set (`train_H`) and the human (`test_H`) and mouse (`test_M`) test sets; `*_PWM` = best single monoPWM (the baseline), `*_PWM_mono` / `*_PWM_di` = best single monoPWM / diPWM |
| `Sup_Tables/` | Supplementary Tables 1-6 (xlsx, and csv per sheet); `headline_numbers.json` and `fig4_refit_numbers.json` hold the numbers quoted in the text |
| `Figures/panels/` | the figure panels as drawn by the scripts (the figures of the manuscript were laid out from them in Adobe Illustrator) |
| `Figures/source_data/` | one csv per figure with exactly the plotted values |
| `scripts/` | the generators of the results table, the supplementary tables, the source data and the panels; `run_local.sh` runs them in order |
| `analysis/` | the per-TF analyses (run on a compute server from the pipeline output of notebooks 0-2) with their outputs, and `inputs/` |

### Figures and tables

| item | panel file (`Figures/panels/`) | drawn by | data |
|---|---|---|---|
| Fig. 2 | `Figure_2_panels_human_test.pdf`, `Figure_2_and_S1_legend.pdf` | `scripts/Figure_2_H_H.R` | `HUMAN_MOUSE_total_100k.csv` |
| Fig. S1 | `Figure_S1_panels_mouse_test.pdf` | `scripts/Figure_S1_H_M.R` | `HUMAN_MOUSE_total_100k.csv` |
| Fig. 3A-C | `Figure_3_panels_ABC_human_test.pdf` | `scripts/Figure_3_H_H.R` | `HUMAN_MOUSE_total_100k.csv` |
| Fig. S3 | `Figure_S3_panels_ABC_mouse_test.pdf` | `scripts/Figure_S3_H_M.R` | `HUMAN_MOUSE_total_100k.csv` |
| Fig. 4 / S4 | `Figure_4_human_test.pdf`, `Figure_S4_mouse_test.pdf` | `scripts/Figure_4_and_S4.R` | `analysis/fig4_refit/` via `assemble_fig4.py` |
| Fig. S2 | `Figure_S2.pdf` | `scripts/make_figures_S2_S5_S6.py` | `HUMAN_MOUSE_total_100k.csv` |
| Fig. S5 | `Figure_S5_saturation.pdf` | `scripts/Figure_S5_saturation.R` | `analysis/saturation/` |
| Fig. S6 | `Figure_S6_cross_species.pdf` | `scripts/Figure_S6_cross_species.R` | Sup. Table 5 |
| Fig. S7 | `Figure_S7_motif_subtypes.pdf` | `analysis/subtypes/Figure_S7_motif_subtypes.R` | `analysis/subtypes/build_source_data.py` |
| Sup. Tables 1, 2 | | `scripts/make_sup_tables_1_2.py` | `analysis/inputs/`, `analysis/gc/`, `analysis/hm_mm/`, `analysis/fig4_refit/rows/` |
| Sup. Tables 3-6 | | `scripts/make_sup_tables.py` | results table, `analysis/operational/`, `analysis/hm_mm/`, `analysis/crossspecies/`, `analysis/saturation/`, `analysis/runtime/` |

Fig. 1 (scheme) and the logo panel of Fig. 3D were drawn in Adobe Illustrator.

### `analysis/`

| folder | what | producer | output |
|---|---|---|---|
| `inputs/` | the results table written by notebook 4 (`HUMAN_MOUSE_total_100k_notebook4.csv`; its human values of the ArChIPelago models are those of the manuscript), the two Random Forest result tables of notebook 2 used for the human values of the single Slim and diChIPMunk models of Fig. 4, the GTRD metadata, the experiment and TF tables of the pipeline, `df_Names_seq_list.csv` (training positives per TF), `per_pwm_with_names.csv` (feature index -> PWM name) | pipeline | |
| `common/` | `archi_data.py`: loader of the per-TF feature matrices of the pipeline output | | |
| `results_table/` | all models and baselines evaluated on both test sets | `run_results_table.py`, `assemble_results_table.py` | `HUMAN_MOUSE_total_100k_recomputed.csv`, `best_pwm_names.csv` |
| `hm_mm/` | human-trained models on the mouse test set and mouse-trained control models (Random Forest, seeds 0 and 1) | `run_hm_mm.py` | `hm_mm_results.csv` (+ per-TF score vectors, not tracked) |
| `operational/` | false positives and precision at fixed recall, true sites in the top of the ranking (Sup. Table 4) | `operational_metrics.py` | `operational_metrics.csv` |
| `saturation/` | Random Forest on random and top-k subsets of k = 1, 2, 4, ..., 128 PWMs (Fig. S5) | `run_saturation.py`, `analyze_saturation.py` | `saturation_results.csv`, `per_tf/`, summaries |
| `fig4_refit/` | the RF2f family of Fig. 4 / S4 (best monoPWM + best diPWM, alone and with Slim / diChIPMunk features) and the single Slim and diChIPMunk models on the mouse test set | `scan_one.sh`, `munk_one.sh`, `train_slim_ANDR_m1.sh`, `run_fig4_refit.py`, `assemble_fig4.py` | `rows/<TF>.json`, `fig4_refit_table.csv`, `slim_ANDR_m1/SlimDimont_1.xml` |
| `crossspecies/` | per-TF cross-species table and the similarity of the baseline human monoPWM to the mouse monoPWMs (Fig. S6, Sup. Table 5) | `motif_similarity.py`, `crossspecies_analysis.py` | `motif_similarity.csv`, `crossspecies_table.csv` |
| `gc/` | GC content of every positive and negative set, pooled peak counts (Sup. Table 2) | `gc_content.py`, `peaks_raw_check.py`, `peaks_raw_check_mouse_allmodels.py` | `gc_content.csv`, `peaks_raw_check*.csv` |
| `runtime/` | raw timings of SPRY-SARUS scanning and Random Forest training (Sup. Table 6) | | `sarus_timings.tsv`, `rf_process_timings.tsv` |
| `subtypes/` | Random Forest feature importances and probability matrices of the top-ranked PWMs (Fig. S7) | `build_source_data.py` | `Figure_S7_source_data.csv`, `matrices/`, `importances/`, `tf_ranking.csv` |

## Rebuilding the tables and figures

```bash
bash Manuscript_analysis/scripts/run_local.sh
```

rebuilds the results table, Sup. Tables 1-6, the source data and all panels from the outputs in `analysis/`
(a few minutes). Requirements: Python 3.9+ with pandas, numpy, scipy, matplotlib and openpyxl; R 4.4 with
ggplot2, dplyr, tidyr, patchwork, cowplot, ggrepel, ggbeeswarm and ggseqlogo; the PWM files of the Zenodo
archive for `motif_similarity.py` and `build_source_data.py` (`ARCHI_ZENODO_DIR`, default `<repository>/14927304`
with `PWMs_mono_HUMAN/`, `PWMs_di_HUMAN/` and `PWMs_mono_MOUSE.tar.gz`). With pandas 2.3, numpy 1.26 and
R 4.4.3 the run reproduces every csv and json file of this folder byte for byte and every panel except the
random horizontal jitter of the dots in Fig. 2 and S1.

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

# results table: 5 algorithms x 3 PWM sets, baselines, both test sets
python results_table/run_results_table.py --workers 4 --threads-per-worker 4 --rows-dir results_table/rows
python results_table/assemble_results_table.py --rows-dir results_table/rows \
       --published inputs/HUMAN_MOUSE_total_100k_notebook4.csv --out-dir results_table

# human-trained models on the mouse test set, mouse-trained control; score vectors for the operational metrics
python hm_mm/run_hm_mm.py --tfs $TFS --workers 2 --rf-jobs 4 --out-csv hm_mm/hm_mm_results.csv \
       --rows-dir hm_mm/rows --scores-dir hm_mm/scores
python operational/operational_metrics.py hm_mm/scores operational/operational_metrics.csv

# saturation curves (Fig. S5)
python saturation/run_saturation.py --workers 2 --rf-jobs 4

# Fig. 4 / S4: Slim (m = 0, 1, -5) and diChIPMunk scans of the mouse test set, then the RF2f family
for tf in $TFS; do for m in 0 1 5; do fig4_refit/scan_one.sh $tf $m MOUSE_10000_control; done; fig4_refit/munk_one.sh $tf; done
python fig4_refit/run_fig4_refit.py --tfs $TFS --rows-dir fig4_refit/rows

# GC content of every set (Sup. Table 2)
python gc/gc_content.py
```

`train_slim_ANDR_m1.sh` trains the Slim m=1 model of ANDR (`fig4_refit/slim_ANDR_m1/SlimDimont_1.xml`), which is
then scanned with `scan_one.sh ANDR 1 <set> fig4_refit/slim_ANDR_m1/SlimDimont_1.xml` for the three sequence sets.
The scripts never overwrite existing outputs; a TF whose output exists is skipped, so runs can be resumed.
