# 3_analysis

Analysis code. The `server/` folder holds the compute-heavy stages and is described in its own README. Every other folder
here runs locally: it reads the results of the compute stages in `results/server/` and the data in `results/`, and writes the
tables behind the figures, the supplementary tables and the numbers quoted in the text.

## Run

From the repository root, one command runs the whole local chain (analysis and visualisation) in dependency order:

```
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
ARCHI_STAGE=H5 \
nice -n 15 python3 3_analysis/integration/run_h3.py all
```

`run_h3.py` first checks every registered input (`preflight_inputs.py`: present, md5 equal to the manifest
`results/server/MIRROR_MANIFEST_md5.txt`, split labels of the GC table), unpacks the PWM archives of the Zenodo deposit
v1.0.0 (`zenodo_v1.py`, md5 checked), and then runs the stages one process at a time. It stops at the first failing
command and writes one log per command to `3_analysis/integration/logs/`. Single stages can be run by name, for example
`python3 3_analysis/integration/run_h3.py sup`. Stage order of `all` (the optional stage `gate` is listed after the table):

| stage | what it runs |
|---|---|
| pkg | `make_pkg_copy.py` (final server tables into `package/analysis/...`), then the package producers: results table, saturation, motif similarity, cross-species table, package Sup. Tables 3-6, figure source data |
| figs | `make_inputs.py`; Fig. 2 / S1, Fig. 3 A-C / S3, Fig. S2, S5, S6, Fig. 4 / S4 (in `4_visualisation/figures`); `fig4_tokens_h3b.py`, `fig4_paired_h4.py` |
| b234 | low FPR metrics (Fig. S8), top peaks (Fig. S9), feature importance (Fig. S10), ablation subtypes against HOCOMOCO v11 |
| b5 | transfer drivers and headroom (Fig. S11) |
| b9 | random split retraining (Fig. S12) |
| sup | Sup. Tables 1-10 |
| s7 | Fig. S7 and Fig. 3D source data, Fig. S7, Fig. 3 composite |
| porig, v1 | PWM origin table, mono vs mono+di, GC tables |
| fig1 | Fig. 1 values, Table 1 rows, Fig. 1 (SVG, PNG 300 dpi, PDF) |
| numbers | `outputs/numbers/NUMBERS_human_35TF*.csv`, `NUMBERS_R2_35TF.csv`, flip report |
| copy | every figure to `outputs/figures/` |

Optional stage `gate` (not part of `all`): check of the low FPR metrics against the 35 score files (`B2_gate_h3_*.csv` in
`low_fpr_metrics/` are tracked). The score files are not in git; they are in the Zenodo v2.0.0 deposit of this release (DOI in
the release notes) as `scores_H5_human.tar.gz`, which unpacks to `scores/<TF>.npz`. Unpack it and run

```
ARCHI_STAGE=H5 ARCHI_H5_SCORES=<folder>/scores python3 3_analysis/integration/run_h3.py gate
```

With `ARCHI_H5_SCORES` unset, or not naming a folder with `.npz` files, the stage stops with a message. The provenance records
of the score files in `low_fpr_metrics/_provenance_h3.yaml` are written `zenodo:scores_H5_human.tar.gz/scores/<TF>.npz`.

`ARCHI_STAGE=H5` selects the final data (training on all autosomes except the test chromosomes). The registry of every
input is `integration/h3_inputs.py`; the value `H4` selects the previous training run, whose inputs are not part of this
repository apart from the test-set tables in `results/server_h4_test_sets/` (its evaluation table is no longer tracked).

## Folders

| folder | contents | outputs (next to the scripts unless stated) |
|---|---|---|
| `integration/` | driver `run_h3.py`, input registry `h3_inputs.py`, pre-flight, shared helpers (`w2.py`: backups, md5, provenance; `p53drop.py`, `p53drop.R`: the 35 reference TFs), and the producers of the values used across tables and text | `top1000_check_h3.csv`, `figure1_values.csv`, `table1_rows.csv`, `sup6_methods_values_h3.csv`, `pwm_origin_h3.csv`, `mono_vs_monodi_h3.csv`, `gc_per_tf_h3.csv`, `gc_vs_gain_h3.csv`, `positive_count_ties_h3.csv`, `inputs/` (copies of two package tables for the human figure scripts); `outputs/numbers/` |
| `package/` | the analysis package of the first submission, patched for 35 TFs and the final server tables; inputs and outputs inside (`scripts/archi_paths.py`) | `results_table.csv`, `Sup_Tables/`, `Figures/source_data/`, `analysis/saturation/`, `analysis/cross_species/` |
| `low_fpr_metrics/` | partial auROC and TPR at low FPR | `B2_*.csv`, Fig. S8 (`B2_low_FPR_metrics.pdf/png`) |
| `top_peaks/` | top-N peaks, peak widths, gain against the number of training positives | `b3_*.csv`, Fig. S9 (`B3_peaks_figure.pdf/png`) |
| `random_split/` | random split retraining (5 TFs) | `b9_*.csv`, Fig. S12 (`B9_random_split_figure.pdf/png`) |
| `feature_importance/` | Gini importance, permutation, PWM subtypes, ablation (`scripts/`) | `tables/B4_*.csv`, Fig. S10 (`B4_importance_figure.pdf/png`) |
| `transfer_drivers/` | drivers of the human to mouse transfer and headroom | `B5_*.csv`, `B5_B_cloud_members.txt`, Fig. S11 (`B5_figure.pdf/png`) |
| `supplementary_tables/` | builders of the supplementary tables | `outputs/supplementary_tables/Sup_Table_*.xlsx` |

## Producers of the supplementary tables and numbers

| output | producer |
|---|---|
| Sup. Tables 1, 2, 3, 4 | `supplementary_tables/make_sup_tables_1_3_R2.py 1 2 3 4` |
| Sup. Table 5 | `supplementary_tables/make_sup_table_5_R2.py` |
| Sup. Table 6 | `supplementary_tables/make_sup_tables_1_3_R2.py 6` (calls `integration/sup_table_6_h3.py`) |
| Sup. Tables 7, 8, 9, 10 | `supplementary_tables/make_sup_tables_7_10.py 7 8 9 10` |
| `outputs/numbers/NUMBERS_human_35TF.csv`, `NUMBERS_human_35TF_unchanged.csv` | `integration/make_numbers.py` |
| `outputs/numbers/NUMBERS_R2_35TF.csv` | `integration/make_numbers_w2.py` |
| `outputs/numbers/flip_report_h5_vs_h4.csv` | `integration/flip_report.py --ref h4` |

The figure producers are listed in `4_visualisation/README.md`.

## Provenance

Every producer appends one block to a `_provenance_h3.yaml` (or `_provenance_p53drop.yaml`) in its output folder: script,
run time, git commit of this repository (or `not a git repository`, with a warning on screen), Python version, and the md5 of
every input and output. Before a producer overwrites an output it keeps one copy named `<stem>_backup_20261107<ext>`; these
copies are not tracked.

## Software

Python 3.9.12 with pandas 2.3.0, numpy 1.26.4, scipy 1.13.1, matplotlib 3.5.2, scikit-learn 1.1.2, openpyxl 3.1.5,
pypdf 6.9.1, PyYAML 6.0, Pillow 11.3.0. R 4.4.3 with dplyr 1.2.0, ggplot2 4.0.2, patchwork 1.3.2, cowplot 1.2.0,
ggrepel 0.9.8, ggExtra 0.11.0, viridis 0.6.5, tidyr 1.3.2, ggbeeswarm 0.7.3, ggseqlogo 0.2, png 0.1.9. Poppler 24.04
(`pdftoppm`, `pdfimages`, `pdftotext`), `rsvg-convert` 2.50.7, macOS `sips`, the Arial font. Tools are taken from `PATH`;
`PDFTOPPM`, `PDFIMAGES` and `PDFTOTEXT` override them.
