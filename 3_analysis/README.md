# 3_analysis

Analysis code. `server/` holds the compute-heavy stages and has its own README. Every other folder here runs locally: it
reads the result tables in `../results/` and writes the tables behind the figures, the supplementary tables and the numbers
of the text.

## Run

From the repository root, one command runs the local chain (analysis and figures) in dependency order:

```
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
ARCHI_STAGE=H5 \
nice -n 15 python3 3_analysis/integration/run_h3.py all
```

`run_h3.py` unpacks the PWM archives in `results/datasets/zenodo_v1.0.0/` (`integration/zenodo_v1.py`, md5 checked),
checks every registered input (`integration/preflight_inputs.py`: present, md5 equal to
`results/server/MIRROR_MANIFEST_md5.txt`), then runs the stages one process at a time. It stops at the first failing
command and writes one log per command to `3_analysis/integration/logs/`. Single stages run by name, for example
`python3 3_analysis/integration/run_h3.py sup`. `ARCHI_STAGE=H5` selects the inputs of the final run; the registry of
every input is `integration/h3_inputs.py`.

Stages of `all`, in order:

| stage | what it runs |
|---|---|
| `pkg` | `integration/make_pkg_copy.py` (server tables into `package/`), then the package producers in `package/`: results table, saturation, motif similarity, cross-species table, package Sup. Tables 3-6, figure source data |
| `figs` | `integration/make_inputs.py`; Fig. 2 and S1, Fig. 3 A-C and S3, Fig. S2, S5, S6, Fig. 4 and S4 (producers in `4_visualisation/`); `integration/fig4_tokens_h3b.py`, `integration/fig4_paired_h4.py` |
| `b234` | low FPR metrics and Fig. S8, top peaks and Fig. S9, feature importance and Fig. S10, ablation subtypes against HOCOMOCO v11 |
| `b5` | transfer drivers and headroom, Fig. S11 |
| `b9` | random split retraining, Fig. S12 |
| `sup` | `integration/top1000_check_h3.py`, Sup. Tables 1-10 |
| `s7` | Fig. S7 source data and figure, Fig. 3D source data and logos, Fig. 3 composite |
| `porig` | `integration/pwm_origin_h3.py` |
| `v1` | `integration/v1_tables_h3.py` (mono against mono+di, GC tables, positive-count ties) |
| `fig1` | `integration/figure1_values.py` (`figure1_values.csv`, the numbers of the hand-drawn Fig. 1), `integration/table1_inputs.py` |
| `numbers` | `integration/make_numbers.py`, `integration/make_numbers_w2.py`, `integration/flip_report.py --ref h4` |

Optional stage `gate` (not part of `all`): checks the low FPR metrics against the 35 test-set score files. Input:
`scores_H5_human.tar.gz` of the Zenodo record named in `../README.md`; it unpacks to `scores/<TF>.npz`. Unpack it and run

```
ARCHI_STAGE=H5 ARCHI_H5_SCORES=<folder>/scores python3 3_analysis/integration/run_h3.py gate
```

Without `ARCHI_H5_SCORES`, or with a folder that holds no `.npz` files, the stage stops with a message. Its outputs
`low_fpr_metrics/B2_gate_h3_per_tf.csv` and `low_fpr_metrics/B2_gate_h3_summary.csv` are tracked.

## Folders

| folder | contents | outputs (next to the scripts unless stated) |
|---|---|---|
| `integration/` | driver `run_h3.py`, input registry `h3_inputs.py`, pre-flight, shared helpers (`w2.py`: backups, md5, provenance; `p53drop.py`, `p53drop.R`: the 35 reference TFs), producers of the values used across tables and text | `top1000_check_h3.csv`, `figure1_values.csv`, `table1_rows.csv`, `sup6_methods_values_h3.csv`, `pwm_origin_h3.csv`, `mono_vs_monodi_h3.csv`, `gc_per_tf_h3.csv`, `gc_vs_gain_h3.csv`, `positive_count_ties_h3.csv`, `inputs/` (tables for the human figure scripts); `../outputs/numbers/` |
| `package/` | analysis package: results table, saturation, cross-species and motif similarity, package supplementary tables, figure source data; inputs in `package/analysis/inputs/`, paths in `package/scripts/archi_paths.py` | `results_table.csv`, `Sup_Tables/`, `Figures/source_data/`, `analysis/` |
| `low_fpr_metrics/` | partial auROC and TPR at low FPR | `B2_*.csv` |
| `top_peaks/` | top-N peaks, peak widths, gain against the number of training positives | `b3_*.csv` |
| `random_split/` | random split retraining (5 TFs) | `b9_*.csv` |
| `feature_importance/` | Gini importance, permutation, PWM subtypes, ablation (`scripts/b4_05_aggregate.py`) | `tables/B4_*.csv` |
| `transfer_drivers/` | drivers of the human to mouse transfer and headroom (`B5_01_drivers_and_ceiling.py`) | `B5_*.csv`, `B5_B_cloud_members.txt` |
| `supplementary_tables/` | builders of the supplementary tables | `../outputs/supplementary_tables/Sup_Table_*.xlsx` |

The figures drawn from these tables are in `../4_visualisation/` (Fig. S8 to S12 in `../4_visualisation/supplementary_figures/`).

## Supplementary tables and numbers

| output | producer |
|---|---|
| Sup. Tables 1, 2, 3, 4 and 6 | `supplementary_tables/make_sup_tables_1_3_R2.py 1 2 3 4 6` (Sup. Table 6 calls `integration/sup_table_6_h3.py`) |
| Sup. Table 5 | `supplementary_tables/make_sup_table_5_R2.py` |
| Sup. Tables 7, 8, 9, 10 | `supplementary_tables/make_sup_tables_7_10.py 7 8 9 10` |
| `outputs/numbers/NUMBERS_human_35TF.csv`, `outputs/numbers/NUMBERS_human_35TF_unchanged.csv` | `integration/make_numbers.py` |
| `outputs/numbers/NUMBERS_R2_35TF.csv` | `integration/make_numbers_w2.py` |
| `outputs/numbers/flip_report_h5_vs_h4.csv` | `integration/flip_report.py --ref h4` |

## Provenance

Python producers of tables and figures append one block to a `_provenance_h3.yaml` (or `_provenance_p53drop.yaml`) in their output folder (in the analysis folder for `Figure_S8.py` and `Figure_S10.py`); `transfer_drivers/B5_01_drivers_and_ceiling.py` and `server/h5_provenance.py` write a `_provenance.yaml`. The R figure producers and `Figure_S11.py` write no record. A block holds
script, run time, git commit of the repository, Python version, and the md5 of every input and output. Before a producer
overwrites an output it keeps one copy named `<stem>_backup_<tag><ext>`; these copies are not tracked.

## Software

Python 3.9.12 with pandas 2.3.0, numpy 1.26.4, scipy 1.13.1, matplotlib 3.5.2, scikit-learn 1.1.2, openpyxl 3.1.5,
pypdf 6.9.1, PyYAML 6.0, Pillow 11.3.0. R 4.4.3 with PRROC 1.4 and rpy2 3.5.17 for auROC and auPRC (rpy2 must load the
R that holds PRROC, set `R_HOME`), and the R packages listed in `../4_visualisation/README.md`. Poppler 24.04
(`pdftoppm`, `pdfimages`, `pdftotext`); the environment variables `PDFTOPPM`, `PDFIMAGES` and `PDFTOTEXT` override the
tools on `PATH`.
