# 1_data_preparation

Builds the training and test data of the models for human and mouse: positive windows, GC-matched negatives, SARUS
feature tables, Slim inputs, and the data checks. CPU only. Later stages import the loader `archi_data_v5.py` from this
folder.

## Files

| file | role |
|---|---|
| `h5_data_server.sh` | driver of the whole build |
| `build_positives_v5.py` | positive windows per TF and pooled |
| `make_negatives_v3.py` | GC-matched negatives (BiasAway) |
| `build_slim_inputs_v5.py` | Slim inputs of the 35 TFs |
| `make_scan_jobs_sym.py` | SARUS scan jobs (feature tables) |
| `check_h5_test_identity.py`, `check_h5_data.py`, `check_C_sym.py`, `check_E.py`, `check_neg_v3_h5.py`, `gc_check_v3_h5.py`, `counts_D_human_h5.py`, `counts_D_mouse_h5.py` | data checks and counts |
| `h5_splits.py` | the chromosome split: human training chr2-7, 9-20, 22, test chr1, 8, 21; mouse training chr2-7, 9-18, test chr1, 8, 19 |
| `archi_data_v5.py` | loader of the data builds |
| `prroc_metrics.py` | auROC and auPRC by R PRROC through rpy2 |

## Environment

Paths outside the repository come from environment variables, read in `../archi_env.py` (Python) and `../archi_env.sh`
(drivers); a variable that is needed and not set stops the script with a message. `ARCHI_WORK` is the working root that
receives `human_sym_h5/` and `mouse_sym_h5/`, `ARCHI_RELEASE` the pipeline output folder `TF-ML`, `ARCHI_PY` the Python
interpreter of the drivers (default `python3`). The inputs are found through `ARCHI_MACS` (MACS peak files,
`macs.tar.gz` on Zenodo), `ARCHI_GENOMES` (hg38 and GRCm38.p6 fasta), `ARCHI_FAMILY` (folder with
`core_annotation_HUMAN_mono.csv`), `ARCHI_HOCOMOCO` (HOCOMOCO v11 models, `hocomoco11.tar.gz` on Zenodo),
`ARCHI_MASK_HUMAN` (`track_out_HUMAN.bed`, from `repeatmasker_hg38_mm10.tar.gz` on Zenodo), `ARCHI_SEQTK` (seqtk,
submodule `../seqtk`), `ARCHI_BIASAWAY` (biasaway executable) and `ARCHI_SARUS_JAR` (`sarus-2.0.2.jar`, submodule
`../sarus`). The Zenodo records are listed in `../README.md`.

Three more inputs are read from fixed places under `$ARCHI_WORK` by `h5_data_server.sh`: `human_sym/inputs/human_sym_peak_ids.csv` (peak-id table of the 35 human TFs, columns `TF,peak_id`), `mouse_sym/inputs/mouse_sym_peak_ids.csv` (same for mouse; the file `../results/datasets/peak_ids/mouse_sym_peak_ids.csv` has this name and format) and `mouse_sym/ref/track_out_MOUSE_mm10.bed` (mouse repeat mask; the Zenodo archive `repeatmasker_hg38_mm10.tar.gz` holds the RepeatMasker masks). The human table comes from an earlier server run and is not published under that name; `../results/datasets/peak_ids/human_used_peak_ids.csv` lists the used human peak sets (the 35 TFs and P53). Copy or link these files to the paths above before the run.

## Run

1. From the repository root, `bash 1_data_preparation/h5_data_server.sh` builds everything in one driver (about 1.5 h;
   every producer refuses to overwrite): `build_positives_v5.py` (human, mouse), `check_h5_test_identity.py positives`,
   `make_negatives_v3.py` (140 jobs, `xargs -P12`), `check_h5_test_identity.py negatives`, `build_slim_inputs_v5.py`
   for 35 TFs, `make_scan_jobs_sym.py` for human and mouse, then 15,984 SARUS jobs (`xargs -P24`) and
   `check_h5_test_identity.py sarus`. It stops on any difference in a test file and writes the end marker
   `H5_DATA_DONE` to `$ARCHI_WORK/human_sym_h5/logs/data_driver.log`.
2. Checks after the driver, run in `1_data_preparation/` with `PY=${ARCHI_PY:-python3}` and
   `DH=$ARCHI_WORK/human_sym_h5`:
   - `$PY check_h5_data.py $DH/checks/data`
   - `$PY check_C_sym.py --tfml <TF-ML> --peaks <species>_sym_peak_ids.csv --mask <bed> --tmp <dir> --species HUMAN|MOUSE`
   - `$PY check_E.py <TF-ML> <peaks csv> HUMAN|MOUSE <out csv>`
   - `$PY counts_D_human_h5.py` and `$PY counts_D_mouse_h5.py`
   - `$PY check_neg_v3_h5.py $DH/checks_neg_v3`
   - `$PY gc_check_v3_h5.py --tfs <35 TFs> --out $DH/checks_neg_v3/gc_check_v3.csv`

   The peak-id tables of the used peak sets are in `../results/datasets/peak_ids/`.

## Loader and file names

`archi_data_v5.py` reads the build folders, where the mouse files named `train` hold the test chromosomes 1, 8 and 19;
its `_file_split()` maps the logical split to those names. In the Zenodo files `train` always marks the training
chromosomes, so to read them `_file_split()` must return the split unchanged for MOUSE as well.

## Software

Versions in `../2_training/VERSIONS.txt`: Python 3.8.18, scikit-learn 1.3.0, pandas 2.0.3, numpy 1.24.3, scipy 1.10.1,
pybedtools 0.9.1, bedtools 2.31.1, seqtk 1.3-r107-dirty, BiasAway 3.3.0, SARUS 2.0.2. `prroc_metrics.py` needs R with
PRROC 1.4 and rpy2 3.5; rpy2 must load the R that holds PRROC (set `R_HOME`).
