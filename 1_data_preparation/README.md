# 1_data_preparation

Builds the training and test data of the RF models for human and mouse: positive windows, GC-matched negatives, SARUS feature tables, Slim inputs, and the data checks. Runs on the build server server (CPU, no GPU).

## Run order

The scripts keep the server paths (`~/archi_work/...`). `setup_server_layout.sh` in the repository root recreates the layout the drivers expect (`~/archi_work/human_sym_h5/code` points to this folder):

```
bash setup_server_layout.sh <repo> ~/archi_work
```

Set `PY=~/miniconda3/envs/rpy2_env/bin/python` and `DH=~/archi_work/human_sym_h5`.

1. `bash 1_data_preparation/h5_data_server.sh` builds everything in one driver (new folders only, every producer refuses to overwrite, about 1.5 h):
   `build_positives_v5.py` per TF and pooled (human, mouse) - `check_h5_test_identity.py positives` - `make_negatives_v3.py` (140 jobs, `xargs -P12`) - `check_h5_test_identity.py negatives` - `build_slim_inputs_v5.py` for 35 TFs and the tarballs for server (in the background) - `make_scan_jobs_sym.py` for human and mouse, then 15,984 SARUS jobs (`xargs -P24`) - `check_h5_test_identity.py sarus`. The driver stops on any difference in a test file and writes the end marker `H5_DATA_DONE` to `$DH/logs/data_driver.log`.
2. Tables and checks, run by hand after the driver. The commands are reconstructed from the usage lines and log names; they were not recorded:
   - `$PY check_h5_data.py $DH/checks/data`
   - `$PY check_C_sym.py --tfml TF-ML --peaks <sp>_sym_peak_ids.csv --mask <bed> --tmp <dir> --species HUMAN|MOUSE > logs/check_C.log`
   - `$PY check_E.py TF-ML <peaks csv> HUMAN|MOUSE checks/check_E.csv`
   - `$PY counts_D_human_h5.py > logs/counts_D_v2.log` and `$PY counts_D_mouse_h5.py > logs/counts_D.log`
   - `$PY check_neg_v3_h5.py $DH/checks_neg_v3`
   - `$PY gc_check_v3_h5.py --tfs <35 TFs> --out $DH/checks_neg_v3/gc_check_v3.csv`

`archi_data_v5.py` is the loader that all later stages import (training and analysis scripts add this folder to `sys.path`). `h5_splits.py` defines the chromosome split once: human training chr2-7, 9-20, 22, test chr1, 8, 21; mouse training chr2-7, 9-18, test chr1, 8, 19. On disk, mouse files named `train` hold the test chromosomes 1, 8 and 19; the loader maps the logical split to those names.

## Software

Versions are listed in `../2_training/VERSIONS.txt`: Python 3.8.18, scikit-learn 1.3.0, pandas 2.0.3, numpy 1.24.3, scipy 1.10.1, pybedtools 0.9.1, bedtools 2.31.1, seqtk 1.3-r107-dirty, BiasAway 3.3.0, SARUS 2.0.2 (`sarus-2.0.2.jar`), conda environment `rpy2_env`. The diChIPMunk scans of `../3_analysis/server/scan_tf_h5.sh` used SARUS 2.0.1.

## Not in the repository

- Peak sets (`~/macs`, `PEAKS041820.interval`), the peak-id tables `human_sym_peak_ids.csv` and `mouse_sym_peak_ids.csv`, the repeat masks (`track_out_HUMAN.bed`, `track_out_MOUSE_mm10.bed`), the genomes (hg38, GRCm38.p6) and `core_annotation_HUMAN_mono.csv`.
- The published release folder `~/archi_release/TF-ML` (PWM files, feature folders, published tables) and the HOCOMOCO v11 models (`~/hocomoco11`).
- Tools: seqtk, BiasAway, bedtools, the SARUS jars, `chipmunk.jar` with `ytilib/` (see `../2_training/README.md`).
- The Zenodo archive holds the resulting data; this folder only holds the code that built it.
