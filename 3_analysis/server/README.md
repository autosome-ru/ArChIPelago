# 3_analysis/server

Compute-heavy analyses that run after the models exist: operational and low FPR metrics, top peaks, random split,
feature importance, PWM saturation, the Fig. 4 and Fig. S4 refits with the Slim and diChIPMunk scans, the runtime
benchmark (Sup. Table 6) and the provenance record. Their result tables are tracked in `../../results/server/`, which the
local stages read.

## Files

| file | role |
|---|---|
| `operational_metrics.py` | operational metrics |
| `b2_compute.py`, `b2_lib.py` | low FPR metrics |
| `b3_topN_h5.py`, `b3_join_sym.py` | top-N peaks; peak widths and join of the test windows with the peak files |
| `feature_importance/b4_01_importance.py` to `feature_importance/b4_07_checks.py`, `feature_importance/b4_common.py` | feature importance, PWM subtypes, ablation |
| `b9_random_split_h5.py`, `b9_summary_h2b.py` | random split retraining (5 TFs) |
| `saturation_HUMAN/run_saturation_h5.py`, `saturation_MOUSE/run_saturation_h5.py` | PWM saturation |
| `h5_fig4_server_v2.sh`, `one_tf_h5.sh`, `scan_tf_h5.sh`, `run_fig4_h5.py`, `assemble_fig4_s4.py`, `s4_numbers.py` | Fig. 4 and Fig. S4: Slim and diChIPMunk scans, RF refit, assembly |
| `runtime/make_h5_runtime_v2.sh`, `runtime/build_h3a_matrices.py`, `runtime/run_runtime_h3a.sh`, `runtime/rf_bench_h3a.py`, `runtime/rf_predict_chunked_h3a.py`, `runtime/assemble_sup_table6_h3a.py` | runtime benchmark |
| `h5_provenance.py` | `_provenance.yaml` of the compute stages |

## Environment

Variables are read in `../../archi_env.py` and `../../archi_env.sh`: `ARCHI_WORK` (the drivers work in
`$ARCHI_WORK/H5_alltrain/human_sym`), `ARCHI_RELEASE`, `ARCHI_PY` (default `python3`), `ARCHI_SARUS_JAR`
(`sarus-2.0.2.jar`, runtime benchmark), `ARCHI_SARUS201_JAR` (`sarus-2.0.1.jar`, diChIPMunk scans), `ARCHI_MASK_HUMAN`
(`h5_provenance.py`), `ARCHI_MACS` (`b3_join_sym.py`). The Slim scans need `Slim/` in the repository root (unpacked
`Slim.tar.gz`, see `../../README.md`).

## Run

1. `../../2_training/h5_server_v2.sh` starts the saturation scripts and the analyses of this folder (see
   `../../2_training/README.md`). `b3_join_sym.py` is run by hand: `python b3_join_sym.py OUTDIR --species H|M`.
2. `bash 3_analysis/server/h5_fig4_server_v2.sh` (about 26 h): per TF, `one_tf_h5.sh` runs `scan_tf_h5.sh` (Slim m0,
   m1, m5 and diChIPMunk scans of the training table, human test and mouse test) and `run_fig4_h5.py` (RF refit with the
   training scaler), 24 TFs in parallel with 1 thread each. A TF starts once
   `slim/models_server/<TF>/COPY_VERIFIED.txt` exists (`../../2_training/README.md`, step 2). The driver then runs
   `assemble_fig4_s4.py` (human and mouse; it also reads the reference tables
   `$ARCHI_WORK/human_sym/H4_trainscale/fig4_refit_trainscale/Figure_{4,S4}_source_data.csv`) and `s4_numbers.py`. End
   marker `FIG4_H5_DONE` in `logs/fig4_driver.log`.
3. Runtime benchmark (about 10 min, on an otherwise idle machine): `bash 3_analysis/server/runtime/make_h5_runtime_v2.sh`
   creates `$ARCHI_WORK/H5_alltrain/human_sym/runtime_v2` (it exits if the folder exists) and runs
   `build_h3a_matrices.py`, `run_runtime_h3a.sh` and `assemble_sup_table6_h3a.py`. Inputs:
   `$ARCHI_WORK/runtime/data/bench_50k_300bp.fa` and the timings in `$ARCHI_WORK/runtime`.
4. `ARCHI_GIT_HASH=<hash> ARCHI_MASK_HUMAN=<track_out_HUMAN.bed> $ARCHI_PY 3_analysis/server/h5_provenance.py` writes
   `_provenance.yaml` (chromosome lists, md5 of scripts, inputs and outputs).

## Software

Versions in `../../2_training/VERSIONS.txt`: Python 3.8.18, scikit-learn 1.3.0, pandas 2.0.3, numpy 1.24.3, scipy 1.10.1,
joblib 1.2.0, xgboost 1.7.3. The diChIPMunk scans use SARUS 2.0.1, the data build SARUS 2.0.2. Slim scans use
`TrainAndApplySlim.jar` with JDK 8u232.
