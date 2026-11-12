# 3_analysis/server

Analyses that run after the models exist: operational and low-FPR metrics, top-peak and random-split analyses, feature importance, PWM saturation, the Fig. 4 and S4 refit with the Slim and diChIPMunk scans, the runtime benchmark (Sup. Table 6), and the provenance file. The result tables that the local analysis stages read are tracked in `../../results/server`.

## Environment

Paths outside the repository come from environment variables, read in `../../archi_env.py` and `../../archi_env.sh`; a variable that is needed and not set stops the script with a message. `ARCHI_WORK` is the working root (the drivers work in `$ARCHI_WORK/H5_alltrain/human_sym`), `ARCHI_RELEASE` the published `TF-ML` folder, `ARCHI_PY` the Python interpreter (default `python3`, packages as in `../../2_training/VERSIONS.txt`), `ARCHI_SARUS_JAR` the `sarus-2.0.2.jar` (runtime benchmark), `ARCHI_SARUS201_JAR` the `sarus-2.0.1.jar` (diChIPMunk scans), `ARCHI_MASK_HUMAN` the human repeat mask (`h5_provenance.py`), `ARCHI_MACS` the MACS peak files (`b3_join_sym.py`). The drivers call the scripts of the repository directly.

## Run order

1. Most analyses of this folder are started by `../../2_training/h5_server_v2.sh` (phase 1: `saturation_*/run_saturation_h5.py`; phase 2: `operational_metrics.py`, `b2_compute.py`, `b3_topN_h5.py`, `feature_importance/b4_01` to `b4_07`, `b9_random_split_h5.py`, `b9_summary_h2b.py`). See `../../2_training/README.md`. `b3_join_sym.py` produces the width and join tables of the top-peak analysis (`python b3_join_sym.py OUTDIR --species H|M`).
2. `bash 3_analysis/server/h5_fig4_server_v2.sh` (working folder `$ARCHI_WORK/H5_alltrain/human_sym`, about 26 h): per TF `one_tf_h5.sh` runs `scan_tf_h5.sh` (Slim m0, m1, m5 scans and diChIPMunk scans of the training table, human test and mouse test) and `run_fig4_h5.py` (RF refit with the training scaler), 24 TFs in parallel with 1 thread each. The driver starts a TF once its models are in `slim/models_server/<TF>/` with `COPY_VERIFIED.txt` (`2_training/README.md`, step 2). It then runs `assemble_fig4_s4.py` (human and mouse) and `s4_numbers.py`. End marker `FIG4_H5_DONE` in `logs/fig4_driver.log`.
3. Runtime benchmark (about 10 min, run when nothing else runs on the machine): `bash 3_analysis/server/runtime/make_h5_runtime_v2.sh` creates `$ARCHI_WORK/H5_alltrain/human_sym/runtime_v2` (the driver exits when it exists) and runs `build_h3a_matrices.py data STAT3 HNF4A CTCF`, `run_runtime_h3a.sh` and `assemble_sup_table6_h3a.py` from `runtime/`. The input of block D is `$ARCHI_WORK/runtime/data/bench_50k_300bp.fa`; the published timings used for the comparison are read from `$ARCHI_WORK/runtime`.
4. `ARCHI_GIT_HASH=<hash> ARCHI_MASK_HUMAN=<track_out_HUMAN.bed> $PY 3_analysis/server/h5_provenance.py` writes `_provenance.yaml` (chromosome lists, md5 of scripts, inputs and outputs).
5. Checks: `$PY tests/server_checks/h5_post.py` (equality and difference checks against the earlier build, headline numbers) and `$PY tests/server_checks/h5_post_v2.py <table1_rows.csv>` (Table 1 definitions).

`$PY` in the commands is the interpreter of `ARCHI_PY`.

## Software

Python 3.8.18, scikit-learn 1.3.0, pandas 2.0.3, numpy 1.24.3, scipy 1.10.1, joblib 1.2.0, xgboost 1.7.3 (versions in `../../2_training/VERSIONS.txt`). The diChIPMunk scans of `scan_tf_h5.sh` used SARUS 2.0.1 (`ARCHI_SARUS201_JAR`), not the 2.0.2 jar of the data build. Slim scans use `TrainAndApplySlim.jar` with JDK 8u232.

## Not in the repository

- The trained Slim and diChIPMunk models (`slim/models_server/<TF>`; on Zenodo, `Slim_diChIPMunk_models_v2.tar.gz`) and the data of `../../1_data_preparation` (on Zenodo). The RF results and score files of `../../2_training` (`hm_mm_results.csv`, `scores/*.npz`) are tracked in `../../results/server/human/`.
- Reference tables of the earlier build, read by `h5_fig4_server_v2.sh`, `h5_post.py` and `h5_post_v2.py` (default paths below `$ARCHI_WORK`; `H5POST_G`, `H5POST_G4` and `H5POST_RED` set the roots and the red control).
- The benchmark input `$ARCHI_WORK/runtime/data/bench_50k_300bp.fa` and the earlier timings in `$ARCHI_WORK/runtime`.
- The published release folder (`$ARCHI_RELEASE`) and `notebook4_results_table.csv`.
