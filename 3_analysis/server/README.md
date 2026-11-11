# 3_analysis/server

Analyses that run on server after the models exist: operational and low-FPR metrics, top-peak and random-split analyses, feature importance, PWM saturation, the Fig. 4 and S4 refit with the Slim and diChIPMunk scans, the runtime benchmark (Sup. Table 6), and the provenance file. The results are mirrored to the Mac (`../../results/server`), where the local analysis stages read them.

## What runs where

Everything in this folder runs on server. The only Mac script is `../../results/server/mirror_h5_to_mac_v2.sh`. The scripts keep the server paths. `setup_server_layout.sh` in the repository root recreates the layout under `~/archi_work/H5_alltrain/human_sym`: `code/` (links to the scripts of this folder and of `../../2_training`), `b2_v2/` (`b2_compute.py`, `b2_lib.py`), `saturation_HUMAN/` and `../mouse_sym/saturation_MOUSE/` (the saturation scripts write `per_tf/` next to themselves, so these are real folders with linked scripts), `B4_importance/scripts` (link to `feature_importance/`) and `runtime_v2/code` (link to `runtime/`).

## Run order

Set `PY=~/miniconda3/envs/rpy2_env/bin/python`.

1. Most analyses of this folder are started by `../../2_training/h5_server_v2.sh` (phase 1: `saturation_*/run_saturation_h5.py`; phase 2: `operational_metrics.py`, `b2_compute.py`, `b3_topN_h5.py`, `feature_importance/b4_01` to `b4_07`, `b9_random_split_h5.py`, `b9_summary_h2b.py`). See `../../2_training/README.md`. `b3_join_sym.py` produces the width and join tables of the top-peak analysis (`python b3_join_sym.py OUTDIR --species H|M`).
2. `bash 3_analysis/server/h5_fig4_server_v2.sh` (working folder `~/archi_work/H5_alltrain/human_sym`, about 26 h): per TF `one_tf_h5.sh` runs `scan_tf_h5.sh` (Slim m0, m1, m5 scans and diChIPMunk scans of the training table, human test and mouse test) and `run_fig4_h5.py` (RF refit with the training scaler), 24 TFs in parallel with 1 thread each. The driver starts a TF once the relay has written `slim/models_server/<TF>/COPY_VERIFIED.txt`. It then runs `assemble_fig4_s4.py` (human and mouse) and `s4_numbers.py`. End marker `FIG4_H5_DONE` in `logs/fig4_driver.log`.
3. Runtime benchmark (about 10 min, run when nothing else of ours runs on the node): `bash 3_analysis/server/runtime/make_h5_runtime_v2.sh` copies the benchmark code of an earlier run of the benchmark (`~/archi_work/human_sym/.../runtime/code`, not in the repository) into `runtime_v2/code`, changes the loader import and the output root with two `sed` edits, and then runs `build_h3a_matrices.py data STAT3 HNF4A CTCF`, `run_runtime_h3a.sh` and `assemble_sup_table6_h3a.py`. The files in `runtime/` are the versions that this driver produced. The driver exits when `runtime_v2` exists, and `setup_server_layout.sh` creates it, so run the driver only on a layout set up without the server mode or remove `runtime_v2` first. To rerun without the copy step, work in `~/archi_work/H5_alltrain/human_sym/runtime_v2` (created by `setup_server_layout.sh`, with `code` linked to `runtime/`; add the folders `scratch` and `data`) and run the three steps by hand: `$PY code/build_h3a_matrices.py data STAT3 HNF4A CTCF`, `bash code/run_runtime_h3a.sh`, `$PY code/assemble_sup_table6_h3a.py $PWD ~/archi_work/runtime`.
4. `ARCHI_GIT_HASH=<hash> $PY 3_analysis/server/h5_provenance.py` writes `_provenance.yaml` (chromosome lists, md5 of scripts, inputs and outputs).
5. Checks: `$PY tests/server_checks/h5_post.py` (equality and difference checks against the earlier build, headline numbers) and `$PY tests/server_checks/h5_post_v2.py <table1_rows.csv>` (Table 1 definitions).
6. Mac: `results/server/mirror_h5_to_mac_v2.sh` copies 822 result files and writes `MIRROR_MANIFEST_md5.txt` (one ssh connection to server; refuses to run when the mirror folders exist).

## Software

Python 3.8.18, scikit-learn 1.3.0, pandas 2.0.3, numpy 1.24.3, scipy 1.10.1, joblib 1.2.0, xgboost 1.7.3 (`~/miniconda3/envs/rpy2_env`, versions in `../../2_training/VERSIONS.txt`). The diChIPMunk scans of `scan_tf_h5.sh` used SARUS 2.0.1 (`~/sarus/releases/sarus-2.0.1.jar`), not the 2.0.2 jar of the data build. Slim scans use `TrainAndApplySlim.jar` with JDK 8u232.

## Not in the repository

- The trained Slim and diChIPMunk models (`slim/models_server/<TF>`), the RF results and score files of `../../2_training` (`hm_mm_results.csv`, `scores/*.npz`), and the data of `../../1_data_preparation`.
- Reference tables of the earlier build, read by `h5_fig4_server_v2.sh`, `h5_post.py` and `h5_post_v2.py` (default paths on server; `H5POST_G`, `H5POST_G4` and `H5POST_RED` set the roots and the red control).
- The earlier benchmark folder used by `make_h5_runtime_v2.sh` and the benchmark input `~/archi_work/runtime/data/bench_50k_300bp.fa`.
- The published release folder `~/archi_release/TF-ML` and `notebook4_results_table.csv`.
