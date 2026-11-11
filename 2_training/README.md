# 2_training

Trains and evaluates the models: random forests (RF) on server, Slim and diChIPMunk models on server, the evaluation tables, and the release models (`Models.tar.gz`).

## What runs where

| step | machine | script |
|---|---|---|
| ship Slim inputs to server | Mac relay | `server/relay/ship_server_to_server_h5.sh` |
| Slim (105 models) and diChIPMunk (35) training | server | `server/h5_server_queue.sh` |
| copy trained models back to server | Mac relay | `server/relay/sync_server_to_server_h5.sh` |
| RF training, evaluation, evaluation tables, and the analyses of `../3_analysis/server` | server | `h5_server_v2.sh` |
| release models | server | `train_models_v2d.py` |

The scripts keep the server paths. `setup_server_layout.sh` in the repository root recreates the layout (server: `~/archi_work/H5_alltrain/human_sym/code` holds links to the scripts of this folder and of `../3_analysis/server`; `--server` links `~/archi_work/H5_alltrain/code` to `server/`; `--relay` adds the link that `ship_server_to_server_h5.sh` reads as `$HERE/server`).

## Run order

Data preparation (`../1_data_preparation`) comes first.

1. Mac: `2_training/slim_chipmunk/relay/ship_server_to_server_h5.sh` waits for the Slim inputs on server, streams the shipment and the server code to `~/archi_work/H5_alltrain` on server, checks the md5 manifest and starts the queue. Run it as `nohup ./ship_server_to_server_h5.sh > ../logs/ship.log 2>&1 &`.
2. server: `bash code/h5_server_queue.sh` (working folder `~/archi_work/H5_alltrain`; started by the relay). It rebuilds the 35 Slim inputs with `build_slim_inputs_from_selection.py`, checks them md5-identical to server, then runs one queue (`xargs -P24`, 4 threads per job, longest first, lists in `server/slim_jobs_lpt.txt` and `server/munk_jobs_lpt.txt`): `train_slim_h5.sh <TF> <0|1|5>` for 105 Slim models and `run_munk_h5.sh <TF>` for 35 diChIPMunk runs (`ruby chipmunk/run_dichiphorde8.rb`). End marker `SERVER_H5_QUEUE_DONE` in `logs/queue.log`.
3. Mac: `2_training/slim_chipmunk/relay/sync_server_to_server_h5.sh` copies the Slim xml files and the diChIPMunk dpwm of each finished TF to server (`slim/models_server/<TF>/`), checks md5 and stops when 35 of 35 TFs are verified.
4. server: `bash 2_training/h5_server_v2.sh` (working folder `~/archi_work/H5_alltrain/human_sym`, about 3 h with 32 threads):
   - phase 1: both `run_saturation_h5.py` (in `../3_analysis/server`) and `run_hm_mm_h5.py --tfs $TFS --workers 4 --rf-jobs 4 --out-csv hm_mm_results.csv --rows-dir rows --scores-dir scores`;
   - phase 2 (the analyses of `../3_analysis/server`): `operational_metrics.py`, `b2_v2/b2_compute.py`, `b3_topN_h5.py`, the feature-importance scripts `b4_01`, `b4_02`, `b4_03`, `b4_04` (`--thr 80` and `--thr 90`) and `b4_07`, `N_JOBS=8 b9_random_split_h5.py b9 PRGR RXRA SPI1 ANDR REST`, `b9_summary_h2b.py`;
   - phase 3: `evaluate_models_h5.py --species HUMAN|MOUSE --workers 4 --threads-per-worker 3`, then `assemble_evaluation_table_h3a.py --mouse-trained ... --notebook4-table <notebook4_results_table.csv> --exclude-tfs P53` for both species.
   End marker `SERVER_H5_DONE` in `logs/driver.log`.
5. Release models, by hand (command not recorded): `$PY 2_training/train_models_v2d.py <out_dir> <threads> <workers>`.

Run order for the Fig. 4 scans and the runtime benchmark is in `../3_analysis/server/README.md`; the Fig. 4 driver can start as soon as step 3 has verified a TF.

## Software

`VERSIONS.txt` lists the versions. Server environment: `~/miniconda3/envs/rpy2_env` with Python 3.8.18, scikit-learn 1.3.0, pandas 2.0.3, numpy 1.24.3, scipy 1.10.1, joblib 1.2.0, xgboost 1.7.3, pybedtools 0.9.1; one BLAS thread per worker, `nice -n 10`. RF settings: `RandomForestClassifier(max_depth=6, max_samples=0.8, n_estimators=100, random_state=0)` with the StandardScaler fitted on the training matrix. The diChIPMunk scans of `../3_analysis/server/scan_tf_h5.sh` used SARUS 2.0.1; the feature tables of the data build used SARUS 2.0.2.

## Not in the repository

- diChIPMunk: `chipmunk.jar` (v8) and `dpmflogo3.rb` from the ChIPMunk v8 supporting scripts (https://autosome.org/ChIPMunk/). `server/chipmunk/run_dichiphorde8.rb` expects `chipmunk.jar` and `ytilib/` next to it. The Slim jar `TrainAndApplySlim.jar`, the JDK 8u232 and `ytilib/` are in `../Slim/` of this repository; on the servers they were kept in `~/Slim` and in a tools folder under `~/archi_work/human_sym_slim` on server.
- `notebook4_results_table.csv` (schema of the published results table), the data of `../1_data_preparation`, and the published release folder `~/archi_release/TF-ML`.
- The relay scripts write their logs to `2_training/slim_chipmunk/logs` (next to the `relay` folder), which is not tracked.
