# 2_training

Trains and evaluates the models: random forests (RF), Slim and diChIPMunk models, the evaluation tables, and the release models (`Models.tar.gz`).

## Scripts

| step | script |
|---|---|
| Slim (105 models) and diChIPMunk (35) training | `slim_chipmunk/h5_server_queue.sh` |
| collecting the trained models for the Fig. 4 scans | copy step, see step 2 below |
| RF training, evaluation, evaluation tables, and the analyses of `../3_analysis/server` | `h5_server_v2.sh` |
| release models | `train_models_v2d.py` |

Paths outside the repository come from environment variables, read in `../archi_env.py` and `../archi_env.sh`; a variable that is needed and not set stops the script with a message. `ARCHI_WORK` is the working root (`H5_alltrain/`, `human_sym_h5/`, `mouse_sym_h5/`), `ARCHI_RELEASE` the published `TF-ML` folder, `ARCHI_PY` the Python interpreter (default `python3`, packages as in `VERSIONS.txt`), `ARCHI_CHIPMUNK` the folder with `run_dichiphorde8.rb`, `chipmunk.jar` and `ytilib/`, `ARCHI_NB4` the file `notebook4_results_table.csv`, `ARCHI_MASK_HUMAN`, `ARCHI_MACS`, `ARCHI_GENOMES`, `ARCHI_FAMILY`, `ARCHI_HOCOMOCO`, `ARCHI_SEQTK`, `ARCHI_BIASAWAY` and `ARCHI_SARUS_JAR` as in `../1_data_preparation/README.md`. The Slim jar, the JDK and `ytilib/` are taken from `../Slim/`.

## Run order

Data preparation (`../1_data_preparation`) comes first.

1. `bash 2_training/slim_chipmunk/h5_server_queue.sh` (working folder `$ARCHI_WORK/H5_alltrain`). It unpacks the Slim inputs that `../1_data_preparation/h5_data_server.sh` wrote to `$ARCHI_WORK/H5_alltrain/human_sym/slim` (`ship/`, `inputs_md5.txt`, `munk_src_md5.txt`), rebuilds the 35 Slim inputs with `build_slim_inputs_from_selection.py`, checks them against `inputs_md5.txt`, then runs one queue (`xargs -P24`, 4 threads per job, longest first, lists in `slim_chipmunk/slim_jobs_lpt.txt` and `slim_chipmunk/munk_jobs_lpt.txt`): `train_slim_h5.sh <TF> <0|1|5>` for 105 Slim models and `run_munk_h5.sh <TF>` for 35 diChIPMunk runs (`ruby run_dichiphorde8.rb`). End marker `SERVER_H5_QUEUE_DONE` in `logs/queue.log`.
2. For each TF whose three Slim models (log line `exit code 0`) and diChIPMunk dpwm exist, the files the scans need are
   copied to `H5_alltrain/human_sym/slim/models_server/<TF>/` without overwriting, checked by md5, and the md5 list is
   written to `COPY_VERIFIED.txt`, which `h5_fig4_server_v2.sh` waits for (35 TFs):
   ```bash
   t=<TF>; cd $ARCHI_WORK/H5_alltrain/models
   F="$t/ChIPMunk.txt $t/train_pos_id_HUMAN_p.fasta"
   for o in 0 1 5; do F="$F $t/${t}_SlimModel_$o/Motif_1/SlimDimont_1.xml $t/${t}_SlimModel_$o/protocol_slimdimont.txt"; done
   for e in dpwm dpcm dat; do F="$F $t/HUMAN_seq_HUMAN_pwm_mono/train_H_1_ChIPMunk_no_repeats_0.$e"; done
   D=$ARCHI_WORK/H5_alltrain/human_sym/slim/models_server; L=$ARCHI_WORK/H5_alltrain/logs/md5_$t.txt
   md5sum $F > $L && mkdir -p $D && tar cf - $F | (cd $D && tar xkf -) && (cd $D && md5sum -c --quiet $L) && cp -n $L $D/$t/COPY_VERIFIED.txt
   ```
3. `bash 2_training/h5_server_v2.sh` (working folder `$ARCHI_WORK/H5_alltrain/human_sym`, about 3 h with 32 threads):
   - phase 1: both `run_saturation_h5.py` (in `../3_analysis/server`) and `run_hm_mm_h5.py --tfs $TFS --workers 4 --rf-jobs 4 --out-csv hm_mm_results.csv --rows-dir rows --scores-dir scores`;
   - phase 2 (the analyses of `../3_analysis/server`): `operational_metrics.py`, `b2_compute.py`, `b3_topN_h5.py`, the feature-importance scripts `b4_01`, `b4_02`, `b4_03`, `b4_04` (`--thr 80` and `--thr 90`) and `b4_07`, `N_JOBS=8 b9_random_split_h5.py b9 PRGR RXRA SPI1 ANDR REST`, `b9_summary_h2b.py`;
   - phase 3: `evaluate_models_h5.py --species HUMAN|MOUSE --workers 4 --threads-per-worker 3`, then `assemble_evaluation_table_h3a.py --mouse-trained ... --notebook4-table <notebook4_results_table.csv> --exclude-tfs P53` for both species.
   End marker `SERVER_H5_DONE` in `logs/driver.log`.
4. Release models, by hand (command not recorded): `$PY 2_training/train_models_v2d.py <out_dir> <threads> <workers>`.

Run order for the Fig. 4 scans and the runtime benchmark is in `../3_analysis/server/README.md`; the Fig. 4 driver can start as soon as step 2 has verified a TF.

## Software

`VERSIONS.txt` lists the versions. Environment: Python 3.8.18, scikit-learn 1.3.0, pandas 2.0.3, numpy 1.24.3, scipy 1.10.1, joblib 1.2.0, xgboost 1.7.3, pybedtools 0.9.1; one BLAS thread per worker, `nice -n 10`. RF settings: `RandomForestClassifier(max_depth=6, max_samples=0.8, n_estimators=100, random_state=0)` with the StandardScaler fitted on the training matrix. The diChIPMunk scans of `../3_analysis/server/scan_tf_h5.sh` used SARUS 2.0.1; the feature tables of the data build used SARUS 2.0.2.

## Not in the repository

- diChIPMunk: `chipmunk.jar` (v8) and `dpmflogo3.rb` from the ChIPMunk v8 supporting scripts (https://autosome.org/ChIPMunk/). `slim_chipmunk/chipmunk/run_dichiphorde8.rb` expects `chipmunk.jar` and `ytilib/` next to it. The Slim jar `TrainAndApplySlim.jar`, the JDK 8u232 and `ytilib/` are in `../Slim/` of this repository.
- `notebook4_results_table.csv` (schema of the published results table), the data of `../1_data_preparation`, and the published release folder (`$ARCHI_RELEASE`).
