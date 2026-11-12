# 2_training

Trains and evaluates the models: Random Forest (RF) and the other algorithms, Slim and diChIPMunk models, the evaluation
tables, and the release models of `Models.tar.gz`. `1_data_preparation/` runs first.

## Files

| file | role |
|---|---|
| `slim_chipmunk/h5_server_queue.sh` | driver: 105 Slim models and 35 diChIPMunk models |
| `slim_chipmunk/build_slim_inputs_from_selection.py`, `slim_chipmunk/train_slim_h5.sh`, `slim_chipmunk/run_munk_h5.sh`, `slim_chipmunk/chipmunk/run_dichiphorde8.rb` | Slim inputs, Slim training, diChIPMunk training |
| `slim_chipmunk/slim_jobs_lpt.txt`, `slim_chipmunk/munk_jobs_lpt.txt` | job lists of the queue, longest first |
| `h5_server_v2.sh` | driver: RF training and evaluation, the analyses of `../3_analysis/server/`, evaluation tables |
| `run_hm_mm_h5.py`, `evaluate_models_h5.py`, `assemble_evaluation_table_h3a.py` | RF training, evaluation, evaluation tables |
| `train_models_v2d.py` | release models (`Models/<TF>/ArChIPelago_<TF>[_MOUSE]_<mode>.{sav,json}`) |
| `VERSIONS.txt` | software versions of steps 1-3 |

## Environment

Paths come from environment variables, read in `../archi_env.py` and `../archi_env.sh`; a variable that is needed and
not set stops the script with a message. `ARCHI_WORK` is the working root (`H5_alltrain/`, `human_sym_h5/`,
`mouse_sym_h5/`), `ARCHI_RELEASE` the pipeline output folder `TF-ML`, `ARCHI_PY` the Python interpreter (default
`python3`), `ARCHI_CHIPMUNK` the folder with `run_dichiphorde8.rb`, `chipmunk.jar` and `ytilib/`, `ARCHI_NB4` the
results-table schema (`../3_analysis/package/analysis/inputs/notebook4_results_table.csv`). The other variables are as in
`../1_data_preparation/README.md`.

Inputs: `Slim.tar.gz` of Zenodo v1.0.0, unpacked in the repository root (`Slim/`: the Slim jar
`TrainAndApplySlim.jar`, JDK 8u232, `ytilib/`; record and md5 in `../README.md`); `chipmunk.jar` (v8) and `dpmflogo3.rb` from the
ChIPMunk v8 supporting scripts (https://autosome.org/ChIPMunk/), placed in `$ARCHI_CHIPMUNK`; the peak-id table `$ARCHI_WORK/human_sym/inputs/human_sym_peak_ids.csv` (from the data build, see `../1_data_preparation/README.md`; not published under that name).

## Run

1. `bash 2_training/slim_chipmunk/h5_server_queue.sh` (works in `$ARCHI_WORK/H5_alltrain`). It unpacks the Slim inputs
   that `../1_data_preparation/h5_data_server.sh` wrote to `$ARCHI_WORK/H5_alltrain/human_sym/slim`, rebuilds the 35
   Slim inputs, checks them against `inputs_md5.txt`, then runs one queue (`xargs -P24`, 4 threads per job):
   `train_slim_h5.sh <TF> <0|1|5>` for 105 Slim models and `run_munk_h5.sh <TF>` for 35 diChIPMunk runs. End marker
   `SERVER_H5_QUEUE_DONE` in `logs/queue.log`.
2. For each TF whose three Slim models (log line `exit code 0`) and diChIPMunk dpwm exist, copy the files the scans need
   to `H5_alltrain/human_sym/slim/models_server/<TF>/` without overwriting, check them by md5 and write the md5 list to
   `COPY_VERIFIED.txt`, which `../3_analysis/server/h5_fig4_server_v2.sh` waits for:
   ```bash
   t=<TF>; cd $ARCHI_WORK/H5_alltrain/models
   F="$t/ChIPMunk.txt $t/train_pos_id_HUMAN_p.fasta"
   for o in 0 1 5; do F="$F $t/${t}_SlimModel_$o/Motif_1/SlimDimont_1.xml $t/${t}_SlimModel_$o/protocol_slimdimont.txt"; done
   for e in dpwm dpcm dat; do F="$F $t/HUMAN_seq_HUMAN_pwm_mono/train_H_1_ChIPMunk_no_repeats_0.$e"; done
   D=$ARCHI_WORK/H5_alltrain/human_sym/slim/models_server; L=$ARCHI_WORK/H5_alltrain/logs/md5_$t.txt
   md5sum $F > $L && mkdir -p $D && tar cf - $F | (cd $D && tar xkf -) && (cd $D && md5sum -c --quiet $L) && cp -n $L $D/$t/COPY_VERIFIED.txt
   ```
3. `bash 2_training/h5_server_v2.sh` (works in `$ARCHI_WORK/H5_alltrain/human_sym`, about 3 h with 32 threads):
   - phase 1: `../3_analysis/server/saturation_HUMAN/run_saturation_h5.py`,
     `../3_analysis/server/saturation_MOUSE/run_saturation_h5.py` and `run_hm_mm_h5.py`;
   - phase 2: the analyses of `../3_analysis/server/` (operational metrics, low FPR metrics, top peaks, feature
     importance, random split);
   - phase 3: `evaluate_models_h5.py --species HUMAN|MOUSE`, then `assemble_evaluation_table_h3a.py` for both species.

   End marker `SERVER_H5_DONE` in `logs/driver.log`.
4. Release models: `$ARCHI_PY 2_training/train_models_v2d.py <out_dir> <threads_per_worker> <workers>`.

The Fig. 4 scans and the runtime benchmark are in `../3_analysis/server/README.md`.

## Software

Versions in `VERSIONS.txt`: Python 3.8.18, scikit-learn 1.3.0, pandas 2.0.3, numpy 1.24.3, scipy 1.10.1, joblib 1.2.0,
xgboost 1.7.3, pybedtools 0.9.1; R with PRROC and rpy2 for auROC and auPRC (rpy2 must load the R that holds PRROC, set
`R_HOME`). RF settings: `RandomForestClassifier(max_depth=6, max_samples=0.8, n_estimators=100, random_state=0)` with
the StandardScaler fitted on the training matrix.
