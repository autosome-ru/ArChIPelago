#!/bin/bash
# Server driver: RF training and evaluation (human->mouse and mouse->mouse), saturation, operational metrics, low-FPR
# metrics, top peaks, feature importance, random split and the evaluation tables. The Slim scans / Fig. 4 refit and the
# runtime table are in 3_analysis/server (h5_fig4_server_v2.sh, runtime/make_h5_runtime_v2.sh). Run from
# ~/archi_work/H5_alltrain/human_sym with: bash h5_server_v2.sh
# Thread budget 32:
#  phase 1: run_hm_mm 4 x 4 = 16 | saturation HUMAN 2 x 4 = 8 | saturation MOUSE 2 x 4 = 8
#  phase 2 (after run_hm_mm): chain A (operational, B2 v2, B3 1 thread; B4 8 threads) | B9 8 threads  (+ saturations)
#  phase 3 (after phases 1-2): evaluate HM 4 x 3 | evaluate MM 4 x 3, then both evaluation tables
# Evaluation runs 4 x 3 threads so that the per-TF Slim scans / Fig. 4 refits (h5_fig4_server_v2.sh, 4 x 1 thread while
# this driver runs) fit in the 32-thread budget.
# Producers refuse to overwrite; a failed step stops its chain. Log G5/logs/driver.log, end marker SERVER_H5_DONE.
R=~/archi_work; G=$R/H5_alltrain/human_sym; GM=$R/H5_alltrain/mouse_sym; DH=$R/human_sym_h5; DM=$R/mouse_sym_h5
cd $G || exit 1
PY=$HOME/miniconda3/envs/rpy2_env/bin/python; P="nice -n 10 $PY"
TFS=$(tail -n +2 $R/human_sym/inputs/human_sym_peak_ids.csv | cut -d, -f1 | sort -u | tr "\n" " ")
log() { echo "$(date +%F_%T) $*" >> $G/logs/driver.log; }
die() { log "FAILED $*"; exit 1; }
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
[ $(echo $TFS | wc -w) -eq 35 ] || die "TF list"
mkdir -p $GM/logs && mkdir models || die mkdir
log start
sat() {  # $1 dir with run_saturation_h5.py, $2 log
  (cd $1 && OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 $P run_saturation_h5.py --tfs $TFS --workers 2 --rf-jobs 4; echo SAT_EXIT $?) > $2 2>&1
  grep -q "SAT_EXIT 0" $2
}
(sat $G/saturation_HUMAN $G/logs/saturation_HUMAN.log && log SAT_HUMAN_DONE || log "FAILED saturation_HUMAN") &
(sat $GM/saturation_MOUSE $GM/logs/saturation_MOUSE.log && log SAT_MOUSE_DONE || log "FAILED saturation_MOUSE") &
(time $P code/run_hm_mm_h5.py --tfs $TFS --workers 4 --rf-jobs 4 --out-csv hm_mm_results.csv --rows-dir rows --scores-dir scores; echo HMMM_EXIT $?) > logs/run_hm_mm.log 2>&1
grep -q "HMMM_EXIT 0" logs/run_hm_mm.log || die run_hm_mm
log run_hm_mm_done
chainA() {
  mkdir operational b2_v2/tables b3 || die mkdir_A
  ln -s $DH/TF-ML/out_tab_HUMAN_10000_control.tab b2_v2/tables/ && ln -s $DM/TF-ML/out_tab_MOUSE_10000_train.tab b2_v2/tables/ || die ln
  $P code/operational_metrics.py scores operational/operational_metrics.csv > logs/operational.log 2>&1 || die operational
  $P b2_v2/b2_compute.py --scores scores --tables b2_v2/tables --out b2_v2 > logs/b2_v2_compute.log 2>&1 || die b2
  $P code/b3_topN_h5.py b3 > logs/b3_topN.log 2>&1 || die b3
  log downstream_done
  mkdir B4_importance/cache && cd B4_importance/scripts || die mkdir_B4
  OMP_NUM_THREADS=8 $P b4_01_importance.py --tfs $TFS --n-jobs 8 > ../../logs/b4_01.log 2>&1 || die b4_01
  $P b4_02_subtypes.py --tfs $TFS --workers 8 > ../../logs/b4_02.log 2>&1 || die b4_02
  OMP_NUM_THREADS=8 $P b4_03_usage.py --tfs $TFS --n-jobs 8 > ../../logs/b4_03.log 2>&1 || die b4_03
  ABL=$($PY -c "
import pandas as pd, glob
s = pd.concat([pd.read_csv(f) for f in glob.glob('../subtypes/*_sim_summary.csv')]).set_index('TF')
new = set(s.index[s.n_subtypes_ge10pct_80 >= 2]); old = {'IRF4','TAL1','E2F4','TFE2','SRF','PRGR','SOX2','STA5A','RXRA'}
print(' '.join(sorted(t for t in new | old if s.loc[t, 'n_subtypes_80'] >= 2)))")
  ABL9=$($PY -c "
import pandas as pd, glob
s = pd.concat([pd.read_csv(f) for f in glob.glob('../subtypes/*_sim_summary.csv')]).set_index('TF')
print(' '.join(t for t in ('MAFK', 'SPI1', 'IRF4','TAL1','E2F4','TFE2','SRF','PRGR','SOX2','STA5A','RXRA') if s.loc[t, 'n_subtypes_80'] < 2 and s.loc[t, 'n_subtypes_90'] >= 2))")
  echo "ablation thr 0.8: $ABL | thr 0.9: $ABL9" > ../../logs/b4_04_tf_lists.txt
  OMP_NUM_THREADS=8 $P b4_04_ablation.py --tfs $ABL --n-jobs 8 --thr 80 > ../../logs/b4_04.log 2>&1 || die b4_04
  if [ -n "$ABL9" ]; then OMP_NUM_THREADS=8 $P b4_04_ablation.py --tfs $ABL9 --n-jobs 8 --thr 90 > ../../logs/b4_04_thr90.log 2>&1 || die b4_04_thr90; fi
  $P b4_07_checks.py > ../../logs/b4_07.log 2>&1 || die b4_07
  cd $G; log CHAIN_A_DONE
}
chainB9() {
  mkdir b9 || die mkdir_b9
  N_JOBS=8 OMP_NUM_THREADS=8 $P code/b9_random_split_h5.py b9 PRGR RXRA SPI1 ANDR REST > logs/b9_h5.log 2>&1 || die b9
  $P code/b9_summary_h2b.py b9/b9_random_split_per_fit.csv b9/b9_summary.csv > logs/b9_summary.log 2>&1 || die b9_summary
  log CHAIN_B9_DONE
}
chainA & A=$!
chainB9 & B=$!
wait $A; wait $B
until grep -q "SAT_HUMAN_DONE\|FAILED saturation_HUMAN" logs/driver.log && grep -q "SAT_MOUSE_DONE\|FAILED saturation_MOUSE" logs/driver.log; do sleep 60; done
grep -q "CHAIN_A_DONE" logs/driver.log && grep -q CHAIN_B9_DONE logs/driver.log && grep -q SAT_HUMAN_DONE logs/driver.log && grep -q SAT_MOUSE_DONE logs/driver.log || die "phase 1-2 incomplete"
log phase12_done
mkdir eval_models $GM/eval_models || die mkdir_eval
(time $P code/evaluate_models_h5.py --species HUMAN --workers 4 --threads-per-worker 3 --tfs $TFS --rows-dir eval_models/rows_HM; echo EVAL_HM_EXIT $?) > logs/eval_models_HM.log 2>&1 &
(time $P code/evaluate_models_h5.py --species MOUSE --workers 4 --threads-per-worker 3 --tfs $TFS --rows-dir $GM/eval_models/rows_MM; echo EVAL_MM_EXIT $?) > $GM/logs/eval_models_MM.log 2>&1 &
wait
grep -q "EVAL_HM_EXIT 0" logs/eval_models_HM.log && grep -q "EVAL_MM_EXIT 0" $GM/logs/eval_models_MM.log || die evaluate
mkdir eval_models/HM $GM/eval_models/MM || die mkdir_tables
NB4=$R/mouse_sym/w1g1/orig/notebook4_results_table.csv
$P code/assemble_evaluation_table_h3a.py --mouse-trained --rows-dir eval_models/rows_HM --notebook4-table $NB4 --out-dir eval_models/HM --exclude-tfs P53 > logs/assemble_HM.log 2>&1 || die assemble_HM
$P code/assemble_evaluation_table_h3a.py --mouse-trained --rows-dir $GM/eval_models/rows_MM --notebook4-table $NB4 --out-dir $GM/eval_models/MM --exclude-tfs P53 > $GM/logs/assemble_MM.log 2>&1 || die assemble_MM
log SERVER_H5_DONE
