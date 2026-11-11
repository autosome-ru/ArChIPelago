#!/bin/bash
# Fig. 4 / S4 on server: per TF (once the relay wrote slim/models_server/<TF>/COPY_VERIFIED.txt) Slim m=0/1/-5 and
# diChIPMunk scans of the training table, human test and mouse test (code/scan_tf_h5.sh), then the RF2f-family refit with
# the training scaler (code/run_fig4_h5.py) -> fig4_refit_trainscale/rows; 24 TFs in parallel, 1 thread each. Then the assembly
# with the hm_mm results (assemble_fig4_s4.py, s4_numbers.py). Log logs/fig4_driver.log, end marker FIG4_H5_DONE.
# Starts at once (does not wait for SERVER_H5_DONE); per-TF jobs start as soon as a TF is verified, at most 4 in parallel
# while the RF driver runs (32-thread budget), 20 after SERVER_H5_DONE; the assembly waits for SERVER_H5_DONE (needs hm_mm).
# Usage: bash h5_fig4_server_v2.sh
R=~/archi_work; G=$R/H5_alltrain/human_sym; G4=$R/human_sym/H4_trainscale; cd $G || exit 1
P="nice -n 10 $HOME/miniconda3/envs/rpy2_env/bin/python"; HM=$G/hm_mm_results.csv
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
log() { echo "$(date +%F_%T) $*" >> logs/fig4_driver.log; }
TFS=$(tail -n +2 $R/human_sym/inputs/human_sym_peak_ids.csv | cut -d, -f1 | sort -u)
mkdir -p fig4_refit_trainscale/rows slim/scans
log start
until [ $(ls slim/models_server/*/COPY_VERIFIED.txt 2>/dev/null | wc -l) -ge 35 ]; do
  for t in $TFS; do                       # start every verified TF that has no row and no running job
    [ -e slim/models_server/$t/COPY_VERIFIED.txt ] && [ ! -e logs/fig4_started_$t ] && echo $t
  done > logs/fig4_ready.txt
  NP=4; grep -q SERVER_H5_DONE logs/driver.log && NP=20
  NR=$(ps -u $USER -o args | grep -c "[o]ne_tf_h5.sh")
  for t in $(head -n $(( NP > NR ? NP - NR : 0 )) logs/fig4_ready.txt); do
    touch logs/fig4_started_$t; nohup code/one_tf_h5.sh $t >> logs/fig4_one_tf.log 2>&1 &
  done
  sleep 300
done
until grep -q "SERVER_H5_DONE\|FAILED" logs/driver.log; do sleep 300; done
for t in $TFS; do [ -e logs/fig4_started_$t ] || { touch logs/fig4_started_$t; echo $t; }; done | xargs -r -P 20 -L 1 code/one_tf_h5.sh >> logs/fig4_one_tf.log 2>&1
while ps -u $USER -o args | grep -q "[o]ne_tf_h5.sh"; do sleep 60; done
n=$(ls fig4_refit_trainscale/rows/*.json 2>/dev/null | wc -l); log "rows $n/35"
[ $n -eq 35 ] || { log "FAILED: rows $n/35 (see logs/fig4_one_tf.log, logs/v2_<TF>_*.log)"; exit 1; }
REF="h4_d10a_human_35TF=$G4/fig4_refit_trainscale/Figure_4_source_data.csv h4_d10a_mouse_35TF=$G4/fig4_refit_trainscale/Figure_S4_source_data.csv"
$P code/assemble_fig4_s4.py fig4_refit_trainscale/rows $HM fig4_refit_trainscale/Figure_4_source_data.csv 35 H > logs/fig4_assemble.log 2>&1 && \
$P code/assemble_fig4_s4.py fig4_refit_trainscale/rows $HM fig4_refit_trainscale/Figure_S4_source_data.csv 35 M >> logs/fig4_assemble.log 2>&1 && \
$P code/s4_numbers.py fig4_refit_trainscale/Figure_4_and_S4_numbers.csv new_human_35TF=fig4_refit_trainscale/Figure_4_source_data.csv new_mouse_35TF=fig4_refit_trainscale/Figure_S4_source_data.csv $REF >> logs/fig4_assemble.log 2>&1 \
  && log FIG4_H5_DONE || log "FAILED assembly"
