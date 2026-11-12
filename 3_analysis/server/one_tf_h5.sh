#!/bin/bash
# One TF: scans (scan_tf_h5.sh) then refit (run_fig4_h5.py, 1 thread); log per TF in logs/v2_<TF>.log. Usage: one_tf_h5.sh <TF>
. "$(dirname "$(readlink -f "$0")")/../../archi_env.sh"
t=$1; H2=$ARCHI_WORK/H5_alltrain/human_sym; P=$PY; A=$ARCHI_CODE/3_analysis/server; cd $H2 || exit 1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
[ -s fig4_refit_trainscale/rows/$t.json ] && { echo "$t: row exists, skip"; exit 0; }
L=logs/v2_${t}_$(date +%s).log
{ t0=$(date +%s); echo "$(date +%F_%T) $t: scans"; $A/scan_tf_h5.sh $t || { echo "ERROR scans $t"; exit 1; }
  echo "$(date +%F_%T) $t: refit"; nice -n 10 $P $A/run_fig4_h5.py --tfs $t --rows-dir fig4_refit_trainscale/rows || { echo "ERROR refit $t"; exit 1; }
  echo "$(date +%F_%T) $t: done $(( $(date +%s)-t0 ))s"; } > $L 2>&1
tail -n 1 $L
