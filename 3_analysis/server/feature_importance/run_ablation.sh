#!/bin/bash
# Ablation of the secondary motif subtype (Fig. S10C, Sup. Table 9 sheet 'ablation'): one driver for one similarity cut.
# usage: ARCHI_WORK=<work root> bash run_ablation.sh <out root> <cut x100> <workers> <TF> [<TF> ...]
#   <out root> holds the existing folders ablation/ and logs/; each TF is one process with 8 threads, niced,
#   all four conditions with seeds 0-9 (--n-seeds 10). An existing <TF>_ablation.csv is never overwritten.
# The paper: cut 80 for the 20 TFs with at least two subtypes at 0.8, cut 90 for the 9 TFs with one subtype at 0.8 and
# at least two at 0.9 (the selection is made by 2_training/h5_server_v2.sh from the *_sim_summary.csv files).
T=$1; THR=$2; W=$3; shift 3
. "$(dirname "$(readlink -f "$0")")/../../../archi_env.sh"
B4=$ARCHI_CODE/3_analysis/server/feature_importance
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
echo "$(date +%F_%T) start cut=$THR workers=$W n_seeds=10 script_md5=$(md5sum $B4/b4_04_ablation.py | cut -c1-32) TFs: $*" >> $T/logs/driver.log
printf "%s\n" "$@" | xargs -P $W -I{} sh -c "nice -n 10 $PY $B4/b4_04_ablation.py --tfs {} --n-jobs 8 --thr $THR --n-seeds 10 --out-dir $T/ablation > $T/logs/{}_cut$THR.log 2>&1 || echo \"\$(date +%F_%T) FAILED {}\" >> $T/logs/driver.log"
echo "$(date +%F_%T) cut${THR}_DONE files=$(ls $T/ablation/*_ablation.csv 2>/dev/null | wc -l)" >> $T/logs/driver.log
