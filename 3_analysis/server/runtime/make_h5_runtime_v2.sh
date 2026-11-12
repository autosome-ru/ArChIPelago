#!/bin/bash
# Sup. Table 6 runtime benchmark on the current training matrices: builds the three matrices (build_h3a_matrices.py, loader
# archi_data_v5), runs the benchmark (run_runtime_h3a.sh) and assembles the table (assemble_sup_table6_h3a.py) in
# $ARCHI_WORK/H5_alltrain/human_sym/runtime_v2. Run when nothing else runs on the machine (timing). Log runtime/run_runtime.log ends
# RUNTIME_H3A_DONE. Usage: bash make_h5_runtime_v2.sh
set -e
. "$(dirname "$(readlink -f "$0")")/../../../archi_env.sh"
CODE=$ARCHI_CODE/3_analysis/server/runtime; H=$ARCHI_WORK/H5_alltrain/human_sym/runtime_v2
[ ! -e $H ] || { echo "exists $H"; exit 1; }
mkdir -p $H/scratch $H/data
cd $H
nice -n 10 $PY $CODE/build_h3a_matrices.py data STAT3 HNF4A CTCF > build_matrices.log 2>&1
bash $CODE/run_runtime_h3a.sh > run_runtime.log 2>&1
$PY $CODE/assemble_sup_table6_h3a.py $H $ARCHI_WORK/runtime > assemble.log 2>&1
echo RUNTIME_H5_DONE >> run_runtime.log
