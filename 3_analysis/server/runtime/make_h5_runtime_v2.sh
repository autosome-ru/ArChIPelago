#!/bin/bash
# Sup. Table 6 runtime benchmark on the current training matrices: copies the code of the benchmark (build_h3a_matrices.py,
# run_runtime_h3a.sh, ...) with only the loader (archi_data_v5) and the output root (H5_alltrain/human_sym/runtime_v2) changed,
# runs it, and assembles the table (assemble_sup_table6_h3a.py). The sed edits use "#" or "|" as delimiter because the
# replacement text holds "#". Run on server when nothing else of ours runs (timing). Log runtime/run_runtime.log ends
# RUNTIME_H3A_DONE. Usage: bash make_h5_runtime_v2.sh
set -e
R=~/archi_work; S=$R/human_sym/H3a/runtime/code; H=$R/H5_alltrain/human_sym/runtime_v2; PY=~/miniconda3/envs/rpy2_env/bin/python
[ ! -e $H ] || { echo "exists $H"; exit 1; }
mkdir -p $H/code $H/scratch $H/data
sed -e 's#os.path.expanduser("~/archi_work/human_sym/code")#os.path.expanduser("~/archi_work/human_sym_h5/code")#' \
    -e 's#import archi_data_v3 as ad#import archi_data_v5 as ad#' $S/build_h3a_matrices.py > $H/code/build_h3a_matrices.py
sed -e 's|^H=~/archi_work/human_sym/H3a/runtime$|H=~/archi_work/H5_alltrain/human_sym/runtime_v2   # H5 (D16): only change|' $S/run_runtime_h3a.sh > $H/code/run_runtime_h3a.sh
cp $S/rf_bench_h3a.py $S/rf_predict_chunked_h3a.py $S/assemble_sup_table6_h3a.py $H/code/
for f in build_h3a_matrices.py run_runtime_h3a.sh; do diff $S/$f $H/code/$f > $H/code/$f.diff || true; done
[ $(grep -c "^[<>]" $H/code/build_h3a_matrices.py.diff) -eq 4 ] && [ $(grep -c "^[<>]" $H/code/run_runtime_h3a.sh.diff) -eq 2 ] || { echo "unexpected diff"; exit 1; }
cd $H
nice -n 10 $PY code/build_h3a_matrices.py data STAT3 HNF4A CTCF > build_matrices.log 2>&1
bash code/run_runtime_h3a.sh > run_runtime.log 2>&1
$PY code/assemble_sup_table6_h3a.py $H $R/runtime > assemble.log 2>&1
echo RUNTIME_H5_DONE >> run_runtime.log
