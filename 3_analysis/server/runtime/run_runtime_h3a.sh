#!/bin/bash
# Sup. Table 6 runtime benchmark on the human training matrices (negatives v3).
# Settings: OMP/OPENBLAS/MKL = 4, n_jobs 1 and 4, nice -n 10, /usr/bin/time -f "%e %U %S %M", SARUS JVM -XX:ActiveProcessorCount=1.
# Blocks: A load / train / predict, B chunked inference, C published matrices retrained the same day as a machine-state control,
# D single-PWM scans. Called by make_h5_runtime_v2.sh.
set -uo pipefail
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
H=~/archi_work/H5_alltrain/human_sym/runtime_v2
PY=~/miniconda3/envs/rpy2_env/bin/python
TF_FMT="%e %U %S %M"
( while true; do echo "$(date +%T) $(cut -d' ' -f1-4 /proc/loadavg)"; sleep 15; done ) > $H/load_during_run.log &
LOGGER=$!
snap() { echo "--- $(date '+%F %T') $1"; uptime; ps -eo user,pid,pcpu,etime,comm --sort=-pcpu | head -12; }
snap START > $H/machine_state.txt
{ hostname; lscpu | grep -E "Model name|Socket|Core|Thread|^CPU\(s\)"; free -g | head -2; uname -r; } > $H/machine.txt

# A. load / train / predict (rf_bench_h3a.py)
export BENCH_R=$H BENCH_CSV="$H/data/{TF}_H3a_human_train_mono.csv"
TSV=$H/rf_process_timings.tsv
echo -e "TF\tstep\tn_jobs\twall_s\tuser_s\tsys_s\tmax_rss_kb" > $TSV
for T in STAT3 HNF4A CTCF; do
  for spec in "load 1" "train 1" "train 4" "predict 1" "predict 4"; do
    set -- $spec; echo "$(date +%T) A $T $1 n_jobs=$2"
    /usr/bin/time -f "$TF_FMT" -o $H/scratch/rf_t.txt nice -n 10 $PY $H/code/rf_bench_h3a.py $T $1 $2
    echo -e "$T\t$1\t$2\t$(tail -1 $H/scratch/rf_t.txt | tr ' ' '\t')" >> $TSV
  done
done
# B. chunked inference
TSV=$H/rf_predict_chunked_process.tsv
echo -e "TF\tn_jobs\twall_s\tuser_s\tsys_s\tmax_rss_kb" > $TSV
for T in STAT3 HNF4A CTCF; do for J in 1 4; do
  echo "$(date +%T) B $T chunked n_jobs=$J"
  /usr/bin/time -f "$TF_FMT" -o $H/scratch/rf_t.txt nice -n 10 $PY $H/code/rf_predict_chunked_h3a.py $T $J
  echo -e "$T\t$J\t$(tail -1 $H/scratch/rf_t.txt | tr ' ' '\t')" >> $TSV
done; done
snap AFTER_H3A >> $H/machine_state.txt
# C. same-day control: published matrices, train n_jobs 1 and 4 (models to calib_round1/, never to ~/archi_work/runtime)
mkdir -p $H/calib_round1/models $H/calib_round1/scratch
export BENCH_R=$H/calib_round1 BENCH_CSV="~/archi_release/TF-ML/outputdir/{TF}/result_H_train_full.csv"
TSV=$H/calib_round1/rf_process_timings.tsv
echo -e "TF\tstep\tn_jobs\twall_s\tuser_s\tsys_s\tmax_rss_kb" > $TSV
for T in STAT3 HNF4A CTCF; do for J in 1 4; do
  echo "$(date +%T) C $T train n_jobs=$J (round 1)"
  /usr/bin/time -f "$TF_FMT" -o $H/scratch/rf_t.txt nice -n 10 $PY $H/code/rf_bench_h3a.py $T train $J
  echo -e "$T\ttrain\t$J\t$(tail -1 $H/scratch/rf_t.txt | tr ' ' '\t')" >> $TSV
done; done
# D. SARUS single-PWM scans on the published benchmark FASTA (input independent of the training data)
JAR=~/sarus-2.0.2.jar; BENCH=~/archi_work/runtime/data/bench_50k_300bp.fa
JVM="java -XX:ActiveProcessorCount=1 -Xmx2G -cp $JAR"; OPTS="--skipn --show-non-matching --output-scoring-mode score besthit"
TSV=$H/sarus_timings.tsv; S=$H/scratch
echo -e "task\tTF\tmode\tk_pwms\tn_seqs\tMb\twall_s\tuser_s\tsys_s\tmax_rss_kb\tjvm_flags" > $TSV
for T in STAT3 HNF4A CTCF; do
  M1=$(ls ~/archi_release/TF-ML/outputdir/$T/HUMAN_seq_HUMAN_pwm_mono/*.pwm | sort -V | head -1)
  D1=$(ls ~/archi_release/TF-ML/outputdir/$T/HUMAN_seq_HUMAN_pwm_di/*.dpwm | sort -V | head -1)
  echo "$(date +%T) D $T $(basename $M1) $(basename $D1)"
  for i in 1 2 3; do
    /usr/bin/time -f "$TF_FMT" -o $S/t.txt nice -n 10 $JVM ru.autosome.SARUS $BENCH $M1 $OPTS 2>/dev/null | grep -v ">" > $S/one_mono_$T.tab
    echo -e "scan_1pwm\t$T\tmono\t1\t50000\t15.00\t$(tail -1 $S/t.txt | tr ' ' '\t')\tAP1" >> $TSV
  done
  for i in 1 2 3; do
    /usr/bin/time -f "$TF_FMT" -o $S/t.txt nice -n 10 $JVM ru.autosome.di.SARUS $BENCH $D1 $OPTS 2>/dev/null | grep -v ">" > $S/one_di_$T.tab
    echo -e "scan_1pwm\t$T\tdi\t1\t50000\t15.00\t$(tail -1 $S/t.txt | tr ' ' '\t')\tAP1" >> $TSV
  done
  cmp -s $S/one_mono_$T.tab ~/archi_work/runtime/scratch/one_mono_$T.tab && echo "D $T mono output identical to round 1" || echo "D $T mono output DIFFERS from round 1"
  cmp -s $S/one_di_$T.tab ~/archi_work/runtime/scratch/one_di_$T.tab && echo "D $T di output identical to round 1" || echo "D $T di output DIFFERS from round 1"
done
snap END >> $H/machine_state.txt
kill $LOGGER
echo "$(date +%T) RUNTIME_H3A_DONE"
