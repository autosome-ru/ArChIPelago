#!/bin/bash
# Data build on server, human and mouse (run from anywhere; new folders only, every producer refuses to overwrite; nice 10).
# Steps:
#  1 positives per_tf + pooled (build_positives_v5.py: chromosome lists from h5_splits.py)
#  2 test-set identity of the positives vs the earlier build (check_h5_test_identity.py) -> stop on any difference
#  3 negatives v3 (make_negatives_v3.py; 140 jobs, 12 in parallel)
#  4 test-set identity of the negatives -> stop on any difference
#  5 SARUS feature tables (make_scan_jobs_sym.py; 2 x 7,992 jobs, 24 in parallel)
#  6 test-set identity of the feature tables
# Usage: bash h5_data_server.sh. Log: $DH/logs/data_driver.log, end marker H5_DATA_DONE.
R=~/archi_work; DH=$R/human_sym_h5; DM=$R/mouse_sym_h5; C=$DH/code
PY=$HOME/miniconda3/envs/rpy2_env/bin/python; P="nice -n 10 $PY"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
log() { echo "$(date +%F_%T) $*" >> $DH/logs/data_driver.log; }
die() { log "FAILED $*"; exit 1; }
TFS=$(tail -n +2 $R/human_sym/inputs/human_sym_peak_ids.csv | cut -d, -f1 | sort -u | tr "\n" " ")
[ $(echo $TFS | wc -w) -eq 35 ] || die "TF list"
log "start; TFs: $(echo $TFS | wc -w)"

# 1 positives
build() {  # $1 species, $2 root, $3 peaks, $4 mask
  cd $2 || return 1
  $P $C/build_positives_v5.py per_tf --peaks $3 --mask $4 --out TF-ML --species $1 > logs/build_positives_per_tf.log 2>&1 || return 1
  $P $C/build_positives_v5.py pooled --mfa-dir TF-ML --out TF-ML --species $1 > logs/build_positives_pooled.log 2>&1 || return 1
}
build HUMAN $DH $R/human_sym/inputs/human_sym_peak_ids.csv $HOME/Cleaning_data/track_out_HUMAN.bed & PH=$!
build MOUSE $DM $R/mouse_sym/inputs/mouse_sym_peak_ids.csv $R/mouse_sym/ref/track_out_MOUSE_mm10.bed & PM=$!
wait $PH || die build_HUMAN
wait $PM || die build_MOUSE
grep -q "POOLED DONE" $DH/logs/build_positives_pooled.log && grep -q "POOLED DONE" $DM/logs/build_positives_pooled.log || die pooled
log positives_done
$P $C/check_h5_test_identity.py positives $DH/checks > $DH/logs/test_identity_positives.log 2>&1 || die "STOP test positives differ from H4"
log test_identity_positives_ok

# 3 negatives v3
for t in $TFS; do for s in train control; do echo "HUMAN $DH $t $s"; echo "MOUSE $DM $t $s"; done; done > $DH/logs/neg_v3_jobs.txt
echo "$(date) start" > $DH/logs/negatives_v3.log
xargs -P 12 -L 1 sh -c 'cd $1 && o=$(nice -n 10 '$PY' '$C'/make_negatives_v3.py --tfml TF-ML --outdir TF-ML/outputdir --tf $2 --split $3 --species $0) && echo "$0 $o" || echo "NEG_FAILED $0 $2 $3"' \
  < $DH/logs/neg_v3_jobs.txt >> $DH/logs/negatives_v3.log 2>&1
grep -q NEG_FAILED $DH/logs/negatives_v3.log && die negatives
[ $(grep -c "neg_out=" $DH/logs/negatives_v3.log) -eq 140 ] || die "negatives: not 140 jobs"
log negatives_done
$P $C/check_h5_test_identity.py negatives $DH/checks > $DH/logs/test_identity_negatives.log 2>&1 || die "STOP test negatives differ from H4"
log test_identity_negatives_ok

# 3b Slim inputs (loader v5 selection, build_slim_inputs_v5.py --export) and the server shipment, in the background
slim_inputs() {
  S=$R/H5_alltrain/human_sym/slim; mkdir -p $S/inputs $S/ship || return 1
  echo $TFS | tr " " "\n" | grep -v '^$' | xargs -P 4 -I{} sh -c "$P $C/build_slim_inputs_v5.py $C {} $S/inputs/{} --export || echo SLIM_IN_FAILED {}" > $S/build_slim_inputs.log 2>&1
  grep -q SLIM_IN_FAILED $S/build_slim_inputs.log && return 1
  [ $(grep -c " positives " $S/build_slim_inputs.log) -eq 35 ] || return 1
  cd $S/ship || return 1
  gzip -c $DH/TF-ML/all_mfa_file_HUMAN_SLIM_10000_train_no_NF_no_N.fasta > all_mfa_file_HUMAN_SLIM_10000_train_no_NF_no_N.fasta.gz || return 1
  gzip -c $DH/TF-ML/out_tab_HUMAN_10000_train.tab > out_tab_HUMAN_10000_train.tab.gz || return 1
  (cd $S/inputs && tar czf $S/ship/selections.tar.gz */selection.npz) || return 1
  (cd $DH/TF-ML/outputdir && tar czf $S/ship/munk_src.tar.gz */neg_v3_work_train/train_pos_id_HUMAN.fasta) || return 1
  (cd $S/inputs && md5sum */train_H_1_SLIM.fasta */train_H_0_SLIM.fasta */selection.npz) > $S/inputs_md5.txt || return 1
  (cd $DH/TF-ML/outputdir && md5sum */neg_v3_work_train/train_pos_id_HUMAN.fasta) > $S/munk_src_md5.txt || return 1
  md5sum *.gz > $S/ship_MANIFEST_md5.txt
}
(slim_inputs && log SLIM_INPUTS_DONE || log "FAILED slim_inputs") &

# 5 SARUS
(cd $DH && $PY $C/make_scan_jobs_sym.py --tfml TF-ML --tfs $TFS --species HUMAN > logs/scan_jobs.txt 2> logs/scan_jobs.err) || die scan_jobs_H
(cd $DM && $PY $C/make_scan_jobs_sym.py --tfml TF-ML --tfs $TFS --species MOUSE > logs/scan_jobs.txt 2> logs/scan_jobs.err) || die scan_jobs_M
{ sed "s|^|cd $DH \&\& nice -n 10 |" $DH/logs/scan_jobs.txt; sed "s|^|cd $DM \&\& nice -n 10 |" $DM/logs/scan_jobs.txt; } > $DH/logs/scan_jobs_all.txt
echo "$(date) start $(wc -l < $DH/logs/scan_jobs_all.txt) jobs" > $DH/logs/scan_run.log
xargs -P 24 -d '\n' -n 1 bash -c 'set -o pipefail; eval "$0" || echo "SCAN_FAILED $0"' < $DH/logs/scan_jobs_all.txt >> $DH/logs/scan_run.log 2>&1
grep -q SCAN_FAILED $DH/logs/scan_run.log && die sarus
echo "$(date) end" >> $DH/logs/scan_run.log
log sarus_done
$P $C/check_h5_test_identity.py sarus $DH/checks > $DH/logs/test_identity_sarus.log 2>&1 || die "STOP test feature tables differ from H4"
log test_identity_sarus_ok
wait
grep -q SLIM_INPUTS_DONE $DH/logs/data_driver.log || die slim_inputs
log H5_DATA_DONE
