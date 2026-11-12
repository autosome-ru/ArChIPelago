#!/bin/bash
# Slim and diChIPMunk queue (working folder $ARCHI_WORK/H5_alltrain: bash h5_server_queue.sh).
# (1) unpack the packaged inputs written by 1_data_preparation/h5_data_server.sh (H5_alltrain/human_sym/slim/ship), rebuild the 35 Slim
# inputs from the selections (build_slim_inputs_from_selection.py) and check them against inputs_md5.txt;
# (2) one queue, 24 slots x 4 threads = 96: 105 Slim models longest first, then 35 diChIPMunk.
# Log logs/queue.log, end marker SERVER_H5_QUEUE_DONE (or SERVER_H5_FAILED).
. "$(dirname "$(readlink -f "$0")")/../../archi_env.sh"
P5=$ARCHI_WORK/H5_alltrain; SRC=$P5/human_sym/slim; K=$ARCHI_CODE/2_training/slim_chipmunk; export K; cd $P5 || exit 1
mkdir -p logs
log() { echo "$(date +%F_%T) $*" >> logs/queue.log; }
die() { log "SERVER_H5_FAILED $*"; exit 1; }
mkdir data munk_src inputs models || die mkdir
gunzip -c $SRC/ship/out_tab_HUMAN_10000_train.tab.gz > data/out_tab_HUMAN_10000_train.tab || die gunzip1
gunzip -c $SRC/ship/all_mfa_file_HUMAN_SLIM_10000_train_no_NF_no_N.fasta.gz > data/all_mfa_file_HUMAN_SLIM_10000_train_no_NF_no_N.fasta || die gunzip2
(cd inputs && tar xzf $SRC/ship/selections.tar.gz) || die sel
(cd munk_src && tar xzf $SRC/ship/munk_src.tar.gz) || die munk_src
(cd munk_src && md5sum -c --quiet $SRC/munk_src_md5.txt) || die munk_src_md5
ls inputs | xargs -P 8 -I{} sh -c "nice -n 10 $PY $K/build_slim_inputs_from_selection.py data/out_tab_HUMAN_10000_train.tab data/all_mfa_file_HUMAN_SLIM_10000_train_no_NF_no_N.fasta inputs/{}/selection.npz inputs/{} || echo BUILD_FAILED {}" > logs/build_inputs.log 2>&1
grep -q BUILD_FAILED logs/build_inputs.log && die build_inputs
(cd inputs && md5sum -c --quiet $SRC/inputs_md5.txt) > logs/inputs_md5_check.log 2>&1 || die "inputs md5 differ from inputs_md5.txt"
log "inputs rebuilt, $(grep -c ': OK' <(cd inputs && md5sum -c $SRC/inputs_md5.txt)) files md5-identical to inputs_md5.txt"
{ sed 's/^/slim /' $K/slim_jobs_lpt.txt; sed 's/^/munk /' $K/munk_jobs_lpt.txt; } > logs/jobs.txt
log "queue start: $(wc -l < logs/jobs.txt) jobs, 24 slots"
xargs -P 24 -L 1 sh -c 'if [ "$0" = slim ]; then $K/train_slim_h5.sh $1 $2; else $K/run_munk_h5.sh $1; fi' < logs/jobs.txt >> logs/queue.log 2>&1
ns=$(grep -c "^[A-Z0-9]* m[015] rc=0" logs/queue.log); nm=$(grep -c "^munk .* rc=0 .* 1 dpwm" logs/queue.log)
log "slim rc=0: $ns/105, munk rc=0 with dpwm: $nm/35"
[ $ns -eq 105 ] && [ $nm -eq 35 ] || die "incomplete"
log SERVER_H5_QUEUE_DONE
