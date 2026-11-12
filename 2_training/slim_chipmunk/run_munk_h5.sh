#!/bin/bash
# diChIPMunk training as notebook 2 l.2982-2987
# (run_dichiphorde8.rb 15:15 filter y 1.0 m:<pos> 200 20 1 4 random auto flat); ARCHI_CHIPMUNK = folder with run_dichiphorde8.rb, chipmunk.jar and ytilib/.
# Input = the human training positives P5/munk_src/<TF>/neg_v3_work_train/train_pos_id_HUMAN.fasta, headers '>150'.
# usage: run_munk_h5.sh <TF>
. "$(dirname "$(readlink -f "$0")")/../../archi_env.sh"
T=${ARCHI_CHIPMUNK:?CHECK FAILED: ARCHI_CHIPMUNK (folder with run_dichiphorde8.rb, chipmunk.jar, ytilib/) is not set}
TF=$1; P5=$ARCHI_WORK/H5_alltrain
D=$P5/models/$TF; P=$D/HUMAN_seq_HUMAN_pwm_mono; SRC=$P5/munk_src/$TF/neg_v3_work_train/train_pos_id_HUMAN.fasta
if [ -s $P/train_H_1_ChIPMunk_no_repeats_0.dpwm ]; then echo "skip munk $TF"; exit 0; fi
[ -s $SRC ] || { echo "ERROR munk input missing $TF"; exit 1; }
mkdir -p $D; mkdir $P 2>/dev/null || { echo "busy munk $TF (folder exists)"; exit 0; }
cd $D || exit 1; t0=$(date +%s)
awk '/^>/ {FILENAME=150;printf(">%s\n",FILENAME);next;} {print}' $SRC > train_pos_id_HUMAN_p.fasta
nice -n 10 ruby $T/run_dichiphorde8.rb $P/train_H_1_ChIPMunk_no_repeats 15:15 filter y 1.0 m:$D/train_pos_id_HUMAN_p.fasta 200 20 1 4 random auto flat > ChIPMunk.txt 2>&1
rc=$?; echo "munk $TF rc=$rc $(( $(date +%s)-t0 ))s $(date +%F_%T) $(ls $P/*.dpwm 2>/dev/null | wc -l) dpwm"
