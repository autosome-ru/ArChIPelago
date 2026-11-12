#!/bin/bash
# Slim Dimont training as notebook 2 SlimModel.fit
# (p=peak v=signal m=<o> moobm=-1 w=0.2 Starts=20 threads=4, -XX:ParallelGCThreads=1), JDK 8u232 + TrainAndApplySlim.jar
# (md5 5186b158) from the Slim/ folder of the repository root (unpacked from Slim.tar.gz of Zenodo v1.0.0, not tracked). Inputs P5/inputs/<TF>/train_H_{1,0}_SLIM.fasta. Atomic mkdir claim.
# usage: train_slim_h5.sh <TF> <0|1|5>
. "$(dirname "$(readlink -f "$0")")/../../archi_env.sh"
TF=$1; D=$2; case $D in 0) M=0;; 1) M=1;; 5) M=-5;; *) echo "bad order $D"; exit 1;; esac
P5=$ARCHI_WORK/H5_alltrain; SL=$ARCHI_CODE/Slim
JAVA=$SL/jdk8u232-b09/bin/java; JAR=$SL/TrainAndApplySlim.jar
IN=$P5/inputs/$TF; OUT=$P5/models/$TF/${TF}_SlimModel_$D
if [ -s $OUT/Motif_1/SlimDimont_1.xml ]; then echo "skip $TF m$D"; exit 0; fi
[ -s $IN/train_H_1_SLIM.fasta ] && [ -s $IN/train_H_0_SLIM.fasta ] || { echo "ERROR inputs missing for $TF"; exit 1; }
mkdir -p $(dirname $OUT); mkdir $OUT 2>/dev/null || { echo "busy $TF m$D (folder exists)"; exit 0; }
LOG=$P5/logs/train_${TF}_m${D}.log; t0=$(date +%s)
echo "start $TF m=$M $(date)" > $LOG
nice -n 10 $JAVA -XX:ParallelGCThreads=1 -jar $JAR slimdimont i=$IN/train_H_1_SLIM.fasta b=$IN/train_H_0_SLIM.fasta p=peak v=signal m=$M moobm=-1 w=0.2 Starts=20 outdir=$OUT threads=4 >> $LOG 2>&1
rc=$?; echo "exit code $rc end $(date)" >> $LOG
echo "$TF m$D rc=$rc $(( $(date +%s)-t0 ))s $(date +%F_%T)"
