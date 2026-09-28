#!/bin/bash
# Slim m=1 model of ANDR (slim_ANDR_m1/SlimDimont_1.xml), trained with the protocol of notebook 2 on the
# human training sequences of the TF; then scanned with scan_one.sh (sets HUMAN_10000_train,
# HUMAN_10000_control and the mouse test set, 4th argument = this model).
# Environment: ARCHI_RELEASE_DIR (default ~/Release/TF-ML), SLIM_DIR (default ../../../Slim), ARCHI_FIG4_SCANS_DIR.
HERE=$(cd "$(dirname "$0")" && pwd)
TFML=${ARCHI_RELEASE_DIR:-~/Release/TF-ML}
SLIM=${SLIM_DIR:-$HERE/../../../Slim}
IN=$TFML/outputdir/ANDR
OUT=${ARCHI_FIG4_SCANS_DIR:-$HERE}/slim_ANDR_m1_training
mkdir -p $OUT
echo "start $(date)"
nice -n 10 $SLIM/jdk8u232-b09/bin/java -jar $SLIM/TrainAndApplySlim.jar slimdimont i=$IN/train_H_1_SLIM.fasta b=$IN/train_H_0_SLIM.fasta p=peak v=signal m=1 moobm=-1 w=0.2 Starts=20 outdir=$OUT threads=4
echo "exit code $? end $(date)"
