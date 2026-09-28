#!/bin/bash
# Slim scan of one sequence set with one Slim model of a TF (inputs of run_fig4.py).
# usage: scan_one.sh <TF> <m> <set> [model_xml]
#   <m>   Slim model order as in the pipeline folder name (0, 1, 5 for LSlim m=-5)
#   <set> global sequence set: HUMAN_10000_train, HUMAN_10000_control, or the mouse test set
#         (MOUSE_10000_control; MOUSE_10000_train in a pipeline directory with ARCHI_MOUSE_FILES_SWAPPED=1)
# -> $ARCHI_FIG4_SCANS_DIR/slim_scans/<TF>/m<m>_<set>/model1_predictions.txt
# Environment: ARCHI_RELEASE_DIR (pipeline output, default ~/Release/TF-ML), SLIM_DIR (folder with
# TrainAndApplySlim.jar and jdk8u232-b09/, default ../../../Slim of this repository), ARCHI_FIG4_SCANS_DIR (default this folder).
HERE=$(cd "$(dirname "$0")" && pwd)
TFML=${ARCHI_RELEASE_DIR:-~/Release/TF-ML}
SLIM=${SLIM_DIR:-$HERE/../../../Slim}
SCANS=${ARCHI_FIG4_SCANS_DIR:-$HERE}
TF=$1; M=$2; SET=$3; XML=${4:-$TFML/outputdir/$TF/${TF}_SlimModel_$M/Motif_1/SlimDimont_1.xml}
out=$SCANS/slim_scans/$TF/m${M}_$SET
if [ -s $out/model1_predictions.txt ]; then echo "skip $TF m$M $SET"; exit 0; fi
if [ -e $out/predictions.txt ]; then echo "ERROR $out/predictions.txt exists but no model1_predictions.txt"; exit 1; fi
[ -s "$XML" ] || { echo "ERROR no model $XML"; exit 1; }
mkdir -p $out
t0=$(date +%s)
nice -n 10 $SLIM/jdk8u232-b09/bin/java -Xmx4G -XX:ParallelGCThreads=1 -jar $SLIM/TrainAndApplySlim.jar scan i=$TFML/all_mfa_file_${SET}_no_NF_no_N.fasta m=$XML outdir=$out > $out/log.txt 2>&1
rc=$?
[ $rc -eq 0 ] && [ -s $out/predictions.txt ] && ln $out/predictions.txt $out/model1_predictions.txt
echo "$TF m$M $SET rc=$rc $(($(date +%s)-t0))s $(wc -l < $out/predictions.txt) rows"
