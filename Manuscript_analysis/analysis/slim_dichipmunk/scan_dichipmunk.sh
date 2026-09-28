#!/bin/bash
# SARUS scan of the mouse test set with the diChIPMunk diPWM of a TF (built by notebook 2 from the human
# training positives); input of fit_RF2f_models.py.
# usage: scan_dichipmunk.sh <TF> [mouse_test_set]
#   mouse_test_set: MOUSE_10000_control (default; MOUSE_10000_train in a pipeline directory with ARCHI_MOUSE_FILES_SWAPPED=1)
# -> $ARCHI_SLIM_SCANS_DIR/munk/<TF>_full_train_M_1_ChIPMunk_no_repeats_0.tab
# Environment: ARCHI_RELEASE_DIR (default ~/Release/TF-ML), SARUS_JAR (default ../../../sarus/releases/sarus-2.0.1.jar), ARCHI_SLIM_SCANS_DIR.
HERE=$(cd "$(dirname "$0")" && pwd)
TFML=${ARCHI_RELEASE_DIR:-~/Release/TF-ML}
JAR=${SARUS_JAR:-$HERE/../../../sarus/releases/sarus-2.0.1.jar}
SCANS=${ARCHI_SLIM_SCANS_DIR:-$HERE}
TF=$1; SET=${2:-MOUSE_10000_control}
mkdir -p $SCANS/munk
out=$SCANS/munk/${TF}_full_train_M_1_ChIPMunk_no_repeats_0.tab
if [ -s $out ]; then echo "skip $TF"; exit 0; fi
t0=$(date +%s)
nice -n 10 java -Xmx2G -XX:ParallelGCThreads=1 -cp $JAR ru.autosome.di.SARUS $TFML/all_mfa_file_${SET}_no_NF_no_N.fasta $TFML/outputdir/$TF/HUMAN_seq_HUMAN_pwm_mono/train_H_1_ChIPMunk_no_repeats_0.dpwm --skipn --show-non-matching --output-scoring-mode --transpose score besthit | grep -v \> > $out
rc=${PIPESTATUS[0]}
echo "$TF rc=$rc $(($(date +%s)-t0))s $(wc -l < $out) rows"
