#!/bin/bash
# Needs Slim/ in the repository root (unpacked from Slim.tar.gz of Zenodo v1.0.0, not tracked).
# Slim m=0/1/5 scans of one TF on the human training, human test and mouse test fasta (java, Slim jar),
# then the diChIPMunk (SARUS 2.0.1) scans. Usage: scan_tf_h5.sh <TF>. ARCHI_SARUS201_JAR is the SARUS 2.0.1 jar; the java log is written to log_v2_<epoch>.txt (never overwritten).
. "$(dirname "$(readlink -f "$0")")/../../archi_env.sh"
: "${ARCHI_SARUS201_JAR:?CHECK FAILED: ARCHI_SARUS201_JAR (sarus-2.0.1.jar) is not set}"
TF=$1; H2=$ARCHI_WORK/H5_alltrain/human_sym/slim; MD=$H2/models_server/$TF; SD=$H2/scans/$TF
JAVA=$ARCHI_CODE/Slim/jdk8u232-b09/bin/java; JAR=$ARCHI_CODE/Slim/TrainAndApplySlim.jar
declare -A FA=( [HUMAN_10000_train]=$ARCHI_WORK/human_sym_h5/TF-ML/all_mfa_file_HUMAN_10000_train_no_NF_no_N.fasta
                [HUMAN_10000_control]=$ARCHI_WORK/human_sym_h5/TF-ML/all_mfa_file_HUMAN_10000_control_no_NF_no_N.fasta
                [MOUSE_10000_train]=$ARCHI_WORK/mouse_sym_h5/TF-ML/all_mfa_file_MOUSE_10000_train_no_NF_no_N.fasta )
for S in "${!FA[@]}"; do [ -s "${FA[$S]}" ] || { echo "ERROR fasta missing ${FA[$S]}"; exit 1; }; done
mkdir -p $SD
for S in HUMAN_10000_train HUMAN_10000_control MOUSE_10000_train; do
  for o in 0 1 5; do
    out=$SD/m${o}_$S; XML=$MD/${TF}_SlimModel_$o/Motif_1/SlimDimont_1.xml
    [ -s $out/model1_predictions.txt ] && continue
    [ -e $out/predictions.txt ] && { echo "ERROR $out/predictions.txt without model1_predictions.txt"; exit 1; }
    [ -s $XML ] || { echo "ERROR no model $XML"; exit 1; }
    mkdir -p $out
    nice -n 10 $JAVA -Xmx4G -XX:ParallelGCThreads=1 -jar $JAR scan i=${FA[$S]} m=$XML outdir=$out > $out/log_v2_$(date +%s).txt 2>&1 || { echo "ERROR scan $TF m$o $S"; exit 1; }
    ln $out/predictions.txt $out/model1_predictions.txt
  done
  out=$SD/munk_$S.tab
  if [ ! -s $out ]; then
    [ -e $out.part ] && { echo "ERROR $out.part exists from an earlier run"; exit 1; }
    nice -n 10 java -Xmx2G -XX:ParallelGCThreads=1 -cp $ARCHI_SARUS201_JAR ru.autosome.di.SARUS ${FA[$S]} $MD/HUMAN_seq_HUMAN_pwm_mono/train_H_1_ChIPMunk_no_repeats_0.dpwm --skipn --show-non-matching --output-scoring-mode --transpose score besthit | grep -v \> > $out.part
    [ ${PIPESTATUS[0]} -eq 0 ] && ln $out.part $out || { echo "ERROR sarus $TF $S"; exit 1; }
  fi
done
echo "scans $TF done $(date +%T)"
