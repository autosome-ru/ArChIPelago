"""Locations of the data that are not part of the repository, read from environment variables in this one place.

A variable that is needed and not set stops the script with a message; there is no default path.

ARCHI_WORK      root of the working folders: human_sym_h5/ and mouse_sym_h5/ (data builds), H5_alltrain/ (model runs),
                mouse_swap/, human_sym_slim/, runtime/ (see the README of each stage)
ARCHI_RELEASE   the published pipeline output, folder TF-ML (outputdir/<TF>/..., train/, test/)
ARCHI_HOCOMOCO  HOCOMOCO v11 models folder (models/pwm/{mono,di}/all/<TF_SP>/)
ARCHI_SARUS_JAR SARUS 2.0.2 jar (make_scan_jobs_sym.py, runtime benchmark)
ARCHI_SARUS201_JAR SARUS 2.0.1 jar (diChIPMunk scans, 3_analysis/server/scan_tf_h5.sh)
ARCHI_CHIPMUNK  folder with run_dichiphorde8.rb, chipmunk.jar and ytilib/ (2_training/slim_chipmunk/run_munk_h5.sh)
ARCHI_MASK_HUMAN track_out_HUMAN.bed (repeat mask of the human build)
ARCHI_NB4       notebook4_results_table.csv (2_training/h5_server_v2.sh)
ARCHI_MACS      folder of the MACS peak files (<peak id>.interval, PEAKS041820.interval)
ARCHI_GENOMES   folder of the genome fasta files (GRCm38.p6.genome.fa, hg38_genome.fasta)
ARCHI_FAMILY    folder of core_annotation_HUMAN_mono.csv
ARCHI_SEQTK     seqtk executable
ARCHI_BIASAWAY  biasaway executable
"""
import os
import sys


def var(name, what):
    v = os.environ.get(name)
    if not v:
        sys.exit("CHECK FAILED: environment variable %s is not set (%s)" % (name, what))
    return os.path.expanduser(v)


def work(*parts):
    return os.path.join(var("ARCHI_WORK", "working root with human_sym_h5/, mouse_sym_h5/, H5_alltrain/"), *parts)


def release(*parts):
    return os.path.join(var("ARCHI_RELEASE", "published pipeline output folder TF-ML"), *parts)


def hocomoco(*parts):
    return os.path.join(var("ARCHI_HOCOMOCO", "HOCOMOCO v11 folder"), *parts)


def macs(*parts):
    return os.path.join(var("ARCHI_MACS", "folder of the MACS peak files"), *parts)


def genomes(*parts):
    return os.path.join(var("ARCHI_GENOMES", "folder of the genome fasta files"), *parts)


def family(*parts):
    return os.path.join(var("ARCHI_FAMILY", "folder with core_annotation_HUMAN_mono.csv"), *parts)


def root():
    """Repository root."""
    return os.path.dirname(os.path.realpath(__file__))


def loader_dir():
    """Folder with archi_data_v5.py and h5_splits.py (the data loader)."""
    return os.path.join(root(), "1_data_preparation")
