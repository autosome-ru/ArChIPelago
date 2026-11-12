# Source this file from a driver: it checks the environment variables and sets PY.
# ARCHI_WORK    working root (human_sym_h5/, mouse_sym_h5/, H5_alltrain/, ...)
# ARCHI_RELEASE published pipeline output folder TF-ML (only where a driver reads it)
# ARCHI_PY      Python interpreter with the packages of 2_training/VERSIONS.txt (default python3)
: "${ARCHI_WORK:?CHECK FAILED: environment variable ARCHI_WORK is not set}"
PY=${ARCHI_PY:-python3}
ARCHI_CODE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
