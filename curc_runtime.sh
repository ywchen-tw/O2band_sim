# Runtime setup shared by the CURC batch scripts.
#
# This file is meant to be sourced, not submitted directly.  Keep the compiler
# first: CURC's hierarchical module system may unload modules tied to a different
# compiler family when the compiler changes.

module purge
module load intel/2022.1.2
module load zlib/1.2.11 hdf5/1.10.1 netcdf/4.8.1 swig/4.1.1 gsl/2.7

# Use the environment's interpreter directly.  This does not depend on the
# Anaconda module leaving a `conda` shell function on PATH, and it works in
# non-interactive Slurm shells without `conda init`.
export O2BAND_CONDA_ENV="${O2BAND_CONDA_ENV:-/projects/yuch8913/software/anaconda/envs/er3t}"
export O2BAND_PYTHON="${O2BAND_PYTHON:-${O2BAND_CONDA_ENV}/bin/python}"

if [[ ! -x "${O2BAND_PYTHON}" ]]; then
    echo "[curc_runtime] ERROR: Python is not executable: ${O2BAND_PYTHON}" >&2
    echo "[curc_runtime] Set O2BAND_CONDA_ENV or O2BAND_PYTHON to the er3t environment." >&2
    return 1
fi

export PATH="${O2BAND_CONDA_ENV}/bin:${PATH}"
export CONDA_PREFIX="${O2BAND_CONDA_ENV}"
hash -r

echo "[curc_runtime] python=${O2BAND_PYTHON}"
