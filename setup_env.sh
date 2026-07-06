# Environment setup for the O2-band RT benchmark on CURC (Alpine).
#
#   source setup_env.sh
#
# Sets the er3t location, MCARaTS v0.10.4 executable, and the scratch output
# base.  er3t is editable-installed in the `er3t` conda env, so activating that
# env is usually enough; ER3T_HOME/PYTHONPATH here make a plain `python3` work too.

# --- er3t (directory CONTAINING the `er3t` package) ------------------------- #
export ER3T_HOME="/projects/yuch8913/wen_soft/er3t"
export PYTHONPATH="${ER3T_HOME}:${PYTHONPATH}"

# --- MCARaTS v0.10.4 executable (read by O2BandConfig) ---------------------- #
export MCARATS_V010_EXE="${MCARATS_V010_EXE:-/projects/yuch8913/wen_soft/mcarats/v0.10.4/src/mcarats}"

# --- data / output on CURC Alpine scratch ----------------------------------- #
# Prescribed inputs (HITRAN/AFGL/solar) and the TIPS-2021 QTpy tables live under
# the scratch data dir; results go to the scratch base.  The driver already
# defaults to these on Linux; exported here so they are explicit + overridable.
export O2BAND_DATA_DIR="${O2BAND_DATA_DIR:-/scratch/alpine/yuch8913/O2band_sim/data}"
export O2BAND_QTPY_DIR="${O2BAND_QTPY_DIR:-${O2BAND_DATA_DIR}/TIPS_2021_PYTHON/QTpy}"
export O2BAND_OUT_DIR="${O2BAND_OUT_DIR:-/scratch/alpine/yuch8913/O2band_sim}"

echo "[setup_env] ER3T_HOME=${ER3T_HOME}"
echo "[setup_env] MCARATS_V010_EXE=${MCARATS_V010_EXE}"
echo "[setup_env] O2BAND_DATA_DIR=${O2BAND_DATA_DIR}"
echo "[setup_env] O2BAND_QTPY_DIR=${O2BAND_QTPY_DIR}"
echo "[setup_env] O2BAND_OUT_DIR=${O2BAND_OUT_DIR}"
