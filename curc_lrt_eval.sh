#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=00:30:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=Yu-Wen.Chen@colorado.edu
#SBATCH --output=sbatch-output_%x_%j.txt
#SBATCH --job-name=o2band_lrt
#SBATCH --account=blanca-airs
#SBATCH --qos=preemptable
#SBATCH --requeue
# Independent RT-solver cross-check (EVAL_PLAN.md E2): our MCARaTS reflectance vs
# libRadtran/DISORT at the band window (near-Rayleigh + Lambertian). uvspec needs
# GSL/NetCDF runtime libs, so this uses the SAME module set as the arcsix
# libRadtran workflow. Run on a node where uvspec works:
#     bash curc_lrt_eval.sh   (or sbatch)

set -euo pipefail

PROJECT_ROOT="/projects/yuch8913/O2band_sim"
cd "$PROJECT_ROOT"
source curc_runtime.sh
source setup_env.sh
export PYTHONPATH="$PROJECT_ROOT/src:${PYTHONPATH:-}"
# libRadtran install used by the arcsix workflow (already the env default here)
export LIBRADTRAN_V2_DIR="${LIBRADTRAN_V2_DIR:-/projects/yuch8913/wen_soft/libRadtran-2.0.6}"

# Band to check (o2a|o2b), and its produced per-band file. OUR_H5 defaults to the
# z120 production file for the chosen BAND, so `BAND=o2b` uses o2b.h5.
BAND="${BAND:-o2a}"
# Run directory is built from the SAME stamp curc_stage_blanca_o2band.sh uses, so
# this tracks the current production run instead of naming one that goes stale
# (it previously hardcoded z120_p1e6_n3, two runs behind).  Override OUR_H5 to
# evaluate a different run.
ZT="${Z_TOP:-120}"; PH="${PHOTONS:-1e7}"; NR="${NRUN:-3}"
GRID="${GRID:-vac}"; CUTOFF="${CUTOFF_CM:-50}"; CIA="${CIA:-none}"
STAMP="${GRID}_c$(printf '%.0f' "$CUTOFF")"
[ "$CIA" = "none" ] || STAMP="${STAMP}_cia"
RUN_DIR="${O2BAND_OUT_DIR:-/scratch/alpine/yuch8913/O2band_sim}/z$(printf '%.0f' "$ZT")_p${PH}_n${NR}_${STAMP}"
OUR_H5="${OUR_H5:-${RUN_DIR}/${BAND}.h5}"
# INBAND=1 runs the in-band check (injects our per-layer gas OT into DISORT,
# spanning column gas OD ~0.05-3); default runs the window (pure-Rayleigh) check.
if [ "${INBAND:-0}" = "1" ]; then
    echo "[lrt] INBAND  uvspec=$LIBRADTRAN_V2_DIR/bin/uvspec  our=$OUR_H5  band=$BAND"
    "$O2BAND_PYTHON" src/eval_lrt_inband.py "$OUR_H5" --band "$BAND" --streams "${STREAMS:-16}"
else
    echo "[lrt] window  uvspec=$LIBRADTRAN_V2_DIR/bin/uvspec  our=$OUR_H5  band=$BAND  n_wvl=${N_WVL:-6}"
    "$O2BAND_PYTHON" src/eval_lrt.py "$OUR_H5" --band "$BAND" --streams "${STREAMS:-16}" --n-wvl "${N_WVL:-6}"
fi
