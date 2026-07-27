#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=08:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=Yu-Wen.Chen@colorado.edu
#SBATCH --output=sbatch-output_%x_%j.txt
#SBATCH --account=blanca-airs
#SBATCH --qos=preemptable
#SBATCH --requeue
# Parametrized stage runner for the parallel O2-band pipeline. Do NOT sbatch this
# directly -- use submit_o2band_array.sh, which sets --ntasks / --array /
# --dependency / --job-name per stage and inherits the config via --export.
# First arg selects the stage: prep | run | assemble.
#
#   prep     : build + cache the absorption object(s) once (no RT), so the array
#              tasks load the cache instead of racing to write it.
#   run      : execute this task's shard of work units (SLURM array job). The
#              shard is a stride over the digest-ordered unit list -> size-balanced
#              and geometry-interleaved. Skip-if-done + --requeue = resumable.
#   assemble : stitch all chunk files into per-band + merged HDF5, then report MC
#              noise per band (non-fatal gate).

set -euo pipefail

PROJECT_ROOT="/projects/yuch8913/O2band_sim"
cd "$PROJECT_ROOT"
source curc_runtime.sh
source setup_env.sh
export PYTHONPATH="$PROJECT_ROOT/src:${PYTHONPATH:-}"

STAGE="${1:?usage: curc_stage_blanca_o2band.sh prep|run|assemble}"

# Config (inherited from the wrapper via --export=ALL; defaults match production).
BANDS="${BANDS:-o2a o2b}"
NTASKS="${NTASKS:-48}"
ZT="${Z_TOP:-120}"; PH="${PHOTONS:-1e6}"; NR="${NRUN:-3}"
GRID="${GRID:-vac}"; CUTOFF="${CUTOFF_CM:-50}"; CIA="${CIA:-none}"
# 1001 divides 15001 into 15 chunks (last 987); 1000 leaves a 1-point chunk
# whose 6 geometry work units each spin up MCARaTS for a single g-point.
CHUNK="${CHUNK:-1001}"
BASE_OUT="${O2BAND_OUT_DIR:-/scratch/alpine/yuch8913/O2band_sim}"
# Config-stamped subdir, so chunk files compose across nodes AND runs with
# different physics can never land in the same directory.  Chunk files are keyed
# by (band, sza, albedo, index) only and `run` is skip-if-done, so mixing a new
# wing cutoff / grid convention / CIA setting into an existing directory would
# silently return the old numbers.  (sim_o2band.py also enforces this via
# _physics_config.json; the name is belt-and-braces so it is visible in `ls`.)
STAMP="${GRID}_c$(printf '%.0f' "$CUTOFF")"
[ "$CIA" = "none" ] || STAMP="${STAMP}_cia"
OUT="${BASE_OUT}/z$(printf '%.0f' "$ZT")_p${PH}_n${NR}_${STAMP}"
OVERWRITE_FLAG=""; [ "${OVERWRITE:-0}" = "1" ] && OVERWRITE_FLAG="--overwrite"

COMMON="--bands ${BANDS} --z-top ${ZT} --photons ${PH} --nrun ${NR} --out-dir ${OUT}"
COMMON="${COMMON} --grid ${GRID} --cutoff-cm ${CUTOFF} --cia ${CIA}"
COMMON="${COMMON} --chunk-size ${CHUNK}"

case "$STAGE" in
  prep)
    echo "[prep] bands=[${BANDS}] out=${OUT}"
    "$O2BAND_PYTHON" src/sim_o2band.py --stage prep $COMMON
    ;;
  run)
    T="${SLURM_ARRAY_TASK_ID:?run stage must be launched as a SLURM array job}"
    NCPU="${SLURM_NTASKS:-8}"
    echo "[run] shard ${T}/${NTASKS}  ncpu=${NCPU}  bands=[${BANDS}]  out=${OUT}"
    "$O2BAND_PYTHON" src/sim_o2band.py --stage run --shard "$T" "$NTASKS" \
        --ncpu "$NCPU" $COMMON $OVERWRITE_FLAG
    ;;
  assemble)
    echo "[assemble] bands=[${BANDS}] out=${OUT}"
    "$O2BAND_PYTHON" src/sim_o2band.py --stage assemble $COMMON
    # Non-fatal noise gate: report worst-case MC noise per band.
    for b in $BANDS; do
        "$O2BAND_PYTHON" src/noise_report.py "${OUT}/${b}.h5" --threshold "${NOISE_THRESHOLD:-0.01}" \
            || echo "[assemble] WARNING: ${b} exceeds noise threshold ${NOISE_THRESHOLD:-0.01}"
    done
    ;;
  *)
    echo "unknown STAGE: ${STAGE}" >&2; exit 2 ;;
esac
