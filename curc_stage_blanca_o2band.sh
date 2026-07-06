#!/bin/env bash

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

module load anaconda intel/2022.1.2 hdf5/1.10.1 zlib/1.2.11 netcdf/4.8.1 swig/4.1.1 gsl/2.7
conda activate er3t

PROJECT_ROOT="/projects/yuch8913/O2band_sim"
cd "$PROJECT_ROOT"
source setup_env.sh
export PYTHONPATH="$PROJECT_ROOT/src:$PYTHONPATH"

STAGE="${1:?usage: curc_stage_blanca_o2band.sh prep|run|assemble}"

# Config (inherited from the wrapper via --export=ALL; defaults match production).
BANDS="${BANDS:-o2a o2b}"
NTASKS="${NTASKS:-48}"
ZT="${Z_TOP:-120}"; PH="${PHOTONS:-1e6}"; NR="${NRUN:-3}"
BASE_OUT="${O2BAND_OUT_DIR:-/scratch/alpine/yuch8913/O2band_sim}"
# Same config-stamped subdir as the single-node runner, so chunk files compose.
OUT="${BASE_OUT}/z$(printf '%.0f' "$ZT")_p${PH}_n${NR}"
OVERWRITE_FLAG=""; [ "${OVERWRITE:-0}" = "1" ] && OVERWRITE_FLAG="--overwrite"

COMMON="--bands ${BANDS} --z-top ${ZT} --photons ${PH} --nrun ${NR} --out-dir ${OUT}"

case "$STAGE" in
  prep)
    echo "[prep] bands=[${BANDS}] out=${OUT}"
    python src/sim_o2band.py --stage prep $COMMON
    ;;
  run)
    T="${SLURM_ARRAY_TASK_ID:?run stage must be launched as a SLURM array job}"
    NCPU="${SLURM_NTASKS:-8}"
    echo "[run] shard ${T}/${NTASKS}  ncpu=${NCPU}  bands=[${BANDS}]  out=${OUT}"
    python src/sim_o2band.py --stage run --shard "$T" "$NTASKS" \
        --ncpu "$NCPU" $COMMON $OVERWRITE_FLAG
    ;;
  assemble)
    echo "[assemble] bands=[${BANDS}] out=${OUT}"
    python src/sim_o2band.py --stage assemble $COMMON
    # Non-fatal noise gate: report worst-case MC noise per band.
    for b in $BANDS; do
        python src/noise_report.py "${OUT}/${b}.h5" --threshold "${NOISE_THRESHOLD:-0.01}" \
            || echo "[assemble] WARNING: ${b} exceeds noise threshold ${NOISE_THRESHOLD:-0.01}"
    done
    ;;
  *)
    echo "unknown STAGE: ${STAGE}" >&2; exit 2 ;;
esac
