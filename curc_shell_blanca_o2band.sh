#!/bin/env bash

#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Yu-Wen.Chen@colorado.edu
#SBATCH --output=sbatch-output_%x_%j.txt
#SBATCH --job-name=o2band_sim
#SBATCH --account=blanca-airs
#### #SBATCH --partition=blanca-airs
#SBATCH --qos=preemptable
# preemptable jobs can be killed mid-run; --requeue + the driver's skip-if-done
# chunk logic let a requeued job pick up where it left off instead of restarting.
#SBATCH --requeue



module load anaconda intel/2022.1.2 hdf5/1.10.1 zlib/1.2.11 netcdf/4.8.1 swig/4.1.1 gsl/2.7
conda activate er3t

PROJECT_ROOT="/projects/yuch8913/O2band_sim"
cd "$PROJECT_ROOT"

# Sets ER3T_HOME, PYTHONPATH (er3t), MCARATS_V010_EXE, and O2BAND_OUT_DIR
# (scratch). er3t is editable-installed in the `er3t` conda env, so this mainly
# pins MCARaTS + the scratch output base and is harmless if already set.
source setup_env.sh
export PYTHONPATH="$PROJECT_ROOT/src:$PYTHONPATH"

# Usage: sbatch curc_shell_blanca_o2band.sh [MODE] [BANDS...]
#   MODE  : test | full  (default test)
#           test -> narrow 763.20-763.30 nm window, few photons (fast sanity run)
#           full -> frozen Phase-1 grid (SZA 0/30/60, albedo 0/0.1); production
#                   defaults z_top=120 km (full AFGL), P=1e6, Nrun=3
#   BANDS : optional band list to restrict a full run, e.g. 'o2a' (default: both)
# Env overrides: Z_TOP, PHOTONS, NRUN, NCPU, OVERWRITE=1
# Examples:
#   sbatch curc_shell_blanca_o2band.sh test
#   sbatch curc_shell_blanca_o2band.sh full
#   sbatch curc_shell_blanca_o2band.sh full o2a
#   Z_TOP=70 sbatch curc_shell_blanca_o2band.sh full o2a   # 70 km scaling test
MODE="${1:-test}"
shift 2>/dev/null
BANDS=("$@")   # remaining args = explicit band list (may be empty)

# MCARaTS parallelises the g-points of one chunk across cores, so give the run
# the cores this allocation holds (override with NCPU). Each g-point is a
# lightweight monochromatic solve, so -- unlike the uvspec/CRE runner -- memory
# is not the binding constraint here; core count is.
NCPU="${NCPU:-$SLURM_NTASKS}"
[ -z "$NCPU" ] && NCPU="auto"
echo "Alloc ${SLURM_NTASKS:-?} cores -> Ncpu=${NCPU}, mode=${MODE}, bands=[${BANDS[*]:-default}]"

# Resumable by default: valid existing chunk HDF5s are skipped and reused, so a
# preempted/requeued job continues instead of restarting. Set OVERWRITE=1 to
# force a full recompute.
OVERWRITE_FLAG=""
[ "${OVERWRITE:-0}" = "1" ] && OVERWRITE_FLAG="--overwrite"

# Optional band restriction (applies to full mode; --test fixes its own band).
BANDS_FLAG=""
[ "${#BANDS[@]}" -gt 0 ] && BANDS_FLAG="--bands ${BANDS[*]}"

# Optional geometry restriction (full mode) -- e.g. shake out ONE SZA/albedo on
# the full band before committing the whole grid:
#   SZAS=30 ALBEDOS=0 sbatch curc_shell_blanca_o2band.sh full o2a
# Space-separated lists are allowed: SZAS="0 30 60", ALBEDOS="0 0.1".
SZAS_FLAG="";    [ -n "${SZAS:-}" ]    && SZAS_FLAG="--szas ${SZAS}"
ALBEDOS_FLAG=""; [ -n "${ALBEDOS:-}" ] && ALBEDOS_FLAG="--albedos ${ALBEDOS}"

# Config-stamped output subdir so runs that differ in physics (z_top / photons /
# Nrun) NEVER share chunk files. This matters because the skip-if-done resume key
# (band, SZA, albedo, grid-index) does NOT encode z_top/photons/Nrun -- without
# the stamp, re-running at a different z_top would silently reuse stale chunks and
# pair them with freshly-computed optical depths in the same output file.
BASE_OUT="${O2BAND_OUT_DIR:-/scratch/alpine/yuch8913/O2band_sim}"

if [ "$MODE" = "test" ]; then
    # Fast end-to-end validation of the RT wiring (closes PLAN.md V7 / step 2f).
    # These match the driver's --test preset; kept in sync so the stamp is honest.
    ZT="${Z_TOP:-70}"; PH="${PHOTONS:-1e4}"; NR="${NRUN:-2}"
    OUT="${BASE_OUT}/test_z$(printf '%.0f' "$ZT")_p${PH}_n${NR}"
    echo "mode=test  z_top=${ZT} photons=${PH} Nrun=${NR}  out=${OUT}"
    python src/sim_o2band.py --test \
        --z-top "$ZT" --photons "$PH" --nrun "$NR" \
        --ncpu "$NCPU" --out-dir "$OUT" $OVERWRITE_FLAG
else
    # Frozen Phase-1 production grid. Defaults are the signed-off production run
    # (z_top=120 km full AFGL, P=1e6, Nrun=3); override via Z_TOP/PHOTONS/NRUN.
    ZT="${Z_TOP:-120}"; PH="${PHOTONS:-1e6}"; NR="${NRUN:-3}"
    OUT="${BASE_OUT}/z$(printf '%.0f' "$ZT")_p${PH}_n${NR}"
    echo "mode=full  z_top=${ZT} photons=${PH} Nrun=${NR} szas=[${SZAS:-all}] albedos=[${ALBEDOS:-all}]  out=${OUT}"
    python src/sim_o2band.py \
        --z-top "$ZT" --photons "$PH" --nrun "$NR" \
        --ncpu "$NCPU" --out-dir "$OUT" \
        $BANDS_FLAG $SZAS_FLAG $ALBEDOS_FLAG $OVERWRITE_FLAG
fi
