#!/bin/env bash
#
# Submit the parallel O2-band pipeline as three dependent SLURM jobs:
#
#     prep (cache absorption)  ->  run (job array of shards)  ->  assemble (+noise)
#          1 small job              NTASKS array tasks             1 small job
#                          afterok                    afterok
#
# Run this on a LOGIN node (it only submits; it is not itself an sbatch job):
#
#     ./submit_o2band_array.sh [BANDS...]        # default: o2a o2b
#
# Tunable via environment (all optional; defaults = signed-off production):
#     NTASKS  (48)   number of array shards (work units are digest-balanced across them)
#     CORES   (8)    cores per array task (MCARaTS spreads a chunk's g-points over these)
#     CAP     (20)   max array tasks running at once (concurrency = CAP*CORES cores)
#     Z_TOP   (120)  km ;  PHOTONS (1e6) ;  NRUN (3)
#     TIME_RUN(08:00:00)   walltime per array task
#     OVERWRITE(0)   1 = recompute chunks even if valid
#     NOISE_THRESHOLD(0.01)  p95 relative-stderr gate reported at assemble
#
# Makespan ~ (total core-hours) / (CAP*CORES). One (band,SZA,albedo) ~= 63 min on
# 32 cores; full grid ~= 400 core-hours -> e.g. CAP=20,CORES=8 (160 cores) ~ 2.5 h.
#
# For ONE band, halve NTASKS (96 units -> keep 48 for 2 units/task, or 24 for 4).

set -euo pipefail
cd "$(dirname "$0")"

export BANDS="${*:-o2a o2b}"
export NTASKS="${NTASKS:-48}"
CORES="${CORES:-8}"
CAP="${CAP:-20}"
export Z_TOP="${Z_TOP:-120}"
export PHOTONS="${PHOTONS:-1e6}"
export NRUN="${NRUN:-3}"
export OVERWRITE="${OVERWRITE:-0}"
export NOISE_THRESHOLD="${NOISE_THRESHOLD:-0.01}"
TIME_RUN="${TIME_RUN:-08:00:00}"

STAGE=curc_stage_blanca_o2band.sh
LAST=$((NTASKS - 1))

echo "Pipeline: bands=[${BANDS}]  NTASKS=${NTASKS} CORES=${CORES} CAP=${CAP}" \
     "z_top=${Z_TOP} P=${PHOTONS} Nrun=${NRUN}"

# 1) prep -- build the absorption cache once (small, fast; no RT)
jid_prep=$(sbatch --parsable --job-name=o2b_prep \
    --ntasks=2 --time=00:20:00 --partition=blanca --export=ALL "$STAGE" prep)
echo "  prep     job ${jid_prep}"

# 2) run -- job array of shards, starts only if prep succeeded
jid_run=$(sbatch --parsable --dependency=afterok:${jid_prep} --job-name=o2b_run \
    --ntasks=${CORES} --array=0-${LAST}%${CAP} --time=${TIME_RUN} \
    --partition=blanca --export=ALL "$STAGE" run)
echo "  run      job ${jid_run}  (array 0-${LAST}%${CAP}, ${CORES} cores/task)"

# 3) assemble -- stitch + noise report, starts only if ALL array tasks succeeded
jid_asm=$(sbatch --parsable --dependency=afterok:${jid_run} --job-name=o2b_asm \
    --ntasks=2 --time=00:30:00 --partition=blanca --export=ALL "$STAGE" assemble)
echo "  assemble job ${jid_asm}"

echo "Submitted. Monitor:  squeue -u \$USER   |   outputs -> ${O2BAND_OUT_DIR:-/scratch/alpine/yuch8913/O2band_sim}/z${Z_TOP%.*}_p${PHOTONS}_n${NRUN}/"
