#!/bin/env bash

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=01:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=Yu-Wen.Chen@colorado.edu
#SBATCH --output=sbatch-output_%x_%j.txt
#SBATCH --job-name=o2band_hapi
#SBATCH --account=blanca-airs
#SBATCH --qos=preemptable
#SBATCH --requeue
# HAPI cross-check (EVAL_PLAN.md #3): install hitran-api and compare our O2
# line-by-line optical depth against HAPI -- an independent implementation of the
# same HITRAN 2020 Voigt LBL.
#
# Needs internet for BOTH `pip install` and HAPI's line fetch. Blanca compute
# nodes usually have no outbound network, so run this on a LOGIN node:
#     bash curc_hapi_eval.sh
# (It also carries SBATCH headers in case a given partition's nodes do have
# network and you prefer to submit it.)

module load anaconda intel/2022.1.2 hdf5/1.10.1 zlib/1.2.11 netcdf/4.8.1 swig/4.1.1 gsl/2.7
conda activate er3t

PROJECT_ROOT="/projects/yuch8913/O2band_sim"
cd "$PROJECT_ROOT"
source setup_env.sh
export PYTHONPATH="$PROJECT_ROOT/src:$PYTHONPATH"

# Install HAPI into the active er3t env (idempotent; skips if already present).
python -c "import hapi" 2>/dev/null || pip install --quiet hitran-api

BANDS="${BANDS:-o2a o2b}"
ZT="${Z_TOP:-120}"
echo "[hapi] bands=[${BANDS}] z_top=${ZT}"
python src/eval_hapi.py --bands ${BANDS} --z-top "${ZT}"
