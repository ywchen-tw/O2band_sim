#!/bin/bash

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

set -euo pipefail

PROJECT_ROOT="/projects/yuch8913/O2band_sim"
cd "$PROJECT_ROOT"
source curc_runtime.sh
source setup_env.sh
export PYTHONPATH="$PROJECT_ROOT/src:${PYTHONPATH:-}"

# Install HAPI into the active er3t env (idempotent; skips if already present).
"$O2BAND_PYTHON" -c "import hapi" 2>/dev/null \
    || "$O2BAND_PYTHON" -m pip install --quiet hitran-api

BANDS="${BANDS:-o2a o2b}"
ZT="${Z_TOP:-120}"
echo "[hapi] bands=[${BANDS}] z_top=${ZT}"
"$O2BAND_PYTHON" src/eval_hapi.py --bands ${BANDS} --z-top "${ZT}"
