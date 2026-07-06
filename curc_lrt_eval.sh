#!/bin/env bash

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

module load anaconda intel/2022.1.2 hdf5/1.10.1 zlib/1.2.11 netcdf/4.8.1 swig/4.1.1 gsl/2.7
conda activate er3t

PROJECT_ROOT="/projects/yuch8913/O2band_sim"
cd "$PROJECT_ROOT"
source setup_env.sh
export PYTHONPATH="$PROJECT_ROOT/src:$PYTHONPATH"
# libRadtran install used by the arcsix workflow (already the env default here)
export LIBRADTRAN_V2_DIR="${LIBRADTRAN_V2_DIR:-/projects/yuch8913/wen_soft/libRadtran-2.0.6}"

# Band to check (o2a|o2b), and its produced per-band file. OUR_H5 defaults to the
# z120 production file for the chosen BAND, so `BAND=o2b` uses o2b.h5.
BAND="${BAND:-o2a}"
OUR_H5="${OUR_H5:-${O2BAND_OUT_DIR:-/scratch/alpine/yuch8913/O2band_sim}/z120_p1e6_n3/${BAND}.h5}"
echo "[lrt] uvspec=$LIBRADTRAN_V2_DIR/bin/uvspec  our=$OUR_H5  band=$BAND  n_wvl=${N_WVL:-6}"
python src/eval_lrt.py "$OUR_H5" --band "$BAND" --streams "${STREAMS:-16}" --n-wvl "${N_WVL:-6}"
