# Running on CURC (Blanca / Alpine)

Deployment notes for CU Research Computing. The physics/usage docs are in
[`README.md`](README.md); this file covers only the CURC-specific environment and
the batch run scripts.

## Scratch layout

Inputs, the TIPS-2021 QTpy tables, and all outputs live under scratch; the driver
defaults to these paths on Linux (override with the `O2BAND_*` env vars):

```
/scratch/alpine/yuch8913/O2band_sim/
├── data/                          # prescribed inputs (staged here, not in the repo)
│   ├── hitran2020_lines.txt
│   ├── afglms.dat
│   ├── CU_composite_solar.dat
│   └── TIPS_2021_PYTHON/QTpy/      # TIPS-2021 partition sums (182 .QTpy files)
└── z<ztop>_p<photons>_n<nrun>/    # config-stamped output tree (per physics config)
    ├── _absorption_cache/          # cached per-band absorption pickles
    ├── <band>/sza<SZA>_alb<ALB>/chunk_<i0>_<i1>.h5     # chunk files (resume unit)
    ├── o2a.h5  o2b.h5              # per-band assembled outputs
    └── o2band_benchmark.h5         # merged (both bands as groups)
```

The output subdir is **stamped with the physics config** (`z120_p1e6_n3`, …) so
runs that differ in `z_top` / `photons` / `Nrun` never share chunk files. The
skip-if-done resume key is `(band, SZA, albedo, grid-index)` and does *not* encode
those, so the stamp is what keeps configs from silently colliding.

## Software environment

`er3t` is editable-installed (`pip install -e`) in the `er3t` conda env, so
activating it usually needs no `sys.path` hacks. The MCARaTS v0.10.4 executable is
`/projects/yuch8913/wen_soft/mcarats/v0.10.4/src/mcarats`.

[`setup_env.sh`](setup_env.sh) sets everything in one place — the er3t location,
the MCARaTS executable, and the scratch data/QTpy/output bases:

```bash
module load anaconda intel/2022.1.2 hdf5/1.10.1 zlib/1.2.11 netcdf/4.8.1 swig/4.1.1 gsl/2.7
conda activate er3t
source setup_env.sh          # ER3T_HOME, MCARATS_V010_EXE, O2BAND_{DATA,QTPY,OUT}_DIR
export PYTHONPATH="$PWD/src:$PYTHONPATH"

python tests/test_absorption.py     # 10/10 physics validation checks (no MCARaTS needed)
```

Use the `er3t` conda env's interpreter
(`/projects/yuch8913/software/anaconda/envs/er3t/bin/python`) — the system
`python3` lacks `h5py`/er3t.

## Single-node runner

[`curc_shell_blanca_o2band.sh`](curc_shell_blanca_o2band.sh) — one SBATCH job on
one node (account `blanca-airs`, `preemptable` QOS, `--requeue`). Loads modules,
activates `er3t`, sources `setup_env.sh`, and calls the driver CLI.

```bash
sbatch curc_shell_blanca_o2band.sh test        # fast RT sanity run (closes PLAN.md V7)
sbatch curc_shell_blanca_o2band.sh full        # frozen Phase-1 grid, both bands
sbatch curc_shell_blanca_o2band.sh full o2a    # restrict a full run to one band
```

Env overrides (defaults = signed-off production for `full`: z_top=120 km, P=1e6,
Nrun=3):

```bash
Z_TOP=70 PHOTONS=1e4 sbatch curc_shell_blanca_o2band.sh full o2a   # 70 km / low-photon test
SZAS=30 ALBEDOS=0    sbatch curc_shell_blanca_o2band.sh full o2a   # one geometry (shakeout)
OVERWRITE=1          sbatch curc_shell_blanca_o2band.sh full       # force recompute
NCPU=16              sbatch curc_shell_blanca_o2band.sh full       # override core count
```

Because runs are **skip-if-done**, a `preemptable` job that gets requeued resumes
from the last completed chunk instead of restarting. MCARaTS parallelizes a
chunk's g-points across `Ncpu` (defaults to the SLURM allocation); memory is not
the binding constraint (each g-point is a lightweight monochromatic solve).

Reference timing: one `(band, SZA, albedo)` full band ≈ **63 min on 32 cores**
(z_top=120, P=1e6, Nrun=3), so the full 12-combo grid is ≈ 400 core-hours.

## Parallel job-array pipeline

For the full grid, split the work across many SLURM array tasks instead of one
serial node. The work units `(band, SZA, albedo, chunk)` are independent and
idempotent (PLAN.md §2h), so they parallelize cleanly and a preempted task only
loses its in-flight chunk.

[`submit_o2band_array.sh`](submit_o2band_array.sh) submits **three dependent
jobs** — run it on a **login node** (it only submits; it is not an sbatch job):

```
prep (cache absorption)  ─afterok→  run (job array of shards)  ─afterok→  assemble (+noise)
```

```bash
./submit_o2band_array.sh              # both bands, medium: 48 tasks × 8 cores, CAP 20
CAP=40 ./submit_o2band_array.sh       # more concurrency (~1.3 h vs ~2.5 h)
./submit_o2band_array.sh o2a          # one band only
```

Tunable via environment (defaults in parentheses):

| var | meaning | default |
|---|---|---|
| `NTASKS` | number of array shards (work units are digest-balanced across them) | 48 |
| `CORES` | cores per array task | 8 |
| `CAP` | max array tasks running at once (concurrency = `CAP*CORES` cores) | 20 |
| `Z_TOP` / `PHOTONS` / `NRUN` | physics config (sets the output stamp) | 120 / 1e6 / 3 |
| `TIME_RUN` | walltime per array task | 08:00:00 |
| `OVERWRITE` | 1 = recompute chunks even if valid | 0 |
| `NOISE_THRESHOLD` | p95 relative-stderr gate reported at assemble | 0.01 |

Makespan ≈ (total core-hours) / (`CAP*CORES`): CAP=20 → ~2.5 h, CAP=40 → ~1.3 h,
subject to Blanca backfill. Each shard is a stride over a stable digest-ordered
unit list, so shards are size-balanced *and* geometry-interleaved.

Under the hood the array uses the parametrized stage runner
[`curc_stage_blanca_o2band.sh`](curc_stage_blanca_o2band.sh) (`prep | run |
assemble`); the wrapper sets `--ntasks` / `--array` / `--dependency` / `--job-name`
per stage. `assemble` runs only if **all** array tasks succeed (`afterok`), so
incomplete data is never assembled.

Monitor: `squeue -u $USER`. Outputs land in
`$O2BAND_OUT_DIR/z<ztop>_p<photons>_n<nrun>/`.

## Noise check

After assembly (the array does this automatically), gate the MC noise:

```bash
python src/noise_report.py /scratch/alpine/yuch8913/O2band_sim/z120_p1e6_n3/o2a.h5 --threshold 0.01
```

Reports per-`(SZA, albedo)` relative + absolute reflectance standard error and
exits non-zero if any p95 relative stderr exceeds the threshold. Reference: at
P=1e6/Nrun=3 the o2a shakeout gave p95 ≈ 0.005 (well under 1%).
