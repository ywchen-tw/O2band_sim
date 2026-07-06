# O2 A/B-band RT Intercomparison — Phase 1

A high-spectral-resolution line-by-line benchmark of **top-of-atmosphere (TOA)
reflectance** in the molecular-oxygen **A-band (~760 nm)** and **B-band
(~688 nm)**, computed under tightly prescribed clear-sky conditions so that
different radiative-transfer (RT) models can be compared line-for-line.

Correct simulation of these bands underpins satellite retrievals of surface
pressure, cloud/aerosol height, and CO₂ light-path correction.

> **Scope.** This repository implements **Phase 1 only**: clear-sky, everything
> prescribed. No clouds, aerosols, polarization, instrument convolution, or
> real-data comparison. See [`PLAN.md`](PLAN.md) for the scientific/architectural
> blueprint and the validation table.

---

## What it computes

For each **band × solar-zenith-angle × surface-albedo** combination, on a uniform
**0.001 nm air-wavelength grid**:

- **TOA reflectance** ρ(λ) and **absolute radiance** I(λ), each with a Monte-Carlo
  standard error.
- Separable, independently inspectable **optical thickness**: O₂ absorption, H₂O
  absorption, and Rayleigh scattering — layer-resolved and column-integrated, on
  the same spectral/vertical grid.

### Prescribed Phase-1 settings (frozen)

| | |
|---|---|
| Bands | B-band **680–695 nm**, A-band **757–772 nm** |
| Spectral resolution | **0.001 nm** (air wavelengths) |
| Line data | **HITRAN 2020** |
| Line shape | **Voigt** (Doppler ⊛ Lorentz, via Faddeeva `wofz`) |
| Partition sums | **TIPS-2021** (consistent with HITRAN 2020) |
| Atmosphere | AFGL **mid-latitude summer** |
| Solar zenith angle | **0°, 30°, 60°** |
| Viewing | **nadir**, relative azimuth 0° |
| Surface | **Lambertian**, albedo **0.0** and **0.1** |
| Rayleigh | Bodhaine et al. (1999) |
| RT engine | **MCARaTS v0.10.4** (1D, IPA solver) via **er3t** |

---

## Physics conventions (stated explicitly)

- **Line intensity T-scaling** uses HITRAN lower-state energy E″, the
  stimulated-emission factor, and the partition ratio Q(296)/Q(T) from TIPS-2021.
  Reference T = 296 K, P = 1013.25 hPa, c₂ = 1.4387769 cm·K.
- **Lorentz HWHM** γ_L = (296/T)^n_air·[γ_air(P−P_self)+γ_self·P_self].
  **Doppler HWHM** α_D = ν₀·√(2 ln2 k_B T / m)/c, per-isotopologue mass.
- **Wavelengths are air**; HITRAN wavenumbers are vacuum. Conversion via the
  Edlén (1966) dispersion formula (λ_air = λ_vac / n(λ)).
- **Line-wing cutoff** ±N·(ν₀/R) per line, default N = 3, R = 20000.
- **H₂O** absorption is optional (`include_h2o`); water vapor is always part of the
  foreign (air) pressure that broadens O₂ (no separate γ_H₂O).
- **O₂ CIA / continuum** excluded in Phase 1 (flagged in metadata).
- **Solar source**: the CU composite solar spectrum is the incident irradiance
  F₀(λ). MCARaTS runs a unit source (`Src_flx = 1`); F₀ is folded in afterward
  (MC transport is linear in incident flux). Reflectance
  **ρ(λ) = π·I/(μ₀·F₀) = π·R_raw/μ₀** is F₀-independent (μ₀ = cos SZA).

Every setting actually used is written into the output HDF5 metadata, and any
deviation from the prescribed settings is flagged there.

---

## Repository layout

```
o2band_sim/
├── PLAN.md                     # scientific + architectural blueprint (read this)
├── README.md
├── CURC_NOTES.md               # CURC/Blanca deployment: env, scratch paths, batch runs
├── setup_env.sh                # env helper (ER3T_HOME, MCARaTS, data/out paths)
├── curc_shell_blanca_o2band.sh # single-node SBATCH runner        (CURC — see CURC_NOTES.md)
├── curc_stage_blanca_o2band.sh # parametrized prep|run|assemble stage runner   (CURC)
├── submit_o2band_array.sh      # submits the parallel prep->array->assemble pipeline (CURC)
├── data/                       # prescribed inputs (not committed; provided separately)
│   ├── hitran2020_lines.txt    # HITRAN 2020 O2 + H2O lines (160-col .par format)
│   ├── afglms.dat              # AFGL mid-latitude-summer profile
│   ├── CU_composite_solar.dat  # CU composite solar reference spectrum
│   └── TIPS_2021_PYTHON/QTpy/  # TIPS-2021 partition sums (provided separately)
├── src/
│   ├── sim_o2band.py           # driver: config, cached absorption, chunked/sharded RT, assemble
│   ├── noise_report.py         # per-(SZA,albedo) MC-noise report / threshold gate
│   └── util/
│       ├── atmosphere.py       # AFGL profile -> layers (p, T, gas columns)
│       ├── absorption.py       # HITRAN parse + Voigt LBL -> per-layer O2/H2O OT; Rayleigh OT
│       ├── tips.py             # TIPS-2021 Q(296)/Q(T)
│       ├── optics.py           # air <-> vacuum wavelength (Edlén 1966)
│       ├── solar.py            # CU solar spectrum reader/interpolator
│       ├── er3t_abs.py         # er3t atm/abs adapters + exact per-0.001nm Rayleigh
│       └── mca_out_lbl.py      # per-g-point radiance reader -> reflectance/radiance spectrum
├── tests/
│   └── test_absorption.py      # physics validation suite (V1–V8)
└── out/                        # simulation output (not committed)
```

---

## Requirements

- Python 3 with **numpy**, **scipy**, **h5py**.
- **er3t** (provides the MCARaTS Python interface). The driver first tries a plain
  `import er3t` (e.g. an editable `pip install -e` in the active conda env); if that
  fails it adds **`$ER3T_HOME`** — the directory *containing* the `er3t` package —
  to `sys.path`:
  ```bash
  export ER3T_HOME=/path/to/er3t          # the dir that holds the er3t/ package
  ```
- **MCARaTS v0.10.4** compiled executable, located via the environment variable
  **`MCARATS_V010_EXE`**:
  ```bash
  export MCARATS_V010_EXE=/path/to/mcarats/v0.10.4/src/mcarats
  ```

The prescribed input data files, the TIPS-2021 partition sums
(`TIPS_2021_PYTHON/QTpy/`), and the compiled RT engine are not part of this
repository and must be provided separately. Their locations are configurable via
environment variables:

| variable | what | default |
|---|---|---|
| `O2BAND_DATA_DIR` | prescribed inputs (HITRAN/AFGL/solar) | in-repo `data/` |
| `O2BAND_QTPY_DIR` | TIPS-2021 QTpy tables | `$O2BAND_DATA_DIR/TIPS_2021_PYTHON/QTpy`, else in-repo `src/TIPS_2021_PYTHON/QTpy/` |
| `O2BAND_OUT_DIR` | simulation + result files | in-repo `out/` |
| `ER3T_HOME` | dir containing the `er3t` package | (unset — uses an installed er3t) |
| `MCARATS_V010_EXE` | MCARaTS v0.10.4 executable | (unset — required) |

On Linux the data/output defaults instead point at a scratch base rather than the
in-repo dirs. For the CURC/Blanca deployment (concrete scratch paths, module
loads, conda env, and batch/array run scripts), see
[`CURC_NOTES.md`](CURC_NOTES.md).

---

## Usage

### Run the validation suite

```bash
python tests/test_absorption.py
```

Checks HITRAN parsing, the AFGL profile, TIPS-2021, air↔vacuum conversion, Voigt
width & area normalization, Rayleigh OT vs. er3t, O₂ band OT, the sub-range guard,
and reproducibility.

### Run a simulation

```python
from sim_o2band import O2BandConfig, O2BandSim

cfg = O2BandConfig(
    bands=('o2a', 'o2b'),
    szas=(0.0, 30.0, 60.0),
    albedos=(0.0, 0.1),
    z_top=120.0,            # full AFGL profile (default 70 km)
    photons=1e6, Nrun=3,    # per g-point photons; floored by min_photons_per_g
    wvl_range=None,         # or (a, b) in air nm for a test sub-window
)
sim = O2BandSim(cfg)
sim.run()                  # chunked, skip-if-done, resumable
sim.assemble()             # per-band o2a.h5 / o2b.h5 + merged o2band_benchmark.h5
```

> **macOS note.** `sim.run()` uses `multiprocessing.Pool`, so it must be called
> under an `if __name__ == '__main__':` guard.

#### Command-line interface

`sim_o2band.py` also has a CLI (used by the CURC batch runner):

```bash
python src/sim_o2band.py --test                    # fast prototype (narrow window, low photons)
python src/sim_o2band.py                            # frozen Phase-1 grid (both bands)
python src/sim_o2band.py --bands o2a --wvl-range 763.20 763.30 --photons 1e4 --nrun 2
python src/sim_o2band.py --ncpu 16 --overwrite      # force recompute on 16 cores
python src/sim_o2band.py --help                     # all flags
```

Flags: `--bands`, `--szas`, `--albedos`, `--wvl-range A B`, `--photons`, `--nrun`,
`--chunk-size`, `--z-top`, `--ncpu`, `--out-dir`, `--test`, `--overwrite`,
`--no-assemble`. Explicit flags override the `--test` presets and frozen defaults.

Key knobs:

- **`wvl_range`** — restrict to a contiguous in-band sub-window (validated against
  the band edges). The 0.001 nm grid is lattice-anchored, so a sub-window's points
  are an exact subset of the full-band grid and its chunk files compose with a full
  run. A partial band is flagged as a deviation in metadata.
- **`chunk_size`** — g-points per MCARaTS call (checkpoint granularity only; not a
  physics knob — Rayleigh is exact per 0.001 nm via `set_per_g_rayleigh`).
- **`photons` / `min_photons_per_g`** — each g-point is an independent
  monochromatic solve and receives at least `min_photons_per_g` photons.

### Resume / re-run

Every work unit `(band, SZA, albedo, chunk)` writes a deterministic, lattice-
anchored chunk file with an atomic (tmp+rename) write. A valid existing chunk is
skipped, so an interrupted run resumes by simply rerunning; `assemble()` stitches
the chunk files into the final HDF5s.

> **Running on a cluster.** For the CURC/Blanca deployment — environment setup,
> scratch paths, the single-node SBATCH runner, and the parallel job-array
> pipeline — see [`CURC_NOTES.md`](CURC_NOTES.md).

---

## Output HDF5 schema

`out/o2band_benchmark.h5` (merged) and `out/o2a.h5` / `out/o2b.h5` (per band):

```
metadata/                     # every setting used, input identities, git commit
<band>/
  wvl                (Nwvl,)              # air wavelength (nm)
  sza                (Nsza,)              # deg
  albedo             (Nalb,)
  f0                 (Nwvl,)              # incident solar irradiance (W m-2 nm-1)
  reflectance        (Nsza, Nalb, Nwvl)   # TOA reflectance
  reflectance_stderr (Nsza, Nalb, Nwvl)   # MC standard error
  radiance           (Nsza, Nalb, Nwvl)   # W m-2 nm-1 sr-1
  radiance_stderr    (Nsza, Nalb, Nwvl)
  optical_thickness/
    o2_layer         (Nlay, Nwvl)         # O2 absorption OT per layer
    h2o_layer        (Nlay, Nwvl)
    rayleigh_layer   (Nlay, Nwvl)
    o2_column        (Nwvl,)              # column-integrated
    h2o_column       (Nwvl,)
    rayleigh_column  (Nwvl,)
  atmosphere/        # z/p/T on layers and levels
```

---

## Reproducibility

All runs are reproducible from committed configuration: pinned constants and
conventions, deterministic outputs, and output metadata that records the bands,
resolution, SZA/albedo, line shape, HITRAN/TIPS versions, the MCARaTS executable
and version, input-file identities, git commit, and any deviations from the
prescribed settings.
