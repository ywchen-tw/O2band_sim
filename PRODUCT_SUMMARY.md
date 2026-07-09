# O2 A/B-band RT benchmark — experiment & product summary

High-spectral-resolution line-by-line benchmark of **top-of-atmosphere (TOA)
reflectance** in the molecular-oxygen **A-band (757–772 nm)** and **B-band
(680–695 nm)** under tightly prescribed clear-sky conditions, produced for the
KNMI-led O2 A/B-band radiative-transfer model **intercomparison** (Wang, Ferlay,
Herbin, Preusker, Wang, Vidot, Duan, Stammes).  

---

## 1. Experiment

A Monte-Carlo radiative-transfer simulation (MCARaTS via er3t) of a horizontally
homogeneous, clear-sky, plane-parallel atmosphere, with line-by-line molecular
absorption computed from HITRAN 2020 and Rayleigh scattering from Bodhaine (1999).

### Prescribed settings (frozen)

| | |
|---|---|
| Bands | B-band **680–695 nm**, A-band **757–772 nm** |
| Spectral grid | **0.001 nm**, **air** wavelengths, 15001 points/band |
| Line data | **HITRAN 2020** |
| Line shape | **Voigt** (Doppler ⊛ Lorentz, Faddeeva `wofz`) |
| Partition sums | **TIPS-2021** (consistent with HITRAN 2020) |
| Atmosphere | AFGL **mid-latitude summer**, surface–120 km (49 layers) |
| Solar zenith angle | **0°, 30°, 60°** |
| Viewing | **nadir**, relative azimuth 0° |
| Surface | **Lambertian**, albedo **0.0** and **0.1** |
| Rayleigh | **Bodhaine et al. (1999)** |
| RT engine | **MCARaTS v0.10.4**, 1-D IPA, 10⁷ photons/wavelength, 3 runs |
| Cloud / aerosol / polarization | **none** (clear-sky) |
| O2 collision-induced absorption + line mixing | **excluded** |

### Method & conventions

- **Absorption:** per-line HITRAN intensity T-scaling (lower-state energy,
  stimulated emission, TIPS-2021 partition ratio); Voigt shape with Lorentz
  γ_L = (296/T)^n_air·[γ_air(P−P_self)+γ_self·P_self] and Doppler α_D; per-line wing
  cutoff ±3·ν₀/R (R=20000); column integration exact for exponential-in-z density.
- **Air wavelengths:** HITRAN vacuum wavenumbers converted via Edlén (1966).
- **Reflectance:** MCARaTS is run unit-source (`Src_flx=1`); ρ(λ) = π·I/(μ₀·F₀) =
  **π·R_raw/μ₀** (μ₀ = cos SZA), which is F₀-independent.  The CU composite solar
  spectrum F₀ is folded in only for the absolute-radiance product.
- **Optical thickness** for O2, H2O, and Rayleigh is output separately (layer and
  column) so that any RT-model difference can be attributed to absorption,
  scattering, transport, or radiometric convention.

---

## 2. Product

Delivered as **HDF5**, self-describing (all settings in `metadata/`):

| file | contents |
|---|---|
| `o2band_benchmark.h5` | **merged** — both bands as groups `o2a`, `o2b` + `metadata` (42 MB) |
| `o2a.h5`, `o2b.h5` | per-band (datasets at root) |
| `reflectance_o2ab.png`, `mc_noise_o2ab.png` | quick-look figures: TOA reflectance and relative MC noise, all geometries (`src/plot_reflectance.py`) |

### Layout (per band)

```
metadata/                     # every setting, input identities, git commit
<band>/
  wvl                (15001,)              # air wavelength (nm)
  sza                (3,)   albedo (2,)    # 0/30/60 deg ; 0.0/0.1
  reflectance        (3, 2, 15001)         # TOA reflectance rho (SZA, albedo, wvl)
  reflectance_stderr (3, 2, 15001)         # Monte-Carlo standard error
  optical_thickness/
    o2_layer  o2_column                    # (49,15001) and (15001,)
    h2o_layer h2o_column
    rayleigh_layer rayleigh_column
  atmosphere/                              # z/p/T on 50 levels and 49 layers
```

### Reading it

```python
import h5py
with h5py.File('o2band_benchmark.h5') as f:
    wvl = f['o2a/wvl'][:]                          # nm (air)
    rho = f['o2a/reflectance'][:]                  # (SZA, albedo, wvl)
    o2_od = f['o2a/optical_thickness/o2_column'][:]
    meta = dict(f['metadata'].attrs)               # provenance
```

### Quality

- All reflectance / optical-thickness values finite.
- MC noise (unbiased sample std over the 3 runs, ddof=1, at 10⁷ photons/g):
  median relative reflectance stderr ~2×10⁻⁴; 95th percentile ≤1.8×10⁻³ (O2A) /
  ≤5×10⁻⁴ (O2B) per (SZA, albedo) — well below the 0.01 sign-off gate.  Noise
  scaled by the textbook 1/√N (3.2×) from the 10⁶-photon predecessor run,
  confirming pure photon statistics with no systematic noise floor.


