# O2 A/B-band RT benchmark — experiment & product summary

High-spectral-resolution line-by-line benchmark of **top-of-atmosphere (TOA)
reflectance** in the molecular-oxygen **A-band (757–772 nm)** and **B-band
(680–695 nm)** under tightly prescribed clear-sky conditions, produced for the
KNMI-led O2 A/B-band radiative-transfer model **intercomparison** (Wang, Ferlay,
Herbin, Preusker, Wang, Vidot, Duan, Stammes).  This document summarizes what was
computed, the delivered product, and how it was validated.  See
[PLAN.md](PLAN.md) for the physics blueprint and [EVAL_REPORT.md](EVAL_REPORT.md)
for the full validation.

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
| Cloud / aerosol / polarization | **none** (Phase-1 clear-sky) |

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
  f0                 (15001,)              # CU solar irradiance (W m-2 nm-1)
  reflectance        (3, 2, 15001)         # TOA reflectance rho (SZA, albedo, wvl)
  reflectance_stderr (3, 2, 15001)         # Monte-Carlo standard error
  radiance           (3, 2, 15001)         # W m-2 nm-1 sr-1 (= rho*mu0*F0/pi)
  radiance_stderr    (3, 2, 15001)
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

- All reflectance / radiance / optical-thickness values finite.
- MC noise (unbiased sample std over the 3 runs, ddof=1, at 10⁷ photons/g):
  median relative reflectance stderr ~2×10⁻⁴; 95th percentile ≤1.8×10⁻³ (O2A) /
  ≤5×10⁻⁴ (O2B) per (SZA, albedo) — well below the 0.01 sign-off gate.  Noise
  scaled by the textbook 1/√N (3.2×) from the 10⁶-photon predecessor run,
  confirming pure photon statistics with no systematic noise floor.
- Reflectance is albedo-independent at saturated line cores (surface screened by
  the optically thick atmosphere) — correct physics, exploited by pressure/height
  retrievals.

---

## 3. Validation

Each physics component was cross-checked against an independent public reference
(difference statistics; full detail in [EVAL_REPORT.md](EVAL_REPORT.md)):

| component | reference | agreement |
|---|---|---|
| Rayleigh cross-section | Bucholtz (1995) | ~0.03% |
| Rayleigh column OT | Hansen & Travis (1974) | <0.1% |
| O2 absorber amount | canonical 0.2095 dry-air VMR | −0.17% (H₂O dilution) |
| O2 line-by-line engine | HAPI (independent Voigt LBL, matched HITRAN 2020) | ~0.1–0.5% |
| O2 A-band absorption | OCO ABSCO v5.2 | line cores ~1%; continuum by design (see §4) |
| RT solver + convention | libRadtran/DISORT (window + injected in-band OD, both bands) | <0.4% rel RMS, corr 1.0 |

The absorption engine, scattering, and RT transport/convention all match
independent references to well under 1%.

---

## 4. Known Phase-1 choices (state these to users)

These are deliberate, documented simplifications — not errors:

1. **O2 collision-induced absorption + line mixing: excluded.**  Against the
   operational OCO ABSCO tables this leaves the O2 A-band *continuum* optical
   depth ~0.01 too low (comparable to the Rayleigh OT), affecting window/continuum
   reflectance.  Candidate for a future phase.
2. **HITRAN 2020** (as prescribed).  Current HITRAN 2024 raised O2 **A-band** line
   intensities by ~1.3% (B-band unchanged); migrating editions would shift A-band
   absorption accordingly.
3. **Absolute radiance uses the smooth CU solar spectrum** (no Fraunhofer structure
   at 0.001 nm).  Reflectance is F₀-independent and therefore unaffected; only the
   absolute-radiance product would need a high-resolution solar spectrum for
   realistic Fraunhofer lines.

Out of scope for Phase 1: clouds, aerosols, polarization, instrument convolution,
and real-data comparison.

---

## 5. Reproducibility

Every run is reproducible from committed configuration and the `metadata/` group:
bands, resolution, SZA/albedo, line shape, HITRAN/TIPS versions, wavelength
convention, Rayleigh model, reflectance definition, photons/Nrun, z_top,
input-file identities, MCARaTS executable/version, and git commit.  Code:
`src/sim_o2band.py` (driver) + `src/util/` (absorption, atmosphere, TIPS, optics,
solar, er3t/MCARaTS adapters); evaluation: `src/eval_*.py`.
