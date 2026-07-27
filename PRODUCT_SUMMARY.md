# O2 A/B-band RT benchmark — product summary

**Delivery v2 (2026-07-26)**, run `z120_p1e7_n3_vac_c50`.  Supersedes the first
delivery (air grid, ~2 cm⁻¹ wing cutoff) — see §4.  For the full experiment
record, validation status, and Phase-1 caveats see [SUMMARY.md](SUMMARY.md).

High-spectral-resolution line-by-line benchmark of **top-of-atmosphere (TOA)
reflectance** in the molecular-oxygen **A-band (757–772 nm)** and **B-band
(680–695 nm)** under tightly prescribed clear-sky conditions, produced for the
KNMI-led O2 A/B-band radiative-transfer model **intercomparison** (Wang, Ferlay,
Herbin, Preusker, Wang, Vidot, Duan, Stammes).

---

## 1. Experiment

A Monte-Carlo radiative-transfer simulation (MCARaTS; Iwabuchi, 2006) through
the er3t / EaR3T interface (Chen et al., 2023) of a horizontally homogeneous,
clear-sky, plane-parallel atmosphere, with line-by-line molecular absorption
computed from HITRAN 2020 and Rayleigh scattering from Bodhaine (1999).

### Prescribed settings (frozen)

| | |
|---|---|
| Bands | B-band **680–695 nm**, A-band **757–772 nm** |
| Spectral grid | **0.001 nm**, **vacuum** wavelengths, 15001 points/band |
| Spectral sampling | **point samples**, not 0.001 nm bin averages |
| Line data | **HITRAN 2020** |
| Line shape | **Voigt** (Doppler ⊛ Lorentz, Faddeeva `wofz`) |
| Line-wing cutoff | **±50 cm⁻¹**, plain truncation |
| Partition sums | **TIPS-2021** (consistent with HITRAN 2020) |
| Atmosphere | AFGL **mid-latitude summer**, surface–120 km (49 layers) |
| Solar zenith angle | **0°, 30°, 60°** |
| Viewing | **nadir**, relative azimuth 0° |
| Surface | **Lambertian**, albedo **0.0** and **0.1** |
| Rayleigh | **Bodhaine et al. (1999)**, at vacuum wavelength |
| RT engine | **MCARaTS v0.10.4**, 1-D IPA, 10⁷ photons/wavelength, 3 runs |
| Cloud / aerosol / polarization | **none** (clear-sky) |
| O2 collision-induced absorption + line mixing | **excluded** |

### Method & conventions

- **Absorption:** per-line HITRAN intensity T-scaling (lower-state energy,
  stimulated emission, TIPS-2021 partition ratio); Voigt shape with Lorentz
  γ_L = (296/T)^n_air·[γ_air(P−P_self)+γ_self·P_self] and Doppler α_D; per-line
  wing cutoff **±50 cm⁻¹** with the line-selection margin derived from it;
  column integration exact for exponential-in-z density.
- **Vacuum wavelengths:** HITRAN ν is vacuum and stays vacuum; air wavelengths
  (Edlén 1966) are written alongside, and every wavelength dataset carries a
  `convention` attribute.
- **Point sampling:** τ(λ) is evaluated *at* each grid wavelength — no
  integration over the ±0.0005 nm cell.
- **Reflectance:** MCARaTS is run unit-source (`Src_flx=1`); ρ(λ) = π·I/(μ₀·F₀) =
  **π·R_raw/μ₀** (μ₀ = cos SZA), which is F₀-independent.  The CU composite solar
  spectrum × Toon SPTS transmittance is folded in only for absolute radiance.
- **Optical thickness** for O2, H2O, O2-O2 CIA, and Rayleigh is output separately
  (layer and column) so that any RT-model difference can be attributed to
  absorption, scattering, transport, or radiometric convention.

---

## 2. Product

Delivered as **HDF5**, self-describing (all settings in `metadata/`):

| file | contents |
|---|---|
| `o2band_benchmark_noradiance.h5` | **intercomparison delivery** — merged file minus every solar-dependent dataset (`radiance`, `radiance_stderr`, `f0`).  Everything in it is F₀-independent (51.9 MB) |
| `o2band_benchmark.h5` | **merged** — both bands as groups `o2a`, `o2b` + `metadata` (55 MB) |
| `o2a.h5`, `o2b.h5` | per-band (datasets at root, 27.5 MB each) |
| `reflectance_o2ab.png`, `mc_noise_o2ab.png` | quick-look figures: TOA reflectance and relative MC noise, all geometries (`src/plot_reflectance.py`) |

### Layout (per band)

```
metadata/                     # every setting, input identities, git commit
<band>/                       # attrs: band, band_range_nm, wavelength_convention
  wvl                (15001,)              # primary grid: VACUUM nm (attrs: units, convention)
  wvl_vac  wvl_air   (15001,)              # both conventions, always present
  nu_vac             (15001,)              # vacuum wavenumber (cm-1)
  sza                (3,)   albedo (2,)    # 0/30/60 deg ; 0.0/0.1
  reflectance        (3, 2, 15001)         # TOA reflectance rho (SZA, albedo, wvl)
  reflectance_stderr (3, 2, 15001)         # Monte-Carlo standard error (ddof=1)
  optical_thickness/
    o2_layer  o2_column                    # (49,15001) and (15001,)
    h2o_layer h2o_column
    o2_cia_layer o2_cia_column             # all zero: CIA excluded by design
    rayleigh_layer rayleigh_column
  atmosphere/                              # z/p/T on 50 levels and 49 layers
```

### Reading it

```python
import h5py
with h5py.File('o2band_benchmark_noradiance.h5') as f:
    wvl = f['o2a/wvl'][:]                          # nm, VACUUM
    wvl_air = f['o2a/wvl_air'][:]                  # nm, air (Edlen 1966)
    rho = f['o2a/reflectance'][:]                  # (SZA, albedo, wvl)
    o2_od = f['o2a/optical_thickness/o2_column'][:]
    meta = dict(f['metadata'].attrs)               # provenance
    print(f['o2a/wvl'].attrs['convention'])        # -> 'vacuum'
```

### Quality

- All reflectance / optical-thickness values finite; ρ ∈ [0, 0.109] (O2A) /
  [0, 0.114] (O2B).
- Grid exact: uniform 0.001 nm in vacuum to 2×10⁻¹⁴ nm; `air_to_vac(wvl_air) ==
  wvl_vac` and `nu_vac == 1e7/wvl_vac` to machine precision; layer optical
  thicknesses sum to the column values exactly.
- MC noise (ddof=1 over 3 runs, 10⁷ photons/g, where ρ > 10⁻³): median relative
  stderr 2.3×10⁻⁴ (O2A) / 1.8×10⁻⁴ (O2B); worst per-(SZA, albedo) 95th percentile
  1.9×10⁻³ / 4.9×10⁻⁴ — well below the 0.01 gate.  Absolute stderr ≤ 8.4×10⁻⁵.
- Reflectance is albedo-independent at saturated line cores (column OD > 20:
  |ρ(0.1) − ρ(0.0)| ≤ 1.5×10⁻⁵) — correct physics, the basis of pressure/height
  retrievals.

---

## 3. Conventions to confirm across the ensemble

Each moves the result at the few-percent level, so a like-for-like comparison
needs them pinned down:

| convention | ours | why it matters |
|---|---|---|
| wavelength | **vacuum** | air vs vacuum is a 0.19–0.21 nm displacement, ~210 grid cells |
| wing cutoff | **50 cm⁻¹** | 25→50 cm⁻¹ moves A-band micro-window column OD by +0.0022 |
| truncation | **plain** | LBLRTM-style pedestal subtraction moves it by −0.0039 — larger, opposite sign |
| sampling | **point** | vs 0.001 nm bin averages: up to 3.5% in OD, 6.9% in ρ on steep flanks (band mean unchanged) |

---

## 4. What changed from delivery v1

1. **Wavelength convention** — v1 reported air wavelengths while the ensemble
   used vacuum, displacing our spectrum by −0.19 to −0.21 nm.  Physics unchanged:
   the A-band head is at 759.564 nm in v1 and 759.559 nm in v2, the same feature
   relabelled.
2. **Line-wing cutoff** — v1's ±N·ν₀/R ≈ 2 cm⁻¹ discarded ~0.055 of column
   optical depth in the A-band micro-windows (more than twice the entire Rayleigh
   OT), leaving between-line reflectance too high.  The A-band peak envelope was
   flat at ~0.107 across the whole band in v1; in v2 it dips to 0.0857 at
   763–764 nm and recovers beyond 769 nm.  Band-mean A-band reflectance −3.2%
   (B-band −0.14%).
3. **Rayleigh** — σ now evaluated at vacuum wavelength, as Bodhaine's fit
   requires (~0.12% correction).

---

## References

- Iwabuchi, H.: Efficient Monte Carlo methods for radiative transfer modeling,
  J. Atmos. Sci., 63, 2324-2339, doi:10.1175/JAS3755.1, 2006.
- Iwabuchi, H., and Okamura, R.: Multispectral Monte Carlo radiative transfer
  simulation by using the maximum cross-section method, J. Quant. Spectrosc.
  Radiat. Transfer, 193, 40-46, doi:10.1016/j.jqsrt.2017.01.025, 2017.
- Chen, H., Schmidt, K. S., Massie, S. T., Nataraja, V., Norgren, M. S.,
  Gristey, J. J., Feingold, G., Holz, R. E., and Iwabuchi, H.: The Education
  and Research 3D Radiative Transfer Toolbox (EaR3T) - Towards the Mitigation
  of 3D Bias in Airborne and Spaceborne Passive Imagery Cloud Retrievals,
  Atmos. Meas. Tech., 16, 1971-2000,
  https://doi.org/10.5194/amt-16-1971-2023, 2023.
