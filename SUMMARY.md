# O2 A/B-band RT benchmark — experiment & product summary

**Delivery v2 (2026-07-26)**, run `z120_p1e7_n3_vac_c50`. Supersedes the first
delivery, which used an **air**-wavelength grid and a ~2 cm⁻¹ line-wing cutoff;
both were corrected after participant feedback — see §5.

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
| Line-wing cutoff | **±50 cm⁻¹**, plain truncation (no pedestal subtraction) |
| Partition sums | **TIPS-2021** (consistent with HITRAN 2020) |
| Atmosphere | AFGL **mid-latitude summer**, surface–120 km (49 layers) |
| Solar zenith angle | **0°, 30°, 60°** |
| Viewing | **nadir**, relative azimuth 0° |
| Surface | **Lambertian**, albedo **0.0** and **0.1** |
| Rayleigh | **Bodhaine et al. (1999)**, at vacuum wavelength |
| RT engine | **MCARaTS v0.10.4**, 1-D IPA, 10⁷ photons/wavelength, 3 runs |
| Cloud / aerosol / polarization | **none** (Phase-1 clear-sky) |
| O2 CIA + line mixing | **excluded** |

### Method & conventions

- **Absorption:** per-line HITRAN intensity T-scaling (lower-state energy,
  stimulated emission, TIPS-2021 partition ratio); Voigt shape with Lorentz
  γ_L = (296/T)^n_air·[γ_air(P−P_self)+γ_self·P_self] and Doppler α_D; per-line
  wing cutoff **±50 cm⁻¹** with the line-selection margin derived from it;
  column integration exact for exponential-in-z density.  At a given grid point
  ~900 individual line profiles are summed (A-band: 168 O2 + 730 H2O).
- **Vacuum wavelengths:** HITRAN ν is vacuum and stays vacuum.  Air wavelengths
  (Edlén 1966) are written alongside as `wvl_air`; every wavelength dataset
  carries a `convention` attribute.
- **Point sampling:** τ(λ) is the spectrum evaluated *at* each grid wavelength;
  there is no integration over the ±0.0005 nm cell.  Verified grid-independent —
  the same wavelength on 0.001 / 0.0002 / 0.00005 nm grids gives column OD
  0.985486109 / 0.985486111 / 0.985486111.
- **Reflectance:** MCARaTS is run unit-source (`Src_flx=1`); ρ(λ) = π·I/(μ₀·F₀) =
  **π·R_raw/μ₀** (μ₀ = cos SZA), which is F₀-independent.  F₀ — the CU composite
  continuum times the **Toon (JPL) SPTS** disk-integrated solar transmittance
  (Fraunhofer lines resolved on a 0.01 cm⁻¹ grid) — is folded in only for the
  absolute-radiance product.
- **Optical thickness** for O2, H2O, O2-O2 CIA, and Rayleigh is output separately
  (layer and column) so that any RT-model difference can be attributed to
  absorption, scattering, transport, or radiometric convention.

---

## 2. Product

Delivered as **HDF5**, self-describing (all settings in `metadata/`):

| file | contents |
|---|---|
| `o2band_benchmark_noradiance.h5` | **intercomparison delivery** — merged file minus every solar-dependent dataset (`radiance`, `radiance_stderr`, `f0`) via `src/make_delivery_h5.py`.  Everything it contains is F₀-independent, so no choice of solar model can enter a comparison against it.  51.9 MB |
| `o2band_benchmark.h5` | **merged** — both bands as groups `o2a`, `o2b` + `metadata`.  55 MB |
| `o2a.h5`, `o2b.h5` | per-band (datasets at root).  27.5 MB each |
| `reflectance_o2ab.png`, `mc_noise_o2ab.png` | quick-look figures: TOA reflectance and relative MC noise, all geometries (`src/plot_reflectance.py`) |

### Layout (per band)

```
metadata/                     # every setting, input identities, git commit
<band>/                       # attrs: band, band_range_nm, wavelength_convention
  wvl                (15001,)              # primary grid: VACUUM nm (attrs: units, convention)
  wvl_vac            (15001,)              # vacuum wavelength (nm)
  wvl_air            (15001,)              # air wavelength (nm, Edlen 1966)
  nu_vac             (15001,)              # vacuum wavenumber (cm-1)
  sza                (3,)   albedo (2,)    # 0/30/60 deg ; 0.0/0.1
  f0                 (15001,)              # F0 = CU continuum x SPTS transmittance (W m-2 nm-1)
  reflectance        (3, 2, 15001)         # TOA reflectance rho (SZA, albedo, wvl)
  reflectance_stderr (3, 2, 15001)         # Monte-Carlo standard error (ddof=1)
  radiance           (3, 2, 15001)         # W m-2 nm-1 sr-1 (= rho*mu0*F0/pi)
  radiance_stderr    (3, 2, 15001)
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

- All reflectance / radiance / optical-thickness values finite; ρ ∈ [0, 0.109]
  (O2A) / [0, 0.114] (O2B).
- Grid verified exact: uniform 0.001 nm in vacuum to 2×10⁻¹⁴ nm;
  `air_to_vac(wvl_air) == wvl_vac` and `nu_vac == 1e7/wvl_vac` to machine
  precision; layer optical thicknesses sum to the column values exactly.
- MC noise (unbiased sample std over the 3 runs, ddof=1, at 10⁷ photons/g,
  where ρ > 10⁻³): median relative stderr 2.3×10⁻⁴ (O2A) / 1.8×10⁻⁴ (O2B);
  worst per-(SZA, albedo) 95th percentile 1.9×10⁻³ (O2A) / 4.9×10⁻⁴ (O2B) —
  well below the 0.01 sign-off gate.  Absolute stderr ≤ 8.4×10⁻⁵ everywhere.
- Reflectance is albedo-independent at saturated line cores (column OD > 20:
  |ρ(0.1) − ρ(0.0)| ≤ 1.5×10⁻⁵) — the surface is screened by the optically thick
  atmosphere, correct physics and the basis of pressure/height retrievals.

---

## 3. Validation

Each physics component is cross-checked against an independent public reference
(difference statistics; full detail in [EVAL_REPORT.md](EVAL_REPORT.md)).

**Status:** re-run against v2 on 2026-07-26 except where marked.  The evaluation
scripts were themselves corrected for the new conventions first — the HAPI check
had been comparing an air-window, ±5 cm⁻¹ line subset with an unmatched wing
against our ±51 cm⁻¹ selection, so it could not have seen the v1 cutoff defect.

| component | reference | agreement (v2) | status |
|---|---|---|---|
| Rayleigh cross-section | Bucholtz (1995) | +0.03% (O2A) / +0.02% (O2B), corr 1.00000 | re-run |
| Rayleigh column OT | Hansen & Travis (1974) | +0.04% @550, −0.05% @688, −0.07% @760 nm | re-run |
| O2 absorber amount | canonical 0.2095 dry-air VMR | −0.17% (H₂O dilution) | re-run |
| O2 line-by-line engine | HAPI, matched HITRAN 2020 **and matched 50 cm⁻¹ wing** | peak +0.39% / −0.04% (O2A), −1.05% / −0.22% (O2B); corr 0.9997–0.99997 | re-run |
| O2 A-band continuum | OCO ABSCO v5.2 | — | **pending** (ABSCO tables are on CURC only) |
| RT solver + convention | libRadtran/DISORT 16-stream | window rel RMS 0.44%; in-band (col OD 0.10–2.67) rel RMS 0.19%; corr 1.00000 | re-run |
| Solar transmittance (Toon SPTS 2024) | OCO L2 solar model | bias −7×10⁻⁶, RMS 3×10⁻⁴ | **pending** (v1 value; solar model unchanged in v2) |

Product gates on the v2 run: MC noise worst rel-stderr p95 **0.0019** (O2A) /
**0.0005** (O2B) against a 0.01 threshold — **PASS** both bands;
`tests/test_rt_regression.py` **PASS** both bands (Rayleigh τ within 0.07%,
saturated cores, albedo monotonicity across 3 SZA × 2 albedo).

Additionally, `tests/test_wings.py` (7/7) bounds the wing-truncation error
against a brute-force **no-cutoff** Voigt sum — a reference the HAPI check
structurally cannot provide, since it truncates its own comparison the same way.
`tests/test_absorption.py` (10/10) covers parsing, profile, TIPS, Voigt width and
area, Rayleigh, and reproducibility.

---

## 4. Known Phase-1 choices (state these to users)

These are deliberate, documented simplifications — not errors:

1. **O2 collision-induced absorption + line mixing: excluded.**  A HITRAN CIA
   reader exists (`util/cia.py`, `O2-O2_2024.cia` covers both bands) and measures
   a column CIA OT of 0.0022–0.0034 (A-band) / ~0.0003 (B-band), but it is off by
   design.  Note CIA and wider wings are **not additive corrections**: with
   untruncated wings our between-line σ already exceeds ABSCO — which itself
   includes CIA *and* line mixing — by ~7×, because a pure Voigt far wing
   overshoots the real sub-Lorentzian O2 wing.  Either plain Voigt with a stated
   cutoff and no CIA (this delivery), or line mixing + χ-factor + CIA together.
2. **HITRAN 2020** (as prescribed).  Current HITRAN 2024 raised O2 **A-band** line
   intensities by ~1.3% (B-band unchanged); migrating editions would shift A-band
   absorption accordingly.
3. **Three conventions that differ between models** and should be confirmed
   across the ensemble, since each moves the result at the few-percent level:
   the **wing cutoff value** (ours 50 cm⁻¹), the **truncation convention**
   (ours plain; LBLRTM-style pedestal subtraction shifts A-band micro-window
   column OD by −0.0039 against +0.0022 for 25→50 cm⁻¹, i.e. larger and opposite
   in sign), and **point sampling vs 0.001 nm bin averaging** (ours point; the
   difference reaches 3.5% in column OD and 6.9% in reflectance on steep line
   flanks, though the band mean is unchanged to +0.0001%).
4. **Absolute radiance folds F₀ = CU continuum × Toon SPTS solar transmittance**
   (disk-integrated, 2024-07-31 release), so Fraunhofer lines are resolved on the
   0.001 nm grid; validated against the OCO L2 solar model in the A-band overlap
   (§3).  The *absolute scale* still comes from the smooth CU composite continuum.
   Reflectance is F₀-independent throughout, and the delivery file
   `o2band_benchmark_noradiance.h5` omits the solar-dependent products
   entirely (`radiance`, `radiance_stderr`, `f0`).

Out of scope for Phase 1: clouds, aerosols, polarization, instrument convolution,
and real-data comparison.  Deferred to a future experiment (PLAN.md §8a):
cell-averaged rather than point-sampled output.

---

## 5. What changed from delivery v1

Two defects were reported by the participants and confirmed:

1. **Wavelength convention.**  v1 reported **air** wavelengths while the rest of
   the ensemble used vacuum, so our spectrum arrived displaced by −0.19 nm (O2B)
   to −0.21 nm (O2A), about 210 grid cells.  The protocol never stated a
   convention.  The physics was never affected: the A-band head sits at
   759.564 nm in v1 and 759.559 nm in v2 — the same physical feature, relabelled.
2. **Line-wing cutoff.**  v1 truncated every Voigt at ±N·ν₀/R ≈ 2 cm⁻¹, whose
   justification ("3× the resolution element gives 99.7% coverage") is
   Gaussian-tail reasoning that does not hold for a Lorentzian far wing.  Against
   a brute-force no-cutoff sum this discarded ~0.055 of column optical depth in
   the A-band micro-windows — more than twice the entire Rayleigh OT — leaving
   between-line reflectance too high.  In v1 the A-band peak envelope was flat at
   ~0.107 (99.8% of continuum) across the *whole* band; in v2 it dips to 0.0857
   at 763–764 nm where lines are densest, recovering toward continuum beyond
   769 nm.  Points with ρ > 0.10 in 762–766 nm fall from 24.6% to 3.1%, and
   band-mean A-band reflectance drops 3.2% (B-band 0.14%).

Also corrected: Rayleigh σ is now evaluated at vacuum wavelength (Bodhaine's fit
is parameterised that way; feeding it air λ biased σ high by ~0.12%).

---

## 6. Reproducibility

Every run is reproducible from committed configuration and the `metadata/` group:
bands, resolution, SZA/albedo, line shape, wing cutoff, wavelength convention,
spectral sampling, HITRAN/TIPS versions, Rayleigh model, CIA setting, reflectance
definition, solar source (continuum + SPTS file), photons/Nrun, z_top,
input-file identities, MCARaTS executable/version, and git commit.  The output
directory additionally carries `_physics_config.json`, and the driver refuses to
mix runs of differing physics into one directory.  Code:
`src/sim_o2band.py` (driver) + `src/util/` (absorption, atmosphere, TIPS, optics,
solar, CIA, er3t/MCARaTS adapters); evaluation: `src/eval_*.py`; regression:
`tests/`.

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
