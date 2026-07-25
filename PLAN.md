# PLAN.md — O2 A/B-band RT intercomparison (Phase 1)

Scientific + architectural blueprint. Per CLAUDE.md §4.1 this must be reviewed
and signed off before physics-affecting code is considered landed. Status of
each item is marked **[PROPOSED]**, **[NEEDS SIGN-OFF]**, or **[AGREED]**.

> Process note: an initial *draft* of the absorption code
> (`src/util/atmosphere.py`, `src/util/solar.py`, `src/util/absorption.py`) was
> written before this plan existed. It is to be treated as a concrete proposal
> to review against this plan — **not** as validated code — until the open
> decisions in §7 are resolved and the §6 validation table is green.

---

## 1. Goal

High-spectral-resolution benchmark of TOA reflectance in the O2 A-band
(757–772 nm) and B-band (680–695 nm) under the frozen Phase-1 settings in
CLAUDE.md §3, plus separable O2 and Rayleigh optical thickness on matching grids.

RT engine: **er3t / MCARaTS** 1D (per user instruction), using
`sim_reference.py` + `util_reference/` as the structural template.

---

## 2. Data flow

```
HITRAN 2020 .par ─┐
AFGL mid-lat-sum ─┼─► line-by-line Voigt cross-sections ─► layer O2/H2O abs. OT ─┐
CU solar spectrum ┘                                                             │
                                                                                ▼
                        er3t abs object (Ng = fine-wvl points, coef[...]) ─► mca_atm_1d
                                                                                │
                        AFGL layers ─► atm object (er3t) ──────────────────────┤
                        Rayleigh OT (Bodhaine, computed in mca_atm_1d) ─────────┤
                        Lambertian albedo (0.0, 0.1) ─► mca_sfc ────────────────┤
                                                                                ▼
                                       MCARaTS 1D radiance ─► TOA reflectance
                                                                                │
                        O2 OT + Rayleigh OT + reflectance + metadata ─► output h5
```

## 3. Modules

| module | responsibility | status |
|---|---|---|
| `util/atmosphere.py` | read AFGL profile → levels + layers (p, T, gas columns) | implemented |
| `util/solar.py` | read CU composite solar spectrum, interpolate | implemented |
| `util/tips.py` | TIPS-2021 partition sums Q(296)/Q(T) per isotopologue | implemented |
| `util/optics.py` | vacuum↔air wavelength (Edlén 1966) | implemented |
| `util/cia.py` | HITRAN `.cia` reader + O2-O2 CIA OT (∫n²dz); coverage-based file choice | implemented (off by default, §7.4) |
| `util/absorption.py` | HITRAN parse + Voigt LBL → per-layer O2/H2O OT; Rayleigh OT | implemented |
| `util/er3t_abs.py` | er3t `atm` adapter + `abs` wrapper + `set_per_g_rayleigh` (exact per-0.001nm Rayleigh) | implemented |
| `util/mca_out_lbl.py` | per-g radiance reader (spectrum); absolute radiance + reflectance + MC stderr | implemented |
| `sim_o2band.py` | driver: config, cached absorption, chunked skip-if-done RT, assemble per-band + merged h5 + metadata | implemented (prototype-validated) |

All four §7 physics decisions are locked and reflected in code. Absorption modules
run end-to-end (~1.5 s/band at 0.001 nm). Informal sanity checks pass (see below);
**formal validation suite (§6) and RT wiring still pending review/sign-off.**

## 4. Physics & conventions (state every one — CLAUDE.md §4.2)

- **Line intensity T-scaling**: HITRAN convention with lower-state energy E″,
  stimulated-emission factor, and partition ratio Q(296)/Q(T). Ref T = 296 K,
  ref P = 1 atm = 1013.25 hPa. c₂ = 1.4387769 cm·K.
- **Line shape**: Voigt = Gaussian(Doppler) ⊛ Lorentzian(pressure), via the
  Faddeeva function `scipy.special.wofz`. (Mandated — CLAUDE.md §2.)
- **Lorentz HWHM**: γ_L = (296/T)^n_air · [γ_air(P−P_self) + γ_self·P_self],
  pressures in atm; P_self = P · VMR_gas.
- **Doppler HWHM**: α_D = ν₀·√(2 ln2 kB T / m)/c, per-isotopologue mass m.
- **Isotopic abundance**: HITRAN S already includes terrestrial abundance ⇒
  multiply by *total* species column (AFGL o2/h2o). No double-counting.
- **H2O treatment** (**[AGREED]**): `include_h2o` flag toggles H2O *absorption*
  lines only. Broadening is fixed to the current setup — water vapor is always
  part of the foreign (air) pressure and broadens O2 at `γ_air` (no separate
  `γ_H2O`). Toggling H2O absorption leaves the O2 optical depth bit-identical
  (verified). A differential `γ_H2O` is a possible future data-dependent add-on.
- **Column integration**: gas number density assumed exponential in z within a
  layer ⇒ exact integral N = dz·(n_bot−n_top)/ln(n_bot/n_top).
- **Rayleigh**: Bodhaine et al. (1999) cross-section × air column, evaluated at
  **vacuum** wavelength (Bodhaine's fit is parameterised in vacuum λ; feeding it
  air λ biased σ high by ~0.12%, corrected 2026-07-25). er3t's `mca_atm_1d` also
  computes Bodhaine internally; we must ensure our diagnostic matches what the
  RT actually uses (see §6).
- **Wavenumber↔wavelength**: λ_vac[nm] = 1e7 / ν[cm⁻¹]; λ_air = λ_vac / n(λ)
  (Edlén 1966). Output grid is vacuum — §7 decision 2.

## 5. Discretization

- **Spectral**: uniform 0.001 nm **vacuum**-wavelength grid per band (§7.2).
  O2A: 757–772 nm → 15001 pts; O2B: 680–695 nm → 15001 pts. The air wavelength
  of every grid point is carried alongside (`wvl_air`).
- **Spectral sampling** (**[documented 2026-07-25]**, `src/eval_grid_sampling.py`):
  outputs are **point samples** at each grid wavelength, not averages over the
  0.001 nm interval. The grid is marginal for the narrowest lines — Voigt FWHM
  is 0.0054 nm at the surface but 0.0015 nm near 10 hPa, i.e. **1.5 grid points
  per FWHM** where lines are Doppler-limited. Measured against a 5× finer
  (0.0002 nm) grid over 763.0–763.2 nm: column OD point-vs-cell-average differs
  by median +0.01%, p95 1.2%, max 3.5%; reflectance by p95 0.34%, max 6.9%
  locally on steep flanks — above the ~0.18% MC noise. **The window-mean
  reflectance is unaffected (+0.0001%)**, so this is a point-by-point
  comparison issue, not a bias. Decomposition: the exp-nonlinearity
  (⟨ρ(τ)⟩ vs ρ(⟨τ⟩)) contributes only ~3% of the difference; the rest is plain
  sampling of a curved spectrum. Deep line cores show the largest OD difference
  (−3.5%) but are radiatively irrelevant there (ρ underflows to 0 either way) —
  it matters for the **optical-thickness deliverable**, which the protocol
  explicitly asks for. Recorded as `metadata/spectral_sampling`; a third
  convention to confirm with the participants, alongside §7.3's cutoff value and
  truncation convention.
- **Vertical**: native AFGL mid-lat-summer levels, capped at z_top (proposed
  70 km; absorber amount above is negligible). Layers between adjacent levels;
  representative P = √(P_top·P_bot) (log-mean), T = ½(T_top+T_bot).
- **Line-wing cutoff**: ±50 cm⁻¹ per line, plain truncation
  (**[REVISED 2026-07-25]**, §7.3); the line-selection margin is derived from
  it, never hardcoded.

## 6. Validation table (must be green before "done" — CLAUDE.md §4.3)

| # | check | criterion | status |
|---|---|---|---|
Automated suite: `tests/test_absorption.py` (`python tests/test_absorption.py`).
Status below reflects that run — **9/9 green** (V7 deferred to RT wiring).

| # | check | criterion | status |
|---|---|---|---|
| V1 | HITRAN parse | first record matches raw file; O2=860, H2O=33159; masses mapped | ✅ |
| V2 | AFGL profile | p_sfc=1013 hPa, T=294.2 K, O2 VMR=0.2092, O2 col=4.51e24 | ✅ |
| Vt | TIPS-2021 | Q(296) exact; ratio vs rigid-rotor law dev O2 0.09%, H2O 0.46% | ✅ |
| Vo | air↔vac | shift 0.19–0.21 nm; λ_air = λ_vac/n round-trip < 1e-6 nm | ✅ |
| V3 | Voigt width vs P | FWHM within 1% of Olivero-Longbothum; Lorentz(sfc)/Doppler(top) regimes; monotonic | ✅ |
| V4 | Area normalisation | ∫Voigt dν = 1 to 0.009% (Doppler/Lorentz/mixed) | ✅ |
| V5 | Rayleigh OT | cross-section == er3t Bodhaine (rtol<1e-12); column OD vs er3t <2.3%; λ⁻⁴ | ✅ |
| V6 | O2 band OT | A max 572 ≫ B max 34 (A/B=16.9); line/continuum contrast 7e4 | ✅ |
| V7 | Reflectance limits | saturated core → ρ=0 (both albedos); ρ monotonic in albedo; window ρ matches analytic Rayleigh+Lambertian (A·T↓T↑) to ~2% ⇒ confirms ρ=π·R_raw/μ0 | ✅ (RT wired; prototype window, 1e4 photons) |
| V8 | Reproducibility | two runs bit-identical (deterministic) | ✅ |

Wing-cutoff / grid-convention suite: `tests/test_wings.py` — **7/7 green**.
Added 2026-07-25 because the HAPI cross-check (#4 in EVAL_REPORT.md) is
structurally blind to wing truncation: it builds its comparison table from the
same band ±5 cm⁻¹ line subset, so both codes were truncated identically. The
reference here is a brute-force Voigt sum over ±300 cm⁻¹ with **no cutoff**.

| # | check | criterion | status |
|---|---|---|---|
| W1 | default cutoff accuracy | micro-window column OD within 5% of no-cutoff brute force | ✅ (1.2% low at 50 cm⁻¹) |
| W2 | legacy cutoff is caught | legacy ~2 cm⁻¹ loses >50% of micro-window OD (guards W1's sensitivity) | ✅ (91% lost) |
| W3 | monotonicity | column OD non-decreasing over cutoff 5→10→25→50 cm⁻¹ | ✅ |
| W4 | margin/cutoff coupling | margin < cutoff raises ValueError, never runs silently | ✅ |
| W5 | grid equivalence | air-grid OD interpolated to the vacuum grid's air λ reproduces it | ✅ (0.23% median) |
| W6 | Edlén shift | wvl_vac − wvl_air = 0.21 nm; air→vac round-trip exact | ✅ |
| W7 | CIA separability | off by default (od_cia ≡ 0); on leaves line OD bit-identical | ✅ |

RT-integration validation was run on a prototype sub-window (763.20–763.30 nm,
1e4 photons/g, Nrun=2, MCARaTS v0.10.4 IPA). Full-grid production photon count
(1e6/g) and MC-noise threshold check (2f) still to be benchmarked before sign-off.

## 7. Decisions

1. **Partition function Q(296)/Q(T)** — **[AGREED]** Exact **HITRAN TIPS-2021**.
   Need tabulated total internal partition sums for O2 (66,68,67) and H2O
   (161,181,171,162,172 = isos 1,2,3,4,7 present in the file). Sourcing: embed
   TIPS-2021 tables (Gamache et al. 2021) as a committed data file — see §9.
   Draft power-law `partition_ratio()` to be **replaced**.
2. **Air vs. vacuum wavelength** — **[REVISED 2026-07-25]** Band limits and the
   0.001 nm grid are **vacuum** wavelengths (`grid='vac'`, the default).
   *Superseded:* the original decision was **air** (λ_air = λ_vac / n(λ), Edlén
   1966), and the first delivery used it. The intercomparison protocol
   ([O2_abs_intercomparison.md](O2_abs_intercomparison.md) §"Wavelength range")
   never states a convention, and the other participants (DAK, DISAMAR, MYSTIC,
   UNL-VRTM, GRASP, IAP) all read it as vacuum — so our spectrum arrived
   displaced by −0.19 nm (O2B) to −0.21 nm (O2A), ≈210 grid cells. The physics
   was never affected (HITRAN ν is vacuum throughout, and every line position
   round-trips to <0.001 nm), only the reported axis. Output now carries
   `wvl` (primary, vacuum), plus `wvl_vac`, `wvl_air` and `nu_vac`, each with a
   `convention` attribute on the dataset itself — the root-level
   `metadata/wavelength_convention` alone was too easy for a comparison script
   to miss. `grid='air'` still reproduces the original delivery exactly.
   *Note:* a vacuum 757–772 nm product spans air 756.79–771.79 nm, so it is a
   recomputation, not a reinterpolation of the existing file.
3. **Line-wing cutoff** — **[REVISED 2026-07-25]** Per-line wing cutoff =
   **±50 cm⁻¹**, **plain truncation** (no pedestal subtraction).
   *Superseded:* ±N·(ν₀/R) with N = 3, R = 20000 ⇒ ≈1.97 cm⁻¹ (O2A). The
   justification for that form — "ν₀/R is a characteristic width, 3× ⇒ 99.7%
   coverage" — is **Gaussian-tail reasoning and does not hold for a Voigt**,
   whose far wing is Lorentzian with area converging only as 1/Δν. Against a
   brute-force no-cutoff Voigt sum, the 2 cm⁻¹ cutoff discarded **~0.055 of
   column optical depth in the A-band micro-windows** (0.0009 → 0.0378 over
   762–766 nm, more than twice the whole Rayleigh OT), leaving between-line
   reflectance ~6–10% too high — the second discrepancy the participants
   reported. 25 cm⁻¹ recovers 92% of that deficit, 33 cm⁻¹ 95%, 50 cm⁻¹ 98%.
   **50 cm⁻¹ chosen**: the protocol prescribes "Voigt line shape" and says
   nothing about a cutoff, so the cutoff is a pure numerical artifact that a
   benchmark should minimise, and at 0.001 nm a full band costs 9.9 s at
   50 cm⁻¹ vs 5.5 s at 25 (18.0 s at 100) — cost is not a constraint. The RT
   stage is unaffected.
   The **line-selection margin is now derived from the cutoff**
   (`required_margin_cm`) and a margin narrower than the cutoff raises: the old
   hardcoded 5 cm⁻¹ margin would have silently truncated the wings of
   out-of-window lines and mimicked a too-short cutoff.
   *Open — two questions for the participants, not one:* the cutoff **value**
   and the **truncation convention**. LBLRTM pairs its 25 cm⁻¹ with subtracting
   f(cutoff) inside the cutoff, and that convention change is *larger* than the
   value change and of opposite sign — measured on the A-band micro-windows:
   25→50 cm⁻¹ is **+0.0022** column OD, while plain→pedestal-subtracted at
   25 cm⁻¹ is **−0.0039**. A participant reporting "25 cm⁻¹" is ambiguous until
   that is pinned down. In reflectance these are 0.4–0.8% at albedo 0.1, versus
   ~0.18% p95 MC noise — detectable, but far below the ~6–10% defect being
   fixed.
4. **O2 CIA / continuum** — **[AGREED, with caveat]** **Excluded** by default
   for Phase 1; flagged in metadata. HITRAN line data only. `util/cia.py` can
   now read a HITRAN `.cia` file (`cia='auto'`, selected **by spectral
   coverage** — `O2-O2_exp_2021.cia` is the 1300–1850 cm⁻¹ fundamental and
   covers neither band; `O2-O2_2024.cia` covers both) and adds τ = k(ν,T)·∫n²dz,
   kept separable in `od_cia` / `o2_cia_column`. Measured effect: A-band column
   CIA OT 0.0022–0.0034, B-band ~0.0003.
   *Caveat before switching it on:* with untruncated wings our between-line
   surface σ is already ~1.6e-26 cm², ~7× **above** ABSCO's 2.2e-27 — and ABSCO
   already includes CIA **and** line mixing. A pure Voigt far wing overshoots
   because the real O2 A-band wing is sub-Lorentzian and line mixing suppresses
   the troughs. So CIA and wider wings are **not additive corrections**: adding
   CIA on top of full wings moves further from ABSCO, not closer. Either keep
   plain Voigt + a stated cutoff and no CIA (Phase 1, consistent with the
   ensemble), or do line mixing + χ-factor + CIA together (Phase 2).
5. **Solar spectrum role** — **[AGREED]** The **CU composite solar spectrum is
   the incident source** F₀(λ). MCARaTS runs with `Src_flx = 1.0` (unit source);
   MC transport is linear in incident flux, so the raw per-g radiance R_raw(λ) is
   the per-unit-incident response, and the solar spectrum is folded in post
   (exactly as er3t's `read_radiance_mca_out` uses `coef['solar']`):
   - **Absolute radiance** I(λ) = R_raw(λ)·F₀(λ), F₀ = `solar_cu.interp(λ_air)`
     [W m⁻² nm⁻¹] ⇒ I in [W m⁻² nm⁻¹ sr⁻¹].
   - **Reflectance** ρ(λ) = π·I(λ)/(μ₀·F₀(λ)) = π·R_raw(λ)/μ₀, μ₀ = cos(SZA).
     F₀ cancels ⇒ reflectance is F₀-independent (unchanged by this choice); the
     solar spectrum only adds the absolute-radiance product. The μ₀ in the
     denominator follows from MCARaTS `Src_flx` being the beam-perpendicular
     irradiance (source `fsrc_horiz = Src_flx·fcone` = TOA horizontal-plane flux,
     `mcarSrc.f90:174-181`); re-confirmed empirically by the V7 Lambertian limit
     (ρ→A at a transparent wavelength).
   Output carries ρ(λ), I(λ) and F₀(λ) so downstream users can recompute either.

## 9. TIPS-2021 sourcing (for §7.1) — RESOLVED
Use **TIPS-2021** (`src/TIPS_2021_PYTHON/QTpy/{mol}_{iso}.QTpy`), the partition
sums produced for and consistent with HITRAN2020 (Gamache et al. 2021). Rationale:
S(296) in HITRAN2020 was normalised with TIPS-2021's Q(296), so the Q(296)/Q(T)
ratio must use the same set. TIPS-2024 kept only as a sensitivity check.
QTpy files are per-isotopologue dicts keyed by integer T (K); read Q(296) exactly
and linearly interpolate Q(T). Verified: O2(66) Q(296)=215.73, H2O(161)=174.58.

## 8. Out of scope (CLAUDE.md §4.5)

No clouds, aerosols, polarization, instrument convolution, or real-data
comparison until Phase 1 is validated and signed off.

### 8a. Deferred to a future experiment (not Phase 1)

**Cell-averaged vs point-sampled output.** Phase 1 delivers **point samples**:
τ(λ) is the spectrum evaluated exactly at each grid wavelength, with no
integration over the ±0.0005 nm cell. Verified by construction — the same
wavelength computed on 0.001 / 0.0002 / 0.00005 nm grids returns column OD
0.985486109 / 0.985486111 / 0.985486111, i.e. identical to float noise, so the
result carries no dependence on cell width. The only ± window in the code is the
±50 cm⁻¹ (±2.9 nm at 763 nm) wing cutoff, which selects *which lines* contribute
to a point — a different thing entirely.

A future experiment could deliver cell **averages** instead, ⟨τ⟩ or (better)
⟨ρ⟩ over each 0.001 nm interval. Measured difference vs point sampling
(`src/eval_grid_sampling.py`, §5): column OD p95 1.2% / max 3.5%, reflectance
p95 0.34% / max 6.9% locally on steep line flanks, while the window-mean
reflectance is unchanged (+0.0001%). Note ⟨ρ(τ)⟩ ≠ ρ(⟨τ⟩), so averaging the
optical depth is *not* equivalent to averaging the reflectance; only the latter
is what a bin-averaging model reports. Cost: sub-sampling 5–10× raises the LBL
stage from ~10 s to ~1–2 min per band and leaves the RT stage unchanged, since
the RT would still run on the 0.001 nm points.

Trigger: adopt only if the participants confirm the ensemble reports bin
averages. Until then point sampling is the stated convention
(`metadata/spectral_sampling`).

---

## 10. Step 2 plan — RT driver `sim_o2band.py` (REVIEW BEFORE IMPLEMENTING)

Grounding facts established from er3t source:
- Clear-sky 1D needs **only** `atm_1ds` — no cloud, no `mca_atm_3d`, no `mca_sca`.
  With no 3D atm, MCARaTS uses a single column (Nx=Ny=1). Uniform Lambertian
  albedo is passed as the scalar `surface_albedo`. Incident `Src_flx = 1.0`.
- `mca_out_ng` weight-*sums* g-points into one radiance ⇒ for a spectrum we read
  per-g raw outputs individually (step 2c).
- `mca_atm_1d` takes Rayleigh from a single scalar `abs.wvl` ⇒ chunk finely.

Absorption is geometry/surface independent ⇒ compute once per band, reuse across
all (SZA, albedo).

### Steps
- **2a Config + physics inputs.** Frozen Phase-1 config object (bands, dwvl=0.001,
  R=20000, ncut=3, SZA {0,30,60}, albedo {0.0,0.1}, z_top, file paths, constants),
  plus an optional **`wvl_range`** (air nm), default None = full band, for quick
  testing / prototype sub-windows. Build `afgl_atmosphere`, `mca_atm_lbl`, and
  per-band `o2band_absorption` once; cache absorption to disk (pickle) keyed by a
  config hash.
- **2a′ Sub-range selection.** The 0.001 nm grid is anchored to the band lattice
  (wl = band_start + k·0.001); `wvl_range` selects the subset of indices within
  [a,b] (None → all). It is a *selection on the canonical grid*, not a separate
  grid, so a test sub-range's points are an exact subset of the full-band points
  and its chunk files compose with the full run (nothing recomputed later). The
  actual `wvl_range` is recorded in metadata; a partial band is flagged as a
  documented deviation from the §3 full-band deliverable.
- **2b Chunking (operational only).** Rayleigh is made exact per 0.001 nm by
  `set_per_g_rayleigh` (see 2d), so chunking is **not** a physics constraint —
  it is purely a batch-size choice (g-points per `mcarats_ng` call: process/file
  count, memory). Each chunk → `mca_abs_lbl(absb, idx_chunk)`; Nchunk chosen for
  efficiency (could be the whole band).
- **2c Per-g output reader (`util/mca_out_lbl.py`).** After each `mcarats_ng`
  run, read `fnames_out[ir][ig]` per g (raw radiance R_raw), average over Nrun,
  fold in the per-g solar irradiance F₀(λ) to get absolute radiance
  I(λ)=R_raw·F₀, and normalise to reflectance ρ(λ)=π·I/(μ₀·F₀)=π·R_raw/μ₀.
  Assemble ρ(λ), I(λ), F₀(λ) vs wavelength; return the MC standard error per
  wavelength (from the Nrun spread) for both radiance and reflectance.
- **2d RT invocation.** Per (band, SZA, albedo, chunk): build
  `atm1d = mca_atm_1d(atm, abs_chunk)`, then `set_per_g_rayleigh(atm1d, ...)` to
  make each g-point's Rayleigh exact at its own wavelength. Then
  `mcarats_ng(atm_1ds=[atm1d], target='radiance',
  solar_zenith_angle=SZA, solar_azimuth_angle=0, sensor_zenith_angle=0,
  sensor_azimuth_angle=0, surface_albedo=alb, Ng=chunk, weights=ones,
  photons=P, Nrun=N, solver='IPA')`. IPA = plane-parallel per column, which is
  exact for this horizontally homogeneous clear-sky scene (no 3D transport to
  resolve) and faster / lower-noise than the 3D solver.
- **2e Assemble outputs.** Per (band, SZA, albedo): reflectance ρ(λ) on the
  0.001 nm air grid. Per band (geometry-independent): O2 OT, H2O OT and Rayleigh
  OT — layer-resolved and column — on the same grid (from the absorption object +
  `cal_rayleigh_od`). Write one HDF5 per band with full metadata (§4.4): every
  setting used, TIPS-2021, air convention, constants, photons/Nrun, input file
  identities, git commit, and any deviations from §3.
- **2f Validation (closes V7 + integration).** albedo 0 + saturated core → ρ≈0;
  ρ monotonically increasing in albedo; at a window (near-transparent) wavelength,
  ρ vs an independent Rayleigh+Lambertian estimate within a few %; MC noise
  (ρ_std/ρ) below a set threshold. Add to the test table.
- **2g Scale/runtime.** Full grid = 2 bands × 3 SZA × 2 albedo × 15001 λ × Nrun
  monochromatic solves. Prototype on a narrow sub-window first, benchmark, then
  decide photon count and whether cluster/batch (SBATCH header already present) is
  needed for the full run. Runtime levers: absorption computed once/band (cached);
  Ng spread over Ncpu per chunk; independent work units → cluster job array;
  prototype-to-size photons against a noise target; optional analytic Lambertian
  surface coupling ρ(A)=ρ₀+T↓T↑A/(1−sA) to derive both albedos from atmosphere-only
  runs.
- **2h Robustness & resume (checkpointing).** Work unit = (band, SZA, albedo,
  chunk) with a deterministic output path
  `out/{band}/sza{SZA}_alb{ALB}/chunk_{i0}_{i1}.h5`, where `{i0}_{i1}` are
  band-relative grid indices (lattice-anchored) so `wvl_range` sub-runs and the
  full run share chunk files. Driver is **skip-if-done**:
  a valid existing chunk file is loaded and skipped, so a failed run resumes by
  simply rerunning. **Atomic writes** (tmp file + rename) prevent corrupt "done"
  files; `valid()` (expected point count, no NaNs) guards acceptance. Chunk size
  set for checkpoint granularity (~500–1000 λ ⇒ ≤1 chunk lost on failure). Units
  are idempotent ⇒ safe to parallelise across cluster array tasks. Final
  `assemble()` stitches chunk files into the per-band output HDF5 + metadata.

### Open decisions — LOCKED
1. **Chunk size** — *operational only* (Rayleigh is exact per-0.001nm via
   `set_per_g_rayleigh`). **[AGREED]** ~1000 λ per chunk for checkpoint
   granularity (≤1 chunk lost on failure); no accuracy trade-off.
2. **Photons & Nrun** — **[AGREED]** P=1e6, Nrun=3 (Nrun gives per-λ MC standard
   error). Re-checked against the 2f noise threshold on the prototype before the
   full run; escalate to P=1e7/Nrun=5 only if saturated cores are too noisy.
3. **Reflectance definition** — **[AGREED]** ρ = π I / (μ0 F0), F0 from
   Src_flx=1; exact factor pinned by the 2f limiting cases.
4. **Runtime approach** — **[AGREED]** prototype-then-scale: build/benchmark on a
   narrow `wvl_range` locally, then decide cluster vs local for the full grid.
5. **Output schema** — **[AGREED]** one HDF5 **per band**
   (`o2a.h5`/`o2b.h5`: ρ[SZA,alb,λ] + O2/H2O/Rayleigh OT[λ,layer] + metadata)
   **and** one **merged** file with both bands as groups. Per-band files are the
   assembly unit; the merged file is produced by `assemble()` at the end.
6. **z_top** — **[AGREED]** optional config parameter, **default 70 km**; the
   signed-off production run uses **120 km (full AFGL)**. Recorded in metadata.
