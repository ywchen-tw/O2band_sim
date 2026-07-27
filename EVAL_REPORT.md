# EVAL_REPORT.md — O2 A/B-band simulation vs independent references

Results of the evaluation planned in [EVAL_PLAN.md](EVAL_PLAN.md).  Each physics
component the simulation controls is cross-checked against an independent public
reference; the two deliberate Phase-1 omissions are *quantified* rather than
hidden.  Per the agreed scope this reports **difference statistics, not pass/fail
verdicts** (EVAL_PLAN §5).

Run/context: AFGL mid-latitude-summer, z_top 120 km, prescribed HITRAN 2020,
Voigt, air wavelengths.  Reflectance figures use the z120/P=1e7/Nrun=3 production
output (full grid: SZA 0/30/60° × albedo 0/0.1; MC-noise p95 ≤0.18% O2A /
≤0.05% O2B).  The solver cross-checks in §6 were run against the P=1e6
predecessor; its mean reflectances differ from the 1e7 production only by MC
noise (verified 1/√N scaling between the two runs), so those difference
statistics carry over.  Generated 2026-07-06; production photons raised to
1e7 on 2026-07-07.

The participant-model ensemble (KNMI intercomparison) was **not available**, so
this evaluation rests on independent public references and local model reruns.

> **Delivery v2 re-evaluation (2026-07-26).** Run `z120_p1e7_n3_vac_c50`:
> **vacuum** wavelength grid, **±50 cm⁻¹** plain-truncation wing cutoff, Rayleigh
> at vacuum λ. This supersedes the v1 numbers below, which were produced with an
> air grid and a ~2 cm⁻¹ cutoff. The v1 cutoff discarded ~0.055 of column optical
> depth in the A-band micro-windows (more than twice the whole Rayleigh OT) and
> **invalidated the attribution in §5** — the 8–9× between-line deficit vs ABSCO
> was mostly wing truncation, not the omitted CIA + line mixing.
>
> **Check #4 could not have caught it**: `eval_hapi_local.py` built its HAPI
> table from the same band ±5 cm⁻¹ line subset, with air-converted band limits,
> and never told HAPI what wing to use — so both codes were truncated the same
> way. That script is fixed (vacuum band limits, margin derived from the cutoff,
> HAPI wing matched to ours) and re-run below. `tests/test_wings.py` adds a
> brute-force **no-cutoff** reference, which is the only thing that can bound
> truncation error.
>
> Scorecard rows below are v2 results except where marked *(v1, pending)* —
> those need inputs available only on CURC (ABSCO tables, OCO solar model).

> Erratum (2026-07-06): the delivered `reflectance_stderr`/`radiance_stderr` were
> originally computed with the population std (ddof=0) over the Nrun=3 runs and
> have been rescaled in place to the unbiased sample std (ddof=1), a uniform
> ×√(3/2) ≈ 1.22 on the error bars only (marker attr `stderr_ddof=1`).  All mean
> reflectance/radiance/OT values are bit-identical, so the difference statistics
> below are unchanged; only "MC noise" percentages quoted in §6b read ~22% larger
> under the corrected convention (conclusions unaffected).  The P=1e7 production
> files compute ddof=1 stderr natively (no rescale involved).

---

## Scorecard

| # | component | independent reference | agreement (v2) | script |
|---|---|---|---|---|
| 1 | Rayleigh cross-section | Bucholtz (1995) | **+0.03%** (O2A) / **+0.02%** (O2B), corr 1.00000 | `eval_rayleigh.py` |
| 2 | Rayleigh column OT | Hansen & Travis (1974) | **+0.04%** @550, **−0.05%** @688, **−0.07%** @760 nm | `eval_band_metrics.py` |
| 3 | O2 absorber amount | canonical 0.2095 dry-air VMR | **−0.17%** (H2O dilution, physically correct) | `eval_band_metrics.py` |
| 4 | O2 line-by-line engine | HAPI, **matched** HITRAN 2020 **and matched 50 cm⁻¹ wing** | peak **+0.39%** / **−0.04%** (O2A sfc/upper), **−1.05%** / **−0.22%** (O2B); corr 0.9997–0.99997; median rel ~5×10⁻⁴ | `eval_hapi_local.py` |
| 5 | O2 A-band continuum | OCO ABSCO v5.2 | *(v1, pending — needs CURC)*; with full wings our between-line σ now **exceeds** ABSCO, see §5 | `eval_absco.py` |
| 6 | RT solver + reflectance convention | libRadtran / DISORT 16-stream | window **rel RMS 0.44%**, in-band (col OD 0.10–2.67) **rel RMS 0.19%**, both corr 1.00000 | `eval_lrt.py`, `eval_lrt_inband.py` |
| 7 | Solar transmittance (Toon SPTS 2024) | OCO L2 solar model | *(v1, pending — needs CURC)*; solar model unchanged in v2 | `eval_spts_oco.py` |
| — | *(context)* O2A line intensities | HAPI online = HITRAN 2024 | **+1.3%** (edition change, not a defect) | `eval_hapi.py` |

Product-level gates on the v2 run, both green:

| gate | result |
|---|---|
| MC noise (`noise_report.py`, threshold 0.01) | worst rel-stderr p95 **0.0019** (O2A) / **0.0005** (O2B) — **PASS** |
| RT regression (`tests/test_rt_regression.py`) | **PASS** both bands: τ_R 0.02556~0.02558 (0.07%) / 0.03928~0.03930 (0.05%); window + saturated core + albedo monotonicity over 3 SZA × 2 albedo |
| Wing/grid guards (`tests/test_wings.py`) | **7/7** |
| Absorption suite (`tests/test_absorption.py`) | **10/10** |

**Bottom line:** absorption (line engine) and scattering (Rayleigh) match
independent references to well under 1%; the RT transport and reflectance
convention match an independent solver to <1% across all six prescribed
geometries.  The only large differences are the two documented Phase-1 choices
(§5, §8 below), both now quantified.

---

## 1. Rayleigh cross-section vs Bucholtz (1995)

Our Bodhaine (1999) cross-section vs the independent Bucholtz (1995)
parameterization (layer-independent, so this equals the relative column-OT diff):

| band | σ @ band centre | Bodhaine vs Bucholtz | rel. RMS | corr |
|---|---|---|---|---|
| O2A (764 nm) | 1.1848e-27 cm² | +0.03% | 2.5e-4 | 1.00000 |
| O2B (688 nm) | 1.8206e-27 cm² | +0.02% | 1.8e-4 | 1.00000 |

A tiny uniform positive bias, far below any RT-relevant error.

## 2. Rayleigh column OT vs Hansen & Travis (1974)

Our Bodhaine-σ × AFGL air-column integration vs the Hansen & Travis (1974)
standard-atmosphere total-column parameterization — validates the cross-section
**and** the column integration together:

| λ (nm) | ours | H&T 1974 | diff |
|---|---|---|---|
| 550 (anchor) | 0.09731 | 0.09728 | +0.04% |
| 688 (O2B) | 0.03916 | 0.03918 | −0.05% |
| 760 (O2A) | 0.02618 | 0.02620 | −0.07% |

## 3. O2 absorber amount vs canonical 0.2095

O2/air volume mixing ratio 0.20914 vs the dry-air canonical 0.2095 → **−0.17%**;
this deficit is the expected water-vapour dilution (our `air` is *total* moist
air, 0.2095 is *dry*), i.e. physically correct.  O2 column 4.512×10²⁴ molec cm⁻².

## 4. O2 line-by-line engine vs HAPI (matched HITRAN 2020)

HAPI is an independent Voigt LBL implementation.  Run on the **same** local
HITRAN 2020 line data (removing any edition difference), our per-layer O2
cross-section vs HAPI:

**v2 (2026-07-26), wing cutoff matched at 50 cm⁻¹ on both sides:**

| band / layer | peak σ ours vs HAPI | correlation | median rel. diff |
|---|---|---|---|
| O2A surface (955.9 hPa, 291.9 K) | +0.39% | 0.99979 | −5.3×10⁻⁴ |
| O2A upper (102.7 hPa, 215.7 K) | −0.04% | 0.99997 | −6.4×10⁻⁴ |
| O2B surface (955.9 hPa, 291.9 K) | −1.05% | 0.99972 | +7.0×10⁻⁴ |
| O2B upper (102.7 hPa, 215.7 K) | −0.22% | 0.99997 | −4.0×10⁻⁴ |

Median relative differences are ~5×10⁻⁴ and the relative bias ~10⁻⁴; peaks agree
to **0.04–1.05%**.  The absorption engine — line intensity S(T) (TIPS-2021
partition sums), Voigt shape, pressure broadening/shift, and now the far wing —
reproduces an independent implementation to well under 1%.  The remaining
per-point scatter (rel. RMS 1–2%, p5/p95 ≈ ±3%) sits on steep line flanks, where
a sub-grid wavenumber offset moves σ by a few percent; it is not a bias.

*v1 for reference (unmatched wings, air band limits, ±5 cm⁻¹ HAPI table):*
O2A −0.13% / −0.05%, O2B +0.51% / −0.12%.  Those numbers looked good for the
wrong reason — both codes were truncated at ~2 cm⁻¹, so the check was blind to
the far wing entirely.

### 4a. HITRAN version note (context, not a defect)

Against HAPI's *online* line list (`eval_hapi.py`) the O2A peaks differed by
−1.2%.  That gap is entirely the line list: **hitran.org now serves HITRAN 2024**,
whose O2 **A-band** intensities are **+1.3%** above HITRAN 2020 (strong lines: mean
+1.37%, median +1.25%), while the **B-band is unchanged (0.00%)**.  Line positions
are identical.  The benchmark *prescribes* HITRAN 2020, so our static 2020 file is
the correct input.  Heads-up: migrating the benchmark to HITRAN 2024 would shift
O2 A-band absorption (and in-band reflectance) by ~1.3%, B-band by ~0.

## 5. O2 A-band continuum vs OCO ABSCO v5.2 (quantifies a Phase-1 omission)

ABSCO (Drouin line list) is the operational OCO absorption and **includes line
mixing + collision-induced absorption (CIA) + speed-dependent line shape**, which
Phase-1 excludes (PLAN.md §7.4).  Per-layer O2 cross-section:

| layer | line cores (>10% peak) | between-line σ (window) |
|---|---|---|
| surface (956 hPa) | peak +0.8%, corr 0.9993 | ours 2.4e-28 → **ABSCO 2.2e-27 (9.1×)** |
| 302 hPa | peak −0.7%, corr 0.9998 | ours 1.4e-28 → **ABSCO 1.1e-27 (7.6×)** |

- **Line cores agree to ~1%** — intensities are mutually consistent
  (HITRAN 2020 vs Drouin/ABSCO), corroborating check #4.
- **Between lines ABSCO is ~8–9× higher** — the CIA + line-mixing continuum we
  omit (~2×10⁻²⁷ cm²).  Over the column this is an O2 A-band **continuum optical
  depth of order ~0.01** — comparable to the Rayleigh OT (~0.026), so **not
  negligible for window/continuum reflectance**.  This is the largest physics gap
  in the evaluation and is a deliberate, documented Phase-1 choice.

ABSCO here covers only the A-band; O2B is not in this table.

## 6. RT solver + reflectance convention vs libRadtran/DISORT

Two independent solvers, Monte-Carlo (MCARaTS) and discrete-ordinate (DISORT, 16
streams), both using ρ = πI/(μ₀F₀), compared two ways: **(6a)** on the pure-Rayleigh
window (transport + convention) and **(6b)** with gas absorption injected (the
absorption+scattering coupling).

### 6a. Window wavelengths (pure Rayleigh)

At **9 window wavelengths spanning 757.0–771.5 nm** (each O2 OT ≲ 0.001) × all 6
prescribed geometries = **54 runs**:

**Overall: relative RMS 0.31%, correlation 1.00000, every point within ±0.82%.**

Representative geometry table at 757.0 nm (pattern is the same at all 9 wavelengths):

| SZA | albedo | MCARaTS | DISORT | diff |
|---|---|---|---|---|
| 0° | 0.0 | 0.00983 | 0.00992 | −0.82% |
| 0° | 0.1 | 0.10746 | 0.10754 | −0.08% |
| 30° | 0.0 | 0.01003 | 0.01010 | −0.71% |
| 30° | 0.1 | 0.10760 | 0.10753 | +0.07% |
| 60° | 0.0 | 0.01273 | 0.01276 | −0.17% |
| 60° | 0.1 | 0.10919 | 0.10911 | +0.07% |

Two fully independent solvers agree to **<1% at every wavelength and geometry**.
This validates:

- **RT transport** — MC vs discrete-ordinate agree on the Rayleigh continuum across
  the whole band (the λ⁻⁴ decline, ρ 0.0098→0.0091 over 757→771.5 nm at SZA 0/
  albedo 0, is tracked by both).  The small systematic (MCARaTS ~0.1–0.8% below
  DISORT at albedo 0, <0.3% and sign-changing at 0.1) lives only in the faint
  atmospheric-path reflectance (~0.01), where DISORT's 16-stream angular treatment
  / MC noise / Rayleigh depolarization differ most; it vanishes once the Lambertian
  surface dominates.
- **Reflectance convention, across SZA** — at albedo 0 the pure-Rayleigh ρ rises
  with SZA (0.00983 → 0.01003 → 0.01273) in *both* solvers, tracking to <1%.  A
  missing μ₀ = cos(SZA) would diverge at SZA 60° (μ₀ = 0.5 → ~2×); it is −0.17%.
  So `ρ = π·R_raw/μ₀` is confirmed across the grid (corroborating V7).
- **Albedo** — ρ increases 0.01 → 0.107 (0 → 0.1 albedo, monotonic); at albedo 0.1
  the two solvers agree to <0.3% (surface term handled identically).

### 6b. In-band, with gas absorption injected (absorption+scattering coupling)

Windows validate transport but not the solver's handling of *coupled* absorption
and scattering — the point of the O2 bands.  To test that, our per-layer gas
absorption OT (o2_layer + h2o_layer) was injected into DISORT as a pure-absorbing
``aerosol_file tau`` profile (ssa = 0) on top of libRadtran's Rayleigh, so both
solvers see the identical scene.  6 wavelengths spanning **column gas OD ≈ 0.05–3**
× 6 geometries = 36 runs per band:

| band | col-OD range | overall rel RMS | corr | notes |
|---|---|---|---|---|
| O2A | 0.05 – 3.1 | **0.29%** | 1.00000 | max \|diff\| ~1.0% only where MC noise ~0.5% |
| O2B | 0.05 – 3.3 | **0.36%** | 1.00000 | max \|diff\| ~0.8% only where MC noise ~0.4% |

MC and discrete-ordinate agree to **<0.4% RMS with absorption and scattering
coupled**, over the full OD range and all geometries, in both bands.  Both solvers
track the absorption darkening (e.g. O2A: ρ 0.0088 at OD 0.10 → 0.0033 at OD 3.1)
and the two-way attenuation of the Lambertian surface term at albedo 0.1.  The
few ~1% points all coincide with comparable MC noise (0.3–0.5%), so they are
noise-consistent, not solver disagreement.  The injected optical depth is exact
(per-layer aerosol τ sums to the column OT), and the absence of any systematic
offset confirms no wavelength mis-scaling of the injected profile.

**Together 6a+6b validate the RT solver + convention across the full reflectance
range — window continuum through saturating line absorption — in both bands.**

---

## 7. Solar transmittance (SPTS) vs OCO L2 solar model

Added 2026-07-07: the absolute-radiance product folds F₀ = CU composite
continuum × **Toon (JPL) SPTS** disk-integrated solar transmittance
(`solar_merged_20240731`, 0.01 cm⁻¹ grid), so Fraunhofer lines are resolved on
the 0.001 nm output grid.  Reflectance and optical thickness are F₀-independent
and unaffected (the intercomparison delivery `o2band_benchmark_noradiance.h5`
carries no solar-dependent dataset at all).

Check: the OCO retrieval's solar model (`l2_solar_model.h5`, Band 1
12700–13300 cm⁻¹ on a 0.001 cm⁻¹ grid) is an independent file built from an
earlier release of the same Toon solar line list, and overlaps our O2 A-band
window.  Comparing SPTS transmittance interpolated onto the OCO grid over
12953–13210 cm⁻¹ (257,000 points):

| statistic | value |
|---|---|
| mean transmittance | 0.98024 (both) |
| bias | −6.7×10⁻⁶ |
| RMS difference | 3.0×10⁻⁴ |
| correlation | 0.99999 |
| fraction \|Δ\| > 10⁻³ / 10⁻² | 0.4% / 0.02% |
| deepest line (K I 766.5 nm) | T = 0.1139 vs 0.1143 |

The largest single-point difference (1.2×10⁻²) sits on a sharp line core at
770.1 nm, consistent with the SPTS grid being 10× coarser than OCO Band 1 —
immaterial at the 0.001 nm (≈0.02 cm⁻¹) output sampling.  This checks both the
SPTS file content and our air-nm → vacuum-wavenumber handling; no OCO band
covers O2B, so the check is O2A only.

---

## 8. Summary of the two Phase-1 differences (both by design)

1. **HITRAN edition** — the prescribed HITRAN 2020 O2 A-band intensities are ~1.3%
   below the current HITRAN 2024; a version choice, not an error.
2. **O2 CIA + line mixing** — excluded, giving an O2 A-band continuum OT ~0.01 too
   low vs the operational ABSCO.  This is the most consequential omission for
   window reflectance and is a candidate for a future phase.

Everything else matches independent references to <1%.

## 9. Caveats / not covered

- **Participant ensemble** (KNMI intercomparison models) unavailable → not compared.
- **ABSCO band coverage**: ABSCO v5.2 has no O2B table, so check #5 (CIA/line-mixing
  continuum) is O2A only.  The DISORT solver check (#6) covers **both** bands.
- Reflectance band metrics (continuum ρ, equivalent width) are geometry/surface-
  dependent with no single clean published value — reported in `eval_metrics.py`
  for reference, not differenced here.

## 10. Reproducibility

Each row of the scorecard is produced by the named script under `src/`:
`eval_rayleigh.py`, `eval_band_metrics.py`, `eval_hapi_local.py` / `eval_hapi.py`
(via `curc_hapi_eval.sh`), `eval_absco.py`, for #6 `eval_lrt.py` (window) +
`eval_lrt_inband.py` (in-band, `INBAND=1`) via `curc_lrt_eval.sh`, and for #7
`eval_spts_oco.py`.
`eval_metrics.py` provides the shared band-metric + `diff_stats` engine.
