# EVAL_REPORT.md — O2 A/B-band simulation vs independent references

Results of the evaluation planned in [EVAL_PLAN.md](EVAL_PLAN.md).  Each physics
component the simulation controls is cross-checked against an independent public
reference; the two deliberate Phase-1 omissions are *quantified* rather than
hidden.  Per the agreed scope this reports **difference statistics, not pass/fail
verdicts** (EVAL_PLAN §5).

Run/context: AFGL mid-latitude-summer, z_top 120 km, prescribed HITRAN 2020,
Voigt, air wavelengths.  Reflectance figures use the z120/P=1e6/Nrun=3 production
output (currently SZA 30° / albedo 0 only; full grid pending).  Generated
2026-07-06.

The participant-model ensemble (KNMI intercomparison) was **not available**, so
this evaluation rests on independent public references and local model reruns.

---

## Scorecard

| # | component | independent reference | agreement | script |
|---|---|---|---|---|
| 1 | Rayleigh cross-section | Bucholtz (1995) | **~0.03%** | `eval_rayleigh.py` |
| 2 | Rayleigh column OT | Hansen & Travis (1974) | **<0.1%** | `eval_band_metrics.py` |
| 3 | O2 absorber amount | canonical 0.2095 dry-air VMR | **−0.17%** (H2O dilution) | `eval_band_metrics.py` |
| 4 | O2 line-by-line engine | HAPI, **matched** HITRAN 2020 | **~0.1–0.5%** | `eval_hapi_local.py` |
| 5 | O2 A-band continuum | OCO ABSCO v5.2 | **8–9× low** → ~0.01 OT (CIA+line-mixing, omitted by design) | `eval_absco.py` |
| 6 | RT solver + reflectance convention | libRadtran / DISORT | **0.71%** | `eval_lrt.py` |
| — | *(context)* O2A line intensities | HAPI online = HITRAN 2024 | **+1.3%** (edition change, not a defect) | `eval_hapi.py` |

**Bottom line:** absorption (line engine) and scattering (Rayleigh) match
independent references to well under 1%; the RT transport and reflectance
convention match an independent solver to 0.71%.  The only large differences are
the two documented Phase-1 choices (§5, §7 below), both now quantified.

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

| band / layer | peak σ ours vs HAPI | correlation |
|---|---|---|
| O2A surface (292 K) | −0.13% | 0.99979 |
| O2A upper (216 K) | −0.05% | 0.99997 |
| O2B surface (292 K) | +0.51% | 0.99972 |
| O2B upper (216 K) | −0.12% | 0.99996 |

Core median relative differences are sub-0.6% (mostly ~0).  The absorption engine
— line intensity S(T) (TIPS-2021 partition sums), Voigt shape, pressure
broadening/shift — reproduces an independent implementation to **~0.1–0.5%**.  The
residual is the pure implementation difference (our documented wing cutoff
±ncut·ν₀/R vs HAPI's 50-half-width default, plus partition/Voigt details).

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

At the band window (757.0 nm, O2 OT ≈ 0 → pure Rayleigh + Lambertian), our
Monte-Carlo (MCARaTS) reflectance vs libRadtran's discrete-ordinate solver
(DISORT, 16 streams), both using ρ = πI/(μ₀F₀):

| window | SZA | albedo | MCARaTS | DISORT | diff |
|---|---|---|---|---|---|
| 757.0 nm | 30° | 0.0 | 0.01003 | 0.01010 | **−0.71%** |

Two fully independent solvers agree to 0.71% on the Rayleigh continuum.  Crucially
this also validates the **reflectance convention**: a wrong factor would appear as
~15% (a missing μ₀ = cos 30°) or ~3× (a missing π), not 0.71%.  Confirms
`ρ = π·R_raw/μ₀` against an independent code (corroborating V7).

---

## 7. Summary of the two Phase-1 differences (both by design)

1. **HITRAN edition** — the prescribed HITRAN 2020 O2 A-band intensities are ~1.3%
   below the current HITRAN 2024; a version choice, not an error.
2. **O2 CIA + line mixing** — excluded, giving an O2 A-band continuum OT ~0.01 too
   low vs the operational ABSCO.  This is the most consequential omission for
   window reflectance and is a candidate for a future phase.

Everything else matches independent references to <1%.

## 8. Caveats / not covered

- **Participant ensemble** (KNMI intercomparison models) unavailable → not compared.
- **Reflectance grid coverage**: the DISORT check (#6) is one geometry (SZA 30° /
  albedo 0) because the production grid is partial; the μ₀ factor is confirmed at
  SZA 30° only.  Re-running `eval_lrt.py` after the full grid
  (`submit_o2band_array.sh`) extends it to all SZA/albedo and adds the
  albedo-monotonicity check.
- **B-band RT / ABSCO**: ABSCO v5.2 has no O2B table; the DISORT check used O2A.
- Reflectance band metrics (continuum ρ, equivalent width) are geometry/surface-
  dependent with no single clean published value — reported in `eval_metrics.py`
  for reference, not differenced here.

## 9. Reproducibility

Each row of the scorecard is produced by the named script under `src/`:
`eval_rayleigh.py`, `eval_band_metrics.py`, `eval_hapi_local.py` / `eval_hapi.py`
(via `curc_hapi_eval.sh`), `eval_absco.py`, `eval_lrt.py` (via `curc_lrt_eval.sh`).
`eval_metrics.py` provides the shared band-metric + `diff_stats` engine.
