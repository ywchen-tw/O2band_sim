# EVAL_PLAN.md — Evaluating the simulation against other research

Blueprint for judging whether this simulation's O2 A/B-band reflectance and
optical-thickness products agree with other research. Written in the project's
plan-before-implementing style (cf. [PLAN.md](PLAN.md)); items are **[PROPOSED]**,
**[NEEDS SIGN-OFF]**, or **[AGREED]**. Direction decisions are resolved (§9): no
participant data, **report difference statistics (not pass/fail)**, local reruns
enabled as a secondary cross-check.

> Scope note (2026-07-06): this is an *evaluation* plan. It does not change the
> physics of the simulation; it defines what we compare against, how, and the
> pass/fail criteria.

---

## 1. Context — what "other research's values" means here

Per [O2_abs_intercomparison.md](O2_abs_intercomparison.md), this simulation is
**one participant model in the KNMI-led O2 A/B-band RT model intercomparison**
(P. Wang, N. Ferlay, H. Herbin, R. Preusker, J. Wang, J. Vidot, M. Duan, P.
Stammes), presented at the *2nd Workshop on Remote Sensing in Oxygen Absorption
Bands* (KNMI, 2024). Our er3t/MCARaTS run already uses the prescribed Phase-1
settings (mid-lat summer; SZA 0/30/60; nadir; albedo 0/0.1; 680–695 & 757–772 nm;
HITRAN 2020; Voigt; 0.001 nm; outputs reflectance + O2A/O2B OT + Rayleigh OT).

**Availability (confirmed 2026-07-06):** the other participants' submitted
results are **not available** — [O2_abs_intercomparison.md](O2_abs_intercomparison.md)
(the prescription) is the only document in hand. So the participant ensemble is
*not* a reference we can use now. Evaluation therefore rests on:
1. **Published band-level metrics** and the closest **public reference spectra** —
   the primary comparison.
2. **Independent local model reruns** (HAPI, libRadtran/DISORT, ABSCO) — enabled as
   a secondary, lower-priority cross-check (Appendix A).

**Reality check (from a public-literature search, §3):** *no* public dataset
reproduces the exact prescribed config (0.001 nm, Voigt, no CIA/line-mixing, those
geometries). So comparisons are **descriptive** — we report difference *statistics*
against each reference (not pass/fail), with config differences attributed.

---

## 2. Strategy — separable, attributable comparison

The simulation outputs **O2 OT, H2O OT, Rayleigh OT, *and* reflectance
separately**, so a disagreement can be *localised* rather than merely observed:

| observation | attributed to |
|---|---|
| OT agrees, reflectance differs | RT solver or reflectance-convention difference |
| O2 / H2O OT differs | absorption implementation (partition fn, line shape, wing cutoff, air↔vac, pressure shift) |
| Rayleigh OT differs | scattering formulation (Bodhaine vs Bucholtz/…) |

Every comparison is designed to isolate one component. This decomposition is the
core scientific value of the evaluation.

---

## 3. Reference sources — **[NEEDS SIGN-OFF]**

Tiered by directness. The user has prioritised **published band metrics** and
public reference spectra; local model reruns are optional (Appendix A).

| # | reference | gives | status / notes |
|---|---|---|---|
| R1 | ~~KNMI intercomparison participants~~ | ensemble reflectance + OT | **NOT AVAILABLE** — only the prescription doc was received. Revisit if the group shares results. |
| R2 | **Richter, Emde et al. 2022, AMT 15, 1587** — MYSTIC synthetic dataset, incl. O2A **755–775 nm** | high-res O2A TOA reflectance (public) | cloud/aerosol study; use its **clear-sky Rayleigh** cases; config differs (convolve/attribute). |
| R3 | **Kopparla, Natraj et al. 2017, JQSRT 198:104** — PCA vs line-by-line **DISORT**, O2 A-band | reference TOA reflectance + a stated LBL-DISORT accuracy | validates band-level reflectance magnitude/shape. |
| R4 | **Vidot et al. 2021, JQSRT 275:107847** (H2O near O2A) + HITRAN 2020 O2 line-parameter papers (e.g. Brown & Plymate) | absorption / line-parameter cross-checks | for O2/H2O OT sanity. |
| R5 | **Published band metrics** — O2 A/B-band equivalent width, band-integrated transmittance, continuum reflectance vs SZA/albedo from the literature | scalar sanity numbers robust to resolution | **primary** robust check given resolution mismatches + no participant data. |

Action E0: fetch/catalog R2–R5, record exact configs. R1 shelved (unavailable).

---

## 4. Metrics & methodology — **[PROPOSED]**

**Registration.** All our products are on the prescribed 0.001 nm *air* grid.
Point-wise comparison is valid only where a reference shares that grid/convention;
otherwise convolve both to a common resolution or compare integrated metrics.

**Primary metrics (resolution-robust, per band × SZA × albedo):**
- **Band-integrated reflectance** ∫ρ dλ and **mean continuum reflectance** (window
  points, O2 OT ≈ 0).
- **Equivalent width** of the band and of selected strong lines.
- **Column optical thickness** (O2, H2O, Rayleigh) and A/B band ratio.

**Secondary metrics (where a matched-resolution reference exists):**
- Per-λ relative-difference spectrum; RMS, max, mean bias.
- Line-core depth and line-position registration (air↔vac, pressure shift).

**MC-noise envelope.** Use the simulation's `reflectance_stderr`
([noise_report.py](src/noise_report.py)): a reflectance difference within k·stderr
(k≈2–3) is not significant. Never report a discrepancy smaller than the noise.

**Reflectance-convention alignment.** Verify all compared models share
ρ = π·I/(μ₀·F₀) via the transparent-window analytic Rayleigh+Lambertian limit
(ρ→A). A factor-of-μ₀ or π offset is the classic intercomparison artifact and must
be ruled out before interpreting any difference.

**Attribution.** For every non-trivial difference, use §2 to assign it to
absorption, Rayleigh, RT solver, or convention — and, against operational/measured
references, to the **Phase-1 physical omissions** (§8).

---

## 5. Difference statistics (no pass/fail) — **[AGREED]**

Per user decision, the evaluation **reports statistics, not verdicts** — no
green/amber/red thresholds. For each quantity below, against each available
reference, tabulate the **difference distribution**:

- mean bias, median, RMS, max |Δ| (absolute *and* relative),
- 5th/50th/95th percentiles of the relative difference,
- spectral correlation coefficient (where a per-λ reference exists),
- the value in **MC-noise units** (Δ / reflectance_stderr) so the reader sees
  whether a difference exceeds the simulation's own noise.

Quantities compared (per band × SZA × albedo where applicable):

| # | quantity | against |
|---|---|---|
| Q1 | Rayleigh column OT | Bodhaine literature / R2 |
| Q2 | O2 & H2O column OT (A & B) | R4 / local LBL (HAPI, App. A) |
| Q3 | continuum (window) reflectance vs SZA, albedo | R3/R5 |
| Q4 | band-integrated reflectance | R3/R5 |
| Q5 | O2 A/B-band equivalent width | R5 |
| Q6 | reflectance-convention limit ρ→A (window) | analytic |
| Q7 | per-λ reflectance / OT (if a matched-config reference appears) | R2 / local reruns |

The report presents these as tables + difference plots and *describes* the
differences (with §2 attribution); it does not label them pass/fail.

---

## 6. Phases — **[PROPOSED]**

- **E0 — References & criteria.** Fetch/catalog R2–R5 (record configs), request R1,
  ratify §5 thresholds. *No code beyond a reference catalog.*
- **E1 — Optical-thickness evaluation (no RT).** C1, C2: column + per-λ OT vs
  public references / band metrics; Rayleigh vs Bodhaine literature.
- **E2 — Reflectance band-metrics evaluation.** Q3–Q6: continuum level,
  band-integrated ρ, equivalent width vs R2/R3/R5, all reported in MC-noise units;
  confirm the convention limit.
- **E3 — Participant-ensemble comparison — SHELVED** (R1 unavailable). If the group
  later shares results: per-λ + band-integrated reflectance/OT vs each model,
  ensemble mean/spread, optional convolution to TROPOMI/MetImage.
- **E4 — Report.** Difference-statistics tables (§5) + plots, each difference
  *described* and §2-attributed. No pass/fail labelling.

---

## 7. Deliverables — **[PROPOSED]**

- `src/eval/` harness (loads a band HDF5 + a reference, emits the §4 metrics),
  reusing the file/threshold patterns of [noise_report.py](src/noise_report.py).
- Difference/overlay plots per (band, SZA, albedo) and a band-metrics table.
- `EVAL_REPORT.md` with the filled §5 criteria table.

---

## 8. Known limitations / honest caveats — **[AGREED]**

- **Phase 1 excludes O2 CIA/continuum and line mixing** (PLAN.md §7.4). Against
  *operational/measured* references (ABSCO, GOSAT/OCO, measured spectra) this
  produces **expected systematic differences**, largest in the band wings and
  strong-line cores — a documented modelling choice, *not* a code error. Such
  references are for quantifying completeness, not code correctness.
- **No public spectrum matches the prescribed config exactly** → band-level
  metrics + the participant ensemble carry the evaluation; per-λ public
  comparisons require convolution/attribution.
- Real-data (measurement) comparison stays **out of scope** until a later phase
  (PLAN.md §8).

---

## 9. Decisions — **[RESOLVED 2026-07-06]**

1. **Participant data (R1)** — *not available*; only the prescription doc was
   received. E3 shelved; evaluation uses public references + band metrics.
2. **Criteria** — *none*; report **difference statistics**, not pass/fail (§5).
3. **Local cross-checks** (Appendix A) — *enabled but not priority*: set up HAPI /
   libRadtran / ABSCO as secondary references; the primary path is published band
   metrics + public reference spectra.

---

## Appendix A — optional local cross-check tools (available, not primary)

Deprioritised per the reference choice, but on this system and worth enabling if
the public references prove insufficient:
- **HAPI** (HITRAN API) — Voigt cross-sections from the *same* HITRAN 2020 lines →
  apples-to-apples check of the LBL absorption code (expect ≲0.5 %). Needs install.
- **libRadtran / uvspec (DISORT)** — independent RT solver fed the *same* per-layer
  O2+Rayleigh OT + Lambertian surface → validates transport + convention within MC
  noise. Available via the arcsix `lrt_sim` setup.
- **ABSCO** (OCO O2 A-band tables, `…/oco_retrieval/data/absco`) — operational
  absorption; use only with the §8 CIA/line-mixing caveat.
