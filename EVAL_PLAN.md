# EVAL_PLAN.md — Evaluating the simulation against other research

Blueprint for judging whether this simulation's O2 A/B-band reflectance and
optical-thickness products agree with other research. Written in the project's
plan-before-implementing style (cf. [PLAN.md](PLAN.md)); items are **[PROPOSED]**,
**[NEEDS SIGN-OFF]**, or **[AGREED]**. Nothing here is built until §5 criteria and
§3 references are signed off.

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

So "other research" has two meanings, handled separately:
1. **The other participating RT models** — the authoritative target. Their
   submissions are *not* on disk yet; obtained from the intercomparison group.
2. **Independent public references** — used to self-validate *now*, before/without
   the group ensemble, and to sanity-check band-level behaviour.

**Reality check (from a public-literature search, §3):** *no* public dataset
reproduces the exact prescribed config (0.001 nm, Voigt, no CIA/line-mixing,
those exact geometries). Therefore the evaluation leans on (a) the participant
ensemble once available, and (b) **published band-level metrics** against the
closest public datasets, with config differences explicitly accounted for.

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
| R1 | **KNMI intercomparison participants** (Wang/Stammes et al., 2024) | reflectance + OT spectra from 6+ independent RT models — the authoritative ensemble | **obtain from group**; not public. Primary target. |
| R2 | **Richter, Emde et al. 2022, AMT 15, 1587** — MYSTIC synthetic dataset, incl. O2A **755–775 nm** | high-res O2A TOA reflectance (public) | cloud/aerosol study; use its **clear-sky Rayleigh** cases; config differs (convolve/attribute). |
| R3 | **Kopparla, Natraj et al. 2017, JQSRT 198:104** — PCA vs line-by-line **DISORT**, O2 A-band | reference TOA reflectance + a stated LBL-DISORT accuracy | validates band-level reflectance magnitude/shape. |
| R4 | **Vidot et al. 2021, JQSRT 275:107847** (H2O near O2A) + HITRAN 2020 O2 line-parameter papers (e.g. Brown & Plymate) | absorption / line-parameter cross-checks | for O2/H2O OT sanity. |
| R5 | **Published band metrics** — O2 A/B-band equivalent width, band-integrated transmittance, continuum reflectance vs SZA/albedo from the literature | scalar sanity numbers robust to resolution | primary robust check given resolution mismatches. |

Action E0: fetch/catalog R2–R5, record exact configs, and request R1 from the
intercomparison coordinators.

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

## 5. Agreement-criteria table — **[NEEDS SIGN-OFF]**

Green/amber/red, mirroring PLAN.md §6. Thresholds are proposals to be ratified.

| # | check | reference | green (proposed) |
|---|---|---|---|
| C1 | Rayleigh column OT | R2/literature Bodhaine | < 1 % |
| C2 | O2 column OT (A & B) | R4 / same-config LBL | < 2 % |
| C3 | continuum (window) reflectance vs SZA, albedo | R2/R3/R5 | < 3 % (and within MC noise) |
| C4 | band-integrated reflectance | R3/R5 | < 3 % |
| C5 | O2 A-band equivalent width | R5 | < 5 % |
| C6 | reflectance convention (ρ→A window limit) | analytic | < 1 % |
| C7 | per-λ reflectance vs participant models (R1) | R1 ensemble | within ensemble spread |
| C8 | O2/Rayleigh OT vs participant models (R1) | R1 ensemble | within ensemble spread |

Amber = 1–3× green; red = worse, requiring an attributed explanation.

---

## 6. Phases — **[PROPOSED]**

- **E0 — References & criteria.** Fetch/catalog R2–R5 (record configs), request R1,
  ratify §5 thresholds. *No code beyond a reference catalog.*
- **E1 — Optical-thickness evaluation (no RT).** C1, C2: column + per-λ OT vs
  public references / band metrics; Rayleigh vs Bodhaine literature.
- **E2 — Reflectance band-metrics evaluation.** C3–C6: continuum level,
  band-integrated ρ, equivalent width vs R2/R3/R5, all against the MC-noise
  envelope; confirm the convention limit.
- **E3 — Participant-ensemble comparison** *(when R1 arrives)*. C7, C8: per-λ +
  band-integrated reflectance/OT vs each model; ensemble mean/spread; flag
  outliers; optional convolution to TROPOMI/MetImage for the workshop's channel view.
- **E4 — Report + criteria table.** Green/amber/red with an attributed explanation
  for every non-green; plots; a short summary suitable to feed back to the group.

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

## 9. Open decisions — **[NEEDS SIGN-OFF]**

1. **R1 timeline / access** — who provides the participant submissions, in what
   format, and when? Drives whether E3 leads or trails E1/E2.
2. **Criteria thresholds** (§5) — ratify or adjust the green bands.
3. **Optional local cross-checks** (Appendix A) — enable HAPI (same-HITRAN LBL,
   would decisively validate C2) and/or libRadtran/DISORT (independent RT for
   C3–C6)? Both are available on this system; deprioritised for now.

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
