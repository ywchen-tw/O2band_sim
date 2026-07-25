#!/usr/bin/env python
"""
Regression guard for the line-wing cutoff, the selection margin, and the
air/vacuum grid convention.

Why this exists
---------------
The 2026-07-25 intercomparison feedback exposed two defects the existing suite
could not see:

  (a) the production wing cutoff (N*nu0/R ~ 2 cm-1) discarded ~0.055 of column
      optical depth in the A-band micro-windows, leaving between-line
      reflectance ~10% too high;
  (b) the HAPI cross-check (eval_hapi_local.py) could not detect it, because it
      builds its comparison table from the SAME band +/-5 cm-1 line subset --
      both codes were truncated the same way.

So the reference here is a *brute-force* Voigt sum with NO cutoff over a wide
line set, computed independently of o2band_absorption's windowed summation.
That is the only reference that can bound the truncation error.

    python tests/test_wings.py
"""

import os
import sys
import traceback
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
sys.path.insert(0, os.path.join(_REPO, 'src'))
DATA = os.environ.get('O2BAND_DATA_DIR', os.path.join(_REPO, 'data'))

from util.atmosphere import afgl_atmosphere
from util.tips import tips2021
from util.optics import air_to_vac_nm, vac_to_air_nm
from util.absorption import (hitran_lines, o2band_absorption, voigt_profile,
                             doppler_hwhm, line_cutoff_cm, required_margin_cm,
                             C2, T_REF, P_REF_HPA, MOL_O2, DEFAULT_CUTOFF_CM)

# a 2 nm A-band window: dense lines, deep micro-windows between them
WIN = (763.0, 765.0)          # vacuum nm
DWVL = 0.005                  # coarse enough to keep the brute force quick
Z_TOP = 120.0

_RESULTS = []


def check(name):
    def deco(fn):
        def run(*a, **kw):
            try:
                ok, detail = fn(*a, **kw)
            except Exception:
                traceback.print_exc()
                ok, detail = False, 'raised'
            _RESULTS.append((name, ok, detail))
            print('%-4s %-46s %s' % ('PASS' if ok else 'FAIL', name, detail))
            return ok
        return run
    return deco


# ---------------------------------------------------------------------------- #
def brute_force_column_od(atm, tips, nu_grid, wl_pad=(755.0, 775.0),
                          margin_cm=300.0):
    """Column O2 optical depth by summing EVERY line within +/-300 cm-1 at every
    grid point, with no wing cutoff at all.  Deliberately independent of
    o2band_absorption's windowed/sorted summation path."""
    lines = hitran_lines(os.path.join(DATA, 'hitran2020_lines.txt'),
                         wl_range=wl_pad, grid='vac', margin_cm=margin_cm)
    idx = np.where(lines.subset(MOL_O2))[0]
    nu0, S296 = lines.nu[idx], lines.S[idx]
    gair, gself = lines.gamma_air[idx], lines.gamma_self[idx]
    Epp, nair, dair = lines.Epp[idx], lines.n_air[idx], lines.delta_air[idx]
    mass, iso = lines.mass[idx], lines.iso[idx]

    od = np.zeros(nu_grid.size)
    lay = atm.lay
    for iz in range(lay['z'].size):
        T = lay['temperature'][iz]
        p_atm = lay['p'][iz] / P_REF_HPA
        p_self = p_atm * lay['o2_vmr'][iz]
        qr = np.array([tips.ratio(MOL_O2, int(i), T) for i in iso])
        S_T = (S296 * qr
               * np.exp(-C2 * Epp / T) / np.exp(-C2 * Epp / T_REF)
               * (1 - np.exp(-C2 * nu0 / T)) / (1 - np.exp(-C2 * nu0 / T_REF)))
        nu_c = nu0 + dair * p_atm
        gL = (T_REF / T) ** nair * (gair * (p_atm - p_self) + gself * p_self)
        aD = doppler_hwhm(nu_c, T, mass)
        col = lay['o2'][iz]
        for k in range(nu_c.size):
            od += col * S_T[k] * voigt_profile(nu_grid, nu_c[k], aD[k], gL[k])
    return od


def build(atm, tips, cutoff_cm, grid='vac', **kw):
    margin = required_margin_cm(WIN, grid=grid, cutoff_cm=cutoff_cm)
    lines = hitran_lines(os.path.join(DATA, 'hitran2020_lines.txt'),
                         wl_range=WIN, grid=grid, margin_cm=margin)
    return o2band_absorption(atm, lines, WIN, dwvl=DWVL, cutoff_cm=cutoff_cm,
                             grid=grid, include_h2o=False, gases=('o2',),
                             tips=tips, **kw)


# ---------------------------------------------------------------------------- #
@check('W1 default cutoff ~ untruncated in micro-windows')
def w1_cutoff_accuracy(atm, tips, ref, absb_def):
    """The delivered cutoff must not lose meaningful OD where the spectrum is
    most transparent -- micro-windows set the continuum reflectance."""
    od = absb_def.od['o2'].sum(axis=0)
    k = np.argsort(ref)[:40]                      # 40 most transparent points
    rel = abs(od[k].mean() - ref[k].mean()) / ref[k].mean()
    return rel < 0.05, '%g cm-1: micro-window OD %.5f vs %.5f (%.2f%% low)' % (
        DEFAULT_CUTOFF_CM, od[k].mean(), ref[k].mean(), 100 * rel)


@check('W2 legacy ~2 cm-1 cutoff is detectably deficient')
def w2_legacy_is_bad(atm, tips, ref):
    """Guard the guard: if this ever 'passes', the brute-force reference has
    stopped being sensitive to truncation and W1 means nothing."""
    absb = build(atm, tips, cutoff_cm=None)       # legacy 3*nu0/R
    od = absb.od['o2'].sum(axis=0)
    k = np.argsort(ref)[:40]
    rel = (ref[k].mean() - od[k].mean()) / ref[k].mean()
    return rel > 0.5, 'legacy loses %.0f%% of micro-window OD (%.5f vs %.5f)' % (
        100 * rel, od[k].mean(), ref[k].mean())


@check('W3 column OD increases monotonically with cutoff')
def w3_monotonic(atm, tips):
    ods = []
    for c in (5.0, 10.0, 25.0, 50.0):
        ods.append(build(atm, tips, cutoff_cm=c).od['o2'].sum(axis=0).mean())
    ok = all(b >= a for a, b in zip(ods, ods[1:]))
    return ok, 'mean column OD %s' % ' -> '.join('%.4f' % o for o in ods)


@check('W4 margin < cutoff is rejected, not silently wrong')
def w4_margin_guard(atm, tips):
    """A margin narrower than the cutoff drops the wings of out-of-window lines
    and mimics a too-short cutoff; it must raise rather than run."""
    lines = hitran_lines(os.path.join(DATA, 'hitran2020_lines.txt'),
                         wl_range=WIN, grid='vac', margin_cm=5.0)
    try:
        o2band_absorption(atm, lines, WIN, dwvl=DWVL, cutoff_cm=25.0,
                          grid='vac', include_h2o=False, gases=('o2',), tips=tips)
    except ValueError as e:
        return 'margin' in str(e), 'raised ValueError: %s' % str(e)[:60]
    return False, 'no error raised for margin 5 cm-1 < cutoff 25 cm-1'


@check('W5 vac and air grids describe the same spectrum')
def w5_grid_equivalence(atm, tips, absb_def):
    """Same physics, two axes: the air-grid OD interpolated onto the vacuum
    grid's air wavelengths must reproduce the vacuum-grid OD."""
    win_air = tuple(vac_to_air_nm(np.array(WIN)))
    absb_air = build(atm, tips, cutoff_cm=DEFAULT_CUTOFF_CM, grid='air')
    od_air = absb_air.od['o2'].sum(axis=0)
    od_vac = absb_def.od['o2'].sum(axis=0)
    m = ((absb_def.wvl_air > absb_air.wvl_air[2]) &
         (absb_def.wvl_air < absb_air.wvl_air[-3]))
    interp = np.interp(absb_def.wvl_air[m], absb_air.wvl_air, od_air)
    rel = np.abs(interp - od_vac[m]) / np.maximum(od_vac[m], 1e-6)
    # tolerance is set by DWVL interpolation error on a steep line, not physics
    return np.median(rel) < 0.02, 'air window %.4f-%.4f nm, median rel diff %.3f%%' % (
        win_air[0], win_air[1], 100 * np.median(rel))


@check('W6 grid conventions differ by the Edlen shift')
def w6_shift(atm, tips, absb_def):
    d = absb_def.wvl_vac - absb_def.wvl_air
    return (0.15 < d.min() and d.max() < 0.25 and
            np.allclose(air_to_vac_nm(absb_def.wvl_air), absb_def.wvl_vac, atol=1e-9)), \
        'vac - air = %.4f-%.4f nm over the window' % (d.min(), d.max())


@check('W7 CIA is separable and off by default')
def w7_cia(atm, tips, absb_def):
    fcia = os.path.join(DATA, 'O2-O2_2024.cia')
    if not os.path.isfile(fcia):
        return True, 'SKIP: no O2-O2_2024.cia in %s' % DATA
    from util.cia import hitran_cia
    absb = build(atm, tips, cutoff_cm=DEFAULT_CUTOFF_CM, cia=hitran_cia(fcia))
    off = absb_def.od_cia.sum()
    on = absb.od_cia.sum(axis=0)
    lines_same = np.allclose(absb.od['o2'], absb_def.od['o2'])
    return (off == 0.0 and on.min() > 0.0 and lines_same and
            np.allclose(absb.od_total, absb.od['o2'] + absb.od_cia)), \
        'default 0; with CIA column OT %.5f-%.5f, line OT unchanged=%s' % (
            on.min(), on.max(), lines_same)


# ---------------------------------------------------------------------------- #
if __name__ == '__main__':

    print('Wing-cutoff / grid-convention regression suite')
    print('window %.1f-%.1f nm vacuum @ %.3f nm, z_top %.0f km, %s\n'
          % (WIN[0], WIN[1], DWVL, Z_TOP, DATA))

    atm = afgl_atmosphere(os.path.join(DATA, 'afglms.dat'), z_top=Z_TOP)
    tips = tips2021()

    absb_def = build(atm, tips, cutoff_cm=DEFAULT_CUTOFF_CM)
    print('building brute-force no-cutoff reference (this is the slow part)...')
    ref = brute_force_column_od(atm, tips, absb_def.nu_vac)
    print('reference column OD: min %.5f  median %.4f\n' % (ref.min(), np.median(ref)))

    w1_cutoff_accuracy(atm, tips, ref, absb_def)
    w2_legacy_is_bad(atm, tips, ref)
    w3_monotonic(atm, tips)
    w4_margin_guard(atm, tips)
    w5_grid_equivalence(atm, tips, absb_def)
    w6_shift(atm, tips, absb_def)
    w7_cia(atm, tips, absb_def)

    nfail = sum(1 for _, ok, _ in _RESULTS if not ok)
    print('\n%d/%d checks passed' % (len(_RESULTS) - nfail, len(_RESULTS)))
    sys.exit(1 if nfail else 0)
