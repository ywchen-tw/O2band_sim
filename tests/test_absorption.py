#!/usr/bin/env python
"""
Phase-1 validation suite for the O2-band line-by-line absorption code.

Implements the checkable acceptance criteria from PLAN.md §6 that do not require
the radiative-transfer engine (V7 reflectance limits is deferred to the RT
wiring).  Run directly:

    python tests/test_absorption.py

Prints a PASS/FAIL table and exits non-zero if any check fails, so it doubles as
a regression guard.  Each check maps to one physics rule (CLAUDE.md §4.3).
"""

import os
import sys
import traceback
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
sys.path.insert(0, os.path.join(_REPO, 'src'))
DATA = os.path.join(_REPO, 'data')

from util.atmosphere import afgl_atmosphere
from util.tips import tips2021
from util.optics import air_to_vac_nm, refractive_index_air
from util.absorption import (hitran_lines, o2band_absorption, cal_rayleigh_od,
                             voigt_profile, doppler_hwhm, BANDS)

# er3t's own Rayleigh cross-section, to confirm we feed the RT a matching value
sys.path.insert(0, '/Users/yuch8913/programming/er3t/er3t')
from er3t.util.util import mol_ext_wvl, cal_mol_ext


# ---------------------------------------------------------------------------- #
# small helpers
# ---------------------------------------------------------------------------- #
def _level_crossings(x, f, level):
    s = np.sign(f - level)
    idx = np.where(np.diff(s) != 0)[0]
    xs = []
    for i in idx:
        xs.append(x[i] + (level - f[i]) * (x[i + 1] - x[i]) / (f[i + 1] - f[i]))
    return np.array(xs)


def measured_fwhm(nu0, alpha_D, gamma_L, npts=800001):
    half = 60.0 * max(alpha_D, gamma_L)
    x = np.linspace(nu0 - half, nu0 + half, npts)
    f = voigt_profile(x, nu0, alpha_D, gamma_L)
    xs = _level_crossings(x, f, 0.5 * f.max())
    return xs[-1] - xs[0]


def olivero_longbothum_fwhm(alpha_D, gamma_L):
    """Voigt FWHM approximation (Olivero & Longbothum 1977), accurate to ~0.02%."""
    fL, fG = 2.0 * gamma_L, 2.0 * alpha_D
    return 0.5346 * fL + np.sqrt(0.2166 * fL ** 2 + fG ** 2)


def voigt_area(alpha_D, gamma_L):
    # Lorentz wings decay as 1/nu^2, so the domain must reach ~7000*gamma_L to
    # capture >99.99% of the area; the profile normalisation itself is exact.
    half = max(30.0 * alpha_D, 7000.0 * gamma_L)
    step = min(alpha_D, gamma_L) / 20.0
    x = np.arange(-half, half + step, step)
    f = voigt_profile(x, 0.0, alpha_D, gamma_L)
    return np.trapz(f, x)


# ---------------------------------------------------------------------------- #
# checks (each returns a (detail string); raise AssertionError on failure)
# ---------------------------------------------------------------------------- #
def V1_hitran_parse():
    lines = hitran_lines(os.path.join(DATA, 'hitran2020_lines.txt'))
    # first record, verified against the raw file
    assert lines.mol[0] == 1 and lines.iso[0] == 1
    assert abs(lines.nu[0] - 12658.287072) < 1e-6
    assert abs(lines.S[0] - 2.432e-29) < 1e-34
    assert abs(lines.gamma_air[0] - 0.0501) < 1e-6
    assert abs(lines.Epp[0] - 1525.136) < 1e-3
    n_o2 = int(np.sum(lines.mol == 7))
    n_h2o = int(np.sum(lines.mol == 1))
    assert n_o2 == 860, n_o2
    assert n_h2o == 33159, n_h2o
    assert np.all(np.isfinite(lines.mass)), 'unmapped isotopologue mass'
    return 'first line OK; O2=%d H2O=%d; all masses mapped' % (n_o2, n_h2o)


def V2_afgl_profile():
    atm = afgl_atmosphere(os.path.join(DATA, 'afglms.dat'), z_top=70.0)
    p_sfc = atm.lev['p'][0]
    t_sfc = atm.lev['temperature'][0]
    o2_vmr = atm.lay['o2_vmr'][0]
    o2_col = atm.lay['o2'].sum()
    air_col = atm.lay['air'].sum()
    assert abs(p_sfc - 1013.0) < 1.0, p_sfc
    assert abs(t_sfc - 294.2) < 0.5, t_sfc
    assert 0.207 < o2_vmr < 0.2095, o2_vmr      # slightly < 0.2095 (moist air)
    assert 4.0e24 < o2_col < 5.0e24, o2_col
    assert 2.0e25 < air_col < 2.3e25, air_col
    return ('p_sfc=%.1f hPa T_sfc=%.1f K O2vmr=%.4f O2col=%.3e air=%.3e'
            % (p_sfc, t_sfc, o2_vmr, o2_col, air_col))


def Vtips_partition():
    tips = tips2021()
    assert abs(tips.Q(7, 1, 296.0) - 215.7345) < 1e-3
    assert abs(tips.Q(1, 1, 296.0) - 174.5813) < 1e-3
    worst_o2 = worst_h2o = 0.0
    for T in (200.0, 230.0, 260.0, 296.0):
        r_o2 = tips.ratio(7, 1, T)
        r_h2o = tips.ratio(1, 1, T)
        worst_o2 = max(worst_o2, abs(r_o2 / (296.0 / T) ** 1.0 - 1.0))
        worst_h2o = max(worst_h2o, abs(r_h2o / (296.0 / T) ** 1.5 - 1.0))
    assert worst_o2 < 0.03, worst_o2         # linear rotor T^1
    assert worst_h2o < 0.05, worst_h2o       # asymmetric top T^1.5
    return ('Q296 O2/H2O exact; max dev vs power law: O2 %.2f%%, H2O %.2f%%'
            % (100 * worst_o2, 100 * worst_h2o))


def V_optics_airvac():
    for wl in (688.0, 765.0):
        wv = float(air_to_vac_nm(wl))
        # round-trip and expected NIR shift ~0.19-0.22 nm
        assert 0.15 < (wv - wl) < 0.25, wv - wl
        n = float(refractive_index_air(wv))
        assert abs(wl - wv / n) < 1e-6         # air = vac / n
    return 'air->vac shift 0.19-0.21 nm; round-trip < 1e-6 nm'


def V3_voigt_width_vs_pressure():
    nu0, mass, gair, nair = 13100.0, 31.9898, 0.05, 0.7
    cases = [('surface', 1.000, 294.0), ('mid', 0.100, 230.0), ('high', 0.003, 250.0)]
    fwhms = []
    for _, p_atm, T in cases:
        gamma_L = (296.0 / T) ** nair * gair * p_atm
        alpha_D = doppler_hwhm(nu0, T, mass)
        meas = measured_fwhm(nu0, alpha_D, gamma_L)
        approx = olivero_longbothum_fwhm(alpha_D, gamma_L)
        assert abs(meas / approx - 1.0) < 0.01, (p_atm, meas, approx)
        fwhms.append((p_atm, alpha_D, gamma_L, meas))
    # regime + monotonicity
    (_, aD_s, gL_s, f_s) = fwhms[0]
    (_, aD_h, gL_h, f_h) = fwhms[2]
    assert gL_s > 3 * aD_s and abs(f_s / (2 * gL_s) - 1.0) < 0.15   # Lorentz-limited
    assert gL_h < 0.1 * aD_h and abs(f_h / (2 * aD_h) - 1.0) < 0.05  # Doppler-limited
    assert fwhms[0][3] > fwhms[1][3] > fwhms[2][3]                   # decreases w/ alt
    return ('FWHM(cm-1): sfc %.4f (Lorentz), mid %.4f, high %.4f (Doppler); '
            'all within 1%% of Olivero-Longbothum'
            % (fwhms[0][3], fwhms[1][3], fwhms[2][3]))


def V4_area_normalisation():
    cases = [(0.013, 0.0002), (0.012, 0.050), (0.015, 0.015)]
    worst = 0.0
    for aD, gL in cases:
        area = voigt_area(aD, gL)
        worst = max(worst, abs(area - 1.0))
        assert abs(area - 1.0) < 1.0e-3, (aD, gL, area)
    return 'integral(Voigt)=1 to within %.3f%% (Doppler/Lorentz/mixed)' % (100 * worst)


def V5_rayleigh():
    atm = afgl_atmosphere(os.path.join(DATA, 'afglms.dat'), z_top=70.0)
    wl = np.array([688.0, 765.0])
    od = cal_rayleigh_od(atm, wl)
    col = od.sum(axis=0)
    # (a) our cross-section == er3t's Bodhaine cross-section (bit-level formula)
    our_xsec = col / atm.lay['air'].sum()
    er3t_xsec = 1.0e-28 * mol_ext_wvl(wl * 1e-3)
    assert np.allclose(our_xsec, er3t_xsec, rtol=1e-12), (our_xsec, er3t_xsec)
    # (b) column OD vs er3t pressure-based Bodhaine (<3%: density- vs p-integration)
    p_sfc, p_top = atm.lev['p'][0], atm.lev['p'][-1]
    er3t_col = cal_mol_ext(wl * 1e-3, p_sfc, p_top)
    rel = np.abs(col / er3t_col - 1.0)
    assert np.all(rel < 0.03), (col, er3t_col, rel)
    # (c) literature lambda^-4 scaling
    lit = 0.098 * (550.0 / wl) ** 4
    assert np.all(np.abs(col / lit - 1.0) < 0.06), (col, lit)
    return ('xsec == er3t (rtol<1e-12); column OD @688/765 = %.4f/%.4f, '
            'vs er3t %.4f/%.4f (<%.1f%%)'
            % (col[0], col[1], er3t_col[0], er3t_col[1], 100 * rel.max()))


def V6_o2_band_ot():
    atm = afgl_atmosphere(os.path.join(DATA, 'afglms.dat'), z_top=70.0)
    tips = tips2021()
    res = {}
    for band in ('o2a', 'o2b'):
        lines = hitran_lines(os.path.join(DATA, 'hitran2020_lines.txt'),
                             wl_range=BANDS[band])
        absb = o2band_absorption(atm, lines, band, dwvl=0.001, tips=tips)
        res[band] = absb.od['o2'].sum(axis=0)   # column O2 OT
    a, b = res['o2a'], res['o2b']
    assert a.max() > 50.0, a.max()                       # saturated A-band cores
    assert a.max() > 5.0 * b.max(), (a.max(), b.max())   # A >> B
    assert a.min() >= 0.0 and b.min() >= 0.0
    assert a.max() / np.median(a) > 100.0                # strong line/continuum contrast
    return ('column O2 OT: A max %.1f, B max %.1f (A/B=%.1f); '
            'A max/median=%.0f' % (a.max(), b.max(), a.max() / b.max(),
                                   a.max() / np.median(a)))


def Vsub_subrange_guard():
    atm = afgl_atmosphere(os.path.join(DATA, 'afglms.dat'), z_top=70.0)
    tips = tips2021()
    lines = hitran_lines(os.path.join(DATA, 'hitran2020_lines.txt'), wl_range=(763.0, 764.0))
    absb = o2band_absorption(atm, lines, (763.0, 764.0), dwvl=0.001, tips=tips)

    # None -> full grid
    assert absb.subrange_indices(None).size == absb.wvl.size
    # valid in-band sub-range -> correct subset, points inside [a,b]
    idx = absb.subrange_indices((763.20, 763.50))
    assert idx.size == 301, idx.size
    assert absb.wvl[idx].min() >= 763.20 - 1e-9 and absb.wvl[idx].max() <= 763.50 + 1e-9
    assert np.array_equal(idx, np.arange(idx[0], idx[-1] + 1))  # contiguous
    # out-of-band (below and above) must raise
    for bad in [(762.5, 763.5), (763.5, 764.5), (760.0, 761.0)]:
        try:
            absb.subrange_indices(bad)
            raise AssertionError('expected ValueError for out-of-band %s' % (bad,))
        except ValueError:
            pass
    return 'None→full; in-band→301 pts contiguous; out-of-band raises'


def V8_reproducibility():
    atm = afgl_atmosphere(os.path.join(DATA, 'afglms.dat'), z_top=70.0)
    tips = tips2021()
    lines = hitran_lines(os.path.join(DATA, 'hitran2020_lines.txt'),
                         wl_range=(763.0, 764.0))
    a1 = o2band_absorption(atm, lines, (763.0, 764.0), dwvl=0.001, tips=tips)
    a2 = o2band_absorption(atm, lines, (763.0, 764.0), dwvl=0.001, tips=tips)
    assert np.array_equal(a1.od_total, a2.od_total), 'non-deterministic output'
    assert np.array_equal(a1.wvl, a2.wvl)
    return 'two runs bit-identical (deterministic)'


CHECKS = [
    ('V1  HITRAN parse',            V1_hitran_parse),
    ('V2  AFGL profile',            V2_afgl_profile),
    ('Vt  TIPS-2021 partition',     Vtips_partition),
    ('Vo  air<->vac wavelength',    V_optics_airvac),
    ('V3  Voigt width vs P',        V3_voigt_width_vs_pressure),
    ('V4  Voigt area norm',         V4_area_normalisation),
    ('V5  Rayleigh OT',             V5_rayleigh),
    ('V6  O2 band OT',              V6_o2_band_ot),
    ('Vs  sub-range guard',         Vsub_subrange_guard),
    ('V8  reproducibility',         V8_reproducibility),
]


def main():
    print('=' * 78)
    print('O2-band absorption validation  (PLAN.md §6)')
    print('=' * 78)
    n_pass = 0
    for name, fn in CHECKS:
        try:
            detail = fn()
            print('[PASS] %-26s %s' % (name, detail))
            n_pass += 1
        except Exception as e:
            print('[FAIL] %-26s %s' % (name, e))
            traceback.print_exc()
    print('-' * 78)
    print('%d/%d checks passed  (V7 reflectance limits deferred to RT wiring)'
          % (n_pass, len(CHECKS)))
    return 0 if n_pass == len(CHECKS) else 1


if __name__ == '__main__':
    sys.exit(main())
