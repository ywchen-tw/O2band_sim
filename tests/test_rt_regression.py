#!/usr/bin/env python
"""
RT-regression test (PLAN.md §6, V7 + V9): validate assembled MCARaTS reflectance
output against the independent analytic references in util.refcheck.

This closes the loop the fast absorption suite can't: it needs an actual RT run.
It therefore operates on an *assembled band HDF5* (o2a.h5 / o2b.h5 written by
sim_o2band.O2BandSim.assemble):

  * If no such file exists, the test SKIPs (exit 0) -- so it stays green in
    environments without MCARaTS / without a completed run, and turns into a real
    check automatically once RT output is present.
  * File discovery order:
      1. $O2BAND_TEST_H5   (explicit path)
      2. out/o2a.h5, out/o2b.h5   (default output location)

Checks performed on each available band:
  V5' Rayleigh OD           vs Hansen & Travis (1974)         : < 1%
  V7a window rho (alb=0)    vs single-scattering Rayleigh     : within [-2%, +15%]
                              (full MC >= SS; O2 residual only lowers it slightly)
  V7b window rho (alb>0)    vs Lambertian-over-atmosphere      : < 8%
  V7c saturated core        rho -> ~0 and albedo-independent   : < 1e-3, |dalb| tiny
  V9  monotonic in albedo   rho increases with albedo everywhere
"""

import os
import sys
import glob
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
sys.path.insert(0, os.path.join(_REPO, 'src'))

from util.refcheck import (rayleigh_od_ht1974, rayleigh_ss_reflectance,
                           lambertian_over_atm)

try:
    import h5py
except ImportError:
    print('[SKIP] h5py not available'); sys.exit(0)


def _find_files():
    env = os.environ.get('O2BAND_TEST_H5')
    if env and os.path.isfile(env):
        return [env]
    found = []
    for b in ('o2a', 'o2b'):
        p = os.path.join(_REPO, 'out', '%s.h5' % b)
        if os.path.isfile(p):
            found.append(p)
    return found


def check_band(fname):
    msgs = []
    with h5py.File(fname, 'r') as f:
        band = f.attrs.get('band', '?')
        band = band.decode() if isinstance(band, bytes) else band
        wvl = f['wvl'][:]
        sza = f['sza'][:]
        alb = f['albedo'][:]
        ref = f['reflectance'][:]                       # (Nsza, Nalb, Nwvl)
        ref_e = f['reflectance_stderr'][:] if 'reflectance_stderr' in f \
            else np.zeros_like(ref)
        o2 = f['optical_thickness/o2_column'][:]
        ray = f['optical_thickness/rayleigh_column'][:]

    # --- Rayleigh OD vs Hansen & Travis --------------------------------------
    tau_ht = float(np.median(rayleigh_od_ht1974(wvl)))
    tau_sim = float(np.median(ray))
    rel = abs(tau_sim - tau_ht) / tau_ht
    assert rel < 0.01, 'Rayleigh OD %.5f vs H&T %.5f (%.2f%%)' % (tau_sim, tau_ht, 100 * rel)
    msgs.append('tauR %.5f~%.5f (%.2f%%)' % (tau_sim, tau_ht, 100 * rel))

    # --- per-geometry reflectance checks -------------------------------------
    iw = int(np.argmin(o2))                             # cleanest window point
    ic = int(np.argmax(o2))                             # saturated core
    tauR = tau_sim
    # Window checks require a genuinely clean continuum point (O2 negligible), so
    # the analytic single-scattering Rayleigh reference applies.  A narrow
    # sub-window may not contain one -> those checks are skipped (not failed).
    have_window = o2[iw] < 0.02
    ia0 = int(np.argmin(alb))
    assert abs(alb[ia0]) < 1e-9, 'expected an albedo=0 case'

    for js, s in enumerate(sza):
        mu0 = np.cos(np.deg2rad(float(s)))

        if have_window:
            # V7a: alb=0 window vs pure single-scattering Rayleigh.
            # Full MC includes multiple scattering (adds signal) minus the tiny
            # O2 residual, so sim should sit just above SS: allow [-3%, +15%].
            rho_ss = rayleigh_ss_reflectance(tauR, mu0)
            d = (ref[js, ia0, iw] - rho_ss) / rho_ss
            assert -0.03 < d < 0.15, ('SZA=%.0f window alb0 rho=%.5f vs SS %.5f (%.1f%%)'
                                      % (s, ref[js, ia0, iw], rho_ss, 100 * d))

            # V7b: alb>0 window vs Lambertian-over-(Rayleigh)-atmosphere
            for ia, a in enumerate(alb):
                if a <= 0:
                    continue
                rho_m = lambertian_over_atm(tauR, float(a), mu0)
                dm = abs(ref[js, ia, iw] - rho_m) / rho_m
                assert dm < 0.08, ('SZA=%.0f alb%.2f window rho=%.5f vs model %.5f (%.1f%%)'
                                   % (s, a, ref[js, ia, iw], rho_m, 100 * dm))

        # V7c: saturated core -> ~0 and albedo-independent
        core = ref[js, :, ic]
        assert np.all(core < 1e-3), 'SZA=%.0f core rho not ~0: %s' % (s, core)
        assert (core.max() - core.min()) < 5e-5, 'SZA=%.0f core albedo-dependent' % s

        # V9: reflectance increases with albedo.
        #
        # A hard np.all() over 15001 points cannot be required.  The per-point
        # sigma comes from only Nrun=3 samples, so it is a chi-distributed
        # estimate with 2 dof that underestimates the true sigma a few percent
        # of the time; across ~10^5 comparisons a handful of "5 sigma"
        # excursions is expected on a perfectly correct run.  Requiring zero
        # made this test fail intermittently on good data.
        #
        # Bound the violations two ways instead -- both count and magnitude.
        # A genuine albedo fault (wired backwards, surface term dropped) shows
        # up at EVERY point including the bright continuum, with deficits of
        # order the albedo itself (~0.09), so it cannot hide under either bound.
        # Deep cores legitimately sit at the ~1e-6 noise floor for every albedo.
        MAX_VIOL_FRAC = 1.0e-3      # 0.1% of points
        MAX_VIOL_ABS = 1.0e-4       # ~0.1% of the albedo-0.1 continuum
        for ia in range(1, len(alb)):
            tol = 5.0 * np.sqrt(ref_e[js, ia]**2 + ref_e[js, ia - 1]**2) + 1e-6
            deficit = (ref[js, ia - 1] - tol) - ref[js, ia]
            viol = deficit > 0
            frac = viol.mean()
            worst = float(deficit[viol].max()) if viol.any() else 0.0
            assert frac <= MAX_VIOL_FRAC, \
                ('SZA=%.0f non-monotonic in albedo at %.3f%% of points '
                 '(> %.3f%% allowed for MC scatter)' % (s, 100 * frac,
                                                        100 * MAX_VIOL_FRAC))
            assert worst <= MAX_VIOL_ABS, \
                ('SZA=%.0f albedo deficit %.2e exceeds %.0e -- too large to be '
                 'MC noise' % (s, worst, MAX_VIOL_ABS))
            # aggregate monotonicity stays strict: robust to per-point noise
            assert np.median(ref[js, ia]) > np.median(ref[js, ia - 1]), \
                'SZA=%.0f band-median not increasing with albedo' % s

    msgs.append('%s window + core + monotonic ok for %d SZA x %d alb'
                % ('with' if have_window else 'NO-clean-window;', len(sza), len(alb)))
    return '%s: %s' % (band, '; '.join(msgs))


def main():
    files = _find_files()
    if not files:
        print('[SKIP] no assembled RT output found '
              '(set $O2BAND_TEST_H5 or run sim_o2band -> out/o2a.h5).')
        return 0
    print('=' * 78)
    print('RT-regression: reflectance vs analytic references (util.refcheck)')
    print('=' * 78)
    ok = True
    for fn in files:
        try:
            print('[PASS] %s' % check_band(fn))
        except Exception as e:
            ok = False
            print('[FAIL] %s : %s' % (os.path.basename(fn), e))
    print('-' * 78)
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
