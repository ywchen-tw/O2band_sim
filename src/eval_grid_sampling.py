#!/usr/bin/env python
"""
Is the prescribed 0.001 nm grid a fair *sample* of the spectrum?

The benchmark grid is 0.001 nm, but the narrowest lines are not much wider than
that: the Voigt FWHM of an O2 A-band line is ~0.0054 nm at the surface but only
~0.0015 nm near 10 hPa, i.e. ~1.5 grid points per FWHM where lines are
Doppler-limited.  A point sample on such a grid is not obviously equal to the
average over the 0.001 nm interval it is supposed to represent.

Two distinct questions, both tested here per wavelength:

  1. Optical depth      -- is tau(grid point) == <tau> over the 0.001 nm cell?
  2. **Reflectance**    -- is rho(tau_point) == <rho> over the cell?

(2) is the one that matters and is NOT implied by (1): rho depends on tau
through an exponential, so <rho(tau)> != rho(<tau>) whenever tau varies across
the cell (Jensen's inequality -- the average reflectance is always >= the
reflectance of the average optical depth).  A model reporting point samples and
a model reporting cell averages will therefore disagree even with identical
physics, most strongly on steep line flanks.

Method: build the absorption on a 5x finer 0.0002 nm grid, which contains the
0.001 nm points exactly as every 5th sample, then compare the point value with
the 5-sample cell mean centred on it.

    python src/eval_grid_sampling.py
    python src/eval_grid_sampling.py --band o2a --sza 0 --albedo 0.1
"""

import os
import sys
import argparse
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from util.inputs import require_inputs, MissingInputs
from util.atmosphere import afgl_atmosphere
from util.tips import tips2021
from util.absorption import (hitran_lines, o2band_absorption, cal_rayleigh_od,
                             required_margin_cm, DEFAULT_CUTOFF_CM)
from eval_cutoff_cia import rho_estimate

DATA = os.environ.get('O2BAND_DATA_DIR',
                      os.path.normpath(os.path.join(_HERE, '..', 'data')))

DWVL = 0.001          # prescribed benchmark grid
NSUB = 5              # sub-samples per cell -> 0.0002 nm


def build(atm, tips, win, dwvl, z_top):
    margin = required_margin_cm(win, grid='vac', cutoff_cm=DEFAULT_CUTOFF_CM)
    lines = hitran_lines(os.path.join(DATA, 'hitran2020_lines.txt'),
                         wl_range=win, grid='vac', margin_cm=margin)
    return o2band_absorption(atm, lines, win, dwvl=dwvl,
                             cutoff_cm=DEFAULT_CUTOFF_CM, grid='vac', tips=tips)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--band', default='o2a', choices=('o2a', 'o2b'))
    p.add_argument('--window', nargs=2, type=float, default=[763.0, 763.2],
                   metavar=('A', 'B'), help='vacuum-nm scan window')
    p.add_argument('--z-top', type=float, default=120.0)
    p.add_argument('--albedo', type=float, default=0.1)
    p.add_argument('--sza', type=float, default=0.0)
    args = p.parse_args()

    win = tuple(args.window)
    try:
        require_inputs('hitran2020_lines.txt', 'afglms.dat', qtpy=True)
    except MissingInputs as e:
        sys.exit('[SETUP] %s' % e)
    atm = afgl_atmosphere(os.path.join(DATA, 'afglms.dat'), z_top=args.z_top)
    tips = tips2021()

    # 0.0002 nm grid; the 0.001 nm points are exactly every 5th sample
    fine = build(atm, tips, win, DWVL / NSUB, args.z_top)
    od_f = fine.od_total.sum(axis=0)
    ray_f = cal_rayleigh_od(atm, fine.wvl_vac).sum(axis=0)
    rho_f = np.array([rho_estimate(r, t, args.albedo, args.sza)
                      for r, t in zip(ray_f, od_f)])

    # coarse grid points = indices 0, 5, 10, ... ; cell = point +/- 2 subsamples
    ic = np.arange(2, od_f.size - 2, NSUB)
    cells = ic[:, None] + np.arange(-2, 3)[None, :]

    od_pt = od_f[ic]
    od_cell = od_f[cells].mean(axis=1)
    rho_pt = rho_f[ic]
    rho_cell = rho_f[cells].mean(axis=1)
    rho_of_mean_od = np.array([rho_estimate(r, t, args.albedo, args.sza)
                               for r, t in zip(ray_f[ic], od_cell)])

    # verify the fine grid really contains the coarse points
    coarse = build(atm, tips, win, DWVL, args.z_top)
    j = np.searchsorted(coarse.wvl, fine.wvl[ic[0]])
    same = np.allclose(coarse.wvl[j:j + ic.size], fine.wvl[ic], atol=1e-9)
    print('=' * 78)
    print('0.001 nm point sample vs 0.0002 nm cell average')
    print('=' * 78)
    print('band %s, window %.3f-%.3f nm vacuum, z_top %.0f km, albedo %.2f, SZA %.0f'
          % (args.band, win[0], win[1], args.z_top, args.albedo, args.sza))
    print('grids aligned (0.001 pts are every 5th 0.0002 pt): %s' % same)

    # Pick by REFLECTANCE impact, not OD: a saturated core has rho = 0 either
    # way, so a large OD difference there is radiatively irrelevant.  Relative
    # differences are only meaningful where rho is not underflowing.
    d_rho_abs = rho_cell - rho_pt
    picks = [('line core    (max OD)', int(np.argmax(od_pt))),
             ('worst case   (max |d rho|)', int(np.argmax(np.abs(d_rho_abs)))),
             ('micro-window (min OD)', int(np.argmin(od_pt)))]

    for name, k in picks:
        sub_od = od_f[cells[k]]
        sub_rho = rho_f[cells[k]]
        print('\n--- %s : %.4f nm vacuum (%.4f nm air) ---'
              % (name, fine.wvl_vac[ic[k]], fine.wvl_air[ic[k]]))
        print('  %-14s %12s %12s' % ('sub-sample nm', 'column OD', 'rho'))
        for w, t, r in zip(fine.wvl[cells[k]], sub_od, sub_rho):
            mark = '  <- 0.001 nm grid point' if abs(w - fine.wvl[ic[k]]) < 1e-9 else ''
            print('  %-14.4f %12.5f %12.5f%s' % (w, t, r, mark))
        d_od = 100 * (od_cell[k] - od_pt[k]) / max(od_pt[k], 1e-12)
        print('  point   : OD %12.5f   rho %.6f' % (od_pt[k], rho_pt[k]))
        print('  cell avg: OD %12.5f   rho %.6f' % (od_cell[k], rho_cell[k]))
        print('  point vs cell average:  OD %+.2f%%   rho %+.2e (absolute)'
              % (d_od, rho_cell[k] - rho_pt[k]))
        if rho_pt[k] > 1e-4:
            print('                          rho %+.3f%% relative; curvature term '
                  '<rho(tau)> - rho(<tau>) = %+.2e'
                  % (100 * (rho_cell[k] - rho_pt[k]) / rho_pt[k],
                     rho_cell[k] - rho_of_mean_od[k]))
        else:
            print('                          rho underflows to 0 either way '
                  '-> radiatively irrelevant despite the OD difference')

    # window-wide context: is the worst case rare or typical?
    dr_abs = rho_cell - rho_pt
    ok = rho_pt > 1e-3                        # where relative error is meaningful
    dr_rel = 100 * dr_abs[ok] / rho_pt[ok]
    dt = 100 * (od_cell - od_pt) / np.maximum(od_pt, 1e-12)
    print('\n--- over the whole %.3f-%.3f nm window (%d grid points) ---'
          % (win[0], win[1], ic.size))
    print('  rho (absolute) : median %+.2e  p95 |diff| %.2e  max |diff| %.2e'
          % (np.median(dr_abs), np.percentile(np.abs(dr_abs), 95), np.abs(dr_abs).max()))
    print('  rho (relative, %d pts with rho>1e-3): median %+.3f%%  p95 %.3f%%  max %.3f%%'
          % (ok.sum(), np.median(dr_rel), np.percentile(np.abs(dr_rel), 95),
             np.abs(dr_rel).max()))
    print('  OD  (relative) : median %+.3f%%  p95 %.3f%%  max %.3f%%'
          % (np.median(dt), np.percentile(np.abs(dt), 95), np.abs(dt).max()))
    print('  window-mean rho: point %.6f  cell avg %.6f  (%+.4f%%)'
          % (rho_pt.mean(), rho_cell.mean(),
             100 * (rho_cell.mean() - rho_pt.mean()) / rho_pt.mean()))
    print('\n  for scale: MC noise on the production run is ~0.18%% p95 (O2A, 1e7 photons)')


if __name__ == '__main__':
    main()
