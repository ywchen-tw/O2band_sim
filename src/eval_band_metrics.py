#!/usr/bin/env python
"""
Band metrics vs published values (EVAL_PLAN.md #2 / Q2, Q5-style).

Compares config-robust band metrics of our simulation against citable published
reference values, and reports the difference (statistics, not pass/fail).

Directly comparable published metrics (atmosphere-derived, geometry-independent):

  Rayleigh column optical depth vs **Hansen & Travis (1974)**, Space Sci. Rev. 16,
      527 -- the widely used total-column parameterization for a 1013.25 mb
      standard atmosphere:
          tau_R = 0.008569 * L^-4 (1 + 0.0113 L^-2 + 0.00013 L^-4),  L = lambda[um].
      Comparing to it validates our Bodhaine cross-section AND the air-column
      integration together (cf. eval_rayleigh.py, which checks the cross-section
      alone against Bucholtz 1995).

  O2 volume mixing ratio / vertical column vs the canonical dry-air value
      **0.2095** (e.g. US Standard Atmosphere) -- validates absorber amount and
      the exponential-in-z column integration.

Reflectance band metrics (continuum reflectance, band-integrated reflectance,
equivalent width) are geometry/surface-dependent and have no single clean
published value for the prescribed config; they are printed for reference only
(compare against HAPI for absorption, or the participant ensemble, when available).

    python src/eval_band_metrics.py            # atmosphere-derived metrics vs published
    python src/eval_band_metrics.py OUT.h5     # also print reflectance metrics from a run
"""

import os
import sys
import argparse
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from util.atmosphere import afgl_atmosphere
from util.absorption import cal_rayleigh_od, BANDS


# canonical dry-air O2 volume mixing ratio (US Standard Atmosphere)
O2_VMR_STANDARD = 0.2095


def rayleigh_od_hansen_travis(wvl_nm):
    """Total-column Rayleigh optical depth, Hansen & Travis (1974), 1013.25 mb."""
    L = np.asarray(wvl_nm, dtype=np.float64) * 1.0e-3         # nm -> um
    return 0.008569 * L ** -4 * (1.0 + 0.0113 * L ** -2 + 0.00013 * L ** -4)


def _reldiff(ours, ref):
    return 100.0 * (ours / ref - 1.0)


def compare_atmosphere(atm):
    """Print atmosphere-derived band metrics vs published references."""
    air_col = float(atm.lay['air'].sum())
    o2_col = float(atm.lay['o2'].sum())
    o2_vmr_col = o2_col / air_col

    print('=== O2 vertical column vs canonical dry-air VMR (0.2095) ===')
    print('  air column      : %.4e molec cm^-2' % air_col)
    print('  O2 column (ours): %.4e molec cm^-2' % o2_col)
    print('  O2/air VMR (ours): %.5f   canonical: %.5f   diff %+.2f%%'
          % (o2_vmr_col, O2_VMR_STANDARD, _reldiff(o2_vmr_col, O2_VMR_STANDARD)))
    print('  O2 column @ 0.2095*air = %.4e   ours %.4e   diff %+.2f%%'
          % (O2_VMR_STANDARD * air_col, o2_col,
             _reldiff(o2_col, O2_VMR_STANDARD * air_col)))

    print('\n=== Rayleigh column OT (Bodhaine*column, ours) vs Hansen & Travis 1974 ===')
    print('  %-8s %14s %14s %9s' % ('wvl[nm]', 'ours', 'H&T 1974', 'diff'))
    # 550 nm is a canonical anchor; 688/760 nm are the band reference wavelengths
    for wl in (550.0, 688.0, 760.0):
        ours = float(cal_rayleigh_od(atm, wl).sum())          # sum over layers = column
        ref = float(rayleigh_od_hansen_travis(wl))
        print('  %-8.1f %14.5f %14.5f %+8.2f%%' % (wl, ours, ref, _reldiff(ours, ref)))


def print_reflectance_metrics(path):
    """Reflectance band metrics from a run (config-dependent; reference only)."""
    from eval_metrics import band_metrics
    print('\n=== reflectance band metrics (config-dependent; no single published '
          'value -- reference only) ===')
    M = band_metrics(path)
    for band, rec in M.items():
        print('[%s]' % band)
        for (s, a), g in sorted(rec['geom'].items()):
            print('  SZA %4.1f alb %.2f: continuum rho=%.5f  band-int=%.4f  EW=%.4f nm'
                  % (s, a, g['continuum_reflectance'],
                     g['band_integrated_reflectance'], g['equivalent_width_nm']))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('run', nargs='?', default=None,
                   help='optional run HDF5 for reflectance band metrics')
    p.add_argument('--z-top', type=float, default=120.0, help='km (default 120)')
    args = p.parse_args()

    data_dir = os.environ.get('O2BAND_DATA_DIR',
                              os.path.normpath(os.path.join(_HERE, '..', 'data')))
    atm = afgl_atmosphere(os.path.join(data_dir, 'afglms.dat'), z_top=args.z_top)
    print('AFGL mid-lat-summer, z_top=%.0f km, P_sfc=%.1f hPa\n'
          % (args.z_top, atm.lev['p'][0]))
    compare_atmosphere(atm)
    if args.run:
        if not os.path.isfile(args.run):
            sys.exit('Error: no such file: %s' % args.run)
        print_reflectance_metrics(args.run)
