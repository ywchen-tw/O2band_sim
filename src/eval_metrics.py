#!/usr/bin/env python
"""
Evaluation metrics + difference statistics for the O2-band benchmark output
(EVAL_PLAN.md).  Two independent pieces:

1. ``band_metrics(path)`` -- resolution-robust scalar metrics of *our* result,
   read from a per-band (o2a.h5) or merged (o2band_benchmark.h5) HDF5:
   continuum reflectance, band-integrated + band-mean reflectance, reflectance
   equivalent width, and column O2/H2O/Rayleigh optical thickness, per
   (SZA, albedo).  These are the quantities compared against published band
   metrics (EVAL_PLAN.md Q1-Q5), robust to a reference's spectral resolution.

2. ``diff_stats(a, b, stderr=)`` -- the difference distribution between our
   spectrum ``a`` and a reference ``b`` on a shared grid: mean bias, median,
   RMS, max of the absolute and relative difference, 5/50/95th percentiles of
   the relative difference, spectral correlation, and -- if ``stderr`` is given
   -- the difference in MC-noise units (|a-b|/stderr).  This is the engine for
   every per-lambda reference comparison; it reports statistics, not verdicts.

Run directly to print the band metrics for an output file:

    python src/eval_metrics.py /scratch/.../z120_p1e6_n3/o2a.h5
"""

import os
import sys
import argparse
import numpy as np
import h5py


# ---------------------------------------------------------------------------- #
def _iter_bands(f):
    """Yield (band_name, group) for a per-band file or a merged file."""
    if 'reflectance' in f:
        yield f.attrs.get('band', 'band'), f
    else:
        for name in f:
            g = f[name]
            if isinstance(g, h5py.Group) and 'reflectance' in g:
                yield name, g


def _continuum(wvl, ref, o2_col):
    """Continuum reflectance rho_c(lambda): linear fit of rho vs wvl over the
    lowest-absorption (window) points, so it tracks the Rayleigh slope instead
    of assuming a flat continuum.  Returns (rho_c_array, window_mask)."""
    thr = o2_col.min() + 0.05 * (o2_col.max() - o2_col.min() + 1e-30)
    win = o2_col <= max(thr, np.percentile(o2_col, 10))
    if win.sum() < 2:
        win = o2_col <= np.percentile(o2_col, 25)
    # robust linear continuum vs wavelength
    c = np.polyfit(wvl[win], ref[win], 1)
    return np.polyval(c, wvl), win


# ---------------------------------------------------------------------------- #
def band_metrics(path):
    """Return {band: {'wvl_nm':(a,b), per (sza,alb): {...}, 'column_ot':{...}}}."""
    out = {}
    with h5py.File(path, 'r') as f:
        for band, g in _iter_bands(f):
            wvl = g['wvl'][:]
            sza = g['sza'][:]
            alb = g['albedo'][:]
            ref = g['reflectance'][:]                 # (Nsza, Nalb, Nwvl)
            o2c = g['optical_thickness/o2_column'][:]
            h2oc = g['optical_thickness/h2o_column'][:]
            rayc = g['optical_thickness/rayleigh_column'][:]

            rec = {'wvl_nm': (float(wvl[0]), float(wvl[-1])),
                   'n_wvl': int(wvl.size),
                   'column_ot': {
                       'o2_max': float(o2c.max()), 'o2_mean': float(o2c.mean()),
                       'h2o_max': float(h2oc.max()), 'h2o_mean': float(h2oc.mean()),
                       'rayleigh_mean': float(rayc.mean())},
                   'geom': {}}
            for i, s in enumerate(sza):
                for j, a in enumerate(alb):
                    r = ref[i, j]
                    finite = np.isfinite(r)
                    rho_c, win = _continuum(wvl[finite], r[finite], o2c[finite])
                    # band-depth / equivalent width (nm): area of absorbed fraction
                    frac = 1.0 - r[finite] / np.where(rho_c > 0, rho_c, np.nan)
                    ew = float(np.trapz(np.nan_to_num(frac), wvl[finite]))
                    rec['geom'][(float(s), float(a))] = {
                        'continuum_reflectance': float(np.mean(rho_c)),
                        'band_integrated_reflectance': float(np.trapz(r[finite], wvl[finite])),
                        'band_mean_reflectance': float(np.mean(r[finite])),
                        'min_reflectance': float(np.min(r[finite])),
                        'equivalent_width_nm': ew,
                        'n_window_pts': int(win.sum())}
            out[band] = rec
    return out


# ---------------------------------------------------------------------------- #
def diff_stats(a, b, stderr=None, eps=1e-12):
    """Difference distribution of a (ours) vs b (reference) on a shared grid.

    Returns a dict of descriptive statistics -- no pass/fail.  Relative stats use
    only points where |b| > eps.  ``stderr`` (per-point) adds |a-b|/stderr.
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    d = a - b
    out = {
        'n': int(a.size),
        'mean_bias': float(np.mean(d)),
        'median_diff': float(np.median(d)),
        'rms_diff': float(np.sqrt(np.mean(d ** 2))),
        'max_abs_diff': float(np.max(np.abs(d))) if a.size else np.nan,
        'corr': float(np.corrcoef(a, b)[0, 1]) if a.size > 1 else np.nan,
    }
    rmask = np.abs(b) > eps
    if rmask.any():
        rd = d[rmask] / b[rmask]
        out.update({
            'rel_mean_bias': float(np.mean(rd)),
            'rel_rms': float(np.sqrt(np.mean(rd ** 2))),
            'rel_max_abs': float(np.max(np.abs(rd))),
            'rel_p5': float(np.percentile(rd, 5)),
            'rel_p50': float(np.percentile(rd, 50)),
            'rel_p95': float(np.percentile(rd, 95)),
        })
    if stderr is not None:
        se = np.asarray(stderr, dtype=np.float64)[m]
        smask = se > eps
        if smask.any():
            ns = np.abs(d[smask]) / se[smask]
            out['diff_in_noise_units_mean'] = float(np.mean(ns))
            out['diff_in_noise_units_max'] = float(np.max(ns))
    return out


def print_diff_stats(name, st):
    print('  %s: n=%d  bias=%.3e  rms=%.3e  max|Δ|=%.3e  corr=%.5f'
          % (name, st['n'], st['mean_bias'], st['rms_diff'],
             st['max_abs_diff'], st['corr']))
    if 'rel_rms' in st:
        print('      rel: bias=%+.2e rms=%.2e max=%.2e  p5/50/95=%+.2e/%+.2e/%+.2e'
              % (st['rel_mean_bias'], st['rel_rms'], st['rel_max_abs'],
                 st['rel_p5'], st['rel_p50'], st['rel_p95']))
    if 'diff_in_noise_units_mean' in st:
        print('      in MC-noise units: mean=%.2f  max=%.2f'
              % (st['diff_in_noise_units_mean'], st['diff_in_noise_units_max']))


# ---------------------------------------------------------------------------- #
def _print_band_metrics(path):
    M = band_metrics(path)
    print('band metrics: %s' % path)
    for band, rec in M.items():
        c = rec['column_ot']
        print('\n[%s] %.3f-%.3f nm (%d pts)  column OT: O2 max=%.1f mean=%.2f | '
              'H2O max=%.3f | Rayleigh mean=%.4f'
              % (band, rec['wvl_nm'][0], rec['wvl_nm'][1], rec['n_wvl'],
                 c['o2_max'], c['o2_mean'], c['h2o_max'], c['rayleigh_mean']))
        print('  %-5s %-6s | %10s %12s %10s %10s' %
              ('SZA', 'alb', 'continuum', 'band_int', 'EW_nm', 'min_rho'))
        for (s, a), g in sorted(rec['geom'].items()):
            print('  %-5.1f %-6.2f | %10.5f %12.4f %10.4f %10.3e'
                  % (s, a, g['continuum_reflectance'], g['band_integrated_reflectance'],
                     g['equivalent_width_nm'], g['min_reflectance']))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('path', help='per-band (o2a.h5) or merged (o2band_benchmark.h5) HDF5')
    args = p.parse_args()
    if not os.path.isfile(args.path):
        sys.exit('Error [eval_metrics]: no such file: %s' % args.path)
    _print_band_metrics(args.path)
