#!/usr/bin/env python
"""
Monte-Carlo noise report for an O2-band benchmark HDF5 (PLAN.md step 2f).

Reads a per-band file (``o2a.h5`` / ``o2b.h5`` -- datasets at the top level) or the
merged file (``o2band_benchmark.h5`` -- one group per band) and, for every
(SZA, albedo), summarises the reflectance MC standard error so the noise sign-off
is a single command:

    python src/noise_report.py /scratch/alpine/yuch8913/O2band_sim/o2a.h5
    python src/noise_report.py .../o2band_benchmark.h5 --threshold 0.01

Two complementary metrics, because reflectance spans transparent windows (rho ~ O(0.1))
and saturated line cores (rho -> 0):

- **relative** stderr  rho_stderr / rho, evaluated only where rho exceeds a floor
  (``--ref-floor``, default 1e-3) so near-zero cores do not produce a spurious
  divide-by-zero blow-up.  This is the metric PLAN.md 2f thresholds.
- **absolute** stderr  rho_stderr, whose worst case lives in the saturated cores
  where the relative metric is undefined; reported so core noise is still visible.

With ``--threshold T`` the script exits non-zero if any (SZA, albedo)'s p95
relative stderr exceeds T, so it doubles as a pass/fail gate in a batch script.
"""

import os
import sys
import argparse
import numpy as np
import h5py


def _iter_bands(f):
    """Yield (band_name, group) for each band in either file layout."""
    if 'reflectance' in f:                       # per-band file: datasets at root
        yield f.attrs.get('band', os.path.basename(getattr(f, 'filename', 'band'))), f
    else:                                        # merged file: one group per band
        for name in f:
            g = f[name]
            if isinstance(g, h5py.Group) and 'reflectance' in g:
                yield name, g


def _summarise(ref, stderr, ref_floor):
    """Relative + absolute stderr summary for one (SZA, albedo) slice (1-D)."""
    ref = np.asarray(ref, dtype=np.float64)
    stderr = np.asarray(stderr, dtype=np.float64)

    finite = np.isfinite(ref) & np.isfinite(stderr)
    # absolute noise: worst case anywhere (this is where saturated cores show up)
    abs_max = float(np.nanmax(stderr[finite])) if finite.any() else np.nan
    abs_arg = int(np.nanargmax(np.where(finite, stderr, -np.inf))) if finite.any() else -1

    # relative noise: only where signal exceeds the floor
    sig = finite & (ref > ref_floor)
    if sig.any():
        rel = stderr[sig] / ref[sig]
        med = float(np.median(rel))
        p95 = float(np.percentile(rel, 95))
        rmax = float(np.max(rel))
    else:
        med = p95 = rmax = np.nan
    return dict(n_signal=int(sig.sum()), rel_med=med, rel_p95=p95, rel_max=rmax,
                abs_max=abs_max, abs_arg=abs_arg)


def noise_report(path, ref_floor=1.0e-3, threshold=None, verbose=True):
    """Print the per-(SZA, albedo) noise table; return (worst_p95, ok)."""
    worst_p95 = 0.0
    ok = True
    with h5py.File(path, 'r') as f:
        if verbose:
            print('MC noise report: %s' % path)
            print('  (relative stderr computed where reflectance > %.1e)' % ref_floor)
        for band, g in _iter_bands(f):
            wvl = g['wvl'][:]
            sza = g['sza'][:]
            alb = g['albedo'][:]
            ref = g['reflectance'][:]              # (Nsza, Nalb, Nwvl)
            err = g['reflectance_stderr'][:]
            if verbose:
                print('\n[%s]  %d wvl %.3f-%.3f nm' % (band, wvl.size, wvl[0], wvl[-1]))
                print('  %-5s %-6s | %8s %8s %8s | %10s @ %-9s'
                      % ('SZA', 'alb', 'rel_med', 'rel_p95', 'rel_max',
                         'abs_max', 'wvl_nm'))
            for i, s in enumerate(sza):
                for j, a in enumerate(alb):
                    r = _summarise(ref[i, j], err[i, j], ref_floor)
                    if np.isfinite(r['rel_p95']):
                        worst_p95 = max(worst_p95, r['rel_p95'])
                    flag = ''
                    if threshold is not None and np.isfinite(r['rel_p95']) \
                            and r['rel_p95'] > threshold:
                        ok = False
                        flag = '  <-- over %.3g' % threshold
                    wl_at = wvl[r['abs_arg']] if r['abs_arg'] >= 0 else np.nan
                    if verbose:
                        print('  %-5.1f %-6.2f | %8.4f %8.4f %8.4f | %10.2e @ %-9.3f%s'
                              % (s, a, r['rel_med'], r['rel_p95'], r['rel_max'],
                                 r['abs_max'], wl_at, flag))
    if verbose:
        print('\nworst relative-stderr p95 across all (band, SZA, albedo): %.4f' % worst_p95)
        if threshold is not None:
            print('threshold %.3g -> %s' % (threshold, 'PASS' if ok else 'FAIL'))
    return worst_p95, ok


def _cli():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('path', help='per-band (o2a.h5) or merged (o2band_benchmark.h5) HDF5')
    p.add_argument('--ref-floor', type=float, default=1.0e-3,
                   help='min reflectance for the relative metric (default 1e-3)')
    p.add_argument('--threshold', type=float, default=None,
                   help='fail (exit 1) if any p95 relative stderr exceeds this')
    return p


if __name__ == '__main__':
    args = _cli().parse_args()
    if not os.path.isfile(args.path):
        sys.exit('Error [noise_report]: no such file: %s' % args.path)
    _, ok = noise_report(args.path, ref_floor=args.ref_floor, threshold=args.threshold)
    sys.exit(0 if ok else 1)
