#!/usr/bin/env python
"""
Independent RT-solver cross-check: MCARaTS (ours) vs libRadtran/DISORT (EVAL_PLAN.md E2).

Our TOA reflectance comes from a Monte-Carlo solver (MCARaTS).  This runs the same
clear-sky scene through libRadtran's discrete-ordinate solver (DISORT) and compares
the reflectance -- validating the RT *transport* and the reflectance convention
(rho = pi*I/(mu0*F0)) independently of the absorption physics.

To keep it a solver check (not an absorption-model check), the comparison is done
at several **window wavelengths spanning the band** -- the lowest **total gas
(O2+H2O) OT** grid point in each of n segments (both bands have thousands of
points with total OT ~ 0) -- with libRadtran molecular absorption switched OFF
(``no_absorption mol``): a pure Rayleigh + Lambertian scene our MCARaTS window
points also approximate.  Filtering on total gas OT (not O2 alone) matters in the
B-band, where strong H2O lines (e.g. 693.8 nm) have near-zero O2 but absorb in our
run while the pure-Rayleigh DISORT run does not.  uvspec is run with
``output_quantity reflectivity`` so its ``uu`` output IS rho, matching ours.
Multiple wavelengths check the solver + convention across the Rayleigh-OD range.

uvspec input is written directly (no wrapper), then run via subprocess.  It needs
its runtime libs (GSL/NetCDF) + $LIBRADTRAN_V2_DIR, so run through
``curc_lrt_eval.sh`` (module loads + conda), not a bare shell.

    python src/eval_lrt.py OUR_OUTPUT.h5 [--band o2a] [--streams 16]
"""

import os
import sys
import argparse
import subprocess
import numpy as np
import h5py

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
from eval_metrics import diff_stats

LRT = os.environ.get('LIBRADTRAN_V2_DIR', '/projects/yuch8913/wen_soft/libRadtran-2.0.6')


def _uvspec_input(wvl, sza, alb, atm_file, streams):
    """Monochromatic pure-Rayleigh + Lambertian reflectivity, nadir (umu=+1)."""
    return '\n'.join([
        'data_files_path %s/data' % LRT,
        'atmosphere_file %s' % atm_file,
        'source solar %s/data/solar_flux/kurudz_0.1nm.dat' % LRT,
        'wavelength %.4f %.4f' % (wvl, wvl),
        'rte_solver disort',
        'number_of_streams %d' % streams,
        'no_absorption mol',                   # molecular absorption off; Rayleigh kept
        'albedo %.4f' % alb,
        'sza %.4f' % sza,
        'phi0 0.0', 'umu 1.0', 'phi 0.0',      # nadir view (umu = cos(0) = 1)
        'zout toa',
        'output_user lambda uu',
        'output_quantity reflectivity',        # uu -> reflectance rho
        'quiet', ''])


def disort_reflectivity(wvl, sza, alb, atm_file, streams, workdir):
    inp = _uvspec_input(wvl, sza, alb, atm_file, streams)
    with open(os.path.join(workdir, 'lrt_in_s%02.0f_a%.2f.inp' % (sza, alb)), 'w') as f:
        f.write(inp)
    r = subprocess.run(['%s/bin/uvspec' % LRT], input=inp, capture_output=True,
                       text=True, cwd=workdir)
    if r.returncode != 0:
        raise RuntimeError('uvspec failed (sza=%s alb=%s):\n%s' % (sza, alb, r.stderr[-800:]))
    vals = r.stdout.split()
    if len(vals) < 2:
        raise RuntimeError('unexpected uvspec output: %r (stderr: %s)' % (r.stdout, r.stderr[-400:]))
    return float(vals[1])                       # uu = reflectivity


def select_windows(wvl, gas_ot, n, ot_max, solar_nm=0.1):
    """Indices of ~n window wavelengths spanning the band: the lowest-TOTAL-gas-OT
    grid point in each of n equal segments, kept only if total gas OT < ot_max.

    ``gas_ot`` is the total column absorption OT (O2 + H2O) -- filtering on O2
    alone would admit strong H2O lines (e.g. O2B 693.8 nm, H2O OT ~0.04) where our
    MCARaTS absorbs but the pure-Rayleigh DISORT run does not, corrupting the
    solver comparison.  Candidates are restricted to grid points aligned to the
    solar-flux file resolution (``solar_nm``, kurudz_0.1nm.dat -> 0.1 nm) because
    uvspec requires the requested wavelength to exist in the solar file; our
    0.001 nm grid contains those 0.1 nm points exactly, so MCARaTS and DISORT are
    compared at the identical wavelength."""
    aligned = np.abs(wvl - np.round(wvl / solar_nm) * solar_nm) < 1e-4
    edges = np.linspace(wvl[0], wvl[-1], n + 1)
    idx = []
    for a, b in zip(edges[:-1], edges[1:]):
        seg = np.where((wvl >= a) & (wvl < b + 1e-9) & aligned)[0]
        if seg.size:
            k = seg[int(np.argmin(gas_ot[seg]))]
            if gas_ot[k] < ot_max:
                idx.append(int(k))
    return sorted(set(idx))


def run(our_h5, band, streams, workdir, n_wvl, ot_max):
    os.makedirs(workdir, exist_ok=True)
    atm_file = '%s/data/atmmod/afglms.dat' % LRT

    with h5py.File(our_h5, 'r') as f:
        if band in f and isinstance(f[band], h5py.Group):
            g = f[band]                                  # merged file: pick the band group
        else:
            g = f                                        # per-band file: datasets at root
            file_band = f.attrs.get('band')              # guard against band/file mismatch
            if file_band is not None and str(file_band) != band:
                raise ValueError(
                    "band mismatch: --band %s but %s is band '%s'. Point --band/OUR_H5 "
                    "at the right file (e.g. the o2b.h5 for --band o2b)."
                    % (band, os.path.basename(our_h5), file_band))
        wvl = g['wvl'][:]; sza = g['sza'][:]; alb = g['albedo'][:]
        ref = g['reflectance'][:]
        gas_ot = g['optical_thickness/o2_column'][:] + g['optical_thickness/h2o_column'][:]

    iws = select_windows(wvl, gas_ot, n_wvl, ot_max)
    print('RT-solver check: MCARaTS vs libRadtran/DISORT (%d streams)' % streams)
    print('  %d window wavelengths across %s (total gas OT < %.3g -> near pure-Rayleigh)\n'
          % (len(iws), band, ot_max))
    print('  %-9s %-8s %-5s %-6s | %12s %12s %9s'
          % ('wvl[nm]', 'gas_OT', 'SZA', 'alb', 'MCARaTS', 'DISORT', 'diff'))

    a_ours, a_do = [], []
    for iw in iws:
        wl = float(wvl[iw])
        for i, s in enumerate(sza):
            for j, a in enumerate(alb):
                rho_mc = float(ref[i, j, iw])
                rho_do = disort_reflectivity(wl, float(s), float(a), atm_file, streams, workdir)
                a_ours.append(rho_mc); a_do.append(rho_do)
                print('  %-9.3f %-8.4f %-5.1f %-6.2f | %12.5f %12.5f %+8.2f%%'
                      % (wl, gas_ot[iw], s, a, rho_mc, rho_do,
                         100 * (rho_mc / rho_do - 1) if rho_do else np.nan))
        print()

    st = diff_stats(np.array(a_ours), np.array(a_do))
    # albedo-0 subset (faint atmospheric path, where solver differences show most)
    print('  overall (%d runs): mean bias=%+.3e  rms=%.3e  rel rms=%.2f%%  corr=%.5f'
          % (st['n'], st['mean_bias'], st['rms_diff'],
             100 * st.get('rel_rms', float('nan')), st['corr']))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('our_h5')
    p.add_argument('--band', default='o2a')
    p.add_argument('--streams', type=int, default=16)
    p.add_argument('--n-wvl', type=int, default=6,
                   help='number of window wavelengths across the band (default 6)')
    p.add_argument('--ot-max', type=float, default=0.001,
                   help='max TOTAL gas (O2+H2O) column OT for a window (default 0.001)')
    p.add_argument('--workdir', default=None)
    args = p.parse_args()
    if not os.path.isfile(args.our_h5):
        sys.exit('Error [eval_lrt]: no such file: %s' % args.our_h5)
    if not os.path.isfile('%s/bin/uvspec' % LRT):
        sys.exit('Error [eval_lrt]: uvspec not found under LIBRADTRAN_V2_DIR=%s' % LRT)
    work = args.workdir or os.path.join(
        os.environ.get('O2BAND_OUT_DIR', os.path.join(_HERE, '..', 'out')), 'lrt_eval')
    run(args.our_h5, args.band, args.streams, work, args.n_wvl, args.ot_max)
