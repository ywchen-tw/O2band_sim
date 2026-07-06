#!/usr/bin/env python
"""
In-band RT-solver cross-check: MCARaTS vs libRadtran/DISORT WITH gas absorption
(EVAL_PLAN.md E2, companion to eval_lrt.py).

eval_lrt.py tests window wavelengths (pure Rayleigh) -- that validates transport +
convention but not the solver's handling of absorption+scattering *coupling*,
which is the whole point of the O2 bands.  This runs the comparison at wavelengths
spanning a range of gas optical depth by injecting OUR per-layer gas absorption
OT into DISORT, so both solvers see the identical optical properties.

Method
------
For a chosen wavelength, our per-layer gas absorption OT (o2_layer + h2o_layer,
from the output HDF5) is written as a libRadtran ``aerosol_file tau`` profile and
made a **pure absorber** with ``aerosol_modify ssa set 0`` (asymmetry 0).
libRadtran's own molecular absorption is switched off (``no_absorption mol``) but
its Rayleigh scattering is kept (Bodhaine, matches ours to <0.1%, per eval #1).
So DISORT solves Rayleigh + our gas absorption -- the same scene as MCARaTS -- and
``output_quantity reflectivity`` makes its ``uu`` the reflectance rho.

The atmosphere_file is OUR afglms.dat so the vertical grid matches the injected
per-layer OT.  Wavelengths are chosen (0.1 nm-aligned, for the solar file) to span
target column gas ODs; note that above OD ~ few the reflectance is small and
MC-noise dominated, so the informative range is OD ~ 0.05-3.

Run via curc_lrt_eval.sh with INBAND=1 (needs uvspec's module environment).

    python src/eval_lrt_inband.py OUR.h5 --band o2a
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
DEFAULT_ODS = [0.05, 0.1, 0.3, 0.7, 1.5, 3.0]


def select_by_od(wvl, col_ot, targets, solar_nm=0.1):
    """0.1 nm-aligned grid indices whose column gas OT is closest (in log) to each
    target OD -- a spread of absorption levels."""
    aligned = np.where(np.abs(wvl - np.round(wvl / solar_nm) * solar_nm) < 1e-4)[0]
    logot = np.log(np.clip(col_ot[aligned], 1e-6, None))
    idx = []
    for t in targets:
        k = aligned[int(np.argmin(np.abs(logot - np.log(t))))]
        idx.append(int(k))
    return sorted(set(idx))


def write_aerosol_tau(path, z_lev_km, layer_ot):
    """aerosol_file tau: decreasing altitude; tau at a level = OD of the layer
    above it; top level = 0.  z_lev ascending (surface..top), layer_ot has
    len(z_lev)-1 entries (layer j between lev[j] and lev[j+1])."""
    nlev = z_lev_km.size
    with open(path, 'w') as f:
        f.write('# z[km]  tau\n')
        for m in range(nlev - 1, -1, -1):           # top -> surface
            tau = 0.0 if m == nlev - 1 else float(layer_ot[m])
            f.write('%10.4f  %.6e\n' % (z_lev_km[m], tau))


def disort_inband(wl, sza, alb, atm_file, z_lev, layer_ot, streams, workdir, tag):
    taufile = os.path.join(workdir, 'aer_tau_%s.dat' % tag)
    write_aerosol_tau(taufile, z_lev, layer_ot)
    inp = '\n'.join([
        'data_files_path %s/data' % LRT,
        'atmosphere_file %s' % atm_file,
        'source solar %s/data/solar_flux/kurudz_0.1nm.dat' % LRT,
        'wavelength %.4f %.4f' % (wl, wl),
        'rte_solver disort',
        'number_of_streams %d' % streams,
        'no_absorption mol',                  # libRadtran gas off; ours injected below
        'aerosol_default',
        'aerosol_file tau %s' % taufile,      # our per-layer gas absorption OT
        'aerosol_modify ssa set 0.0',         # pure absorber
        'aerosol_modify gg set 0.0',
        'albedo %.4f' % alb,
        'sza %.4f' % sza,
        'phi0 0.0', 'umu 1.0', 'phi 0.0',
        'zout toa', 'output_user lambda uu', 'output_quantity reflectivity',
        'quiet', ''])
    with open(os.path.join(workdir, 'lrt_in_%s.inp' % tag), 'w') as f:
        f.write(inp)
    r = subprocess.run(['%s/bin/uvspec' % LRT], input=inp, capture_output=True,
                       text=True, cwd=workdir)
    if r.returncode != 0:
        raise RuntimeError('uvspec failed (%s):\n%s' % (tag, r.stderr[-800:]))
    v = r.stdout.split()
    if len(v) < 2:
        raise RuntimeError('bad uvspec output (%s): %r\n%s' % (tag, r.stdout, r.stderr[-400:]))
    return float(v[1])


def run(our_h5, band, streams, workdir, targets):
    os.makedirs(workdir, exist_ok=True)
    atm_file = os.path.join(os.environ.get('O2BAND_DATA_DIR',
                            os.path.join(_HERE, '..', 'data')), 'afglms.dat')

    with h5py.File(our_h5, 'r') as f:
        g = f[band] if (band in f and isinstance(f[band], h5py.Group)) else f
        if g is f and str(f.attrs.get('band', band)) != band:
            sys.exit('band mismatch: --band %s vs file band %s' % (band, f.attrs.get('band')))
        wvl = g['wvl'][:]; sza = g['sza'][:]; alb = g['albedo'][:]
        ref = g['reflectance'][:]; rerr = g['reflectance_stderr'][:]
        o2l = g['optical_thickness/o2_layer'][:]; h2ol = g['optical_thickness/h2o_layer'][:]
        col = o2l.sum(0) + h2ol.sum(0)
        z_lev = g['atmosphere/z_lev_km'][:]

    iws = select_by_od(wvl, col, targets)
    print('In-band RT-solver check: MCARaTS vs DISORT with injected gas OT (%d streams)' % streams)
    print('  %d wavelengths spanning column gas OD, band %s\n' % (len(iws), band))
    print('  %-9s %-9s %-5s %-6s | %11s %11s %9s %9s'
          % ('wvl[nm]', 'col_OT', 'SZA', 'alb', 'MCARaTS', 'DISORT', 'diff', 'MCnoise'))

    a_o, a_d = [], []
    for iw in iws:
        wl = float(wvl[iw]); layer_ot = o2l[:, iw] + h2ol[:, iw]
        for i, s in enumerate(sza):
            for j, a in enumerate(alb):
                rmc = float(ref[i, j, iw]); se = float(rerr[i, j, iw])
                tag = 'w%.1f_s%02.0f_a%.2f' % (wl, s, a)
                rdo = disort_inband(wl, float(s), float(a), atm_file, z_lev, layer_ot,
                                    streams, workdir, tag)
                a_o.append(rmc); a_d.append(rdo)
                print('  %-9.3f %-9.4f %-5.1f %-6.2f | %11.5f %11.5f %+8.2f%% %8.1f%%'
                      % (wl, col[iw], s, a, rmc, rdo,
                         100 * (rmc / rdo - 1) if rdo else np.nan,
                         100 * se / rmc if rmc else np.nan))
        print()

    st = diff_stats(np.array(a_o), np.array(a_d))
    print('  overall (%d runs): mean bias=%+.3e  rms=%.3e  rel rms=%.2f%%  corr=%.5f'
          % (st['n'], st['mean_bias'], st['rms_diff'],
             100 * st.get('rel_rms', float('nan')), st['corr']))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('our_h5')
    p.add_argument('--band', default='o2a')
    p.add_argument('--streams', type=int, default=16)
    p.add_argument('--ods', type=float, nargs='+', default=DEFAULT_ODS,
                   help='target column gas ODs to sample (default %s)' % DEFAULT_ODS)
    p.add_argument('--workdir', default=None)
    args = p.parse_args()
    if not os.path.isfile(args.our_h5):
        sys.exit('Error: no such file: %s' % args.our_h5)
    if not os.path.isfile('%s/bin/uvspec' % LRT):
        sys.exit('Error: uvspec not found under LIBRADTRAN_V2_DIR=%s' % LRT)
    work = args.workdir or os.path.join(
        os.environ.get('O2BAND_OUT_DIR', os.path.join(_HERE, '..', 'out')), 'lrt_eval_inband')
    run(args.our_h5, args.band, args.streams, work, args.ods)
