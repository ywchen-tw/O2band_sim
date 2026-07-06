#!/usr/bin/env python
"""
Independent RT-solver cross-check: MCARaTS (ours) vs libRadtran/DISORT (EVAL_PLAN.md E2).

Our TOA reflectance comes from a Monte-Carlo solver (MCARaTS).  This runs the same
clear-sky scene through libRadtran's discrete-ordinate solver (DISORT) and compares
the reflectance -- validating the RT *transport* and the reflectance convention
(rho = pi*I/(mu0*F0)) independently of the absorption physics.

To keep it a solver check (not an absorption-model check), the comparison is done
at the band's **window wavelength** (grid point of minimum O2 optical depth) with
libRadtran molecular absorption switched OFF (``no_absorption``): a pure
Rayleigh + Lambertian scene, which our MCARaTS window point also approximates
(O2 OT ~ 0).  uvspec is run with ``output_quantity reflectivity`` so its ``uu``
output IS rho, matching ours.

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


def run(our_h5, band, streams, workdir):
    os.makedirs(workdir, exist_ok=True)
    atm_file = '%s/data/atmmod/afglms.dat' % LRT

    with h5py.File(our_h5, 'r') as f:
        g = f[band] if band in f else f
        wvl = g['wvl'][:]; sza = g['sza'][:]; alb = g['albedo'][:]
        ref = g['reflectance'][:]
        o2c = g['optical_thickness/o2_column'][:]

    iw = int(np.argmin(o2c)); wl = float(wvl[iw])
    print('RT-solver check: MCARaTS vs libRadtran/DISORT (%d streams)' % streams)
    print('  window %.4f nm  (O2 column OT = %.4f -> near pure-Rayleigh)\n' % (wl, o2c[iw]))
    print('  %-5s %-6s | %12s %12s %10s' % ('SZA', 'alb', 'MCARaTS', 'DISORT', 'diff'))

    a_ours, a_do = [], []
    for i, s in enumerate(sza):
        for j, a in enumerate(alb):
            rho_mc = float(ref[i, j, iw])
            rho_do = disort_reflectivity(wl, float(s), float(a), atm_file, streams, workdir)
            a_ours.append(rho_mc); a_do.append(rho_do)
            print('  %-5.1f %-6.2f | %12.5f %12.5f %+9.2f%%'
                  % (s, a, rho_mc, rho_do, 100 * (rho_mc / rho_do - 1) if rho_do else np.nan))

    st = diff_stats(np.array(a_ours), np.array(a_do))
    print('\n  overall: n=%d  mean bias=%+.3e  rms=%.3e  rel rms=%.2f%%  corr=%.5f'
          % (st['n'], st['mean_bias'], st['rms_diff'],
             100 * st.get('rel_rms', float('nan')), st['corr']))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('our_h5')
    p.add_argument('--band', default='o2a')
    p.add_argument('--streams', type=int, default=16)
    p.add_argument('--workdir', default=None)
    args = p.parse_args()
    if not os.path.isfile(args.our_h5):
        sys.exit('Error [eval_lrt]: no such file: %s' % args.our_h5)
    if not os.path.isfile('%s/bin/uvspec' % LRT):
        sys.exit('Error [eval_lrt]: uvspec not found under LIBRADTRAN_V2_DIR=%s' % LRT)
    work = args.workdir or os.path.join(
        os.environ.get('O2BAND_OUT_DIR', os.path.join(_HERE, '..', 'out')), 'lrt_eval')
    run(args.our_h5, args.band, args.streams, work)
