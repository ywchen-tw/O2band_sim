"""
Cross-check the Toon SPTS solar transmittance against the OCO L2 solar model
(EVAL_REPORT sec.7).

The OCO retrieval's solar model (``l2_solar_model.h5``) stores a solar
absorption (transmittance-like) spectrum per band, built from an earlier
release of the same Toon solar line list that generates the SPTS.  Its Band 1
(12700-13300 cm-1, 0.001 cm-1 grid) overlaps our O2 A-band window, giving an
independent-file check of the SPTS content and of our vacuum-wavenumber
handling.  There is no OCO band covering O2B.

Usage:
    python src/eval_spts_oco.py [--oco FILE] [--spts FILE]
"""

import os
import argparse

import numpy as np
import h5py

_HERE = os.path.dirname(os.path.abspath(__file__))
import sys
sys.path.insert(0, _HERE)

from util.solar import solar_spts                             # noqa: E402

_OCO_DEFAULT = '/scratch/alpine/yuch8913/oco/data/absco/v5.2_final/l2_solar_model.h5'

# O2 A-band air window 757-772 nm in vacuum wavenumber
_WIN = (12953.0, 13210.0)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    data_dir = os.environ.get('O2BAND_DATA_DIR',
                              os.path.join(_HERE, '..', 'data'))
    p.add_argument('--oco', default=_OCO_DEFAULT, help='OCO l2_solar_model.h5')
    p.add_argument('--spts', default=None,
                   help='SPTS solar_merged_*.out (default: newest in data dir)')
    args = p.parse_args(argv)

    fname_spts = args.spts
    if fname_spts is None:
        import glob
        hits = sorted(glob.glob(os.path.join(data_dir, 'solar_merged_*.out')))
        if not hits:
            raise OSError('no solar_merged_*.out under %s' % data_dir)
        fname_spts = hits[-1]

    spts = solar_spts(os.path.join(data_dir, 'CU_composite_solar.dat'),
                      fname_spts)

    with h5py.File(args.oco, 'r') as f:
        nu_oco = f['Solar/Absorption/Absorption_1/wavenumber'][:]
        T_oco = f['Solar/Absorption/Absorption_1/spectrum'][:]

    m = (nu_oco >= _WIN[0]) & (nu_oco <= _WIN[1])
    nu, To = nu_oco[m], T_oco[m]
    Ts = np.interp(nu, spts.nu, spts.trans)     # both vacuum wavenumber

    d = Ts - To
    imax = int(np.abs(d).argmax())
    print('SPTS: %s' % os.path.basename(fname_spts))
    print('OCO : %s (Band 1, %.3f-%.3f cm-1)' % (args.oco, nu_oco[0], nu_oco[-1]))
    print('O2A overlap window %.1f-%.1f cm-1: %d points (OCO 0.001 cm-1 grid)'
          % (_WIN[0], _WIN[1], nu.size))
    print('  mean T   : SPTS %.5f | OCO %.5f' % (Ts.mean(), To.mean()))
    print('  bias     : %+.2e' % d.mean())
    print('  RMS diff : %.3e' % np.sqrt((d ** 2).mean()))
    print('  max|diff|: %.3e at %.3f cm-1 (%.3f nm vac; sharp line core --'
          % (np.abs(d).max(), nu[imax], 1.0e7 / nu[imax]))
    print('             SPTS grid is 10x coarser than OCO Band 1)')
    print('  corr     : %.6f' % np.corrcoef(Ts, To)[0, 1])
    for thr in (1e-3, 5e-3, 1e-2):
        print('  fraction |diff| > %.0e : %.4f' % (thr, (np.abs(d) > thr).mean()))
    k = int(np.argmin(To))
    print('  deepest line (K I 766.5 nm): nu=%.3f cm-1  T_oco=%.4f  T_spts=%.4f'
          % (nu[k], To[k], Ts[k]))


if __name__ == '__main__':
    main()
