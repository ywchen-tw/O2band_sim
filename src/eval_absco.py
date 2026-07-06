#!/usr/bin/env python
"""
O2 A-band absorption vs the OCO **ABSCO** tables (EVAL_PLAN.md #, tier R2).

ABSCO (OCO-2/3 Absorption Coefficient tables, v5.2, Drouin O2 line list) is the
*operational* O2 A-band cross-section.  Crucially it includes **line mixing,
collision-induced absorption (CIA), and a speed-dependent line shape** -- all of
which Phase-1 deliberately EXCLUDES (PLAN.md sec.7.4).  So this is NOT a code
check (unlike HAPI); it *quantifies our deliberate physical omissions*.  Expect:

  - a nonzero **between-line continuum** in ABSCO (CIA + line-mixing far wings)
    where our plain-Voigt sigma is near zero -> ABSCO >> ours between lines;
  - **line-core / wing** redistribution from line mixing;
  - a ~1% intensity offset (Drouin vs our HITRAN 2020).

ABSCO here covers only the A-band (12745-13245 cm-1); O2B is not in this table.

    python src/eval_absco.py --z-top 120
"""

import os
import sys
import argparse
import numpy as np
import h5py

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from util.atmosphere import afgl_atmosphere
from util.tips import tips2021
from util.absorption import hitran_lines, o2band_absorption, BANDS
from util.optics import air_to_vac_nm
from eval_metrics import diff_stats, print_diff_stats

ABSCO_O2 = '/scratch/alpine/yuch8913/oco/data/absco/v5.2_final/o2_v52.hdf'


class Absco:
    """Lazy reader + (P,T,H2O)-interpolator for the ABSCO O2 cross-section."""

    def __init__(self, path=ABSCO_O2):
        self.f = h5py.File(path, 'r')
        self.P = self.f['Pressure'][:]                 # (64,) Pa, ascending
        self.T = self.f['Temperature'][:]              # (64,17) K per P level
        self.B = self.f['Broadener_01_VMR'][:]         # (3,) H2O VMR
        self.wn = self.f['Wavenumber'][:]              # (50001,) cm-1
        self.sig = self.f['Gas_07_Absorption']         # (64,17,3,50001) cm^2, lazy

    def sigma(self, p_hpa, T_K, h2o_vmr, nu_lo, nu_hi):
        """Interpolated O2 cross-section sigma(nu) over [nu_lo, nu_hi] cm-1."""
        p_pa = p_hpa * 100.0
        w0, w1 = np.searchsorted(self.wn, [nu_lo, nu_hi])
        w0, w1 = max(0, w0 - 1), min(self.wn.size, w1 + 1)
        wn = self.wn[w0:w1]

        ip = int(np.clip(np.searchsorted(self.P, p_pa) - 1, 0, self.P.size - 2))
        # interpolate T within each bracketing P level, then log-P, then VMR
        planes = []
        for k in (ip, ip + 1):
            s = self.sig[k, :, :, w0:w1].astype(np.float64)     # (17,3,nw)
            Tg = self.T[k]                                      # (17,)
            o = np.argsort(Tg)
            s = np.stack([np.array([np.interp(T_K, Tg[o], s[o, b, j])
                                    for j in range(s.shape[2])])
                          for b in range(s.shape[1])])          # (3, nw)
            planes.append(s)
        lp = np.log(self.P[ip:ip + 2])
        wgt = 0.0 if lp[1] == lp[0] else (np.log(p_pa) - lp[0]) / (lp[1] - lp[0])
        wgt = float(np.clip(wgt, 0.0, 1.0))
        s_pt = (1 - wgt) * planes[0] + wgt * planes[1]          # (3, nw) over VMR
        sig = np.array([np.interp(h2o_vmr, self.B, s_pt[:, j]) for j in range(s_pt.shape[1])])
        return wn, sig


def run(z_top, layers=None):
    data_dir = os.environ.get('O2BAND_DATA_DIR',
                              os.path.normpath(os.path.join(_HERE, '..', 'data')))
    atm = afgl_atmosphere(os.path.join(data_dir, 'afglms.dat'), z_top=z_top)
    tips = tips2021()
    lines = hitran_lines(os.path.join(data_dir, 'hitran2020_lines.txt'),
                         wl_range=BANDS['o2a'], margin_cm=5.0)
    absb = o2band_absorption(atm, lines, band='o2a', include_h2o=False, tips=tips)

    nu = absb.nu_vac
    order = np.argsort(nu)
    nu_asc = nu[order]

    absco = Absco()
    if layers is None:
        p = atm.lay['p']
        layers = [int(np.argmax(p)), int(np.argmin(np.abs(p - 300.0)))]

    print('O2 A-band: ours (Voigt, HITRAN2020, no CIA/line-mixing) vs ABSCO v5.2 '
          '(Drouin, +CIA +line-mixing)')
    for iz in layers:
        p, T = atm.lay['p'][iz], atm.lay['temperature'][iz]
        vmr = float(atm.lay['h2o_vmr'][iz])
        sig_ours = (absb.od['o2'][iz] / atm.lay['o2'][iz])[order]
        wn_a, sig_a = absco.sigma(p, T, vmr, nu_asc.min(), nu_asc.max())
        sig_absco = np.interp(nu_asc, wn_a, sig_a)             # onto our grid

        # line mask (near a core) vs between-line, from OUR sigma
        pk = sig_ours.max()
        core = sig_ours > 0.1 * pk
        between = sig_ours < 1e-3 * pk
        print('\n layer %d: p=%.1f hPa T=%.1f K H2O_vmr=%.4f' % (iz, p, T, vmr))
        print('  peak sigma : ours=%.4e  ABSCO=%.4e  (ABSCO %+.1f%%)'
              % (pk, sig_absco.max(), 100 * (sig_absco.max() / pk - 1)))
        print('  between-line median sigma: ours=%.3e  ABSCO=%.3e  (ABSCO/ours = %.1fx)  '
              '<- CIA + line-mixing continuum'
              % (np.median(sig_ours[between]), np.median(sig_absco[between]),
                 np.median(sig_absco[between]) / max(np.median(sig_ours[between]), 1e-30)))
        print_diff_stats('sigma at line cores (>10%% peak)',
                         diff_stats(sig_ours[core], sig_absco[core]))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--z-top', type=float, default=120.0)
    args = p.parse_args()
    if not os.path.isfile(ABSCO_O2):
        sys.exit('Error [eval_absco]: ABSCO O2 table not found: %s' % ABSCO_O2)
    run(args.z_top)
