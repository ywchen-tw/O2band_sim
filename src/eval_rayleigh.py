#!/usr/bin/env python
"""
Rayleigh optical-thickness cross-check: our Bodhaine (1999) implementation vs the
independent Bucholtz (1995) parameterization (EVAL_PLAN.md Q1).

The Rayleigh cross-section is layer-independent (a function of wavelength only), so
column Rayleigh OT = total_air_column * sigma(lambda).  The *relative* difference in
column OT therefore equals the relative difference in the cross-section exactly;
this script reports both the cross-section difference (the fundamental quantity)
and the absolute column OT (for context, using the AFGL atmosphere).

References
- Bodhaine, B.A., et al. (1999), J. Atmos. Oceanic Technol. 16, 1854 -- the
  cross-section our cal_rayleigh_od uses (360 ppm CO2 air).
- Bucholtz, A. (1995), Appl. Opt. 34(15), 2765 -- independent parameterization
  sigma = A * lambda^-(B + C*lambda + D/lambda), lambda in um, sigma in cm^2.

    python src/eval_rayleigh.py           # both bands, table of difference stats
"""

import os
import sys
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from util.absorption import BANDS
from eval_metrics import diff_stats, print_diff_stats


# ---------------------------------------------------------------------------- #
def rayleigh_xsec_bodhaine(wvl_nm):
    """Bodhaine et al. (1999) Rayleigh cross-section (cm^2), 360 ppm CO2 air.

    Identical closed form to util.absorption.cal_rayleigh_od (which multiplies
    this by the layer air column)."""
    wl_um = np.asarray(wvl_nm, dtype=np.float64) * 1.0e-3
    inv2 = 1.0 / wl_um ** 2
    num = 1.0455996 - 341.29061 * inv2 - 0.90230850 * wl_um ** 2
    den = 1.0 + 0.0027059889 * inv2 - 85.968563 * wl_um ** 2
    return 1.0e-28 * num / den


def rayleigh_xsec_bucholtz(wvl_nm):
    """Bucholtz (1995) Rayleigh cross-section (cm^2).

    sigma = A * lambda^-(B + C*lambda + D/lambda), lambda in um.  Two coefficient
    sets split at 0.5 um; the O2 A/B bands (0.68-0.77 um) use the >0.5 um set, but
    both are applied per-point so the function is general."""
    wl_um = np.asarray(wvl_nm, dtype=np.float64) * 1.0e-3
    # (A, B, C, D) for lambda <= 0.5 um and lambda > 0.5 um  (Bucholtz 1995, Table 3)
    lo = (3.01577e-28, 3.55212, 1.35579, 0.11563)
    hi = (4.01061e-28, 3.99668, 1.10298e-3, 2.71393e-2)
    A = np.where(wl_um <= 0.5, lo[0], hi[0])
    B = np.where(wl_um <= 0.5, lo[1], hi[1])
    C = np.where(wl_um <= 0.5, lo[2], hi[2])
    D = np.where(wl_um <= 0.5, lo[3], hi[3])
    return A * wl_um ** (-(B + C * wl_um + D / wl_um))


# ---------------------------------------------------------------------------- #
def compare(bands=('o2a', 'o2b'), dwvl=0.001, air_column=None):
    """Print Bodhaine-vs-Bucholtz difference statistics per band."""
    print('Rayleigh cross-section: Bodhaine (1999, ours) vs Bucholtz (1995)')
    if air_column is not None:
        print('  absolute column OT uses total air column = %.4e molec cm^-2' % air_column)
    for band in bands:
        wl0, wl1 = BANDS[band]
        wvl = np.arange(wl0, wl1 + 0.5 * dwvl, dwvl)
        xb = rayleigh_xsec_bodhaine(wvl)     # ours
        xu = rayleigh_xsec_bucholtz(wvl)     # reference
        st = diff_stats(xb, xu)
        mid = 0.5 * (wl0 + wl1)
        print('\n[%s] %.0f-%.0f nm  (mid %.0f nm)' % (band, wl0, wl1, mid))
        print('  sigma @ %.0f nm: Bodhaine=%.4e  Bucholtz=%.4e cm^2  (Bodhaine %+.2f%%)'
              % (mid, rayleigh_xsec_bodhaine(mid), rayleigh_xsec_bucholtz(mid),
                 100.0 * (rayleigh_xsec_bodhaine(mid) / rayleigh_xsec_bucholtz(mid) - 1.0)))
        print_diff_stats('xsec (Bodhaine-Bucholtz)', st)
        if air_column is not None:
            colb, colu = xb * air_column, xu * air_column
            print('  column OT @ band mean: Bodhaine=%.4f  Bucholtz=%.4f'
                  % (colb.mean(), colu.mean()))


if __name__ == '__main__':
    air_col = None
    # use the real AFGL air column for absolute column OT context, if available
    try:
        from util.atmosphere import afgl_atmosphere
        data_dir = os.environ.get('O2BAND_DATA_DIR',
                                  os.path.normpath(os.path.join(_HERE, '..', 'data')))
        fafgl = os.path.join(data_dir, 'afglms.dat')
        if os.path.isfile(fafgl):
            z_top = float(os.environ.get('O2BAND_ZTOP', '120'))
            atm = afgl_atmosphere(fafgl, z_top=z_top)
            air_col = float(atm.lay['air'].sum())
    except Exception as e:
        print('(no AFGL atmosphere for absolute column OT: %s)' % e)
    compare(air_column=air_col)
