#!/usr/bin/env python
"""
O2 (and H2O) line-by-line optical-depth cross-check vs HAPI (EVAL_PLAN.md #3).

HAPI (the HITRAN Application Programming Interface) is an *independent*
implementation of Voigt line-by-line absorption from the same HITRAN 2020 line
list, so it is the apples-to-apples check of our absorption engine: line
intensity S(T) (partition sums), the Voigt line shape, pressure broadening/shift,
and the wavenumber grid.  It validates the one physics component the Rayleigh
(Bucholtz/Hansen&Travis) and column (0.2095) checks do not touch.

Method
------
For each band and a few representative layers (surface: pressure-broadened;
upper: Doppler-dominated), compute the O2 absorption *cross-section* sigma(nu)
[cm^2/molecule] with HAPI at the layer's (p, T) and self/air mixing, on the SAME
vacuum-wavenumber grid our air-wavelength grid maps to, and diff_stats it against
our per-layer sigma = od['o2'][iz] / o2_column[iz].  Matching the cross-section
per layer validates the column OT (a linear sum over layers).

HAPI is fetched from HITRAN online (needs internet), so run this via
``curc_hapi_eval.sh`` on a login node.  Line-list version differences show up as
isolated single-line residuals, distinguishable from a systematic shape/intensity
bias; the documented wing-cutoff (ncut*nu0/R) makes our far wings slightly lower.

    python src/eval_hapi.py --bands o2a o2b --z-top 120
"""

import os
import sys
import argparse
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from util.atmosphere import afgl_atmosphere
from util.tips import tips2021
from util.absorption import hitran_lines, o2band_absorption, BANDS, MOL_O2
from util.optics import air_to_vac_nm
from eval_metrics import diff_stats, print_diff_stats


# HITRAN global isotopologue IDs for O2 (66, 68, 67) -- the isos in our file
O2_GLOBAL_IDS = [36, 37, 38]


def our_sigma_layer(absb, atm, iz, gas='o2'):
    """Our per-layer cross-section sigma(nu) on an ascending vacuum-wavenumber grid.

    Returns (nu_asc, sigma_asc) with sigma = od[gas][iz] / column[gas][iz]."""
    col = float(atm.lay[gas][iz])
    sigma = absb.od[gas][iz] / col                      # (nwvl,), air-wvl order
    nu = absb.nu_vac                                    # cm-1, air-wvl order
    order = np.argsort(nu)
    return nu[order], sigma[order]


def hapi_sigma_layer(hapi, table, nu_grid, p_hpa, T, vmr):
    """HAPI cross-section sigma(nu) [cm^2/molecule] at (p, T, self-fraction vmr)."""
    p_atm = p_hpa / 1013.25
    nu_out, coef = hapi.absorptionCoefficient_Voigt(
        SourceTables=table,
        Environment={'p': p_atm, 'T': float(T)},
        Diluent={'air': 1.0 - vmr, 'self': vmr},
        WavenumberGrid=np.ascontiguousarray(nu_grid, dtype=np.float64),
        HITRAN_units=True)                              # cm^2/molecule
    return np.asarray(nu_out), np.asarray(coef)


def run(bands, z_top, cache_dir, layers=None):
    import hapi

    data_dir = os.environ.get('O2BAND_DATA_DIR',
                              os.path.normpath(os.path.join(_HERE, '..', 'data')))
    atm = afgl_atmosphere(os.path.join(data_dir, 'afglms.dat'), z_top=z_top)
    tips = tips2021()
    fhit = os.path.join(data_dir, 'hitran2020_lines.txt')

    os.makedirs(cache_dir, exist_ok=True)
    hapi.db_begin(cache_dir)

    nlay = atm.lay['z'].size
    if layers is None:
        # surface (Lorentz) and a mid/upper layer (Doppler), by pressure
        p = atm.lay['p']
        layers = [int(np.argmax(p)), int(np.argmin(np.abs(p - 100.0)))]

    for band in bands:
        wl0, wl1 = BANDS[band]
        # our absorption (O2 only here; matches HAPI SourceTables='O2')
        lines = hitran_lines(fhit, wl_range=(wl0, wl1), margin_cm=5.0)
        absb = o2band_absorption(atm, lines, band=band, include_h2o=False, tips=tips)

        # vacuum-wavenumber span for the HAPI fetch (with margin)
        wl_vac = air_to_vac_nm(np.array([wl0, wl1]))
        nu_lo, nu_hi = 1.0e7 / wl_vac.max() - 5.0, 1.0e7 / wl_vac.min() + 5.0
        table = 'O2_%s' % band
        if not os.path.isfile(os.path.join(cache_dir, table + '.data')):
            print('[hapi] fetching %s O2 lines %.1f-%.1f cm-1 ...' % (band, nu_lo, nu_hi))
            hapi.fetch_by_ids(table, O2_GLOBAL_IDS, nu_lo, nu_hi)

        print('\n================ %s : O2 cross-section, ours vs HAPI ================' % band)
        for iz in layers:
            nu, sig_ours = our_sigma_layer(absb, atm, iz, 'o2')
            _, sig_hapi = hapi_sigma_layer(hapi, table, nu,
                                           atm.lay['p'][iz], atm.lay['temperature'][iz],
                                           float(atm.lay['o2_vmr'][iz]))
            st = diff_stats(sig_ours, sig_hapi)
            print('\n layer %d: p=%.1f hPa T=%.1f K  (sigma peak ours=%.3e hapi=%.3e cm^2)'
                  % (iz, atm.lay['p'][iz], atm.lay['temperature'][iz],
                     sig_ours.max(), sig_hapi.max()))
            print_diff_stats('O2 sigma (ours-HAPI)', st)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--bands', nargs='+', default=['o2a', 'o2b'])
    p.add_argument('--z-top', type=float, default=120.0)
    p.add_argument('--cache-dir', default=None,
                   help='HAPI line cache dir (default: $O2BAND_OUT_DIR/hapi_cache)')
    args = p.parse_args()

    try:
        import hapi  # noqa: F401
    except ImportError:
        sys.exit('Error [eval_hapi]: HAPI not installed. Run curc_hapi_eval.sh '
                 '(pip install hitran-api) on a node with internet.')

    out = os.environ.get('O2BAND_OUT_DIR', os.path.join(_HERE, '..', 'out'))
    cache = args.cache_dir or os.path.join(out, 'hapi_cache')
    run(args.bands, args.z_top, cache)
