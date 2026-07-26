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
from util.absorption import (hitran_lines, o2band_absorption, required_margin_cm,
                             band_nu_range, BANDS, MOL_O2, DEFAULT_CUTOFF_CM)
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


def hapi_sigma_layer(hapi, table, nu_grid, p_hpa, T, vmr,
                     wing_cm=DEFAULT_CUTOFF_CM):
    """HAPI cross-section sigma(nu) [cm^2/molecule] at (p, T, self-fraction vmr).

    The wing cutoff MUST match ours or this is not a like-for-like comparison.
    HAPI's default is 50 half-widths -- about 2 cm-1 at surface pressure --
    against our 50 cm-1, and the far wing is precisely what the 2026-07-25
    finding was about, so leaving it at the default would manufacture a
    disagreement and blame it on our code.

    HAPI renamed OmegaWing/OmegaWingHW to WavenumberWing/WavenumberWingHW; try
    the modern spelling first and fall back, since which one a given hitran-api
    release accepts is not something we can assume.  WingHW=0 makes the absolute
    wing govern rather than a multiple of the half-width.
    """
    p_atm = p_hpa / 1013.25
    kw = dict(SourceTables=table,
              Environment={'p': p_atm, 'T': float(T)},
              Diluent={'air': 1.0 - vmr, 'self': vmr},
              WavenumberGrid=np.ascontiguousarray(nu_grid, dtype=np.float64),
              HITRAN_units=True)                        # cm^2/molecule
    for wing in ({'WavenumberWing': wing_cm, 'WavenumberWingHW': 0.0},
                 {'OmegaWing': wing_cm, 'OmegaWingHW': 0.0}):
        try:
            nu_out, coef = hapi.absorptionCoefficient_Voigt(**dict(kw, **wing))
        except TypeError:
            continue
        return np.asarray(nu_out), np.asarray(coef)
    raise TypeError('hapi.absorptionCoefficient_Voigt accepts neither '
                    'WavenumberWing nor OmegaWing, so our %g cm-1 wing cutoff '
                    'cannot be matched; comparison would be invalid' % wing_cm)


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
        # margin must track the wing cutoff (PLAN.md §7.3), not be hardcoded
        lines = hitran_lines(fhit, wl_range=(wl0, wl1),
                             margin_cm=required_margin_cm((wl0, wl1)))
        absb = o2band_absorption(atm, lines, band=band, include_h2o=False, tips=tips)

        # Vacuum-wavenumber span for the HAPI fetch.  The margin must cover our
        # wing cutoff, else HAPI's line list omits the out-of-window lines whose
        # wings reach into the band and the two codes are not comparable.
        margin = required_margin_cm((wl0, wl1))
        nu_lo, nu_hi = band_nu_range((wl0, wl1), grid='vac')
        nu_lo, nu_hi = nu_lo - margin, nu_hi + margin
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
