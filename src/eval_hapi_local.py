#!/usr/bin/env python
"""
Version-matched HAPI cross-check (EVAL_PLAN.md #3, companion to eval_hapi.py).

eval_hapi.py compares against HAPI's *online* line list, which is now HITRAN2024
-- so its O2 A-band intensities are ~1.3% above our prescribed HITRAN 2020 file,
and that line-list difference dominates the ~1% discrepancy seen there.

This test removes the version difference: it registers *our local HITRAN 2020*
``.par`` as a HAPI table (no internet) and runs HAPI's Voigt LBL on the SAME line
data.  What remains is the pure implementation difference (Voigt evaluation,
partition function, broadening) -- expected to be well below the 1.3% line-list
gap, cleanly validating our absorption engine independent of HITRAN edition.

    python src/eval_hapi_local.py --bands o2a o2b --z-top 120
"""

import os
import sys
import json
import argparse
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from util.atmosphere import afgl_atmosphere
from util.tips import tips2021
from util.absorption import hitran_lines, o2band_absorption, BANDS, MOL_O2
from util.optics import air_to_vac_nm
from eval_metrics import diff_stats, print_diff_stats
from eval_hapi import our_sigma_layer, hapi_sigma_layer

# Standard HITRAN 160-char .par HAPI header (field layout is content-independent,
# so it applies verbatim to our 2020 file).
_PAR_HEADER = {
    "table_type": "column-fixed", "size_in_bytes": -1,
    "extra": [], "extra_format": {}, "extra_separator": ",",
    "order": ["molec_id", "local_iso_id", "nu", "sw", "a", "gamma_air",
              "gamma_self", "elower", "n_air", "delta_air", "global_upper_quanta",
              "global_lower_quanta", "local_upper_quanta", "local_lower_quanta",
              "ierr", "iref", "line_mixing_flag", "gp", "gpp"],
    "format": {"molec_id": "%2d", "local_iso_id": "%1d", "nu": "%12.6f",
               "sw": "%10.3E", "a": "%10.3E", "gamma_air": "%5.4f",
               "gamma_self": "%5.3f", "elower": "%10.4f", "n_air": "%4.2f",
               "delta_air": "%8.6f", "global_upper_quanta": "%15s",
               "global_lower_quanta": "%15s", "local_upper_quanta": "%15s",
               "local_lower_quanta": "%15s", "ierr": "%6s", "iref": "%12s",
               "line_mixing_flag": "%1s", "gp": "%7.1f", "gpp": "%7.1f"},
    "default": {"molec_id": 0, "local_iso_id": 0, "nu": 0.0, "sw": 0.0, "a": 0.0,
                "gamma_air": 0.0, "gamma_self": 0.0, "elower": 0.0, "n_air": 0.0,
                "delta_air": 0.0, "global_upper_quanta": "000",
                "global_lower_quanta": "000", "local_upper_quanta": "000",
                "local_lower_quanta": "000", "ierr": "EEE", "iref": "EEE",
                "line_mixing_flag": "EEE", "gp": "FFF", "gpp": "FFF"},
    "position": {"molec_id": 0, "local_iso_id": 2, "nu": 3, "sw": 15, "a": 25,
                 "gamma_air": 35, "gamma_self": 40, "elower": 45, "n_air": 55,
                 "delta_air": 59, "global_upper_quanta": 67,
                 "global_lower_quanta": 82, "local_upper_quanta": 97,
                 "local_lower_quanta": 112, "ierr": 127, "iref": 133,
                 "line_mixing_flag": 145, "gp": 146, "gpp": 153},
}


def register_local_table(fhitran, band, cache_dir):
    """Write the O2 lines of `band` from our HITRAN 2020 .par into a HAPI table
    (verbatim 160-char lines, so intensities/positions are exactly ours)."""
    wl0, wl1 = BANDS[band]
    wl_vac = air_to_vac_nm(np.array([wl0, wl1]))
    nu_lo, nu_hi = 1.0e7 / wl_vac.max() - 5.0, 1.0e7 / wl_vac.min() + 5.0
    table = 'O2_%s_2020' % band

    kept = []
    with open(fhitran) as f:
        for ln in f:
            if len(ln) < 15 or ln[:2].strip() != str(MOL_O2):
                continue
            try:
                nu = float(ln[3:15])
            except ValueError:
                continue
            if nu_lo <= nu <= nu_hi:
                kept.append(ln.rstrip('\n').ljust(160) + '\n')

    with open(os.path.join(cache_dir, table + '.data'), 'w') as f:
        f.writelines(kept)
    hdr = dict(_PAR_HEADER, table_name=table, number_of_rows=len(kept))
    with open(os.path.join(cache_dir, table + '.header'), 'w') as f:
        json.dump(hdr, f, indent=2)
    return table, len(kept)


def run(bands, z_top, cache_dir, layers=None):
    import hapi

    data_dir = os.environ.get('O2BAND_DATA_DIR',
                              os.path.normpath(os.path.join(_HERE, '..', 'data')))
    atm = afgl_atmosphere(os.path.join(data_dir, 'afglms.dat'), z_top=z_top)
    tips = tips2021()
    fhit = os.path.join(data_dir, 'hitran2020_lines.txt')

    os.makedirs(cache_dir, exist_ok=True)

    if layers is None:
        p = atm.lay['p']
        layers = [int(np.argmax(p)), int(np.argmin(np.abs(p - 100.0)))]

    for band in bands:
        table, nrow = register_local_table(fhit, band, cache_dir)
        hapi.db_begin(cache_dir)                             # (re)load tables
        lines = hitran_lines(fhit, wl_range=BANDS[band], margin_cm=5.0)
        absb = o2band_absorption(atm, lines, band=band, include_h2o=False, tips=tips)

        print('\n=========== %s : O2 sigma, ours vs HAPI on IDENTICAL 2020 lines '
              '(%d lines) ===========' % (band, nrow))
        for iz in layers:
            nu, sig_ours = our_sigma_layer(absb, atm, iz, 'o2')
            _, sig_hapi = hapi_sigma_layer(hapi, table, nu, atm.lay['p'][iz],
                                           atm.lay['temperature'][iz],
                                           float(atm.lay['o2_vmr'][iz]))
            peak_ratio = sig_ours.max() / sig_hapi.max() - 1.0
            st = diff_stats(sig_ours, sig_hapi, ref_floor=0.01 * sig_hapi.max())
            print('\n layer %d: p=%.1f hPa T=%.1f K  peak ours=%.4e hapi=%.4e '
                  '(ours %+.3f%%)' % (iz, atm.lay['p'][iz], atm.lay['temperature'][iz],
                                      sig_ours.max(), sig_hapi.max(), 100 * peak_ratio))
            print_diff_stats('O2 sigma (rel stats at cores >1%% of peak)', st)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--bands', nargs='+', default=['o2a', 'o2b'])
    p.add_argument('--z-top', type=float, default=120.0)
    p.add_argument('--cache-dir', default=None)
    args = p.parse_args()

    try:
        import hapi  # noqa: F401
    except ImportError:
        sys.exit('Error [eval_hapi_local]: HAPI not installed (pip install hitran-api).')

    out = os.environ.get('O2BAND_OUT_DIR', os.path.join(_HERE, '..', 'out'))
    cache = args.cache_dir or os.path.join(out, 'hapi_cache')
    run(args.bands, args.z_top, cache)
