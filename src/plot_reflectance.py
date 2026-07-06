#!/usr/bin/env python
"""
Two-panel TOA-reflectance figure for the O2 A/B-band benchmark: one panel per band
(O2A, O2B), reflectance vs air wavelength for all prescribed geometries.

Encoding: colour = solar zenith angle (Okabe-Ito colourblind-safe, fixed order),
linestyle = surface albedo (the two albedo groups also separate by magnitude --
albedo 0.1 ~0.10 surface-dominated, albedo 0.0 ~0.01 Rayleigh).

    python src/plot_reflectance.py [MERGED_OR_BAND.h5] [-o out.png]
"""

import os
import sys
import argparse
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Okabe-Ito (colourblind-safe), fixed order by SZA
SZA_COLOR = {0.0: '#0072B2', 30.0: '#E69F00', 60.0: '#009E73'}
ALB_STYLE = {0.0: (0, (4, 2)), 0.1: '-'}       # dashed = albedo 0, solid = albedo 0.1
BANDS = [('o2a', 'O$_2$ A-band'), ('o2b', 'O$_2$ B-band')]


def _group(f, band):
    return f[band] if (band in f and isinstance(f[band], h5py.Group)) else f


def plot(h5path, out_png):
    plt.rcParams.update({'font.size': 11, 'axes.titlesize': 12,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.edgecolor': '0.5', 'axes.linewidth': 0.8})
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.0))

    with h5py.File(h5path, 'r') as f:
        for ax, (band, title) in zip(axes, BANDS):
            g = _group(f, band)
            wvl = g['wvl'][:]; sza = g['sza'][:]; alb = g['albedo'][:]
            ref = g['reflectance'][:]
            for i, s in enumerate(sza):
                for j, a in enumerate(alb):
                    ax.plot(wvl, ref[i, j], color=SZA_COLOR.get(float(s), '0.3'),
                            ls=ALB_STYLE.get(float(a), '-'), lw=0.6, alpha=0.9)
            ax.set_title(title)
            ax.set_xlabel('Wavelength (nm, air)')
            ax.set_xlim(wvl[0], wvl[-1])
            ax.set_ylim(0, None)
            ax.margins(x=0)
            ax.grid(True, color='0.85', lw=0.4)
            ax.set_axisbelow(True)

    axes[0].set_ylabel('TOA reflectance')

    # single figure-level legend along the BOTTOM, clear of the data
    handles = [Line2D([0], [0], color=SZA_COLOR[s], lw=2.2, label='SZA %.0f°' % s)
               for s in (0.0, 30.0, 60.0)]
    handles += [Line2D([0], [0], color='0.35', ls='-', lw=2.2, label='albedo 0.1'),
                Line2D([0], [0], color='0.35', ls=(0, (4, 2)), lw=2.2, label='albedo 0.0')]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, 0.0),
               ncol=5, frameon=False, fontsize=10, columnspacing=1.8, handlelength=2.2)

    fig.suptitle('TOA reflectance — O$_2$ A/B-band Phase-1 clear-sky benchmark '
                 '(AFGL mid-lat summer, 0.001 nm)', fontsize=12)
    fig.tight_layout(rect=[0, 0.06, 1, 0.95])       # reserve bottom (legend) + top (suptitle)
    fig.savefig(out_png, dpi=200)
    print('wrote %s' % out_png)


if __name__ == '__main__':
    _HERE = os.path.dirname(os.path.abspath(__file__))
    default_h5 = os.path.join(
        os.environ.get('O2BAND_OUT_DIR', os.path.join(_HERE, '..', 'out')),
        'z120_p1e6_n3', 'o2band_benchmark.h5')
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('h5', nargs='?', default=default_h5)
    p.add_argument('-o', '--out', default=None)
    args = p.parse_args()
    if not os.path.isfile(args.h5):
        sys.exit('Error: no such file: %s' % args.h5)
    out = args.out or os.path.join(os.path.dirname(args.h5), 'reflectance_o2ab.png')
    plot(args.h5, out)
