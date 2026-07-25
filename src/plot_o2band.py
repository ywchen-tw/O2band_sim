"""
Plot O2-band benchmark results from an assembled per-band HDF5
(o2a.h5 / o2b.h5, schema written by sim_o2band.O2BandSim.assemble).

Produces a stacked figure:
  1. TOA reflectance ρ(λ) — one curve per (SZA, albedo)
  2. Absolute radiance I(λ) — one curve per (SZA, albedo)
  3. Optical thickness: O2 (+ H2O if present) column OD on a log axis, and
     Rayleigh column OD on a twin linear axis

Usage:
    python plot_o2band.py <band.h5> [out.png]
"""

import os
import sys
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from util.refcheck import rayleigh_ss_reflectance, lambertian_over_atm


def plot_band(fname, fout=None):

    with h5py.File(fname, 'r') as f:
        band = f.attrs.get('band', os.path.basename(fname))
        if isinstance(band, bytes):
            band = band.decode()
        wvl = f['wvl'][:]
        sza = f['sza'][:]
        alb = f['albedo'][:]
        ref = f['reflectance'][:]           # (Nsza, Nalb, Nwvl)
        ref_e = f['reflectance_stderr'][:]
        rad = f['radiance'][:]
        o2 = f['optical_thickness/o2_column'][:]
        ray = f['optical_thickness/rayleigh_column'][:]
        h2o = f['optical_thickness/h2o_column'][:] if 'optical_thickness/h2o_column' in f else None

    fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(11, 10), sharex=True)

    # colour by albedo, linestyle by SZA
    styles = ['-', '--', ':', '-.']
    cmap = plt.get_cmap('viridis')
    ncol = max(len(alb), 2)

    # analytic reference at the cleanest window point (O2 negligible):
    # single-scattering Rayleigh (+ Lambertian surface).  Overlaid as markers so
    # the MC continuum can be eye-checked against theory.
    iw = int(np.argmin(o2))
    tauR = float(np.median(ray))
    have_window = o2[iw] < 0.02

    for ia, a in enumerate(alb):
        color = cmap(ia / (ncol - 1) if ncol > 1 else 0.0)
        for js, s in enumerate(sza):
            ls = styles[js % len(styles)]
            lab = 'α=%.2f, SZA=%.0f°' % (a, s)
            ax0.plot(wvl, ref[js, ia], ls, color=color, lw=0.8, label=lab)
            ax1.plot(wvl, rad[js, ia], ls, color=color, lw=0.8, label=lab)
            if have_window:
                mu0 = np.cos(np.deg2rad(float(s)))
                rho_ref = (rayleigh_ss_reflectance(tauR, mu0) if a <= 0
                           else lambertian_over_atm(tauR, float(a), mu0))
                ax0.plot(wvl[iw], rho_ref, marker='*', ms=11, mfc='none',
                         mec=color, mew=1.4, zorder=5,
                         label=('analytic ref' if (ia == 0 and js == 0) else None))

    ax0.set_ylabel('TOA reflectance')
    ax0.set_title('O2 %s-band benchmark  (%.3f–%.3f nm, %d pts @ %.4f nm)'
                  % (str(band).upper().replace('O2', ''), wvl[0], wvl[-1], wvl.size,
                     np.median(np.diff(wvl))))
    ax0.legend(fontsize=8, ncol=2, loc='upper right')
    ax0.margins(x=0)
    ax0.grid(alpha=0.25)

    ax1.set_ylabel('Radiance\n(W m$^{-2}$ nm$^{-1}$ sr$^{-1}$)')
    ax1.margins(x=0)
    ax1.grid(alpha=0.25)

    # optical thickness
    ax2.semilogy(wvl, o2, color='C3', lw=0.8, label='O$_2$ absorption')
    if h2o is not None and np.nanmax(h2o) > 0:
        ax2.semilogy(wvl, h2o, color='C0', lw=0.8, label='H$_2$O absorption')
    ax2.set_ylabel('Column optical\nthickness (O$_2$/H$_2$O)')
    ax2.set_xlabel('Wavelength (nm, air)')
    ax2.grid(alpha=0.25, which='both')

    axr = ax2.twinx()
    axr.plot(wvl, ray, color='0.5', lw=0.8, label='Rayleigh')
    axr.set_ylabel('Rayleigh OD', color='0.4')
    axr.tick_params(axis='y', colors='0.4')

    # merged legend for OD panel
    h_a, l_a = ax2.get_legend_handles_labels()
    h_b, l_b = axr.get_legend_handles_labels()
    ax2.legend(h_a + h_b, l_a + l_b, fontsize=8, loc='upper right')
    ax2.margins(x=0)

    fig.tight_layout()

    if fout is None:
        fout = os.path.splitext(fname)[0] + '_spectrum.png'
    fig.savefig(fout, dpi=140)
    plt.close(fig)
    print('wrote %s' % fout)
    return fout


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('usage: python plot_o2band.py <band.h5> [out.png]')
        sys.exit(1)
    plot_band(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
