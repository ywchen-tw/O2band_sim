#!/usr/bin/env python
"""
Local pre-flight for the 2026-07-25 corrections, before committing CURC time.

Answers three questions on a few wavelengths, with no RT solver involved:

  1. grid convention  -- what the vacuum grid does to the wavelength registration
  2. wing cutoff      -- how column O2 optical depth and the implied continuum
                         reflectance respond to cutoff = 25 / 33 / 50 cm-1,
                         referenced to a brute-force NO-cutoff Voigt sum
  3. O2-O2 CIA        -- how much a HITRAN CIA file adds on top of (2)

Only the optical depths here are exact.  Reflectance is a single-scattering
estimate (util.refcheck) used to translate OD differences into the quantity the
intercomparison actually plots; absolute reflectance still comes from MCARaTS.
The estimate is anchored by printing the pure-Rayleigh (zero-gas) value, which
should sit near the 0.1075 that MCARaTS and DISORT both give at alb 0.1, SZA 0.

    python src/eval_cutoff_cia.py                       # A-band, defaults
    python src/eval_cutoff_cia.py --band o2b --cutoffs 25 33 50
"""

import os
import sys
import argparse
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from util.inputs import require_inputs, MissingInputs
from util.atmosphere import afgl_atmosphere
from util.tips import tips2021
from util.cia import hitran_cia, pick_cia_file
from util.optics import vac_to_air_nm
from util.absorption import (hitran_lines, o2band_absorption, cal_rayleigh_od,
                             voigt_profile, doppler_hwhm, required_margin_cm,
                             C2, T_REF, P_REF_HPA, MOL_O2)
from util.refcheck import rayleigh_ss_reflectance

DATA = os.environ.get('O2BAND_DATA_DIR',
                      os.path.normpath(os.path.join(_HERE, '..', 'data')))

# default sub-windows (vacuum nm): dense-line regions where the collaborator's
# plot shows the disagreement
WINDOWS = {'o2a': (762.0, 766.0), 'o2b': (686.0, 690.0)}


def rho_estimate(tau_ray, tau_gas, albedo, sza_deg, mu=1.0):
    """
    TOA reflectance estimate: single-scattering Rayleigh path radiance plus a
    Lambertian surface seen through both the Rayleigh and the gas optical depth.

        rho = rho_ss(tau_ray) * exp(-tau_gas * <path>) + T_dn T_up A / (1 - A s)

    with T = exp(-(tau_ray + tau_gas)/mu).  Approximate (no multiple
    scattering); used for *differences between configurations*, not as a
    benchmark value.
    """
    mu0 = np.cos(np.deg2rad(sza_deg))
    tot = tau_ray + tau_gas
    rho_atm = (rayleigh_ss_reflectance(tau_ray, mu0, mu)
               * np.exp(-tau_gas * 0.5 * (1.0 / mu0 + 1.0 / mu)))
    s = 0.5 * tau_ray
    rho_surf = albedo * np.exp(-tot / mu0) * np.exp(-tot / mu) / (1.0 - albedo * s)
    return rho_atm + rho_surf


def brute_force_column_od(atm, tips, nu_grid, pad_nm, margin_cm=300.0):
    """Column O2 OD summing every line within +/-margin_cm, no wing cutoff."""
    lines = hitran_lines(os.path.join(DATA, 'hitran2020_lines.txt'),
                         wl_range=pad_nm, grid='vac', margin_cm=margin_cm)
    idx = np.where(lines.subset(MOL_O2))[0]
    nu0, S296 = lines.nu[idx], lines.S[idx]
    gair, gself = lines.gamma_air[idx], lines.gamma_self[idx]
    Epp, nair, dair = lines.Epp[idx], lines.n_air[idx], lines.delta_air[idx]
    mass, iso = lines.mass[idx], lines.iso[idx]

    od = np.zeros(nu_grid.size)
    lay = atm.lay
    for iz in range(lay['z'].size):
        T = lay['temperature'][iz]
        p_atm = lay['p'][iz] / P_REF_HPA
        p_self = p_atm * lay['o2_vmr'][iz]
        qr = np.array([tips.ratio(MOL_O2, int(i), T) for i in iso])
        S_T = (S296 * qr
               * np.exp(-C2 * Epp / T) / np.exp(-C2 * Epp / T_REF)
               * (1 - np.exp(-C2 * nu0 / T)) / (1 - np.exp(-C2 * nu0 / T_REF)))
        nu_c = nu0 + dair * p_atm
        gL = (T_REF / T) ** nair * (gair * (p_atm - p_self) + gself * p_self)
        aD = doppler_hwhm(nu_c, T, mass)
        for k in range(nu_c.size):
            od += lay['o2'][iz] * S_T[k] * voigt_profile(nu_grid, nu_c[k], aD[k], gL[k])
    return od


def build(atm, tips, win, dwvl, cutoff_cm, grid='vac', cia=None):
    margin = required_margin_cm(win, grid=grid, cutoff_cm=cutoff_cm)
    lines = hitran_lines(os.path.join(DATA, 'hitran2020_lines.txt'),
                         wl_range=win, grid=grid, margin_cm=margin)
    return o2band_absorption(atm, lines, win, dwvl=dwvl, cutoff_cm=cutoff_cm,
                             grid=grid, cia=cia, include_h2o=False,
                             gases=('o2',), tips=tips)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--band', default='o2a', choices=('o2a', 'o2b'))
    p.add_argument('--window', nargs=2, type=float, default=None,
                   metavar=('A', 'B'), help='vacuum-nm sub-window')
    p.add_argument('--dwvl', type=float, default=0.01,
                   help='grid spacing nm (0.01 keeps the brute force quick)')
    p.add_argument('--cutoffs', nargs='+', type=float, default=[25.0, 33.0, 50.0])
    p.add_argument('--cia', default='auto', metavar='FILE|auto|none')
    p.add_argument('--z-top', type=float, default=120.0)
    p.add_argument('--albedo', type=float, default=0.1)
    p.add_argument('--sza', type=float, default=0.0)
    p.add_argument('--npts', type=int, default=40,
                   help='number of most-transparent points defining "micro-window"')
    args = p.parse_args()

    win = tuple(args.window) if args.window else WINDOWS[args.band]
    try:
        require_inputs('hitran2020_lines.txt', 'afglms.dat', qtpy=True)
    except MissingInputs as e:
        sys.exit('[SETUP] %s' % e)
    atm = afgl_atmosphere(os.path.join(DATA, 'afglms.dat'), z_top=args.z_top)
    tips = tips2021()

    fcia = args.cia
    if fcia in ('none', 'off'):
        fcia = None
    elif fcia == 'auto':
        nu_win = (1.0e7 / max(win), 1.0e7 / min(win))
        try:
            fcia = pick_cia_file(DATA, [nu_win])
        except OSError as e:
            print('  [cia] auto-selection failed: %s' % e)
            fcia = None
    cia = hitran_cia(fcia) if fcia else None

    print('=' * 78)
    print('Pre-flight: grid convention, wing cutoff, O2-O2 CIA')
    print('=' * 78)
    print('band %s, window %.1f-%.1f nm VACUUM @ %.3f nm, z_top %.0f km'
          % (args.band, win[0], win[1], args.dwvl, args.z_top))
    print('reflectance estimate: albedo %.2f, SZA %.0f deg, nadir view'
          % (args.albedo, args.sza))
    print('CIA file: %s' % (os.path.basename(fcia) if fcia else 'none'))

    # ---- 1. grid convention -------------------------------------------------
    ref_cut = max(args.cutoffs)
    absb = build(atm, tips, win, args.dwvl, ref_cut)
    print('\n--- 1. grid convention -------------------------------------------')
    print('  this window in AIR wavelengths: %.4f - %.4f nm'
          % (absb.wvl_air[0], absb.wvl_air[-1]))
    print('  vac - air = %.4f nm  => the old air-grid delivery plotted every'
          % (absb.wvl_vac - absb.wvl_air).mean())
    print('  feature %.0f grid cells (of 0.001 nm) SHORT of a vacuum-grid model.'
          % ((absb.wvl_vac - absb.wvl_air).mean() / 0.001))
    if cia is not None:
        nu_lo, nu_hi = absb.nu_vac.min(), absb.nu_vac.max()
        print('  CIA coverage over %.1f-%.1f cm-1: %s, T tables %s'
              % (nu_lo, nu_hi, cia.covers(nu_lo, nu_hi),
                 cia.temperature_tables(nu_lo, nu_hi)))

    # ---- reference: no cutoff at all ---------------------------------------
    print('\nbuilding brute-force no-cutoff reference...')
    pad = (win[0] - 5.0, win[1] + 5.0)
    ref = brute_force_column_od(atm, tips, absb.nu_vac, pad)

    od_ray = cal_rayleigh_od(atm, absb.wvl_vac).sum(axis=0)
    kmw = np.argsort(ref)[:args.npts]              # micro-window points

    # ---- 2. wing cutoff ----------------------------------------------------
    print('\n--- 2. wing cutoff (column O2 optical depth) ----------------------')
    print('%-22s %10s %10s %10s %10s' % ('config', 'min OD', 'median', 'micro-win',
                                         'rho(est)'))
    rows = []

    def row(tag, od, od_extra=None):
        tot = od if od_extra is None else od + od_extra
        rho = rho_estimate(od_ray[kmw].mean(), tot[kmw].mean(),
                           args.albedo, args.sza)
        print('%-22s %10.5f %10.4f %10.5f %10.5f'
              % (tag, tot.min(), np.median(tot), tot[kmw].mean(), rho))
        rows.append((tag, tot[kmw].mean(), rho))
        return rho

    absb_legacy = build(atm, tips, win, args.dwvl, None)
    row('legacy ~2 cm-1', absb_legacy.od['o2'].sum(axis=0))
    cut_od = {}
    for c in args.cutoffs:
        a = build(atm, tips, win, args.dwvl, c)
        cut_od[c] = a.od['o2'].sum(axis=0)
        row('cutoff %g cm-1' % c, cut_od[c])
    row('no cutoff (reference)', ref)
    print('  Rayleigh column OD over the window: %.5f' % od_ray[kmw].mean())
    print('  pure-Rayleigh rho (zero gas)      : %.5f   <- anchor, cf. 0.1075 MC/DISORT'
          % rho_estimate(od_ray[kmw].mean(), 0.0, args.albedo, args.sza))

    print('\n  convergence of the micro-window OD toward the untruncated value:')
    for c in args.cutoffs:
        d = ref[kmw].mean() - cut_od[c][kmw].mean()
        print('    %5g cm-1: %.5f  (%.2f%% below no-cutoff, %+.0f%% of the '
              'legacy deficit recovered)'
              % (c, cut_od[c][kmw].mean(), 100 * d / ref[kmw].mean(),
                 100 * (cut_od[c][kmw].mean() - absb_legacy.od['o2'].sum(axis=0)[kmw].mean())
                 / (ref[kmw].mean() - absb_legacy.od['o2'].sum(axis=0)[kmw].mean())))

    # ---- 3. CIA on top -----------------------------------------------------
    print('\n--- 3. O2-O2 CIA on top of the 25 cm-1 cutoff --------------------')
    if cia is None:
        print('  no CIA file available -- skipped')
    else:
        c0 = args.cutoffs[0]
        a_cia = build(atm, tips, win, args.dwvl, c0, cia=cia)
        od_cia = a_cia.od_cia.sum(axis=0)
        print('%-22s %10s %10s %10s %10s' % ('config', 'min OD', 'median',
                                             'micro-win', 'rho(est)'))
        r_off = row('cutoff %g, CIA off' % c0, cut_od[c0])
        r_on = row('cutoff %g, CIA on' % c0, cut_od[c0], od_cia)
        print('  CIA column OT alone: min %.5f  mean %.5f  max %.5f'
              % (od_cia.min(), od_cia.mean(), od_cia.max()))
        print('  effect on micro-window rho: %.5f -> %.5f  (%+.2f%%)'
              % (r_off, r_on, 100 * (r_on - r_off) / r_off))

    # ---- a few individual wavelengths --------------------------------------
    print('\n--- representative wavelengths (vacuum / air nm) ------------------')
    order = np.argsort(ref)
    picks = [('most transparent', order[0]),
             ('25th pct', order[int(0.25 * order.size)]),
             ('median', order[int(0.50 * order.size)]),
             ('line core', order[-1])]
    hdr = ['%s' % ('cut %g' % c) for c in args.cutoffs]
    print('%-16s %9s %9s %9s %9s %s'
          % ('point', 'vac nm', 'air nm', 'legacy', 'no-cut', ' '.join('%9s' % h for h in hdr)))
    for name, i in picks:
        vals = ' '.join('%9.4f' % cut_od[c][i] for c in args.cutoffs)
        print('%-16s %9.3f %9.3f %9.4f %9.4f %s'
              % (name, absb.wvl_vac[i], absb.wvl_air[i],
                 absb_legacy.od['o2'].sum(axis=0)[i], ref[i], vals))


if __name__ == '__main__':
    main()
