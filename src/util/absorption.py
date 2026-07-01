"""
Line-by-line absorption for the O2 A- and B-bands from HITRAN 2020 line
parameters, using the Voigt line shape.

Reads a HITRAN 160-column ``.par`` file (``data/hitran2020_lines.txt``,
molecules H2O = 1 and O2 = 7), an AFGL layered atmosphere
(:class:`util.atmosphere.afgl_atmosphere`), and produces the per-layer
absorption optical depth on a fine **air**-wavelength grid (default 0.001 nm),
as needed for the O2-band radiative-transfer intercomparison.

Locked Phase-1 conventions (see PLAN.md §7)
-------------------------------------------
- Partition function Q(296)/Q(T): **HITRAN TIPS-2021** (per isotopologue),
  consistent with the HITRAN2020 line intensities (Gamache et al. 2021).
- Wavelength grid: **air** wavelengths (Edlen 1966); HITRAN nu is vacuum, so
  each air grid point is evaluated at its corresponding vacuum wavenumber.
- Line-wing cutoff: +/- N * (nu0 / R) per line, defaults N = 3, R = 20000.
- O2 collision-induced absorption / continuum: **excluded** for Phase 1.
- H2O: `include_h2o` toggles H2O *absorption* lines only. Broadening is fixed --
  water vapor is always part of the foreign (air) pressure and broadens O2 at the
  air rate `gamma_air`; there is no separate `gamma_H2O`.

Physics
-------
  S(T) = S(296) * Q(296)/Q(T)
                * exp(-c2 E''/T) / exp(-c2 E''/296)
                * (1 - exp(-c2 nu0/T)) / (1 - exp(-c2 nu0/296))
  gamma_L = (296/T)^n_air * [gamma_air (P - P_self) + gamma_self P_self]   (atm)
  alpha_D = nu0 * sqrt(2 ln2 kB T / m) / c
  f(nu)   = sqrt(ln2/pi)/alpha_D * Re[w(x + i y)],
            x = sqrt(ln2)(nu - nu_shift)/alpha_D,  y = sqrt(ln2) gamma_L/alpha_D
  sigma   = S(T) f(nu)                  [cm^2/molecule]
  tau_lay = sigma * N_col               [-]   (N_col in molec cm^-2)

HITRAN S(296) already includes terrestrial isotopic abundance, so it multiplies
the *total* species column (AFGL o2/h2o); no abundance double-counting.
"""

import os
import numpy as np
from scipy.special import wofz

from .atmosphere import afgl_atmosphere
from .tips import tips2021
from .optics import air_to_vac_nm


__all__ = ['hitran_lines', 'o2band_absorption', 'cal_rayleigh_od',
           'voigt_profile', 'doppler_hwhm', 'BANDS']


# ---------------------------------------------------------------------------- #
# physical constants (CODATA / HITRAN)
C2        = 1.4387768775039337  # second radiation constant hc/kB (cm K)
KB        = 1.380649e-23        # Boltzmann constant (J/K)
C_MS      = 2.99792458e8        # speed of light (m/s)
AMU       = 1.66053906660e-27   # atomic mass unit (kg)
T_REF     = 296.0               # HITRAN reference temperature (K)
P_REF_HPA = 1013.25             # 1 atm in hPa
LN2         = np.log(2.0)
SQRT_LN2    = np.sqrt(LN2)
SQRT_LN2_PI = np.sqrt(LN2 / np.pi)

MOL_H2O = 1
MOL_O2  = 7

# molecular mass (g/mol) per (molecule id, isotopologue) for the Doppler width
ISO_MASS = {
    (MOL_H2O, 1): 18.010565, (MOL_H2O, 2): 20.014811, (MOL_H2O, 3): 19.014780,
    (MOL_H2O, 4): 19.016740, (MOL_H2O, 5): 21.020985, (MOL_H2O, 6): 20.020956,
    (MOL_H2O, 7): 20.022915,
    (MOL_O2, 1): 31.989830, (MOL_O2, 2): 33.994076, (MOL_O2, 3): 32.994045,
}

# intercomparison band windows (nm, air wavelengths per PLAN.md §7)
BANDS = {
    'o2b': (680.0, 695.0),   # O2 B-band
    'o2a': (757.0, 772.0),   # O2 A-band
}


# ---------------------------------------------------------------------------- #
def voigt_profile(nu, nu0, alpha_D, gamma_L):
    """
    Area-normalised Voigt profile (units cm, i.e. 1/cm-1) at wavenumbers ``nu``.

        f(nu) = sqrt(ln2/pi)/alpha_D * Re[w(x + i y)]
        x = sqrt(ln2)(nu - nu0)/alpha_D,   y = sqrt(ln2) gamma_L/alpha_D

    alpha_D : Doppler HWHM (cm-1); gamma_L : Lorentz HWHM (cm-1).
    """
    x = SQRT_LN2 * (np.asarray(nu, dtype=np.float64) - nu0) / alpha_D
    y = SQRT_LN2 * gamma_L / alpha_D
    return SQRT_LN2_PI / alpha_D * wofz(x + 1j * y).real


def doppler_hwhm(nu0, temperature, mass_gmol):
    """Doppler half-width at half-maximum (cm-1) for a line at nu0 (cm-1)."""
    return nu0 * np.sqrt(2.0 * LN2 * KB * temperature / (mass_gmol * AMU)) / C_MS


# ---------------------------------------------------------------------------- #
class hitran_lines:

    """
    Parse a HITRAN 160-column ``.par`` file.

    Input:
        fname   : path to the HITRAN file (data/hitran2020_lines.txt)
        wl_range: optional (wl_min, wl_max) *air* nm; lines within a small
                  wavenumber margin of the corresponding vacuum window are kept.
        margin_cm : wavenumber margin (cm-1) kept on each side, default 5.0.

    Attributes (1-D arrays, one entry per line):
        mol, iso, nu (cm-1, vacuum), S (cm-1/(molec cm-2) @296K),
        gamma_air, gamma_self (cm-1/atm @296K), Epp (cm-1), n_air,
        delta_air (cm-1/atm), mass (g/mol)
    """

    _FIELDS = {
        'mol':   (0, 2),   'iso':    (2, 3),   'nu':      (3, 15),
        'S':     (15, 25), 'A':      (25, 35), 'gair':    (35, 40),
        'gself': (40, 45), 'Epp':    (45, 55), 'nair':    (55, 59),
        'dair':  (59, 67),
    }

    def __init__(self, fname, wl_range=None, margin_cm=5.0):
        self.fname = fname
        self._parse(fname)
        if wl_range is not None:
            self._select(wl_range, margin_cm)

    def _parse(self, fname):
        mol, iso, nu, S, gair, gself, Epp, nair, dair = ([] for _ in range(9))
        f = self._FIELDS
        with open(fname) as fo:
            for ln in fo:
                if len(ln) < 67:
                    continue
                sub = ln[f['nu'][0]:f['nu'][1]].strip()
                if not sub.replace('.', '').replace('-', '').isdigit():
                    continue
                try:
                    m = int(ln[f['mol'][0]:f['mol'][1]])
                except ValueError:
                    continue
                mol.append(m)
                iso.append(int(ln[f['iso'][0]:f['iso'][1]]))
                nu.append(float(ln[f['nu'][0]:f['nu'][1]]))
                S.append(float(ln[f['S'][0]:f['S'][1]]))
                gair.append(float(ln[f['gair'][0]:f['gair'][1]]))
                gself.append(float(ln[f['gself'][0]:f['gself'][1]]))
                Epp.append(float(ln[f['Epp'][0]:f['Epp'][1]]))
                nair.append(float(ln[f['nair'][0]:f['nair'][1]]))
                dair.append(float(ln[f['dair'][0]:f['dair'][1]]))

        self.mol        = np.array(mol, dtype=np.int32)
        self.iso        = np.array(iso, dtype=np.int32)
        self.nu         = np.array(nu, dtype=np.float64)
        self.S          = np.array(S, dtype=np.float64)
        self.gamma_air  = np.array(gair, dtype=np.float64)
        self.gamma_self = np.array(gself, dtype=np.float64)
        self.Epp        = np.array(Epp, dtype=np.float64)
        self.n_air      = np.array(nair, dtype=np.float64)
        self.delta_air  = np.array(dair, dtype=np.float64)
        self.mass = np.array(
            [ISO_MASS.get((m, i), np.nan) for m, i in zip(self.mol, self.iso)],
            dtype=np.float64)

    def _select(self, wl_range, margin_cm):
        # air window -> vacuum wavelength -> vacuum wavenumber range
        wl_vac = air_to_vac_nm(np.array([min(wl_range), max(wl_range)]))
        nu_lo = 1.0e7 / wl_vac.max() - margin_cm
        nu_hi = 1.0e7 / wl_vac.min() + margin_cm
        keep = (self.nu >= nu_lo) & (self.nu <= nu_hi)
        for key in ('mol', 'iso', 'nu', 'S', 'gamma_air', 'gamma_self',
                    'Epp', 'n_air', 'delta_air', 'mass'):
            setattr(self, key, getattr(self, key)[keep])

    def subset(self, mol):
        return self.mol == mol

    def __repr__(self):
        return ('hitran_lines(%s): %d lines (O2=%d, H2O=%d), %.2f-%.2f cm-1'
                % (os.path.basename(self.fname), self.nu.size,
                   int(np.sum(self.mol == MOL_O2)), int(np.sum(self.mol == MOL_H2O)),
                   self.nu.min() if self.nu.size else 0,
                   self.nu.max() if self.nu.size else 0))


# ---------------------------------------------------------------------------- #
def _voigt_od_for_gas(nu_grid, idx, lines, tips, p_atm, p_self_atm,
                      temperature, column, R, ncut):
    """
    Optical depth on ``nu_grid`` (ascending vacuum cm-1) for one gas in one layer.

    idx : indices into `lines` of the gas's lines.
    Each line is summed only over grid points within +/- ncut*nu0/R of its
    (pressure-shifted) centre.
    """
    od = np.zeros_like(nu_grid)
    if idx.size == 0:
        return od

    nu0   = lines.nu[idx]
    S296  = lines.S[idx]
    gair  = lines.gamma_air[idx]
    gself = lines.gamma_self[idx]
    Epp   = lines.Epp[idx]
    nair  = lines.n_air[idx]
    dair  = lines.delta_air[idx]
    mass  = lines.mass[idx]
    iso   = lines.iso[idx]
    mol   = int(lines.mol[idx[0]])

    # per-isotopologue partition ratio Q(296)/Q(T) (TIPS-2021)
    qr = np.empty(idx.size)
    for iso0 in np.unique(iso):
        qr[iso == iso0] = tips.ratio(mol, int(iso0), temperature)

    boltz = np.exp(-C2 * Epp / temperature) / np.exp(-C2 * Epp / T_REF)
    stim  = ((1.0 - np.exp(-C2 * nu0 / temperature)) /
             (1.0 - np.exp(-C2 * nu0 / T_REF)))
    S_T = S296 * qr * boltz * stim

    nu_c    = nu0 + dair * p_atm                                   # shifted centre
    gamma_L = (T_REF / temperature) ** nair * (gair * (p_atm - p_self_atm)
                                               + gself * p_self_atm)
    alpha_D = doppler_hwhm(nu_c, temperature, mass)
    amp     = S_T * column                                        # cm
    cutoff  = ncut * nu0 / R                                      # cm-1, per line

    for k in range(nu_c.size):
        lo = np.searchsorted(nu_grid, nu_c[k] - cutoff[k], side='left')
        hi = np.searchsorted(nu_grid, nu_c[k] + cutoff[k], side='right')
        if hi <= lo:
            continue
        od[lo:hi] += amp[k] * voigt_profile(nu_grid[lo:hi], nu_c[k],
                                            alpha_D[k], gamma_L[k])

    return od


# ---------------------------------------------------------------------------- #
class o2band_absorption:

    """
    Line-by-line O2/H2O absorption optical depth for one band.

    Input:
        atm         : afgl_atmosphere object
        lines       : hitran_lines object (may span both bands)
        band        : key in BANDS ('o2a'/'o2b') or an (air wl_min, wl_max) tuple (nm)
        dwvl        : air-wavelength grid spacing (nm), default 0.001
        R           : resolving power defining the wing cutoff, default 20000
        ncut        : cutoff = ncut * nu0/R (line half-widths of the resolution
                      element), default 3
        include_h2o : include H2O *absorption* lines, default True. This ONLY
                      controls whether H2O lines contribute to the optical depth.
                      Broadening is unchanged either way: O2 lines are always
                      broadened with foreign pressure P - P_O2, which includes
                      water vapor at the air-broadening rate (the current setup).
        gases       : explicit molecule list; if None, derived from include_h2o
                      (('o2','h2o') or ('o2',))
        tips        : tips2021 instance (created if None)

    Output attributes:
        self.wvl       : air wavelength grid (nm, ascending), (Nwvl,)
        self.nu_vac    : vacuum wavenumber of each grid point (cm-1), (Nwvl,)
        self.od        : dict of per-gas layer OT (Nlay, Nwvl): od['o2'], od['h2o']
        self.od_total  : total gas OT (Nlay, Nwvl)
        self.od_column : total-column gas OT (Nwvl,)
    """

    _GAS_MOL = {'o2': MOL_O2, 'h2o': MOL_H2O}

    def __init__(self, atm, lines, band, dwvl=0.001, R=20000.0, ncut=3.0,
                 include_h2o=True, gases=None, tips=None):

        if isinstance(band, str):
            wl_min, wl_max = BANDS[band]
            self.band = band
        else:
            wl_min, wl_max = band
            self.band = 'custom'

        if gases is None:
            gases = ('o2', 'h2o') if include_h2o else ('o2',)

        self.atm = atm
        self.dwvl = dwvl
        self.R = R
        self.ncut = ncut
        self.include_h2o = include_h2o
        self.gases = tuple(gases)
        self.tips = tips if tips is not None else tips2021()

        # air-wavelength grid, and the vacuum wavenumber sampled at each point
        self.wvl = np.arange(wl_min, wl_max + 0.5 * dwvl, dwvl)
        wl_vac = air_to_vac_nm(self.wvl)
        self.nu_vac = 1.0e7 / wl_vac
        # ascending vacuum-wavenumber grid for windowed line summation
        self._order = np.argsort(self.nu_vac)          # wvl-order -> nu-asc order
        self._nu_sorted = self.nu_vac[self._order]
        self._inv = np.argsort(self._order)            # nu-asc order -> wvl order

        self._run(lines)

    def _run(self, lines):
        lay = self.atm.lay
        nlay = lay['z'].size
        nwvl = self.wvl.size

        self.od = {g: np.zeros((nlay, nwvl)) for g in self.gases}
        gas_idx = {g: np.where(lines.subset(self._GAS_MOL[g]))[0] for g in self.gases}

        for iz in range(nlay):
            p_atm = lay['p'][iz] / P_REF_HPA
            temperature = lay['temperature'][iz]
            for gas in self.gases:
                p_self_atm = p_atm * lay['%s_vmr' % gas][iz]
                od_nu = _voigt_od_for_gas(
                    self._nu_sorted, gas_idx[gas], lines, self.tips,
                    p_atm, p_self_atm, temperature, lay[gas][iz],
                    self.R, self.ncut)
                self.od[gas][iz, :] = od_nu[self._inv]  # back to wvl order

        self.od_total = sum(self.od.values())
        self.od_column = self.od_total.sum(axis=0)

    def subrange_indices(self, wvl_range=None):
        """
        Indices of the grid falling within `wvl_range` (air nm); None → all.

        Enforces that the requested range lies within the band's grid extent:
        an out-of-band request raises ValueError rather than silently running a
        different range. Returns a contiguous 1-D index array.
        """
        if wvl_range is None:
            return np.arange(self.wvl.size)

        a, b = float(min(wvl_range)), float(max(wvl_range))
        wl0, wl1 = float(self.wvl[0]), float(self.wvl[-1])
        tol = 0.5 * self.dwvl  # half a grid step of slack at the edges

        if (a < wl0 - tol) or (b > wl1 + tol):
            raise ValueError(
                'wvl_range (%.4f, %.4f) nm is outside band %s grid [%.4f, %.4f] nm'
                % (a, b, self.band, wl0, wl1))

        idx = np.where((self.wvl >= a - tol) & (self.wvl <= b + tol))[0]
        if idx.size == 0:
            raise ValueError(
                'wvl_range (%.4f, %.4f) nm selects no grid points on the %.4f nm grid'
                % (a, b, self.dwvl))
        return idx

    def __repr__(self):
        return ('o2band_absorption(band=%s): %d air-wvl (%.3f-%.3f nm @ %.4f nm), '
                '%d layers, gases=%s, cutoff=%g*nu0/%g'
                % (self.band, self.wvl.size, self.wvl[0], self.wvl[-1], self.dwvl,
                   self.atm.lay['z'].size, ','.join(self.gases), self.ncut, self.R))


# ---------------------------------------------------------------------------- #
def cal_rayleigh_od(atm, wvl_nm):
    """
    Rayleigh-scattering optical depth per layer, Bodhaine et al. (1999).

    Input:
        atm    : afgl_atmosphere object (uses layer air column, molec cm-2)
        wvl_nm : wavelength(s) in nm (scalar or array)
    Output:
        od_ray : Rayleigh optical depth, shape (Nlay, Nwvl)
    """
    wl = np.atleast_1d(np.asarray(wvl_nm, dtype=np.float64))
    wl_um = wl * 1.0e-3
    inv2 = 1.0 / wl_um ** 2

    # Bodhaine (1999) Rayleigh cross-section (cm^2/molecule), 360 ppm CO2 air
    num = 1.0455996 - 341.29061 * inv2 - 0.90230850 * wl_um ** 2
    den = 1.0 + 0.0027059889 * inv2 - 85.968563 * wl_um ** 2
    xsec = 1.0e-28 * num / den

    air_col = atm.lay['air']  # molec cm-2
    return air_col[:, None] * xsec[None, :]


# ---------------------------------------------------------------------------- #
if __name__ == '__main__':

    import time
    here = os.path.dirname(os.path.abspath(__file__))
    fdir_data = os.path.join(here, '..', '..', 'data')

    atm = afgl_atmosphere(os.path.join(fdir_data, 'afglms.dat'), z_top=70.0)
    print(atm)
    tips = tips2021()

    for band in ('o2b', 'o2a'):
        lines = hitran_lines(os.path.join(fdir_data, 'hitran2020_lines.txt'),
                             wl_range=BANDS[band])
        print(lines)
        t0 = time.time()
        absb = o2band_absorption(atm, lines, band, dwvl=0.001, tips=tips)
        print(absb, '(%.1f s)' % (time.time() - t0))
        od_ray = cal_rayleigh_od(atm, absb.wvl)
        print('  %s: max column O2+H2O OT = %.3f | column O2 OT max = %.3f | '
              'mean Rayleigh OT = %.4f'
              % (band, absb.od_column.max(), absb.od['o2'].sum(axis=0).max(),
                 od_ray.sum(axis=0).mean()))
