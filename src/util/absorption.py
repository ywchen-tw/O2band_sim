"""
Line-by-line absorption for the O2 A- and B-bands from HITRAN 2020 line
parameters, using the Voigt line shape.

Reads a HITRAN 160-column ``.par`` file (``data/hitran2020_lines.txt``,
molecules H2O = 1 and O2 = 7), an AFGL layered atmosphere
(:class:`util.atmosphere.afgl_atmosphere`), and produces the per-layer
absorption optical depth on a fine wavelength grid (default 0.001 nm, vacuum),
as needed for the O2-band radiative-transfer intercomparison.

Locked Phase-1 conventions (see PLAN.md §7)
-------------------------------------------
- Partition function Q(296)/Q(T): **HITRAN TIPS-2021** (per isotopologue),
  consistent with the HITRAN2020 line intensities (Gamache et al. 2021).
- Wavelength grid: selectable via ``grid=``.  ``'vac'`` (default) puts the
  uniform output grid on **vacuum** wavelengths -- the convention used by the
  other intercomparison participants (PLAN.md §7.2, revised 2026-07-25);
  ``'air'`` reproduces the original air-wavelength delivery (Edlen 1966).
  Either way HITRAN nu is vacuum, so every grid point is evaluated at its own
  vacuum wavenumber, and BOTH ``wvl_air`` and ``wvl_vac`` are exposed.
- Line-wing cutoff: ``cutoff_cm`` cm-1 per line, plain truncation (default 50;
  see DEFAULT_CUTOFF_CM for why not the more common 25).  The legacy
  ``N * nu0 / R`` form is kept for reproducing earlier runs and is used only
  when ``cutoff_cm=None``.  Note the far Voigt wing is
  Lorentzian, whose area converges as 1/dnu -- a cutoff of a few cm-1 discards a
  physically significant amount of between-line absorption (PLAN.md §7.3).
- O2 collision-induced absorption / continuum: excluded by default; pass a
  ``cia=`` object (util.cia) to include it.  Kept separable in ``od_cia``.
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
from .optics import air_to_vac_nm, vac_to_air_nm


__all__ = ['hitran_lines', 'o2band_absorption', 'cal_rayleigh_od',
           'voigt_profile', 'doppler_hwhm', 'line_cutoff_cm', 'required_margin_cm',
           'band_nu_range', 'BANDS', 'DEFAULT_CUTOFF_CM']


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

# intercomparison band windows (nm).  Interpreted in whichever convention the
# grid uses (PLAN.md §7.2): vacuum by default, air for the legacy delivery.
BANDS = {
    'o2b': (680.0, 695.0),   # O2 B-band
    'o2a': (757.0, 772.0),   # O2 A-band
}

# Default line-wing cutoff (cm-1), with PLAIN truncation (no pedestal
# subtraction).  The Voigt far wing is Lorentzian and its area converges only as
# 1/dnu, so the cutoff is a pure numerical artifact -- and the benchmark
# prescribes "Voigt", not a cutoff.  50 cm-1 sits 1.7% below an untruncated sum
# in the A-band micro-windows (25 cm-1: 7.4%) and costs ~10 s for a full band at
# 0.001 nm, so there is no reason to accept the larger artifact.
#
# Two caveats, both documented in PLAN.md §7.3:
#   - 25 cm-1 is the more common LBL convention (LBLRTM), and LBLRTM pairs it
#     with subtracting f(cutoff) inside the cutoff -- a LARGER effect (-0.0039
#     column OD) than the 25->50 change (+0.0022), and of opposite sign.  A
#     participant reporting "25 cm-1" is ambiguous until that is pinned down.
#   - converging onto the untruncated Voigt is fidelity to the *prescribed
#     model*, not to nature: the real O2 far wing is sub-Lorentzian.
DEFAULT_CUTOFF_CM = 50.0


def line_cutoff_cm(nu0, R=20000.0, ncut=3.0, cutoff_cm=DEFAULT_CUTOFF_CM):
    """
    Per-line wing cutoff (cm-1) for line centres ``nu0``.

    ``cutoff_cm`` (absolute, the conventional form) wins when given; passing
    None falls back to the legacy resolution-element form ``ncut * nu0 / R``,
    which is retained only to reproduce pre-2026-07-25 runs.
    """
    nu0 = np.asarray(nu0, dtype=np.float64)
    if cutoff_cm is not None:
        return np.full(nu0.shape, float(cutoff_cm))
    return ncut * nu0 / R


def band_nu_range(wl_range, grid='vac'):
    """Vacuum wavenumber range (cm-1) spanned by a (wl_min, wl_max) nm window
    given in the ``grid`` convention ('vac' or 'air')."""
    wl = np.array([min(wl_range), max(wl_range)], dtype=np.float64)
    wl_vac = wl if grid == 'vac' else air_to_vac_nm(wl)
    return 1.0e7 / wl_vac.max(), 1.0e7 / wl_vac.min()


def required_margin_cm(wl_range, grid='vac', R=20000.0, ncut=3.0,
                       cutoff_cm=DEFAULT_CUTOFF_CM, slack=1.0):
    """
    Wavenumber margin (cm-1) that line selection must keep on each side of a
    window so that no line whose wing reaches into the window is dropped.

    Must be >= the largest per-line cutoff in the band: a margin narrower than
    the cutoff silently truncates the wings of out-of-window lines, which looks
    exactly like a too-small cutoff.  ``slack`` adds a small safety band.
    """
    _, nu_hi = band_nu_range(wl_range, grid=grid)
    return float(line_cutoff_cm(nu_hi, R=R, ncut=ncut, cutoff_cm=cutoff_cm)) + slack


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
        wl_range: optional (wl_min, wl_max) nm in the ``grid`` convention; lines
                  within ``margin_cm`` of the corresponding vacuum wavenumber
                  window are kept.
        grid    : 'vac' (default) or 'air' -- how wl_range is to be read.
        margin_cm : wavenumber margin (cm-1) kept on each side.  MUST be >= the
                  wing cutoff used downstream, else out-of-window lines lose
                  their wings; use ``required_margin_cm`` to derive it.

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

    def __init__(self, fname, wl_range=None, margin_cm=DEFAULT_CUTOFF_CM + 1.0,
                 grid='vac'):
        self.fname = fname
        self.grid = grid
        self.margin_cm = float(margin_cm)
        self._parse(fname)
        if wl_range is not None:
            self._select(wl_range, margin_cm, grid)

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

    def _select(self, wl_range, margin_cm, grid='vac'):
        # window (in its own convention) -> vacuum wavenumber range
        nu_lo, nu_hi = band_nu_range(wl_range, grid=grid)
        keep = (self.nu >= nu_lo - margin_cm) & (self.nu <= nu_hi + margin_cm)
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
                      temperature, column, R, ncut, cutoff_cm=DEFAULT_CUTOFF_CM):
    """
    Optical depth on ``nu_grid`` (ascending vacuum cm-1) for one gas in one layer.

    idx : indices into `lines` of the gas's lines.
    Each line is summed only over grid points within +/- its wing cutoff
    (``line_cutoff_cm``) of the (pressure-shifted) centre.
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
    cutoff  = line_cutoff_cm(nu0, R=R, ncut=ncut, cutoff_cm=cutoff_cm)

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
        band        : key in BANDS ('o2a'/'o2b') or a (wl_min, wl_max) tuple (nm)
        grid        : 'vac' (default) or 'air' -- convention of the uniform
                      output grid, and how ``band`` limits are interpreted
        dwvl        : grid spacing (nm) in that convention, default 0.001
        cutoff_cm   : per-line wing cutoff (cm-1), default 50.  None selects the
                      legacy ``ncut * nu0 / R`` form.
        R, ncut     : legacy cutoff parameters, used only when cutoff_cm is None
        cia         : optional util.cia.hitran_cia object; when given, O2-O2
                      collision-induced absorption is added (kept separable in
                      ``od_cia``).  Default None = excluded, as in Phase 1.
        include_h2o : include H2O *absorption* lines, default True. This ONLY
                      controls whether H2O lines contribute to the optical depth.
                      Broadening is unchanged either way: O2 lines are always
                      broadened with foreign pressure P - P_O2, which includes
                      water vapor at the air-broadening rate (the current setup).
        gases       : explicit molecule list; if None, derived from include_h2o
                      (('o2','h2o') or ('o2',))
        tips        : tips2021 instance (created if None)

    Output attributes:
        self.wvl       : output grid (nm, ascending) in the chosen convention
        self.wvl_vac   : vacuum wavelength of each grid point (nm), (Nwvl,)
        self.wvl_air   : air wavelength of each grid point (nm), (Nwvl,)
        self.nu_vac    : vacuum wavenumber of each grid point (cm-1), (Nwvl,)
        self.od        : dict of per-gas layer OT (Nlay, Nwvl): od['o2'], od['h2o']
        self.od_cia    : O2-O2 CIA layer OT (Nlay, Nwvl); zeros when cia is None
        self.od_total  : total gas OT incl. CIA (Nlay, Nwvl)
        self.od_column : total-column gas OT (Nwvl,)
    """

    _GAS_MOL = {'o2': MOL_O2, 'h2o': MOL_H2O}

    def __init__(self, atm, lines, band, dwvl=0.001, R=20000.0, ncut=3.0,
                 cutoff_cm=DEFAULT_CUTOFF_CM, grid='vac', cia=None,
                 include_h2o=True, gases=None, tips=None):

        if grid not in ('vac', 'air'):
            raise ValueError("grid must be 'vac' or 'air', got %r" % (grid,))

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
        self.cutoff_cm = cutoff_cm
        self.grid = grid
        self.cia = cia
        self.include_h2o = include_h2o
        self.gases = tuple(gases)
        self.tips = tips if tips is not None else tips2021()

        # uniform output grid in the requested convention; both conventions and
        # the vacuum wavenumber (what the HITRAN line data actually lives on)
        # are always carried, so downstream users never have to guess.
        self.wvl = np.arange(wl_min, wl_max + 0.5 * dwvl, dwvl)
        if grid == 'vac':
            self.wvl_vac = self.wvl.copy()
            self.wvl_air = vac_to_air_nm(self.wvl_vac)
        else:
            self.wvl_air = self.wvl.copy()
            self.wvl_vac = air_to_vac_nm(self.wvl_air)
        self.nu_vac = 1.0e7 / self.wvl_vac
        # ascending vacuum-wavenumber grid for windowed line summation
        self._order = np.argsort(self.nu_vac)          # wvl-order -> nu-asc order
        self._nu_sorted = self.nu_vac[self._order]
        self._inv = np.argsort(self._order)            # nu-asc order -> wvl order

        self._check_margin(lines)
        self._run(lines)

    def _check_margin(self, lines):
        """Guard the coupling between wing cutoff and line-selection margin: a
        margin narrower than the cutoff silently drops the wings of lines just
        outside the window, which mimics a too-short cutoff (the 2026-07-25
        finding).  Only checkable when the lines were actually windowed."""
        margin = getattr(lines, 'margin_cm', None)
        if margin is None:
            return
        need = float(line_cutoff_cm(self.nu_vac.max(), R=self.R, ncut=self.ncut,
                                    cutoff_cm=self.cutoff_cm))
        if margin < need:
            raise ValueError(
                'line-selection margin %.2f cm-1 is narrower than the wing '
                'cutoff %.2f cm-1: wings of out-of-window lines would be lost. '
                'Use required_margin_cm() to derive the margin.' % (margin, need))

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
                    self.R, self.ncut, self.cutoff_cm)
                self.od[gas][iz, :] = od_nu[self._inv]  # back to wvl order

        # O2-O2 CIA: a smooth continuum, kept separate from the line OT so the
        # two can always be inspected (and differenced) independently.
        if self.cia is None:
            self.od_cia = np.zeros((nlay, nwvl))
        else:
            from .cia import cal_o2_cia_od
            self.od_cia = cal_o2_cia_od(self.atm, self.nu_vac, self.cia)

        self.od_total = sum(self.od.values()) + self.od_cia
        self.od_column = self.od_total.sum(axis=0)

    def subrange_indices(self, wvl_range=None):
        """
        Indices of the grid falling within `wvl_range` (nm, in this object's
        grid convention -- see ``self.grid``); None → all.

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
        cut = ('%g cm-1' % self.cutoff_cm if self.cutoff_cm is not None
               else '%g*nu0/%g' % (self.ncut, self.R))
        return ('o2band_absorption(band=%s): %d %s-wvl (%.3f-%.3f nm @ %.4f nm), '
                '%d layers, gases=%s, cutoff=%s, cia=%s'
                % (self.band, self.wvl.size, self.grid, self.wvl[0], self.wvl[-1],
                   self.dwvl, self.atm.lay['z'].size, ','.join(self.gases), cut,
                   'on' if self.cia is not None else 'off'))


# ---------------------------------------------------------------------------- #
def cal_rayleigh_od(atm, wvl_nm):
    """
    Rayleigh-scattering optical depth per layer, Bodhaine et al. (1999).

    Input:
        atm    : afgl_atmosphere object (uses layer air column, molec cm-2)
        wvl_nm : **vacuum** wavelength(s) in nm (scalar or array).  Bodhaine's
                 fit is parameterised in vacuum wavelength; feeding it air
                 wavelengths biases sigma high by ~0.12% (sigma ~ lam^-4 and
                 lam_air is 0.03% low).  Pass ``absb.wvl_vac``.
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
        od_ray = cal_rayleigh_od(atm, absb.wvl_vac)
        print('  %s: max column O2+H2O OT = %.3f | column O2 OT max = %.3f | '
              'mean Rayleigh OT = %.4f'
              % (band, absb.od_column.max(), absb.od['o2'].sum(axis=0).max(),
                 od_ray.sum(axis=0).mean()))
