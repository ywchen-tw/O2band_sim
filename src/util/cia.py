"""
O2-O2 collision-induced absorption (CIA) from a HITRAN ``.cia`` file.

CIA is a *binary* absorption process: two O2 molecules absorb during a
collision, so the absorption coefficient scales with the SQUARE of the O2
number density rather than linearly like a line spectrum:

    alpha(nu) = k_CIA(nu, T) * n_O2^2        [cm-1]        (k in cm^5 molec^-2)
    tau_layer = integral alpha dz = k_CIA(nu, T) * integral n_O2^2 dz

That n^2 weighting concentrates CIA far more strongly near the surface than the
line optical depth, so the layer column of n^2 must be integrated properly (the
same exponential-in-height treatment used for the linear columns in
:mod:`util.atmosphere`), not approximated by (integral n dz)^2 / dz.

HITRAN .cia format
------------------
A file is a concatenation of segments.  Each begins with a fixed-width header:

    cols  1-20  molecule pair ('O2-O2')
    cols 21-30  nu_min (cm-1)
    cols 31-40  nu_max (cm-1)
    cols 41-47  number of points
    cols 48-54  temperature (K)
    cols 55-64  max CIA value (cm^5 molec^-2)
    cols 65-74  resolution / uncertainty
    cols 75+    comment + reference

followed by that many "nu  k" lines.  One file holds many segments covering
different spectral regions and temperatures.

Phase-1 status: CIA is EXCLUDED from the prescribed benchmark (PLAN.md §7.4).
This module exists so the omission can be quantified and so a later phase can
switch it on explicitly; using it is a documented deviation.

Temperature handling
--------------------
Where a band is covered by several temperatures, k is linearly interpolated in
T between the two bracketing tables and held constant outside their range.
Where a band has only ONE temperature (the case for the O2 A- and B-bands in
HITRAN's O2-O2 files), that single table is used at every layer temperature and
``temperature_tables`` reports it, so the approximation is visible rather than
silent.
"""

import os
import numpy as np

from .atmosphere import afgl_atmosphere


__all__ = ['hitran_cia', 'cal_o2_cia_od', 'o2_squared_column', 'pick_cia_file']


class hitran_cia:

    """
    Reader / evaluator for a HITRAN ``.cia`` file.

    Input:
        fname : path to the .cia file (e.g. data/O2-O2_2024.cia)
        pair  : molecule pair to keep, default 'O2-O2'

    Attributes:
        segments : list of dicts with keys nu_min, nu_max, npts, T, nu, k, ref
        pair     : the molecule pair

    Methods:
        covers(nu_lo, nu_hi)      : True if some segment spans the whole range
        temperature_tables(...)   : temperatures available over a range
        k(nu, T)                  : CIA coefficient (cm^5 molec^-2), 0 outside
                                    any covering segment
    """

    def __init__(self, fname, pair='O2-O2'):
        self.fname = fname
        self.pair = pair
        self.segments = []
        self._read(fname, pair)

    def _read(self, fname, pair):
        with open(fname) as f:
            lines = f.readlines()

        i, n = 0, len(lines)
        while i < n:
            ln = lines[i]
            if ln[:20].strip() != pair:
                i += 1
                continue
            try:
                nu_min = float(ln[20:30])
                nu_max = float(ln[30:40])
                npts = int(float(ln[40:47]))
                temperature = float(ln[47:54])
            except ValueError:
                i += 1
                continue
            body = lines[i + 1:i + 1 + npts]
            if len(body) < npts:
                break                                   # truncated file
            arr = np.array([[float(x) for x in b.split()[:2]] for b in body],
                           dtype=np.float64)
            self.segments.append(dict(
                nu_min=nu_min, nu_max=nu_max, npts=npts, temperature=temperature,
                nu=np.ascontiguousarray(arr[:, 0]), k=np.ascontiguousarray(arr[:, 1]),
                ref=ln[74:].strip()))
            i += 1 + npts

    # ------------------------------------------------------------------ #
    def _covering(self, nu_lo, nu_hi):
        """Segments spanning the whole [nu_lo, nu_hi] range, sorted by T."""
        segs = [s for s in self.segments
                if s['nu_min'] <= nu_lo and s['nu_max'] >= nu_hi]
        return sorted(segs, key=lambda s: s['temperature'])

    def covers(self, nu_lo, nu_hi):
        return len(self._covering(min(nu_lo, nu_hi), max(nu_lo, nu_hi))) > 0

    def temperature_tables(self, nu_lo, nu_hi):
        return [s['temperature']
                for s in self._covering(min(nu_lo, nu_hi), max(nu_lo, nu_hi))]

    def k(self, nu, temperature):
        """
        CIA coefficient k(nu, T) in cm^5 molec^-2 on the vacuum wavenumbers
        ``nu``.  Returns zeros where no segment covers the range.

        Linear in T between bracketing tables, clamped outside their range (the
        tables are sparse in T; extrapolating a binary-collision coefficient
        beyond the measured range is not defensible).
        """
        nu = np.atleast_1d(np.asarray(nu, dtype=np.float64))
        segs = self._covering(nu.min(), nu.max())
        if not segs:
            return np.zeros_like(nu)

        temps = np.array([s['temperature'] for s in segs])
        if temps.size == 1 or temperature <= temps[0]:
            return np.interp(nu, segs[0]['nu'], segs[0]['k'])
        if temperature >= temps[-1]:
            return np.interp(nu, segs[-1]['nu'], segs[-1]['k'])

        j = int(np.searchsorted(temps, temperature))
        lo, hi = segs[j - 1], segs[j]
        w = (temperature - lo['temperature']) / (hi['temperature'] - lo['temperature'])
        k_lo = np.interp(nu, lo['nu'], lo['k'])
        k_hi = np.interp(nu, hi['nu'], hi['k'])
        return (1.0 - w) * k_lo + w * k_hi

    def __repr__(self):
        return ('hitran_cia(%s, %s): %d segments, %.1f-%.1f cm-1, T %s'
                % (os.path.basename(self.fname), self.pair, len(self.segments),
                   min(s['nu_min'] for s in self.segments) if self.segments else 0,
                   max(s['nu_max'] for s in self.segments) if self.segments else 0,
                   sorted({s['temperature'] for s in self.segments})))


# ---------------------------------------------------------------------------- #
def pick_cia_file(data_dir, nu_ranges, pattern='O2-O2_*.cia'):
    """
    Choose the CIA file in ``data_dir`` that actually covers every
    (nu_lo, nu_hi) range in ``nu_ranges``.

    Selection is by COVERAGE, not by filename: HITRAN ships several O2-O2 files
    whose spectral spans differ wildly (``O2-O2_exp_2021.cia`` is the 1300-1850
    cm-1 fundamental and covers neither O2 band), so picking the
    alphabetically-last one silently gets it wrong.  Ties go to the newest file.
    """
    import glob
    hits = sorted(glob.glob(os.path.join(data_dir, pattern)))
    if not hits:
        raise OSError('no %s found in %s' % (pattern, data_dir))
    good = []
    for fname in hits:
        try:
            c = hitran_cia(fname)
        except (OSError, ValueError):
            continue
        if all(c.covers(lo, hi) for lo, hi in nu_ranges):
            good.append(fname)
    if not good:
        raise OSError('none of %s covers %s cm-1 (checked: %s)'
                      % (pattern, nu_ranges, ', '.join(os.path.basename(h)
                                                       for h in hits)))
    return max(good, key=lambda f: os.stat(f).st_mtime)


def o2_squared_column(atm):
    """
    Layer column of the SQUARED O2 number density, integral n_O2^2 dz
    (molec^2 cm^-5), for the binary CIA optical depth.

    n_O2 is taken to vary exponentially between levels -- the same assumption
    :mod:`util.atmosphere` uses for the linear columns -- so n^2 is exponential
    too and the exact integral is available in closed form.  Using the linear
    column instead would misplace CIA in the vertical (its true scale height is
    half that of n).
    """
    lev = atm.lev
    dz_cm = (lev['z'][1:] - lev['z'][:-1]) * 1.0e5
    return afgl_atmosphere._column_exponential(lev['o2'][:-1] ** 2,
                                               lev['o2'][1:] ** 2, dz_cm)


def cal_o2_cia_od(atm, nu_vac, cia, strict=True):
    """
    O2-O2 CIA optical depth per layer on a vacuum-wavenumber grid.

    Input:
        atm    : afgl_atmosphere object
        nu_vac : vacuum wavenumbers (cm-1), shape (Nwvl,)
        cia    : hitran_cia object
        strict : raise if the file has no segment covering the whole grid
                 (default).  False returns zeros instead.
    Output:
        od_cia : (Nlay, Nwvl) optical depth
    """
    nu_vac = np.asarray(nu_vac, dtype=np.float64)
    nu_lo, nu_hi = float(nu_vac.min()), float(nu_vac.max())
    if not cia.covers(nu_lo, nu_hi):
        msg = ('CIA file %s has no O2-O2 segment covering %.1f-%.1f cm-1'
               % (os.path.basename(cia.fname), nu_lo, nu_hi))
        if strict:
            raise ValueError(msg)
        return np.zeros((atm.lay['z'].size, nu_vac.size))

    n2_col = o2_squared_column(atm)                    # (Nlay,) molec^2 cm^-5
    temperature = atm.lay['temperature']
    od = np.empty((n2_col.size, nu_vac.size))
    for iz in range(n2_col.size):
        od[iz, :] = cia.k(nu_vac, temperature[iz]) * n2_col[iz]
    return od


# ---------------------------------------------------------------------------- #
if __name__ == '__main__':

    here = os.path.dirname(os.path.abspath(__file__))
    fdir_data = os.path.normpath(os.path.join(here, '..', '..', 'data'))
    from .absorption import BANDS, band_nu_range

    atm = afgl_atmosphere(os.path.join(fdir_data, 'afglms.dat'), z_top=120.0)
    cia = hitran_cia(os.path.join(fdir_data, 'O2-O2_2024.cia'))
    print(cia)
    for band in ('o2a', 'o2b'):
        nu_lo, nu_hi = band_nu_range(BANDS[band], grid='vac')
        print('%s: %.1f-%.1f cm-1  covered=%s  T tables=%s'
              % (band, nu_lo, nu_hi, cia.covers(nu_lo, nu_hi),
                 cia.temperature_tables(nu_lo, nu_hi)))
        if cia.covers(nu_lo, nu_hi):
            nu = np.linspace(nu_lo, nu_hi, 501)
            od = cal_o2_cia_od(atm, nu, cia)
            print('   column CIA OT: min %.5f  mean %.5f  max %.5f'
                  % (od.sum(axis=0).min(), od.sum(axis=0).mean(),
                     od.sum(axis=0).max()))
