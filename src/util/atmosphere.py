"""
Read an AFGL atmospheric-constituent profile (e.g. mid-latitude summer,
``afglms.dat``) and build a layered atmosphere for line-by-line absorption.

The AFGL file columns are:

    z(km)  p(mb)  T(K)  air(cm-3)  o3(cm-3)  o2(cm-3)  h2o(cm-3)  co2(cm-3)  no2(cm-3)

Number densities are given *at levels* (interfaces).  For a line-by-line
calculation we need, per layer:

    - a representative pressure and temperature (for the line shapes), and
    - the vertical column amount of each gas (molec cm^-2), for the optical depth.

Column amounts are integrated assuming each gas number density varies
exponentially with height between two levels (exact integral of an exponential),
which is far more accurate than a linear layer-mean for the near-surface layers
where the densities change quickly.
"""

import numpy as np


__all__ = ['afgl_atmosphere']


# gas columns present in a standard AFGL *.dat file (after z, p, T)
AFGL_COLUMNS = ['z', 'p', 'temperature', 'air', 'o3', 'o2', 'h2o', 'co2', 'no2']
GASES = ['air', 'o3', 'o2', 'h2o', 'co2', 'no2']


class afgl_atmosphere:

    """
    Input:
        fname   : path to an AFGL profile file (e.g. data/afglms.dat)
        z_top   : keyword, float, only keep levels with z <= z_top (km), default 120
        levels  : keyword, optional 1-D array of level altitudes (km, ascending)
                  to re-grid onto.  If None, the native AFGL levels are used.

    Output (all ascending in altitude):
        self.lev : dict of level quantities
                   'z'  (km), 'p' (hPa), 'temperature' (K),
                   'air','o2','h2o','o3','co2','no2' (cm-3)
        self.lay : dict of layer quantities (size = n_lev - 1)
                   'z'  (km, midpoint), 'dz' (km),
                   'p'  (hPa, log-mean), 'temperature' (K, mean),
                   'air','o2','h2o','o3','co2','no2' : column amount (molec cm-2)
                   '<gas>_vmr' : layer-mean volume mixing ratio (unitless)
    """

    def __init__(self, fname, z_top=120.0, levels=None):

        self.fname = fname
        self._read(fname)

        if levels is not None:
            self._regrid(np.asarray(levels, dtype=np.float64))

        # keep only levels within z_top
        keep = self.lev['z'] <= (z_top + 1e-6)
        for key in self.lev:
            self.lev[key] = self.lev[key][keep]

        self._build_layers()

    # ------------------------------------------------------------------ #
    def _read(self, fname):

        data = np.genfromtxt(fname, comments='#')
        # sort ascending in altitude
        data = data[np.argsort(data[:, 0])]

        self.lev = {}
        for i, name in enumerate(AFGL_COLUMNS):
            self.lev[name] = np.ascontiguousarray(data[:, i], dtype=np.float64)

    # ------------------------------------------------------------------ #
    def _regrid(self, new_z):
        """Interpolate the profile onto a new set of level altitudes (km)."""

        z0 = self.lev['z']
        lev = {'z': new_z}
        for key in self.lev:
            if key == 'z':
                continue
            if key in ('p',) or key in GASES:
                # log-interpolate quantities that vary exponentially
                lev[key] = np.exp(np.interp(new_z, z0, np.log(self.lev[key])))
            else:
                lev[key] = np.interp(new_z, z0, self.lev[key])
        self.lev = lev

    # ------------------------------------------------------------------ #
    @staticmethod
    def _column_exponential(n_bot, n_top, dz_cm):
        """
        Column amount (molec cm^-2) through a layer, assuming the number
        density varies exponentially between the bottom and top levels:

            N = integral n dz = dz * (n_bot - n_top) / ln(n_bot / n_top)

        Falls back to the arithmetic mean where the exponential form is
        ill-conditioned (near-equal or non-positive densities).
        """
        n_bot = np.asarray(n_bot, dtype=np.float64)
        n_top = np.asarray(n_top, dtype=np.float64)
        col = 0.5 * (n_bot + n_top) * dz_cm  # default: layer-mean
        good = (n_bot > 0) & (n_top > 0) & (np.abs(n_bot - n_top) > 1e-30 * n_bot)
        ratio = np.where(good, n_bot / np.where(n_top > 0, n_top, 1.0), 1.0)
        col_exp = dz_cm * (n_bot - n_top) / np.log(np.where(good, ratio, np.e))
        return np.where(good, col_exp, col)

    # ------------------------------------------------------------------ #
    def _build_layers(self):

        lev = self.lev
        z = lev['z']
        dz_km = z[1:] - z[:-1]
        dz_cm = dz_km * 1.0e5

        lay = {}
        lay['z'] = 0.5 * (z[1:] + z[:-1])
        lay['dz'] = dz_km

        # representative p (log-mean) and T (arithmetic mean) for line shapes
        lay['p'] = np.sqrt(lev['p'][1:] * lev['p'][:-1])
        lay['temperature'] = 0.5 * (lev['temperature'][1:] + lev['temperature'][:-1])

        # column amounts (molec cm-2) for each gas
        for gas in GASES:
            lay[gas] = self._column_exponential(lev[gas][:-1], lev[gas][1:], dz_cm)

        # layer-mean volume mixing ratio relative to total air
        for gas in GASES:
            if gas == 'air':
                continue
            lay['%s_vmr' % gas] = lay[gas] / lay['air']

        self.lay = lay

    # ------------------------------------------------------------------ #
    def __repr__(self):
        return ('afgl_atmosphere(%s): %d levels, %d layers, '
                'surface p=%.1f hPa T=%.1f K'
                % (self.fname, self.lev['z'].size, self.lay['z'].size,
                   self.lev['p'][0], self.lev['temperature'][0]))


if __name__ == '__main__':

    import os
    here = os.path.dirname(os.path.abspath(__file__))
    fdir_data = os.path.join(here, '..', '..', 'data')
    atm = afgl_atmosphere(os.path.join(fdir_data, 'afglms.dat'), z_top=70.0)
    print(atm)
    print('O2 total column  : %.4e molec cm-2' % atm.lay['o2'].sum())
    print('H2O total column : %.4e molec cm-2' % atm.lay['h2o'].sum())
    print('surface O2 vmr   : %.5f' % atm.lay['o2_vmr'][0])
