"""
Read the CU composite solar reference spectrum (``CU_composite_solar.dat``)
and interpolate it onto a requested wavelength grid.

File format: two whitespace-separated columns

    wavelength(nm)   irradiance

Comment lines start with ``#``.  The spectrum runs from ~115 nm to ~200 um at
roughly 1 nm spacing, so it is a *smooth* solar continuum at the 0.001 nm scale
of the line-by-line grid -- it does not resolve solar Fraunhofer structure.
"""

import numpy as np


__all__ = ['solar_cu']


class solar_cu:

    """
    Input:
        fname : path to CU_composite_solar.dat

    Attributes:
        self.wvl  : wavelength (nm), ascending
        self.flux : solar irradiance at 1 AU (file units, ~W m^-2 nm^-1)

    Method:
        self.interp(wvl) : linear interpolation onto wvl (nm)
    """

    def __init__(self, fname):

        self.fname = fname
        data = np.genfromtxt(fname, comments='#')
        idx = np.argsort(data[:, 0])
        self.wvl = np.ascontiguousarray(data[idx, 0], dtype=np.float64)
        self.flux = np.ascontiguousarray(data[idx, 1], dtype=np.float64)

    def interp(self, wvl):
        """Interpolate the solar spectrum onto wvl (nm)."""
        return np.interp(wvl, self.wvl, self.flux)

    def __repr__(self):
        return ('solar_cu(%s): %d points, %.1f-%.1f nm'
                % (self.fname, self.wvl.size, self.wvl.min(), self.wvl.max()))


if __name__ == '__main__':

    import os
    here = os.path.dirname(os.path.abspath(__file__))
    fdir_data = os.path.join(here, '..', '..', 'data')
    sol = solar_cu(os.path.join(fdir_data, 'CU_composite_solar.dat'))
    print(sol)
    for wl in [688.0, 765.0]:
        print('  solar @ %.1f nm = %.4f' % (wl, sol.interp(wl)))
