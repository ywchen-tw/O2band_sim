"""
Solar spectrum models.

``solar_cu`` reads the CU composite solar reference spectrum
(``CU_composite_solar.dat``) and interpolates it onto a requested wavelength
grid.  File format: two whitespace-separated columns

    wavelength(nm)   irradiance

Comment lines start with ``#``.  The spectrum runs from ~115 nm to ~200 um at
roughly 1 nm spacing, so it is a *smooth* solar continuum at the 0.001 nm scale
of the line-by-line grid -- it does not resolve solar Fraunhofer structure.

``solar_spts`` multiplies that continuum by the Toon (JPL) disk-integrated
Solar Pseudo-Transmittance Spectrum (SPTS, mark4sun.jpl.nasa.gov/toon/solar),
which resolves the Fraunhofer lines on a 0.01 cm-1 vacuum-wavenumber grid
(~0.0005 nm in the O2 bands).  Same ``.interp(wvl_air_nm)`` interface, so the
two are drop-in interchangeable as the pipeline's F0(lambda) source.
"""

import os

import numpy as np

from .optics import air_to_vac_nm


__all__ = ['solar_cu', 'solar_spts']


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


class solar_spts:

    """
    High-resolution solar irradiance: CU composite continuum x Toon SPTS
    disk-integrated solar transmittance.

    Input:
        fname_continuum : path to CU_composite_solar.dat
        fname_spts      : path to a ``solar_merged_*.out`` SPTS file
                          (3 header lines, then two columns: vacuum wavenumber
                          cm-1 ascending on a uniform grid, transmittance)

    Attributes:
        self.continuum : the wrapped solar_cu object
        self.nu        : SPTS vacuum wavenumber grid (cm-1, ascending)
        self.trans     : solar transmittance on that grid

    Methods:
        self.transmittance(wvl_air_nm) : SPTS transmittance at air wavelengths
        self.interp(wvl_air_nm)        : F0 = continuum x transmittance
    """

    def __init__(self, fname_continuum, fname_spts):

        self.continuum = solar_cu(fname_continuum)
        self.fname = fname_spts
        self.nu, self.trans = self._load(fname_spts)
        if np.any(np.diff(self.nu) <= 0.0):
            raise ValueError('Error [solar_spts]: wavenumber grid in <%s> is '
                             'not strictly ascending.' % fname_spts)

    @staticmethod
    def _load(fname):
        """Parse the SPTS text file, memoized as <fname>.npy (the text file is
        ~3.3M rows; the cache turns a ~30 s parse into a ~0.1 s load)."""
        cache = fname + '.npy'
        if (os.path.isfile(cache)
                and os.path.getmtime(cache) >= os.path.getmtime(fname)):
            data = np.load(cache)
        else:
            data = np.loadtxt(fname, skiprows=3, dtype=np.float64)
            tmp = '%s.tmp%d' % (cache, os.getpid())
            try:
                with open(tmp, 'wb') as f:
                    np.save(f, data)
                os.replace(tmp, cache)
            except OSError:
                pass                      # read-only data dir: just skip caching
        return (np.ascontiguousarray(data[:, 0]),
                np.ascontiguousarray(data[:, 1]))

    def transmittance(self, wvl):
        """SPTS solar transmittance at air wavelength(s) wvl (nm).

        The pipeline grid is air nm; SPTS is vacuum wavenumber, so convert
        with the same Edlen formula used for the HITRAN line positions.
        """
        nu = 1.0e7 / air_to_vac_nm(wvl)
        return np.interp(nu, self.nu, self.trans)

    def interp(self, wvl):
        """Solar irradiance F0 (continuum x transmittance) at air wvl (nm)."""
        return self.continuum.interp(wvl) * self.transmittance(wvl)

    def __repr__(self):
        return ('solar_spts(%s x %s): %d points, %.2f-%.2f cm-1'
                % (self.continuum.fname, self.fname, self.nu.size,
                   self.nu.min(), self.nu.max()))


if __name__ == '__main__':

    here = os.path.dirname(os.path.abspath(__file__))
    fdir_data = os.path.join(here, '..', '..', 'data')
    sol = solar_cu(os.path.join(fdir_data, 'CU_composite_solar.dat'))
    print(sol)
    for wl in [688.0, 765.0]:
        print('  solar @ %.1f nm = %.4f' % (wl, sol.interp(wl)))
