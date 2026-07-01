"""
Vacuum <-> air wavelength conversion.

HITRAN line positions are vacuum wavenumbers (cm-1).  The intercomparison output
grid is specified in **air** wavelengths, so we convert between the two with the
Edlen (1966) dispersion formula for standard air (dry, 15 C, 101.325 kPa,
0.03% CO2) -- the conventional definition of "air wavelength" in spectroscopy.

    n(sigma) - 1 = 1e-8 * (8342.54 + 2406147/(130 - sigma^2) + 15998/(38.9 - sigma^2))

with sigma = 1/lambda_vac in um^-1.  lambda_air = lambda_vac / n.
"""

import numpy as np

__all__ = ['refractive_index_air', 'vac_to_air_nm', 'air_to_vac_nm']


def refractive_index_air(wl_vac_nm):
    """Refractive index of standard air (Edlen 1966) at vacuum wavelength (nm)."""
    wl_um = np.asarray(wl_vac_nm, dtype=np.float64) * 1.0e-3
    sigma2 = (1.0 / wl_um) ** 2
    n = 1.0 + 1.0e-8 * (8342.54
                        + 2406147.0 / (130.0 - sigma2)
                        + 15998.0 / (38.9 - sigma2))
    return n


def vac_to_air_nm(wl_vac_nm):
    """Convert vacuum wavelength (nm) to air wavelength (nm)."""
    wl_vac_nm = np.asarray(wl_vac_nm, dtype=np.float64)
    return wl_vac_nm / refractive_index_air(wl_vac_nm)


def air_to_vac_nm(wl_air_nm, n_iter=3):
    """
    Convert air wavelength (nm) to vacuum wavelength (nm).

    n depends on the (unknown) vacuum wavelength, so iterate; the correction is
    ~3e-4 and converges in 1-2 passes.
    """
    wl_air_nm = np.asarray(wl_air_nm, dtype=np.float64)
    wl_vac = wl_air_nm.copy()  # initial guess
    for _ in range(n_iter):
        wl_vac = wl_air_nm * refractive_index_air(wl_vac)
    return wl_vac


if __name__ == '__main__':
    for wl in [688.0, 765.0]:
        wv = air_to_vac_nm(wl)
        print('air %.3f nm -> vac %.4f nm (shift %.4f nm), n=%.8f'
              % (wl, wv, wv - wl, refractive_index_air(wv)))
