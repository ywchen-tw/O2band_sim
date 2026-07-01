"""
Adapters that expose the line-by-line atmosphere / absorption to er3t's MCARaTS
1D pipeline (`er3t.rtm.mca.mca_atm_1d`).

Two objects are provided:

- :class:`mca_atm_lbl` -- an er3t-`atm`-compatible view of an
  :class:`util.atmosphere.afgl_atmosphere`, exposing the nested
  ``lev[...]['data']`` / ``lay[...]['data']`` dicts that `mca_atm_1d` and
  er3t's Rayleigh routine (`cal_mol_ext_atm`, method='atm') read.

- :class:`mca_abs_lbl` -- an er3t-`abs`-compatible object whose g-points ARE a
  set of fine (air-)wavelengths from an :class:`util.absorption.o2band_absorption`
  result.  Each g-point carries that wavelength's total gas absorption optical
  depth per layer; MCARaTS then solves one monochromatic RT per g-point.

Design notes
------------
- The atm adapter is built from the SAME `afgl_atmosphere` used for the
  absorption, so the RT layer grid, thicknesses and Rayleigh air column are
  mutually consistent with the absorption optical depths (verified in tests:
  RT Rayleigh == our `cal_rayleigh_od`).
- `mca_atm_1d` derives Rayleigh from a single scalar `abs.wvl`.  A `mca_abs_lbl`
  therefore represents a *narrow* wavelength chunk and reports `wvl` = the chunk
  mean; the driver (sim step 2) chunks each band finely (e.g. ~1 nm) so Rayleigh
  is effectively per-wavelength correct.  The per-wavelength O2 / Rayleigh OT for
  the intercomparison output come from the absorption object directly, not the RT.
- `solar` and `weight` are set to 1 for every g-point; per-g radiances are read
  back individually (not weight-summed), and reflectance normalisation is applied
  downstream.  Absolute solar irradiance is not needed for reflectance.
"""

import copy
import numpy as np

from .absorption import cal_rayleigh_od


__all__ = ['mca_atm_lbl', 'mca_abs_lbl', 'set_per_g_rayleigh']


def _d(data, name, units):
    return {'data': np.asarray(data), 'name': name, 'units': units}


class mca_atm_lbl:

    """
    er3t-`atm`-compatible view of an afgl_atmosphere.

    Exposes exactly the fields `mca_atm_1d` and `cal_mol_ext_atm(method='atm')`
    read:
        lev['altitude'](km), lev['pressure'](hPa)
        lay['altitude'](km), lay['thickness'](km), lay['temperature'](K),
        lay['air'](cm-3 number density), lay['co2'](cm-3 number density)

    The afgl_atmosphere stores gas *columns* (molec cm-2); here they are turned
    back into layer-mean number densities (column / dz) so that er3t's
    density*thickness recovers exactly our column amount.
    """

    ID = 'Atmosphere 1D (line-by-line, AFGL)'

    def __init__(self, atm, lat=0.0):

        self.afgl = atm
        self.lat = lat

        lev = atm.lev
        lay = atm.lay
        dz_cm = lay['dz'] * 1.0e5  # km -> cm

        self.lev = {
            'altitude':    _d(lev['z'],           'Altitude', 'km'),
            'pressure':    _d(lev['p'],           'Pressure', 'hPa'),
            'temperature': _d(lev['temperature'], 'Temperature', 'K'),
        }

        self.lay = {
            'altitude':    _d(lay['z'],           'Altitude', 'km'),
            'thickness':   _d(lay['dz'],          'Thickness', 'km'),
            'temperature': _d(lay['temperature'], 'Temperature', 'K'),
            'pressure':    _d(lay['p'],           'Pressure', 'hPa'),
            # columns (cm-2) -> layer-mean number density (cm-3)
            'air':         _d(lay['air'] / dz_cm, 'Air number density', 'cm-3'),
            'co2':         _d(lay['co2'] / dz_cm, 'CO2 number density', 'cm-3'),
        }

    def __repr__(self):
        return ('mca_atm_lbl: %d layers, z %.1f-%.1f km'
                % (self.lay['altitude']['data'].size,
                   self.lev['altitude']['data'][0], self.lev['altitude']['data'][-1]))


class mca_abs_lbl:

    """
    er3t-`abs`-compatible object over a set of fine wavelengths (g-points).

    Input:
        absb  : util.absorption.o2band_absorption object
        idx   : indices into absb.wvl that become g-points (default: all)
        solar : optional util.solar.solar_cu object.  If given, coef['solar'] is
                the incident solar irradiance F0(lambda) [W m-2 nm-1] interpolated
                at each g-point's air wavelength (used to turn the unit-source RT
                radiance into absolute radiance downstream).  If None, coef['solar']
                is unity (reflectance-only; ratio is F0-independent anyway).

    Attributes required by mca_atm_1d / mca_out_ng:
        Ng, wvl (scalar nm, chunk mean, for Rayleigh), wvl_info,
        coef['abso_coef']['data'] (Nz, Ng)  -- total gas absorption OT per layer,
        coef['weight']['data'] (Ng),
        coef['solar']['data']  (Ng)  -- F0(lambda) or ones,
        coef['slit_func']['data'] (Nz, Ng)

    Extra bookkeeping (for output, not used by er3t):
        wvls_air (Ng, nm), nu_vac (Ng, cm-1), idx (Ng,), f0 (Ng, solar irradiance)
    """

    def __init__(self, absb, idx=None, solar=None):

        if idx is None:
            idx = np.arange(absb.wvl.size)
        idx = np.asarray(idx, dtype=int)

        self.absb = absb
        self.idx = idx
        self.Ng = idx.size

        self.wvls_air = absb.wvl[idx]
        self.nu_vac = absb.nu_vac[idx]
        # single scalar wavelength for Rayleigh: chunk mean (must be a narrow chunk)
        self.wvl = float(self.wvls_air.mean())
        self.wvl_info = ('LBL %s: Ng=%d, %.4f-%.4f nm (air)'
                         % (absb.band, self.Ng, self.wvls_air.min(), self.wvls_air.max()))

        nz = absb.od_total.shape[0]
        abso = absb.od_total[:, idx]          # (Nz, Ng), total gas absorption OT

        # incident solar irradiance F0(lambda) at each g-point (air wvl), or unity.
        # MCARaTS runs with Src_flx=1 (unit source); folding F0 in post gives
        # absolute radiance because MC transport is linear in incident flux.
        self.solar = solar
        if solar is None:
            f0 = np.ones(self.Ng, dtype=np.float64)
        else:
            f0 = np.asarray(solar.interp(self.wvls_air), dtype=np.float64)
        self.f0 = f0

        self.coef = {
            'wvl':       {'name': 'Wavelength', 'units': 'nm', 'data': self.wvls_air.copy()},
            'abso_coef': {'name': 'Absorption optical depth (Nz, Ng)',
                          'data': np.ascontiguousarray(abso, dtype=np.float64)},
            'weight':    {'name': 'Weight (Ng)', 'data': np.ones(self.Ng, dtype=np.float64)},
            'solar':     {'name': 'Incident solar irradiance F0 (Ng)', 'units': 'W m-2 nm-1',
                          'data': f0.copy()},
            'slit_func': {'name': 'Slit function (Nz, Ng)',
                          'data': np.ones((nz, self.Ng), dtype=np.float64)},
        }

    def rayleigh_od(self):
        """
        Per-layer Rayleigh OT at *each* g-point's own wavelength, shape (Nz, Ng).
        Exact at the 0.001 nm grid (no chunk averaging).
        """
        return cal_rayleigh_od(self.absb.atm, self.wvls_air)

    def __repr__(self):
        return ('mca_abs_lbl(band=%s): Ng=%d, %.4f-%.4f nm (air), wvl(Rayleigh)=%.4f nm'
                % (self.absb.band, self.Ng, self.wvls_air.min(),
                   self.wvls_air.max(), self.wvl))


def set_per_g_rayleigh(atm1d, mca_atm, mca_abs):
    """
    Override each g-point's Rayleigh extinction with the exact value at that
    g-point's own wavelength.

    `mca_atm_1d` builds `Atm_ext1d` (the Rayleigh scattering field, omega=1,
    Rayleigh phase function) from a single scalar `abs.wvl` shared by all
    g-points -- accurate only over a narrow chunk.  Because Rayleigh is the sole
    wavelength-dependent field shared across g-points (gas absorption is already
    per-g), replacing it per g-point makes the RT Rayleigh exact at every
    0.001 nm grid point, matching the Rayleigh OT we output.  This removes the
    physics reason for narrow chunking; chunk size then becomes purely a
    batching choice (process/file count per MCARaTS call).

    Input:
        atm1d   : an er3t `mca_atm_1d` object built from (mca_atm, mca_abs)
        mca_atm : the mca_atm_lbl used to build it (for layer thickness)
        mca_abs : the mca_abs_lbl (provides per-g wavelengths / Rayleigh OT)
    Returns the same atm1d, modified in place.
    """
    ray_od = mca_abs.rayleigh_od()                              # (Nz, Ng) OT
    thick_m = mca_atm.lay['thickness']['data'] * 1000.0         # km -> m
    for ig in range(mca_abs.Ng):
        atm1d.nml[ig]['Atm_ext1d(1:, 1)']['data'] = ray_od[:, ig] / thick_m
    return atm1d
