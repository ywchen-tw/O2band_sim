"""
Per-g-point MCARaTS radiance reader for the line-by-line O2-band simulation.

Unlike er3t's :class:`mca_out_ng` (which weight-*sums* all g-points into a single
band radiance), here every g-point is an independent monochromatic wavelength, so
we read each ``fnames_out[ir][ig]`` file separately and assemble a *spectrum*.

Radiometry / conventions (see PLAN.md sec.7.5)
----------------------------------------------
MCARaTS is run with ``Src_flx = 1.0`` (unit collimated source, beam-perpendicular
irradiance).  Because MC transport is linear in the incident flux, the raw
per-g-point radiance is the response per unit incident irradiance:

    R_raw(lambda) = I(lambda) / F0            [sr^-1]

where F0 is the beam-perpendicular solar irradiance.  Therefore

    absolute radiance   I(lambda) = R_raw(lambda) * F0(lambda) * sol_fac
    TOA reflectance     rho(lambda) = pi * I / (mu0 * F0 * sol_fac)
                                    = pi * R_raw(lambda) / mu0

with mu0 = cos(SZA) and F0(lambda) the CU composite solar spectrum carried on
``abs.coef['solar']``.  Reflectance is F0-independent (the solar spectrum cancels);
it is retained only to also produce the absolute-radiance product.  ``sol_fac`` is
the sun-earth distance factor; it defaults to 1.0 (the CU spectrum is referenced to
1 AU and the benchmark prescribes no specific date), and cancels in reflectance.

The MC standard deviation across the ``Nrun`` independent runs is propagated to
both radiance and reflectance, and the standard error (std / sqrt(Nrun)) is also
reported.
"""

import os
import numpy as np

from er3t.rtm.mca.mca_out import mca_out_raw


__all__ = ['mca_out_lbl']


class mca_out_lbl:

    """
    Read per-g-point radiance from an MCARaTS ``radiance`` run and assemble the
    reflectance / radiance spectrum.

    Input:
        mca_obj : the mcarats_ng object that was run (provides fnames_out, Ng,
                  Nrun, solar_zenith_angle, target).
        abs_obj : the mca_abs_lbl used for the run (provides per-g air wavelengths
                  and coef['solar'] = F0(lambda)).
        sol_fac : sun-earth distance factor (default 1.0; cancels in reflectance).

    Attributes (all length Ng, ordered as abs_obj.wvls_air):
        wvl        : air wavelength (nm)
        f0         : incident solar irradiance F0 (W m-2 nm-1)
        r_raw      : mean raw radiance per unit incident flux (sr-1)
        r_raw_std  : std of raw radiance across runs (sr-1)
        rad        : absolute radiance I = r_raw * f0 * sol_fac (W m-2 nm-1 sr-1)
        rad_std    : std of absolute radiance across runs
        rad_stderr : standard error of absolute radiance (rad_std / sqrt(Nrun))
        ref        : TOA reflectance rho = pi * r_raw / mu0 (unitless)
        ref_std    : std of reflectance across runs
        ref_stderr : standard error of reflectance
        mu0        : cos(SZA)
    """

    def __init__(self, mca_obj, abs_obj, sol_fac=1.0):

        if mca_obj.target != 'radiance':
            msg = ('Error [mca_out_lbl]: expected target="radiance", got "%s".'
                   % mca_obj.target)
            raise ValueError(msg)

        Ng = mca_obj.Ng
        Nrun = mca_obj.Nrun

        if Ng != abs_obj.Ng:
            msg = ('Error [mca_out_lbl]: mca_obj.Ng (%d) != abs_obj.Ng (%d).'
                   % (Ng, abs_obj.Ng))
            raise ValueError(msg)

        self.mca = mca_obj
        self.abs = abs_obj
        self.Ng = Ng
        self.Nrun = Nrun
        self.sol_fac = float(sol_fac)

        self.wvl = np.asarray(abs_obj.wvls_air, dtype=np.float64)
        self.f0 = np.asarray(abs_obj.coef['solar']['data'], dtype=np.float64)

        self.mu0 = float(np.cos(np.deg2rad(mca_obj.solar_zenith_angle)))

        # raw radiance per (run, g): shape (Nrun, Ng)
        r = np.full((Nrun, Ng), np.nan, dtype=np.float64)
        for ir in range(Nrun):
            for ig in range(Ng):
                r[ir, ig] = self._read_one(mca_obj.fnames_out[ir][ig])

        # statistics across the Nrun independent runs
        self.r_raw = np.mean(r, axis=0)                 # (Ng,)
        self.r_raw_std = np.std(r, axis=0)              # population std across runs
        self._r = r                                     # keep raw (Nrun, Ng) for diagnostics

        # absolute radiance  I = R_raw * F0 * sol_fac
        scale_I = self.f0 * self.sol_fac
        self.rad = self.r_raw * scale_I
        self.rad_std = self.r_raw_std * scale_I
        self.rad_stderr = self.rad_std / np.sqrt(Nrun)

        # reflectance  rho = pi * R_raw / mu0   (F0 cancels)
        scale_ref = np.pi / self.mu0
        self.ref = self.r_raw * scale_ref
        self.ref_std = self.r_raw_std * scale_ref
        self.ref_stderr = self.ref_std / np.sqrt(Nrun)

    @staticmethod
    def _read_one(fname):
        """Read a single MCARaTS radiance binary; return the scalar TOA radiance.

        For the clear-sky single-column nadir satellite geometry (Nx=Ny=1, one
        radiance level), the output collapses to a single value.  We assert that
        so an unexpected shape fails loudly rather than being silently averaged.
        """
        if not os.path.isfile(fname):
            raise OSError('Error [mca_out_lbl]: missing output <%s>.' % fname)
        out = mca_out_raw(fname)
        val = np.squeeze(out.data[0]['data'])
        if val.size != 1:
            raise ValueError(
                'Error [mca_out_lbl]: expected a single radiance value from <%s>, '
                'got shape %s. Non-scalar radiance output is unexpected for the '
                'single-column nadir geometry.' % (fname, out.data[0]['data'].shape))
        return float(val)

    def valid(self):
        """True if every g-point produced a finite radiance (no NaNs)."""
        return bool(np.all(np.isfinite(self.r_raw)))

    def __repr__(self):
        return ('mca_out_lbl: Ng=%d, Nrun=%d, %.4f-%.4f nm, SZA=%.1f (mu0=%.4f), '
                'ref %.3e-%.3e'
                % (self.Ng, self.Nrun, self.wvl.min(), self.wvl.max(),
                   self.mca.solar_zenith_angle, self.mu0,
                   np.nanmin(self.ref), np.nanmax(self.ref)))
