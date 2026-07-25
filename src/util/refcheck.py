"""
Independent analytic reference formulas for cross-checking the MCARaTS
line-by-line reflectance / optical-depth output.

These are deliberately *simple closed-form* expressions (not another full RT
solve) so they provide an orthogonal sanity bound at the few-percent level:

- ``rayleigh_od_ht1974``      : Rayleigh optical depth vs the standard
                                Hansen & Travis (1974) polynomial in 1/lambda.
- ``rayleigh_ss_reflectance`` : single-scattering Rayleigh reflectance over a
                                black surface, in the rho = pi*I/(mu0*F0) sense.
- ``lambertian_over_atm``     : single-scattering atmosphere + Lambertian surface
                                coupling rho = rho_atm + T_dn*T_up*A/(1-A*s).

The single-scattering forms OMIT multiple scattering, so a correct full-MC
reflectance should sit slightly ABOVE the pure-Rayleigh single-scattering value
(multiple scattering only adds signal) and agree to a few percent for the thin
Rayleigh optical depths of the O2 A/B bands.  They are reference *bounds*, not the
cross-model benchmark (that is a second line-by-line RT model).
"""

import numpy as np

__all__ = ['rayleigh_od_ht1974', 'rayleigh_phase', 'scattering_cosine',
           'rayleigh_ss_reflectance', 'lambertian_over_atm']


def rayleigh_od_ht1974(wvl_nm):
    """Total (sea-level) Rayleigh optical depth, Hansen & Travis (1974):

        tau_R(lambda) = 0.008569 lam^-4 (1 + 0.0113 lam^-2 + 0.00013 lam^-4)

    with lam in micrometres.  Standard reference value for the whole atmosphere.
    """
    lam = np.asarray(wvl_nm, dtype=float) / 1000.0        # nm -> um
    return 0.008569 * lam**-4 * (1.0 + 0.0113 * lam**-2 + 0.00013 * lam**-4)


def rayleigh_phase(cos_theta):
    """Rayleigh scattering phase function, normalised so <P>_4pi = 1."""
    return 0.75 * (1.0 + cos_theta**2)


def scattering_cosine(mu, mu0, raa_deg=0.0):
    """cos of the scattering angle for view cosine mu, solar cosine mu0, and
    relative azimuth raa_deg.  Reflection geometry (both beams referenced to the
    upward/downward hemispheres)."""
    raa = np.deg2rad(raa_deg)
    return -mu * mu0 + np.sqrt(max(0.0, 1 - mu**2)) * np.sqrt(max(0.0, 1 - mu0**2)) * np.cos(raa)


def rayleigh_ss_reflectance(tau, mu0, mu=1.0, raa_deg=0.0, ssa=1.0):
    """Single-scattering reflectance of a homogeneous scattering layer over a
    black surface, in the rho = pi*I/(mu0*F0) convention:

        rho_ss = ssa*P(Theta)/(4(mu+mu0)) * [1 - exp(-tau(1/mu + 1/mu0))]

    For Rayleigh use ssa=1 and P = 3/4(1+cos^2 Theta) (the default via
    ``rayleigh_phase``).
    """
    cosT = scattering_cosine(mu, mu0, raa_deg)
    P = rayleigh_phase(cosT)
    return ssa * P / (4.0 * (mu + mu0)) * (1.0 - np.exp(-tau * (1.0 / mu + 1.0 / mu0)))


def lambertian_over_atm(tau, albedo, mu0, mu=1.0, raa_deg=0.0, spherical_albedo=None):
    """Approximate TOA reflectance of a Lambertian surface (albedo A) beneath a
    single-scattering Rayleigh atmosphere of optical depth tau:

        rho = rho_atm(tau) + T_dn(mu0) T_up(mu) A / (1 - A s)

    with direct-beam transmittances T = exp(-tau/mu) and a small atmospheric
    spherical albedo s (default ~tau/2 for a thin conservative layer).  Reproduces
    the correct limits: A=0 -> rho_atm; tau=0 -> A (a bare Lambertian surface).
    """
    if spherical_albedo is None:
        spherical_albedo = 0.5 * tau            # crude thin-layer estimate
    rho_atm = rayleigh_ss_reflectance(tau, mu0, mu, raa_deg)
    Tdn = np.exp(-tau / mu0)
    Tup = np.exp(-tau / mu)
    rho_surf = albedo * Tdn * Tup / (1.0 - albedo * spherical_albedo)
    return rho_atm + rho_surf
