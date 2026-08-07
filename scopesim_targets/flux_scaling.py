# -*- coding: utf-8 -*-
"""Flux-authority scaling: turn a :class:`Brightness` into a spectrum scale.

This is the single place that answers "by what factor must this spectrum be
multiplied so that its synthetic photometry at the locator equals the stated
amount?", for every branch of the brightness grammar (magnitude in any
photometric system, spectral flux density per frequency or wavelength, and
band-integrated energy flux) at any locator (band, wavelength or frequency).

It deliberately depends on synphot + astropy only. The two inputs that need the
network / heavier stack -- the bandpass (a spextra ``Passband``) and the Vega
reference spectrum -- are resolved by the caller
(:meth:`SpectrumTarget._get_spectrum_scale`) and handed in, so this whole
dispatch is unit-testable offline against real synphot.
"""

import astropy.units as u
from synphot import Observation, units
from synphot.units import VEGAMAG

from .brightness import (
    AmountKind,
    LocatorKind,
    PhotometricSystem,
    AmountError,
)

__all__ = ["synphot_flux_scale"]


_SYSTEM_UNIT = {
    PhotometricSystem.VEGA: VEGAMAG,
    PhotometricSystem.AB: u.ABmag,
    PhotometricSystem.ST: u.STmag,
}


def synphot_flux_scale(spectrum, brightness, *, band=None, vegaspec=None):
    """Dimensionless factor so that `spectrum * factor` matches `brightness`.

    Parameters
    ----------
    spectrum : synphot.SourceSpectrum
        The resolved spectrum to be scaled.
    brightness : .brightness.Brightness
        A normalized, **integrated** brightness. Surface brightness must have
        been reduced to an integrated amount by the caller (profiles do this in
        :meth:`_effective_integrated_brightness`); this routine never sees a
        ``solid_angle`` and does not check for one.
    band : synphot.SpectralElement, optional
        Bandpass for a BAND locator. Required for BAND locators and for any
        magnitude amount.
    vegaspec : synphot.SourceSpectrum, optional
        Vega reference spectrum; required only for VEGA-system magnitudes.

    Returns
    -------
    float
        The multiplicative scale factor.
    """
    # -- magnitudes: one identity across all three systems ----------------
    if brightness.amount_kind is AmountKind.MAG:
        if band is None:
            raise ValueError("a magnitude amount requires a resolved band")
        system_unit = _SYSTEM_UNIT[brightness.system]
        extra = (
            {"vegaspec": vegaspec}
            if brightness.system is PhotometricSystem.VEGA
            else {}
        )
        m_actual = Observation(spectrum, band).effstim(system_unit, **extra)
        delta = brightness.value.to_value(u.mag) - m_actual.to_value(system_unit)
        return float(10 ** (-0.4 * delta))

    # -- band locator: photometry through the passband --------------------
    if brightness.locator_kind is LocatorKind.BAND:
        if band is None:
            raise ValueError("a band locator requires a resolved band")
        obs = Observation(spectrum, band)
        if brightness.amount_kind is AmountKind.ENERGY_FLUX:
            # band-integrated energy flux: the integral of f_lambda * P(lambda)
            actual = obs.integrate(flux_unit=units.FLAM)
        else:
            # band-averaged spectral flux density, measured in the amount's own
            # unit so the ratio comes out dimensionless.
            actual = obs.effstim(brightness.value.unit)
        return float((brightness.value / actual).to_value(u.dimensionless_unscaled))

    # -- monochromatic locator: evaluate the spectrum at the point --------
    if brightness.amount_kind is AmountKind.ENERGY_FLUX:
        # An energy flux is band-/line-integrated; at a single point it is a
        # density, so the pairing is undefined.
        raise AmountError(
            "Energy flux at a monochromatic locator is undefined. Use a band "
            "locator, or give a spectral flux density",
        )

    wavelength = brightness.locator.to(u.AA, u.spectral())
    actual = spectrum(wavelength, flux_unit=brightness.value.unit)
    return float((brightness.value / actual).to_value(u.dimensionless_unscaled))
