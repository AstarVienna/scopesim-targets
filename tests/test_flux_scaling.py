# -*- coding: utf-8 -*-
"""Unit tests for flux_scaling.py.

Pure: synphot + astropy only, no spextra/scopesim, no network. The bandpass and
the Vega reference that the production caller resolves over the network are
supplied here as synthetic fixtures -- a boxcar ``SpectralElement`` and a
constant reference spectrum -- so the whole flux-authority dispatch is exercised
offline against real synphot.
"""

import numpy as np
from numpy import testing as npt
import pytest
import astropy.units as u
from synphot import SourceSpectrum, SpectralElement, Observation, units
from synphot.models import Box1D, ConstFlux1D
from synphot.units import VEGAMAG
from spextra import Spextrum

from scopesim_targets.brightness import (
    parse_brightness,
    LocatorKind,
    AmountKind,
    PhotometricSystem,
    BrightnessError,
)
from scopesim_targets.flux_scaling import synphot_flux_scale


@pytest.fixture
def flat_spec():
    """Simple flat spectrum for testing, produces ~10 PHOTLAM @ 550 nm."""
    return Spextrum.flat_spectrum(
        5*u.ABmag,
        waves=np.linspace(3000, 25000, 6000)*u.AA,
    )


@pytest.fixture
def band():
    # Boxcar "V-ish" passband, fully inside the flat_spec's wavelength range.
    return SpectralElement(
        Box1D, amplitude=1, x_0=5500*u.AA, width=1000*u.AA
    )


@pytest.fixture
def vega():
    # A stand-in reference; the *identity* under test is independent of which
    # spectrum plays Vega, as long as the same one sets and reads the zero point.
    return SourceSpectrum(ConstFlux1D, amplitude=1*units.FLAM)


def _measured_amount(spectrum, brightness, band=None, vega=None):
    """Re-measure the stated amount on a spectrum (the inverse of the scale)."""
    if brightness.amount_kind is AmountKind.MAG:
        if brightness.system is PhotometricSystem.VEGA:
            return Observation(spectrum, band).effstim(VEGAMAG, vegaspec=vega)

        unit = {
            PhotometricSystem.AB: u.ABmag,
            PhotometricSystem.ST: u.STmag,
        }[brightness.system]
        return Observation(spectrum, band).effstim(unit)

    if brightness.locator_kind is LocatorKind.BAND:
        obs = Observation(spectrum, band)
        if brightness.amount_kind is AmountKind.ENERGY_FLUX:
            return obs.integrate(flux_unit=units.FLAM)
        return obs.effstim(brightness.value.unit)

    lam = brightness.locator.to(u.AA, u.spectral())
    return spectrum(lam, flux_unit=brightness.value.unit)


class TestScaleReproducesAmount:
    """The scaled spectrum's synthetic photometry equals the stated amount."""

    # The four ESO-371803 5.3.6 canonical pairs, plus AB/ST and band energy flux.
    @pytest.mark.parametrize("spec", [
        ("R", 15*u.mag),  # Vega magnitude
        ("R", 10.5*u.ABmag),  # AB magnitude
        ("R", 18*u.STmag),  # ST magnitude
        ("K", 3.5*u.mJy),  # flux density (per nu)
        ("V", 1.2e-15*u.erg/(u.s*u.cm**2*u.AA)),  # flux density (per lambda)
        ("V", 2e-15*u.W/u.m**2),  # band-integrated energy flux
    ])
    def test_band_forms(self, spec, flat_spec, band, vega):
        b = parse_brightness(spec)
        vg = vega if b.system is PhotometricSystem.VEGA else None
        scale = synphot_flux_scale(flat_spec, b, band=band, vegaspec=vg)
        got = _measured_amount(flat_spec * scale, b, band=band, vega=vega)
        if b.amount_kind is AmountKind.MAG:
            # Can't use npt here because of AB/STmag
            assert u.isclose(got, b.value.value * got.unit, rtol=1e-9)
        else:
            npt.assert_allclose(got, b.value, rtol=1e-9)

    @pytest.mark.parametrize("spec", [
        (230*u.GHz, 5*u.mJy),  # frequency locator
        (656.3*u.nm, 1.2e-16*u.erg/(u.s*u.cm**2*u.AA)),  # wavelength locator
    ])
    def test_monochromatic_forms(self, spec, flat_spec):
        b = parse_brightness(spec)
        got = _measured_amount(flat_spec * synphot_flux_scale(flat_spec, b), b)
        npt.assert_allclose(got, b.value, rtol=1e-9)


class TestEquivalences:
    def test_four_form_scale_equivalence(self, subtests, flat_spec, band, vega):
        # T-B1 at the scale level: take one physical flux, express it as a Vega
        # mag, an AB mag, a band flux density and a band energy flux (all
        # measured from the same reference level), and confirm the four scales
        # agree. Anchor on an arbitrary target scale applied to the flat_spec.
        target = flat_spec * 3.0
        obs = Observation(target, band)
        specs_and_vega = [
            (("V", obs.effstim(VEGAMAG, vegaspec=vega)), vega),
            (("V", obs.effstim(u.ABmag)), None),
            (("V", obs.effstim(u.Jy)), None),
            (("V", obs.integrate(flux_unit=units.FLAM).to(u.W/u.m**2)), None),
        ]
        for spec, vg in specs_and_vega:
            with subtests.test(amount=spec[1]):
                b = parse_brightness(spec)
                s = synphot_flux_scale(flat_spec, b, band=band, vegaspec=vg)
                npt.assert_allclose(s, 3.0, rtol=1e-9)

    def test_mag_arcsec2_equals_jy_sr(self, flat_spec, band, vega):
        # T-B6: a magnitude surface brightness reduced to an integrated Vega mag
        # and the equivalent Jy (from the same reference) give the same scale --
        # conversion to linear happens before area scaling, not in mag space.
        b_mag = parse_brightness(("V", 18.05*u.mag))
        s_mag = synphot_flux_scale(flat_spec, b_mag, band=band, vegaspec=vega)
        implied_jy = Observation(flat_spec * s_mag, band).effstim(u.Jy)
        b_jy = parse_brightness(("V", implied_jy))
        s_jy = synphot_flux_scale(flat_spec, b_jy, band=band)
        npt.assert_allclose(s_mag, s_jy, rtol=1e-9)

    def test_mag_closed_form(self, flat_spec, band):
        # T8-style: an AB magnitude scale is exactly 10**(-0.4*(m - m_actual)).
        b = parse_brightness(("R", 13*u.ABmag))
        m_actual = Observation(flat_spec, band).effstim(u.ABmag)
        expected = 10**(-0.4 * (13 - m_actual.to_value(u.ABmag)))
        npt.assert_allclose(
            synphot_flux_scale(flat_spec, b, band=band), expected, rtol=1e-12
        )


class TestGuards:
    def test_monochromatic_energy_flux_raises_E1(self, flat_spec):
        b = parse_brightness((656.3*u.nm, 2e-15*u.W/u.m**2))
        with pytest.raises(BrightnessError) as exc:
            synphot_flux_scale(flat_spec, b)
        assert exc.value.code == "E1"

    def test_magnitude_without_band_raises(self, flat_spec, vega):
        b = parse_brightness(("R", 15*u.mag))
        with pytest.raises(ValueError):
            synphot_flux_scale(flat_spec, b, band=None, vegaspec=vega)

    def test_band_flux_without_band_raises(self, flat_spec):
        b = parse_brightness(("K", 3.5*u.mJy))
        with pytest.raises(ValueError):
            synphot_flux_scale(flat_spec, b, band=None)
