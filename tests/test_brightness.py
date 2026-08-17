# -*- coding: utf-8 -*-
"""Unit tests for brightness.py.

These exercise the pure parser only: no FilterSystem, no spextra, no network.
The load-time band-vocabulary check and the synphot scaling live on
``SpectrumTarget`` and are tested in ``test_target.py`` instead.
"""

import pytest
import astropy.units as u
from synphot.units import VEGAMAG

from scopesim_targets.brightness import (
    parse_brightness,
    Brightness,
    LocatorKind,
    AmountKind,
    PhotometricSystem,
    AnchorFrame,
    FromSpectralType,
    BrightnessError,
    solid_angle_unit,
)


class TestLocatorDispatch:
    def test_band(self):
        b = parse_brightness(("R", "15 mag"))
        assert b.locator_kind is LocatorKind.BAND
        assert b.locator == "R"

    def test_wavelength(self):
        b = parse_brightness(("656.3 nm", "1.2e-16 erg / (s cm2 Angstrom)"))
        assert b.locator_kind is LocatorKind.WAVELENGTH
        assert b.locator == 656.3*u.nm

    def test_frequency(self):
        b = parse_brightness(("230 GHz", "5 mJy"))
        assert b.locator_kind is LocatorKind.FREQUENCY
        assert b.locator == 230*u.GHz

    def test_quantity_locator_matches_string(self):
        qty = parse_brightness((656.3*u.nm, "5 mJy"))
        assert qty == parse_brightness(("656.3 nm", "5 mJy"))


class TestAmountDispatch:
    @pytest.mark.parametrize("amount, system", [
        ("15 mag", PhotometricSystem.VEGA),
        ("10.5 mag(AB)", PhotometricSystem.AB),
        ("10 mag(ST)", PhotometricSystem.ST),
        ("10 mag(Vega)", PhotometricSystem.VEGA),
    ])
    def test_magnitude_systems(self, amount, system):
        b = parse_brightness(("R", amount))
        assert b.amount_kind is AmountKind.MAG
        assert b.system is system
        assert b.solid_angle is None

    def test_bare_number_is_vega_mag(self):
        b = parse_brightness(("V", 15))
        assert b.amount_kind is AmountKind.MAG
        assert b.system is PhotometricSystem.VEGA
        assert b.value == 15*u.mag

    @pytest.mark.parametrize("locator, amount, kind", [
        ("K", "3.5 mJy", AmountKind.FLUX_DENSITY_NU),
        ("656.3 nm", "1.2e-15 erg / (s cm2 nm)", AmountKind.FLUX_DENSITY_LAM),
        ("K", "4.2e-18 W / m2", AmountKind.ENERGY_FLUX),
    ])
    def test_flux_amounts(self, locator, amount, kind):
        b = parse_brightness((locator, amount))
        assert b.amount_kind is kind
        assert b.solid_angle is None


class TestSurfaceBrightness:
    @pytest.mark.parametrize("amount, system, sa", [
        ("21.5 mag / arcsec2", PhotometricSystem.VEGA, u.arcsec**2),
        ("21.5 mag(Vega) / arcsec2", PhotometricSystem.VEGA, u.arcsec**2),
        (21.5*u.mag/u.arcsec**2, PhotometricSystem.VEGA, u.arcsec**2),
        ("21.5 mag(AB) / arcsec2", PhotometricSystem.AB, u.arcsec**2),
        ("18 mag(ST) / sr", PhotometricSystem.ST, u.sr),
    ])
    def test_surface_brightness_mag(self, amount, system, sa):
        b = parse_brightness(("V", amount))
        assert b.amount_kind is AmountKind.MAG
        assert b.is_surface_brightness
        assert b.solid_angle == sa
        assert b.system is system

    @pytest.mark.parametrize("amount, kind, sa", [
        ("5 MJy / sr", AmountKind.FLUX_DENSITY_NU, u.sr),
        ("3e-8 W / (m2 arcsec2)", AmountKind.ENERGY_FLUX, u.arcsec**2),
        (3e-8*u.W/(u.m**2*u.arcsec**2), AmountKind.ENERGY_FLUX, u.arcsec**2),
    ])
    def test_surface_brightness_physical(self, amount, kind, sa):
        b = parse_brightness(("V", amount))
        assert b.amount_kind is kind
        assert b.is_surface_brightness
        assert b.solid_angle == sa

    def test_integrated_has_no_solid_angle(self):
        assert not parse_brightness(("V", 15*u.mag)).is_surface_brightness


class TestCanonicalMapping:
    def test_system_field(self):
        b = parse_brightness({"band": "V", "value": 21.4, "system": "AB"})
        assert b.amount_kind is AmountKind.MAG
        assert b.system is PhotometricSystem.AB

    def test_frequency_mapping(self):
        b = parse_brightness({"frequency": "230 GHz", "value": "5 mJy"})
        assert b.locator_kind is LocatorKind.FREQUENCY
        assert b.amount_kind is AmountKind.FLUX_DENSITY_NU

    @pytest.mark.parametrize("mapping", [
        {"value": "15 mag"},  # no locator
        {"band": "V", "frequency": "230 GHz", "value": 1},  # two locators
        {"band": "V"},  # no value
    ])
    def test_requires_exactly_one_locator_and_value(self, mapping):
        with pytest.raises(BrightnessError) as exc:
            parse_brightness(mapping)
        assert exc.value.code == "E2"


@pytest.mark.parametrize("quantity, string", [
    (15*u.mag, "15 mag"),
    (3.5*u.mJy, "3.5 mJy"),
    (18*u.ABmag, "18 mag(AB)"),
    (20*u.STmag, "20 mag(ST)"),
    (23*VEGAMAG, "23 mag(Vega)"),
    (23*VEGAMAG, "23 mag(VEGA)"),
    (23*VEGAMAG, "23 mag"),
])
def test_quote_independence(quantity, string):
    """A YAML author quoting or not must not change the parsed result."""
    assert parse_brightness(("R", quantity)) == parse_brightness(("R", string))


@pytest.mark.parametrize("unit, expected", [
    ("Jy / sr", u.sr),
    ("Jy / arcsec2", u.arcsec**2),
    ("W / (m2 arcsec2)", u.arcsec**2),
    ("Jy", None),
    ("W / m2", None),
])
def test_solid_angle_unit(unit, expected):
    assert solid_angle_unit(u.Unit(unit)) == expected


class TestErrorMatrix:
    @pytest.mark.parametrize("spec, code", [
        (("V", "5 kg"), "E1"),  # bad phys type
        (("5 kg", "15 mag"), "E2"),  # bad locator
        (("230 GHz", "15 mag"), "E3"),  # mag needs band
        ({"band": "K", "value": "3.5 mJy", "system": "AB"}, "E4"),  # system on non-mag
        ({"band": "R", "value": "10 mag(ST)", "system": "AB"}, "E5"),  # double system
    ])
    def test_error_codes(self, spec, code):
        with pytest.raises(BrightnessError) as exc:
            parse_brightness(spec)
        assert exc.value.code == code

    @pytest.mark.parametrize("bad", ["bogus", 42, None])
    def test_wrong_shape_is_typeerror(self, bad):
        # Wrong *type* (not a tuple/mapping) is a TypeError, distinct from the
        # content matrix (E1-E5), matching the pre-existing contract.
        with pytest.raises(TypeError):
            parse_brightness(bad)


class TestRoundTrip:
    @pytest.fixture(params=[
        ("R", "15 mag"),
        ("R", "10.5 mag(AB)"),
        ("V", "21.5 mag / arcsec2"),
        ("V", "21.5 mag(AB) / arcsec2"),
        ("K", "3.5 mJy"),
        ("230 GHz", "5 mJy"),
        ("656.3 nm", "1.2e-16 erg / (s cm2 Angstrom)"),
        {"band": "V", "value": 21.4, "system": "AB"},
    ])
    def sample_spec(self, request):
        return request.param

    def test_parses_to_brightness(self, sample_spec):
        assert isinstance(parse_brightness(sample_spec), Brightness)

    def test_idempotent(self, sample_spec):
        # Parsing is deterministic and equality is structural.
        assert parse_brightness(sample_spec) == parse_brightness(sample_spec)


class TestAnchorFrame:
    """Pure tests for the ``anchor`` frame enum (see defining_brightness.md)."""

    def test_default_is_observed(self):
        assert AnchorFrame.coerce(None) is AnchorFrame.OBSERVED

    @pytest.mark.parametrize("value, expected", [
        ("observed", AnchorFrame.OBSERVED),
        ("intrinsic", AnchorFrame.INTRINSIC),
        ("absolute", AnchorFrame.ABSOLUTE),
        (AnchorFrame.ABSOLUTE, AnchorFrame.ABSOLUTE),
    ])
    def test_coerce_accepts_strings_and_enum(self, value, expected):
        assert AnchorFrame.coerce(value) is expected

    def test_coerce_rejects_unknown(self):
        with pytest.raises(ValueError, match="anchor must be one of"):
            AnchorFrame.coerce("apparent")


class TestFromSpectralTypeMarker:
    """The ``{from_spectral_type: ...}`` form parses to a deferred marker.

    Resolution needs the target's spectrum + distance, so the pure parser only
    recognizes the shape and hands back a :class:`FromSpectralType`; the actual
    (network-backed) lookup is exercised in ``test_target.py`` as a webtest.
    """

    def test_returns_marker_with_default_band(self):
        marker = parse_brightness({"from_spectral_type": "mamajek"})
        assert isinstance(marker, FromSpectralType)
        assert marker.table == "mamajek"
        assert marker.band == "V"

    def test_optional_band(self):
        marker = parse_brightness({"from_spectral_type": "mamajek", "band": "Ks"})
        assert marker.band == "Ks"

    def test_table_must_be_string(self):
        with pytest.raises(BrightnessError) as exc:
            parse_brightness({"from_spectral_type": 5})
        assert exc.value.code == "E2"

    def test_unexpected_keys_rejected(self):
        with pytest.raises(BrightnessError) as exc:
            parse_brightness({"from_spectral_type": "mamajek", "value": 10})
        assert exc.value.code == "E2"
