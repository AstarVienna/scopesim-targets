# -*- coding: utf-8 -*-
"""Unit tests for brightness.py.

These exercise the pure parser only: no FilterSystem, no spextra, no network.
The load-time band-vocabulary check and the synphot scaling live on
``SpectrumTarget`` and are tested in ``test_target.py`` instead.
"""

import pytest
import astropy.units as u

from scopesim_targets.brightness import (
    parse_brightness,
    Brightness,
    LocatorKind,
    AmountKind,
    PhotometricSystem,
    BrightnessError,
    solid_angle_unit,
    unit_includes_per_physical_type,
)


ARCSEC2 = u.Unit("arcsec2")


class TestLocatorDispatch:
    def test_band(self):
        b = parse_brightness(["R", "15 mag"])
        assert b.locator_kind is LocatorKind.BAND
        assert b.locator == "R"

    def test_wavelength(self):
        b = parse_brightness(["656.3 nm", "1.2e-16 erg / (s cm2 Angstrom)"])
        assert b.locator_kind is LocatorKind.WAVELENGTH

    def test_frequency(self):
        b = parse_brightness(["230 GHz", "5 mJy"])
        assert b.locator_kind is LocatorKind.FREQUENCY

    def test_quantity_locator_matches_string(self):
        # The YAML resolver may hand us a Quantity for an unquoted locator.
        assert (parse_brightness([656.3 * u.nm, "5 mJy"])
                == parse_brightness(["656.3 nm", "5 mJy"]))


class TestAmountDispatch:
    @pytest.mark.parametrize("amount, system", [
        ("15 mag", PhotometricSystem.VEGA),
        ("10.5 mag(AB)", PhotometricSystem.AB),
        ("10 mag(ST)", PhotometricSystem.ST),
        ("10 mag(Vega)", PhotometricSystem.VEGA),   # astropy can't parse this; we do
    ])
    def test_magnitude_systems(self, amount, system):
        b = parse_brightness(["R", amount])
        assert b.amount_kind is AmountKind.MAG
        assert b.system is system
        assert b.solid_angle is None

    def test_bare_number_is_vega_mag(self):
        b = parse_brightness(["V", 15])
        assert b.amount_kind is AmountKind.MAG
        assert b.system is PhotometricSystem.VEGA
        assert b.value == 15 * u.mag

    @pytest.mark.parametrize("amount, kind", [
        ("3.5 mJy", AmountKind.FLUX_DENSITY_NU),
        ("1.2e-16 erg / (s cm2 Angstrom)", AmountKind.FLUX_DENSITY_LAM),
        ("4.2e-18 W / m2", AmountKind.ENERGY_FLUX),
    ])
    def test_flux_amounts(self, amount, kind):
        locator = "656.3 nm" if "Angstrom" in amount else "K"
        b = parse_brightness([locator, amount])
        assert b.amount_kind is kind
        assert b.solid_angle is None


class TestSurfaceBrightness:
    @pytest.mark.parametrize("amount, kind, system, sa", [
        ("21.5 mag / arcsec2", AmountKind.MAG, PhotometricSystem.VEGA, ARCSEC2),
        ("21.5 mag(AB) / arcsec2", AmountKind.MAG, PhotometricSystem.AB, ARCSEC2),
        ("18 mag(ST) / sr", AmountKind.MAG, PhotometricSystem.ST, u.sr),
        ("5 MJy / sr", AmountKind.FLUX_DENSITY_NU, PhotometricSystem.VEGA, u.sr),
        ("3e-8 W / (m2 arcsec2)", AmountKind.ENERGY_FLUX, PhotometricSystem.VEGA, ARCSEC2),
    ])
    def test_surface_brightness_forms(self, amount, kind, system, sa):
        b = parse_brightness(["V", amount])
        assert b.amount_kind is kind
        assert b.is_surface_brightness
        assert b.solid_angle == sa
        if kind is AmountKind.MAG:
            assert b.system is system

    def test_integrated_has_no_solid_angle(self):
        assert not parse_brightness(["V", "15 mag"]).is_surface_brightness


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
        {"value": "15 mag"},                                 # no locator
        {"band": "V", "frequency": "230 GHz", "value": 1},   # two locators
        {"band": "V"},                                       # no value
    ])
    def test_requires_exactly_one_locator_and_value(self, mapping):
        with pytest.raises(BrightnessError) as exc:
            parse_brightness(mapping)
        assert exc.value.code == "E2"


class TestQuoteIndependence:
    """A YAML author quoting or not must not change the parsed result."""

    @pytest.mark.parametrize("quantity, string", [
        (15 * u.mag, "15 mag"),
        (3.5 * u.mJy, "3.5 mJy"),
    ])
    def test_amount(self, quantity, string):
        assert parse_brightness(["R", quantity]) == parse_brightness(["R", string])

    def test_locator(self):
        assert (parse_brightness([230 * u.GHz, "5 mJy"])
                == parse_brightness(["230 GHz", "5 mJy"]))


class TestErrorMatrix:
    @pytest.mark.parametrize("spec, code", [
        (["V", "5 kg"], "E1"),                                         # bad phys type
        (["5 kg", "15 mag"], "E2"),                                    # bad locator
        (["230 GHz", "15 mag"], "E3"),                                 # mag needs band
        ({"band": "K", "value": "3.5 mJy", "system": "AB"}, "E4"),     # system on non-mag
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


class TestHelpers:
    @pytest.mark.parametrize("unit, expected", [
        ("Jy / sr", u.sr),
        ("Jy / arcsec2", ARCSEC2),
        ("W / (m2 arcsec2)", ARCSEC2),
        ("Jy", None),
        ("W / m2", None),
    ])
    def test_solid_angle_unit(self, unit, expected):
        assert solid_angle_unit(u.Unit(unit)) == expected

    @pytest.mark.parametrize("unit, per", [
        ("Jy / sr", True),
        ("Jy / arcsec2", True),
        ("Jy", False),
    ])
    def test_unit_includes_per_solid_angle(self, unit, per):
        assert unit_includes_per_physical_type(u.Unit(unit), "solid angle") is per

    def test_function_unit_has_no_bases(self):
        # ABmag is a function unit -> the AttributeError path -> False
        assert unit_includes_per_physical_type(u.ABmag, "solid angle") is False


class TestRoundTrip:
    @pytest.fixture(params=[
        ["R", "15 mag"],
        ["R", "10.5 mag(AB)"],
        ["V", "21.5 mag / arcsec2"],
        ["V", "21.5 mag(AB) / arcsec2"],
        ["K", "3.5 mJy"],
        ["230 GHz", "5 mJy"],
        ["656.3 nm", "1.2e-16 erg / (s cm2 Angstrom)"],
        {"band": "V", "value": 21.4, "system": "AB"},
    ])
    def sample_spec(self, request):
        return request.param

    def test_parses_to_brightness(self, sample_spec):
        assert isinstance(parse_brightness(sample_spec), Brightness)

    def test_idempotent(self, sample_spec):
        # Parsing is deterministic and equality is structural.
        assert parse_brightness(sample_spec) == parse_brightness(sample_spec)


class TestCompatShims:
    """TEMPORARY: the .band/.mag shims bridge the legacy to_source path and are
    removed in Phase 2. These guard their interim behaviour until then."""

    def test_band_and_mag(self):
        b = parse_brightness(["R", "12.5 mag"])
        assert b.band == "R"
        assert b.mag == 12.5 * u.mag

    def test_mag_rejects_non_vega(self):
        b = parse_brightness(["R", "10 mag(AB)"])
        with pytest.raises(NotImplementedError):
            _ = b.mag

    def test_band_rejects_non_band_locator(self):
        b = parse_brightness(["230 GHz", "5 mJy"])
        with pytest.raises(AttributeError):
            _ = b.band
