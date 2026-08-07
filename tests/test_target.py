# -*- coding: utf-8 -*-
"""Unit tests for target.py."""

import pytest

import numpy as np
from numpy import testing as npt
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle
from synphot import SourceSpectrum

from astar_utils import SpectralType
from spextra.exceptions import NotInLibraryError

from scopesim_targets.brightness import (
    parse_brightness,
    BrightnessError,
    AnchorFrame,
    FromSpectralType,
)
from scopesim_targets.point_source import Star
from scopesim_targets.target import Target, SpectrumTarget


@pytest.fixture(scope="function")
def target_subcls():
    """Mock subclass of `Target` for testing.

    This is necessary because Target is an abstract base class which cannot be
    instantiated directly. The hack below allows to run a super basic test on
    it, without using any of its actual subcls.
    """
    class MockTargetSubcls(Target):
        def to_source(self, *args, **kwargs):
            pass
    return MockTargetSubcls()


@pytest.fixture(scope="function")
def spectrum_target_subcls():
    """Like ``target_subcls``, but for `SpectrumTarget`."""
    class MockSpectrumTargetSubcls(SpectrumTarget):
        def to_source(self, *args, **kwargs):
            pass
    return MockSpectrumTargetSubcls()


class TestTarget:
    def test_basic(self, target_subcls):
        assert isinstance(target_subcls, Target)

    @pytest.mark.parametrize("position", (
        (-5, 10),
        {"y": 10/3600 * u.deg, "x": -5/60 * u.arcmin},
        SkyCoord(-5, 10, unit="arcsec"),
    ))
    def test_position(self, position, target_subcls):
        target_subcls.position = position
        assert target_subcls.position == SkyCoord(-5*u.arcsec, 10*u.arcsec)

    def test_position_throws(self, target_subcls):
        with pytest.raises(TypeError):
            target_subcls.position = "bogus"

    def test_position_with_distance(self, target_subcls):
        target_subcls.position = {"x": -2, "y": 10.5, "distance": 20*u.pc}
        expected = SkyCoord(-2*u.arcsec, 10.5*u.arcsec, 20*u.pc)
        assert target_subcls.position == expected

    def test_negative_position_throws(self, target_subcls):
        with pytest.raises(ValueError):
            target_subcls.position = {"distance": -20*u.pc}

    def test_offset(self, target_subcls):
        target_subcls.offset = {
            "separation": 2*u.AU,
            "position_angle": 30*u.deg,
        }
        assert target_subcls.offset["separation"] == 2*u.AU
        assert target_subcls.offset["position_angle"] == Angle(30*u.deg)

    def test_offset_with_default(self, target_subcls):
        target_subcls.offset = {"separation": 5*u.arcsec}
        assert target_subcls.offset["separation"] == 5*u.arcsec
        assert target_subcls.offset["position_angle"] == Angle(0*u.deg)

    def test_offset_throws(self, target_subcls):
        with pytest.raises(TypeError):
            target_subcls.offset = "bogus"

    def test_default_offset(self, target_subcls):
        assert target_subcls.resolve_position() == SkyCoord(0, 0, unit=u.arcsec)

    def test_offset_needs_parent(self, target_subcls):
        target_subcls.offset = {"separation": 5*u.arcsec}
        with pytest.raises(ValueError):
            target_subcls.resolve_position()

    def test_offset_length_with_distance(self, target_subcls):
        # target_subcls.position = {"distance": 10*u.pc}
        target_subcls.offset = {"separation": .1*u.arcsec}
        # TODO: Try to do without this...
        parent_position = SkyCoord(0*u.deg, 0*u.deg, 10*u.pc)
        offset = target_subcls.resolve_position(parent_position)
        assert (offset.dec << u.arcsec).round(7) == .1*u.arcsec


class TestSpectrumTarget:
    @pytest.mark.webtest
    def test_spectrum_synphot(self, spectrum_target_subcls):
        spectrum_target_subcls.spectrum = SourceSpectrum.from_vega()
        assert isinstance(spectrum_target_subcls.spectrum, SourceSpectrum)

    def test_spectrum_str(self, spectrum_target_subcls):
        spectrum_target_subcls.spectrum = "G2V"
        assert isinstance(spectrum_target_subcls.spectrum, SpectralType)
        assert spectrum_target_subcls.spectrum == "g2v"

    @pytest.mark.parametrize("spectrum", ("bogus", 42))
    def test_spectrum_throws(self, spectrum, spectrum_target_subcls):
        with pytest.raises((ValueError, TypeError)):
            spectrum_target_subcls.spectrum = spectrum

    def test_spectrum_from_file(self, spectrum_target_subcls):
        # TODO: Add actual test with mock file
        spectrum_target_subcls.spectrum = "file:bogus"
        with pytest.raises(FileNotFoundError):
            spectrum_target_subcls.resolve_spectrum(spectrum_target_subcls.spectrum)

    @pytest.mark.webtest
    def test_resolves_spectrum(self, spectrum_target_subcls):
        spectrum_target_subcls.spectrum = "G2V"
        resolved = spectrum_target_subcls.resolve_spectrum(spectrum_target_subcls.spectrum)
        assert isinstance(resolved, SourceSpectrum)

    @pytest.mark.webtest
    def test_resolves_spectrum_spex(self, spectrum_target_subcls):
        spectrum_target_subcls.spectrum = "spex:kurucz/g2v"
        resolved = spectrum_target_subcls.resolve_spectrum(spectrum_target_subcls.spectrum)
        assert isinstance(resolved, SourceSpectrum)

    @pytest.mark.webtest
    def test_resolves_spectrum_throws(self, spectrum_target_subcls):
        spectrum_target_subcls.spectrum = "G5V"  # not in current default lib
        with pytest.raises(NotInLibraryError):
            spectrum_target_subcls.resolve_spectrum(spectrum_target_subcls.spectrum)

    @pytest.mark.webtest
    def test_brightness_delegates_to_parser(self, spectrum_target_subcls):
        # setter -> _parse_brightness wrapper -> parse_brightness, valid band survives
        spectrum_target_subcls.brightness = ["V", 10]
        assert spectrum_target_subcls.brightness.locator == "V"
        assert spectrum_target_subcls.brightness.value == 10 * u.mag

    @pytest.mark.webtest
    def test_brightness_throws_filter(self, spectrum_target_subcls):
        with pytest.raises(ValueError):
            spectrum_target_subcls.brightness = ["bogus", 10]

    def test_spectrum_scaling(self, spectrum_target_subcls):
        spectrum_target_subcls.spectrum = "G2V"
        spectrum_target_subcls.brightness = ("R", 18)
        scale = spectrum_target_subcls._get_spectrum_scale(
            spectrum_target_subcls.resolve_spectrum(spectrum_target_subcls.spectrum),
            spectrum_target_subcls.brightness)
        npt.assert_allclose(scale, 6.253e-12, rtol=3e-4)  # TODO: CHECK THIS NUMBER!!!

    def test_surface_brightness_on_point_source_raises_E7(self):
        sb = parse_brightness(["V", "21.5 mag / arcsec2"])
        with pytest.raises(BrightnessError) as exc:
            SpectrumTarget._get_spectrum_scale(None, sb)
        assert exc.value.code == "E7"


class TestAnchor:
    """The ``anchor`` frame attribute (lives on the ``Target`` base)."""

    def test_default_is_observed(self, target_subcls):
        assert target_subcls.anchor is AnchorFrame.OBSERVED

    @pytest.mark.parametrize(
        "value, expected",
        [
            ("intrinsic", AnchorFrame.INTRINSIC),
            ("absolute", AnchorFrame.ABSOLUTE),
            (AnchorFrame.ABSOLUTE, AnchorFrame.ABSOLUTE),
        ],
    )
    def test_set_via_string_or_enum(self, value, expected, target_subcls):
        target_subcls.anchor = value
        assert target_subcls.anchor is expected

    def test_reject_unknown(self, target_subcls):
        with pytest.raises(ValueError, match="anchor must be one of"):
            target_subcls.anchor = "apparent"

    def test_distance_or_none_without_position(self, target_subcls):
        assert target_subcls._distance_or_none() is None

    def test_distance_or_none_without_distance(self, target_subcls):
        target_subcls.position = (1, 2)  # no distance component
        assert target_subcls._distance_or_none() is None

    def test_distance_or_none_with_distance(self, target_subcls):
        target_subcls.position = {"distance": 25 * u.pc}
        npt.assert_allclose(
            target_subcls._distance_or_none().to_value(u.pc), 25.0
        )


class TestAnchorScaling:
    """Anchor semantics in the flux-scaling path (E10/E11 + distance modulus).

    These stay offline: the E-guards raise before any photometry, and the
    distance-modulus factor is separable, so it is checked by stubbing the pure
    photometric scale (:meth:`_get_spectrum_scale`) to 1. The end-to-end
    photometric identity is the webtest ``test_absolute_distance_modulus`` below.
    """

    def _bare(self, spectrum_target_subcls, amount, anchor, position=None):
        # Set the parsed Brightness directly to bypass the FILTER_SYSTEM
        # band-membership check (which is network-backed).
        t = spectrum_target_subcls
        t._brightness = parse_brightness(("V", amount))
        t.anchor = anchor
        if position is not None:
            t.position = position
        return t

    def test_absolute_without_distance_raises_E10(self, spectrum_target_subcls):
        t = self._bare(spectrum_target_subcls, "10 mag(AB)", "absolute")
        with pytest.raises(BrightnessError) as exc:
            t._anchored_spectrum_scale(None, t.brightness)
        assert exc.value.code == "E10"

    def test_absolute_surface_brightness_raises_E11(self, spectrum_target_subcls):
        t = spectrum_target_subcls
        t._brightness = parse_brightness(("V", "21 mag(AB) / arcsec2"))
        t.anchor = "absolute"
        t.position = {"distance": 25 * u.pc}  # E11 must fire despite a distance
        with pytest.raises(BrightnessError) as exc:
            t._anchored_spectrum_scale(None, t.brightness)
        assert exc.value.code == "E11"

    def test_absolute_applies_distance_modulus(
        self, spectrum_target_subcls, monkeypatch
    ):
        # Stub the pure photometric scale so only the distance factor remains.
        monkeypatch.setattr(
            type(spectrum_target_subcls),
            "_get_spectrum_scale",
            staticmethod(lambda spectrum, brightness: 1.0),
        )
        t = self._bare(
            spectrum_target_subcls, "10 mag(AB)", "absolute",
            position={"distance": 25 * u.pc},
        )
        # 5 log10(25/10) = 1.9897 mag -> factor (10/25)**2 = 0.16
        npt.assert_allclose(
            t._anchored_spectrum_scale(None, t.brightness), 0.16, rtol=1e-12
        )

    def test_observed_and_intrinsic_scale_identically_for_now(
        self, spectrum_target_subcls, monkeypatch
    ):
        # Until extinction screens are wired in, observed and intrinsic differ
        # only structurally; the scale is the bare photometric scale for both.
        monkeypatch.setattr(
            type(spectrum_target_subcls),
            "_get_spectrum_scale",
            staticmethod(lambda spectrum, brightness: 3.0),
        )
        t = self._bare(spectrum_target_subcls, "10 mag(AB)", "observed")
        assert t._anchored_spectrum_scale(None, t.brightness) == 3.0
        t.anchor = "intrinsic"
        assert t._anchored_spectrum_scale(None, t.brightness) == 3.0


class TestFromSpectralTypeResolver:
    """The ``{from_spectral_type: ...}`` resolver on a target (E12, provenance).

    The pure guards (non-SpectralType spectrum, unknown table) raise before any
    network access; the successful mamajek lookup is the webtest below.
    """

    def test_non_spectraltype_spectrum_raises_E12(self, spectrum_target_subcls):
        t = spectrum_target_subcls
        t.spectrum = "blackbody:5000 K"  # a str, not a SpectralType
        t.brightness = {"from_spectral_type": "mamajek"}
        with pytest.raises(BrightnessError) as exc:
            _ = t.brightness
        assert exc.value.code == "E12"

    def test_unknown_table_raises(self, spectrum_target_subcls):
        t = spectrum_target_subcls
        t.spectrum = "K5V"
        t.brightness = {"from_spectral_type": "bogus"}
        with pytest.raises(ValueError, match="unknown from_spectral_type table"):
            _ = t.brightness

    @pytest.mark.xfail(reason="anchor not resolvnig, deferred")
    @pytest.mark.webtest
    def test_resolver_matches_manual_absolute(self):
        # T-B10: from_spectral_type is exactly a manual absolute-mag lookup.
        from scopesim_targets.spectral_classes import StellarParameters

        row = StellarParameters().closest_spectral_type(SpectralType("K5V"))
        m_abs = u.Quantity(row["M_V"][0])

        resolved = Star.from_spectral_type("K5V", position={"distance": 25 * u.pc})
        manual = Star(
            spectrum="K5V",
            brightness=("V", m_abs),
            position={"distance": 25 * u.pc},
            anchor="absolute",
        )
        # same resolved value, frame, and provenance
        assert resolved.anchor is AnchorFrame.ABSOLUTE
        npt.assert_allclose(resolved.brightness.value.to_value(u.mag), m_abs)
        assert resolved.brightness_provenance["table"] == "mamajek"
        assert resolved.brightness_provenance["band"] == "V"

        # and the same physical spectrum scale (resolver == manual)
        spec = resolved.resolve_spectrum(resolved.spectrum)
        npt.assert_allclose(
            resolved._anchored_spectrum_scale(spec, resolved.brightness),
            manual._anchored_spectrum_scale(spec, manual.brightness),
            rtol=1e-10,
        )


class TestAbsoluteDistanceModulusEndToEnd:
    @pytest.mark.webtest
    def test_absolute_distance_modulus(self, spectrum_target_subcls):
        # T-B9: absolute M at distance d == apparent (M + 5 log10(d/10)) observed.
        t = spectrum_target_subcls
        t.spectrum = "G2V"
        spec = t.resolve_spectrum(t.spectrum)
        t.position = {"distance": 25 * u.pc}
        d_mod = 5 * np.log10(25 / 10)  # 1.9897 mag (plan rounds apparent to 6.82)

        t.brightness = ("V", "4.83 mag")
        t.anchor = "absolute"
        s_abs = t._anchored_spectrum_scale(spec, t.brightness)

        t.brightness = ("V", f"{4.83 + d_mod} mag")
        t.anchor = "observed"
        s_app = t._anchored_spectrum_scale(spec, t.brightness)

        npt.assert_allclose(s_abs, s_app, rtol=1e-10)
