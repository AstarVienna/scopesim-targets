# -*- coding: utf-8 -*-
"""Unit tests for target.py."""

import pytest

from astropy import units as u
from astropy.coordinates import SkyCoord, Angle
from synphot import SourceSpectrum

from astar_utils import SpectralType
from spextra.exceptions import NotInLibraryError

from scopesim_targets.target import Target, SpectrumTarget


@pytest.fixture(scope="function")
def target_subcls():
    """Mock subclass of `Target` for testing.

    This is necessary because Target is an abstract base class which cannot be
    instantiated directly. The hack below allows to run a super basic test on
    it, without using any of its actual subcls.
    """
    class MockTargetSubcls(Target):
        def to_source(self):
            pass
    return MockTargetSubcls()


@pytest.fixture(scope="function")
def spectrum_target_subcls():
    """Like ``target_subcls``, but for `SpectrumTarget`."""
    class MockSpectrumTargetSubcls(SpectrumTarget):
        def to_source(self):
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
        assert target_subcls.resolve_offset() == SkyCoord(0, 0, unit=u.arcsec)

    def test_offset_needs_parent(self, target_subcls):
        target_subcls.offset = {"separation": 5*u.arcsec}
        with pytest.raises(ValueError):
            target_subcls.resolve_offset()

    def test_offset_length_with_distance(self, target_subcls):
        # target_subcls.position = {"distance": 10*u.pc}
        target_subcls.offset = {"separation": .1*u.arcsec}
        # TODO: Try to do without this...
        parent_position = SkyCoord(0*u.deg, 0*u.deg, 10*u.pc)
        offset = target_subcls.resolve_offset(parent_position)
        assert (offset.dec << u.arcsec).round(7) == .1*u.arcsec


class TestSpectrumTarget:
    # @pytest.mark.webtest
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

    # @pytest.mark.webtest
    def test_resolves_spectrum(self, spectrum_target_subcls):
        spectrum_target_subcls.spectrum = "G2V"
        resolved = spectrum_target_subcls.resolve_spectrum(spectrum_target_subcls.spectrum)
        assert isinstance(resolved, SourceSpectrum)

    # @pytest.mark.webtest
    def test_resolves_spectrum_spex(self, spectrum_target_subcls):
        spectrum_target_subcls.spectrum = "spex:kurucz/g2v"
        resolved = spectrum_target_subcls.resolve_spectrum(spectrum_target_subcls.spectrum)
        assert isinstance(resolved, SourceSpectrum)

    # @pytest.mark.webtest
    def test_resolves_spectrum_throws(self, spectrum_target_subcls):
        spectrum_target_subcls.spectrum = "G5V"  # not in current default lib
        with pytest.raises(NotInLibraryError):
            spectrum_target_subcls.resolve_spectrum(spectrum_target_subcls.spectrum)

    # @pytest.mark.webtest
    def test_brightness_number(self, spectrum_target_subcls):
        spectrum_target_subcls.brightness = ["V", 10]
        assert spectrum_target_subcls.brightness.band == "V"
        assert spectrum_target_subcls.brightness.mag == 10*u.mag

    # @pytest.mark.webtest
    def test_brightness_qty(self, spectrum_target_subcls):
        spectrum_target_subcls.brightness = ["R", 12.5*u.mag]
        assert spectrum_target_subcls.brightness.band == "R"
        assert spectrum_target_subcls.brightness.mag == 12.5*u.mag

    # @pytest.mark.webtest
    def test_brightness_throws(self, spectrum_target_subcls):
        with pytest.raises(TypeError):
            spectrum_target_subcls.brightness = "bogus"

    # @pytest.mark.webtest
    def test_brightness_throws_filter(self, spectrum_target_subcls):
        with pytest.raises(ValueError):
            spectrum_target_subcls.brightness = ["bogus", 10]

    def test_spectrum_scaling(self, spectrum_target_subcls):
        spectrum_target_subcls.spectrum = "G2V"
        spectrum_target_subcls.brightness = ("R", 18)
        scale = spectrum_target_subcls._get_spectrum_scale(
            spectrum_target_subcls.resolve_spectrum(spectrum_target_subcls.spectrum),
            spectrum_target_subcls.brightness)
        assert scale == pytest.approx(1.9e-23)  # TODO: CHECK THIS NUMBER!!!
