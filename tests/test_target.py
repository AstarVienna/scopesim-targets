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

    # @pytest.mark.webtest
    def test_resolves_spectrum(self, spectrum_target_subcls):
        spectrum_target_subcls.spectrum = "G2V"
        resolved = spectrum_target_subcls.resolve_spectrum()
        assert isinstance(resolved, SourceSpectrum)

    # @pytest.mark.webtest
    def test_resolves_spectrum_spex(self, spectrum_target_subcls):
        spectrum_target_subcls.spectrum = "spex:kurucz/g2v"
        resolved = spectrum_target_subcls.resolve_spectrum()
        assert isinstance(resolved, SourceSpectrum)

    # @pytest.mark.webtest
    def test_resolves_spectrum_throws(self, spectrum_target_subcls):
        spectrum_target_subcls.spectrum = "G5V"  # not in current default lib
        with pytest.raises(NotInLibraryError):
            spectrum_target_subcls.resolve_spectrum()
