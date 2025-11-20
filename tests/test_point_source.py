# -*- coding: utf-8 -*-
"""Unit tests for point_source.py."""

from unittest.mock import patch

import pytest
import yaml
import numpy as np
from astropy import units as u

from scopesim_targets.target import Brightness
from scopesim_targets.point_source import (
    PointSourceTarget,
    Star,
    Binary,
    Exoplanet,
)


class TestStar:
    def test_basic(self):
        tgt = Star()
        assert isinstance(tgt, PointSourceTarget)

    # Webtest??
    @pytest.mark.parametrize("position", ((0, 0), (10, -4.2)))
    def test_to_source(self, position):
        # Note: Without any additional info, single point source will be placed
        #       in the center of the FOV for ScopeSim.
        src = Star(position=position, spectrum="A0V",
                   brightness=("V", 10)).to_source()
        assert src.fields[0].field["x"] == 0
        assert src.fields[0].field["y"] == 0
        assert src.fields[0].field["ref"][0] in src.fields[0].spectra

    # Webtest??
    def test_loads_yaml(self):
        tgt = yaml.full_load("""
            !Star
            position: [2, 3]
            spectrum: A0V
            brightness: ["R", 15 mag]
        """)
        assert isinstance(tgt, Star)


@pytest.fixture
def basic_binary():
    return Binary(brightness=("R", 10))


class TestBinary:
    def test_two_brightnesses(self):
        tgt = Binary(brightness=(("R", 10), ("V", 15*u.mag)))
        assert tgt.brightness == Brightness("R", 10*u.mag)
        assert tgt.brightness_secondary == Brightness("V", 15*u.mag)
        with pytest.raises(AttributeError):
            tgt.contrast

    def test_brightness_and_contrast(self):
        tgt = Binary(brightness=("R", 10), contrast=100.)
        assert tgt.brightness == Brightness("R", 10*u.mag)
        assert tgt.contrast == 100
        with pytest.raises(AttributeError):
            tgt.brightness_secondary

    def test_invalid_contrast_throws(self):
        with pytest.raises(TypeError):
            Binary(brightness=("R", 10), contrast="bogus")

    def test_two_brightnesses_and_contrast_throws(self):
        with pytest.raises(TypeError):
            Binary(brightness=(("R", 10), ("V", 15*u.mag)), contrast=100.)

    def test_invalid_brightness_and_contrast_throws(self):
        with pytest.raises(TypeError):
            Binary(brightness="bogus", contrast="bogus")

    @pytest.mark.parametrize(
        ("brightness", "contrast"), (
            (None, None),
            (("R", 10), None),
            (None, 10.),
        ))
    def test_other_cases(self, brightness, contrast):
        # TODO: Replace this with more meaningful tests!
        tgt = Binary(brightness=brightness, contrast=contrast)
        assert isinstance(tgt, Binary)

    def test_to_source(self):
        # TODO: cover more possible cases
        tgt = Binary(
            brightness=("R", 10),
            contrast=100.,
            spectra=["A0V", "M2V"],
            offset={"separation": 5*u.arcsec},
        )
        src = tgt.to_source()
        np.testing.assert_array_equal(src.fields[0].field["x"], [0, 0])
        np.testing.assert_array_equal(src.fields[0].field["y"], [0, 5])

    def test_throws_if_no_contras_or_secondary_brightness(self):
        tgt = Binary(
            spectra=("A0V", "M2V"),
            brightness=("R", 10),
            offset={"separation": 2*u.arcsec},
        )
        with pytest.raises(ValueError):
            tgt.to_source()

    @pytest.mark.parametrize(("spectra", "refs", "called"), (
        (None, None, True),
        (None, 42, True),
        ({0: "foo", 1: "bar"}, None, False),
        ({0: "foo", 1: "bar"}, [0, 1], False),
    ))
    def test_resolve_spectra_refs(self, basic_binary, spectra, refs, called):
        start = refs if isinstance(refs, int) else 0
        expected_spec = {start: "foo", start+1: "bar"}
        expected_refs = (start, start+1)

        with patch.object(basic_binary, "source_spectra") as mock_spectra:
            mock_spectra.return_value = expected_spec

            result = basic_binary._resolve_spectra_refs(spectra, refs)

            if called:
                if refs is not None:
                    mock_spectra.assert_called_once_with(refs)
                else:
                    mock_spectra.assert_called_once()
            else:
                mock_spectra.assert_not_called()

        assert result == (expected_spec, expected_refs)

    @pytest.mark.parametrize(("spectra", "refs", "exc", "msg"), (
        (None, [0, 1], ValueError, "refs sequence must have matching spectra"),
        ({0: "foo"}, [0, 1], ValueError, "not all refs found in spectra"),
        ("bogus", 42, TypeError, "refs and spectra not understood"),
    ))
    def test_resolve_spectra_refs_throws(
            self, basic_binary, spectra, refs, exc, msg):
        with pytest.raises(exc, match=msg):
            basic_binary._resolve_spectra_refs(spectra, refs)


class TestExoplanet:
    def test_basic(self):
        tgt = Exoplanet()
        assert isinstance(tgt, Exoplanet)

    # webtest
    def test_all_attributes(self):
        # TODO: Replace this with a more meaningful test
        Exoplanet(position=(0, 1), offset={"separation": 2*u.arcsec},
                  brightness=("K", 23), spectrum="spex:irtf/Jupiter",
                  contrast=1e3)

    # webtest
    def test_default_spectrum(self):
        tgt = Exoplanet(offset={"separation": 2*u.arcsec}, contrast=1e3)
        assert str(tgt.spectrum) == "Spextrum(irtf/Neptune)"
