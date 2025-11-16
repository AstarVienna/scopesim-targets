# -*- coding: utf-8 -*-
"""Unit tests for point_source.py."""

import pytest
import yaml
from astropy import units as u

from scopesim_targets.target import Brightness
from scopesim_targets.point_source import PointSourceTarget, Star, Binary


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
