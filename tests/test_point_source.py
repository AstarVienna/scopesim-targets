# -*- coding: utf-8 -*-
"""Unit tests for point_source.py."""

import pytest

from scopesim_targets.point_source import PointSourceTarget, Star


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
