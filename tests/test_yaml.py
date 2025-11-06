# -*- coding: utf-8 -*-
"""Unit tests for custom YAML representations."""

import pytest

import yaml
from astropy import units as u
from astropy.coordinates import SkyCoord


@pytest.fixture(scope="class")
def basic_qtys():
    qtys = {
        "distance": 10 * u.kpc,
        "flux": 314000. * u.W / u.m**2
    }
    return qtys


@pytest.fixture(scope="class")
def basic_skycoord():
    return SkyCoord(5*u.deg, -2*u.deg)


class TestQty:
    def test_loading(self, basic_qtys):
        qtys = yaml.full_load("""
            distance: 10 kpc
            flux: 3.14e5 W m-2
        """)
        assert "distance" in qtys
        assert isinstance(qtys["distance"], u.Quantity)
        assert qtys["distance"] == basic_qtys["distance"]

        assert "flux" in qtys
        assert isinstance(qtys["flux"], u.Quantity)
        assert qtys["flux"] == basic_qtys["flux"]

    def test_dumping(self, basic_qtys):
        expected = "distance: 10.0 kpc\nflux: 314000.0 W / m2\n"
        assert yaml.dump(basic_qtys) == expected

    def test_roundtrip(self, basic_qtys):
        assert yaml.full_load(yaml.dump(basic_qtys)) == basic_qtys


class TestCoord:
    def test_loading(self, basic_skycoord):
        coord = yaml.full_load("""
            !Coord
            ra: 5 deg
            dec: -2 deg
        """)
        assert isinstance(coord, SkyCoord)
        assert coord == basic_skycoord

    def test_dumping(self, basic_skycoord):
        expected = "!Coord\ndec: '-2d00m00s'\nra: '5d00m00s'\n"
        assert yaml.dump(basic_skycoord) == expected

    def test_roundtrip(self, basic_skycoord):
        assert yaml.full_load(yaml.dump(basic_skycoord)) == basic_skycoord
