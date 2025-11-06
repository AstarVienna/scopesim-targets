# -*- coding: utf-8 -*-
"""Unit tests for custom YAML representations."""

import yaml
from astropy import units as u
from astropy.coordinates import SkyCoord


def test_qty():
    qtys = yaml.full_load("""
    distance: 10 kpc
    flux: 3.14e5 W m-2
    """)
    assert isinstance(qtys["distance"], u.Quantity)
    assert qtys["distance"].unit == u.kpc
    assert qtys["distance"].value == 10

    assert isinstance(qtys["flux"], u.Quantity)
    assert qtys["flux"].unit == (u.W / u.m**2)
    assert qtys["flux"].value == 314000.


def test_coord():
    coord = yaml.full_load("""
    !Coord
    ra: 5 deg
    dec: -2 deg
    """)
    assert isinstance(coord, SkyCoord)
    assert coord.ra == 5 * u.deg
    assert coord.dec == -2 * u.deg
