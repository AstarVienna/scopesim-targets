# -*- coding: utf-8 -*-
"""Tests for nebular_backgrounds.py (Orion- and Ophiuchus-based backgrounds).

Pure tier: ``load_sed`` is replaced by a sentinel, so no speXtra/synphot data
is needed. The webtest loads the real SED files through speXtra.
"""
import numpy as np
import pytest
from numpy import testing as npt
import astropy.units as u
from astropy.coordinates import SkyCoord

from . import nebular_backgrounds as nb

KING = {"r_core": 0.13 * u.pc, "r_tide": 2.0 * u.pc}  # reference-scene cluster
GRID = {"pixel_scale": 1 * u.arcsec / u.pixel, "width": 64, "height": 64}


def _pos(distance):
    return SkyCoord(84.7 * u.deg, -6.9 * u.deg, distance=distance)


def _arcsec(length, distance):
    return (length / distance * u.rad).to_value(u.arcsec)


@pytest.fixture(autouse=True)
def no_sed(monkeypatch):
    """Keep the pure tier offline: the SED is irrelevant to the geometry."""
    monkeypatch.setattr(nb, "load_sed", lambda filename, components=None: filename)


@pytest.mark.parametrize("factory", [nb.orion_background, nb.ophiuchus_background])
class TestGeometry:
    def test_radii_follow_position_distance(self, factory):
        d = 414 * u.pc
        bkg = factory(_pos(d), **KING)
        m = bkg._model
        npt.assert_allclose(m.r_core.value, _arcsec(KING["r_core"], d), rtol=1e-6)
        npt.assert_allclose(m.r_out.value, _arcsec(KING["r_tide"], d), rtol=1e-6)
        npt.assert_allclose(m.r_ref.value, _arcsec(nb.ORION_R_REF, d), rtol=1e-6)
        assert m.index.value == nb.DEFAULT_INDEX

    def test_r_ref_is_the_peeters_aperture_at_414_pc(self, factory):
        # 111.4 arcsec at 414 pc, recovered as an angle at that distance.
        bkg = factory(_pos(414 * u.pc), **KING)
        npt.assert_allclose(bkg._model.r_ref.value, 111.4, rtol=1e-9)

    def test_pushing_away_shrinks_angles_keeps_brightness(self, factory):
        near = factory(_pos(414 * u.pc), **KING)
        far = factory(_pos(4140 * u.pc), **KING)
        for name in ("r_core", "r_ref", "r_out"):
            npt.assert_allclose(
                getattr(near._model, name).value / getattr(far._model, name).value,
                10.0, rtol=1e-9,
            )
        assert near.brightness == far.brightness  # SB is distance-invariant

    def test_total_finite_thanks_to_r_tide(self, factory):
        # Index 2 alone diverges; truncation at r_tide makes it finite and
        # FOV-independent, scaling as d^-2 when the scene is moved.
        near = factory(_pos(414 * u.pc), **KING)
        far = factory(_pos(4140 * u.pc), **KING)
        assert near.has_finite_total and far.has_finite_total
        ratio = near.total_flux_factor() / far.total_flux_factor()
        npt.assert_allclose(ratio.decompose().value, 100.0, rtol=1e-9)

    @pytest.mark.parametrize("bad", ["r_core", "r_tide"])
    def test_angular_radius_rejected(self, factory, bad):
        kwargs = dict(KING) | {bad: 10 * u.arcsec}
        with pytest.raises(ValueError, match="physical length"):
            factory(_pos(414 * u.pc), **kwargs)

    def test_tidal_radius_inside_r_ref_rejected(self, factory):
        # The SB anchor (0.224 pc) must lie inside the truncation radius.
        with pytest.raises(ValueError, match="r_ref"):
            factory(_pos(414 * u.pc), r_core=0.13 * u.pc, r_tide=0.2 * u.pc)


class TestOrionVsOphiuchus:
    def test_same_geometry_different_level(self):
        pos = _pos(414 * u.pc)
        ori = nb.orion_background(pos, **KING)
        oph = nb.ophiuchus_background(pos, **KING)
        npt.assert_array_equal(ori._weightmap(GRID)[0], oph._weightmap(GRID)[0])
        assert ori.brightness != oph.brightness
        assert (ori.anchor, oph.anchor) == ("intrinsic", "observed")
        assert (ori.spectrum, oph.spectrum) == (nb.ORION_SED, nb.OPHIUCHUS_SED)


@pytest.mark.webtest
@pytest.mark.parametrize("filename", [nb.ORION_SED, nb.OPHIUCHUS_SED])
def test_real_sed_loads_as_spextrum(filename, monkeypatch):
    monkeypatch.undo()  # use the real load_sed
    from spextra import Spextrum

    sed = nb.load_sed(filename)
    assert isinstance(sed, Spextrum)
    wave = sed.waveset
    assert wave.min() <= 0.97 * u.um and wave.max() >= 5.39 * u.um
    assert np.all(sed(wave).value >= 0)
