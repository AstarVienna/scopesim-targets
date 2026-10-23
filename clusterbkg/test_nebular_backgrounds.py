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
from . import nebular_textures as nt

KING = {"r_core": 0.13 * u.pc, "r_tide": 2.0 * u.pc}  # reference-scene cluster
GRID = {"pixel_scale": 1 * u.arcsec / u.pixel, "width": 64, "height": 64}


def _pos(distance):
    return SkyCoord(84.7 * u.deg, -6.9 * u.deg, distance=distance)


def _arcsec(length, distance):
    return (length / distance * u.rad).to_value(u.arcsec)


@pytest.fixture(autouse=True)
def no_sed(monkeypatch):
    """Keep the pure tier offline: the SED is irrelevant to the geometry."""
    monkeypatch.setattr(nb, "load_sed", lambda filename, components=None:
                        filename if components is None else (filename, tuple(components)))


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

    def test_texture_passes_through(self, factory):
        tex = nt.orion_texture(0.3 * u.pc, 0, sigma=0.5, n_grid=256)
        smooth = factory(_pos(414 * u.pc), **KING)
        textured = factory(_pos(414 * u.pc), **KING, texture=tex)
        assert textured.texture is tex and smooth.texture is None
        w_s, w_t = smooth._weightmap(GRID)[0], textured._weightmap(GRID)[0]
        assert w_t.shape == w_s.shape and not np.allclose(w_t, w_s)

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


class TestOrionLayers:
    def test_two_layers_same_geometry(self):
        gas, pah = nb.orion_layers(_pos(414 * u.pc), **KING)
        npt.assert_array_equal(gas._weightmap(GRID)[0], pah._weightmap(GRID)[0])
        assert gas.spectrum == (nb.ORION_SED, nb.GAS_COMPONENTS)
        assert pah.spectrum == (nb.ORION_SED, nb.PAH_COMPONENTS)
        assert gas.anchor == pah.anchor == nb.ORION_ANCHOR

    def test_layers_partition_the_full_flux(self):
        from astropy.table import QTable
        tab = QTable.read(nb.SED_DIR / nb.ORION_SED, hdu="SPECTRUM")
        comps = set(tab.meta["FLUXDEF"].split("=", 1)[1].replace("+", " ").split())
        assert comps == set(nb.GAS_COMPONENTS) | set(nb.PAH_COMPONENTS)
        assert not set(nb.GAS_COMPONENTS) & set(nb.PAH_COMPONENTS)

    def test_pah_anchor_reproduces_sed_level(self):
        from astropy.table import QTable
        tab = QTable.read(nb.SED_DIR / nb.ORION_SED, hdu="SPECTRUM")
        locator, amount = nb.sed_level(nb.ORION_SED, nb.PAH_COMPONENTS, nb.ORION_PAH_WAVELENGTH)
        assert locator == "3291.0 nm"
        expected = np.interp(32910.0, tab["WAVELENGTH"].value, tab["PAH"].value) * tab["PAH"].unit
        npt.assert_allclose(u.Quantity(amount).to_value(expected.unit), expected.value, rtol=1e-6)

    def test_pah_negligible_in_k(self):
        # Justifies keeping ORION_BRIGHTNESS for the gas-only layer.
        import importlib.util
        from pathlib import Path
        from astropy.io.votable import parse_single_table
        from astropy.table import QTable
        tab = QTable.read(nb.SED_DIR / nb.ORION_SED, hdu="SPECTRUM")
        lam = tab["WAVELENGTH"].value
        svo = Path(importlib.util.find_spec("spextra").origin).parent / "data" / "svo"
        rsr = parse_single_table(str(svo / "2MASS/2MASS.Ks")).to_table()
        t = np.interp(lam, np.asarray(rsr["Wavelength"], float), np.asarray(rsr["Transmission"], float), left=0, right=0)
        share = np.trapezoid(tab["PAH"].value * lam * t, lam) / np.trapezoid(tab["FLUX"].value * lam * t, lam)
        assert -2.5 * np.log10(1 - share) < 0.005  # below the rounding of ORION_BRIGHTNESS

    def test_layer_textures_pass_through(self):
        gas_tex = nt.orion_texture(0.3 * u.pc, 0, sigma=0.5, n_grid=256)
        pah_tex = nt.cloud_texture(0.4, 0, n_grid=512)
        gas, pah = nb.orion_layers(_pos(414 * u.pc), **KING, gas_texture=gas_tex, pah_texture=pah_tex)
        assert gas.texture is gas_tex and pah.texture is pah_tex


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
