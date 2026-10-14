# -*- coding: utf-8 -*-
"""Unit tests for texture.py."""

import numpy as np
import pytest
from numpy import testing as npt
import astropy.units as u

from scopesim_targets.texture import LognormalTexture


def _tex(**kwargs):
    base = dict(sigma=0.5, beta=3.5, outer_scale=0.05 * u.pc,
                inner_scale=0.02 * u.pc, beta_inner=5.0, rng_seed=7, n_grid=256)
    return LognormalTexture(**(base | kwargs))


def _slit_slope(tex, k_lo, k_hi):
    """Slope of the mean 1D power spectrum along rows of ln f (a 'slit')."""
    lf = tex.log_field
    dx = tex.spacing.to_value(tex.unit)
    ps = np.mean(np.abs(np.fft.rfft(lf, axis=1)) ** 2, axis=0)
    k = np.fft.rfftfreq(lf.shape[1], dx)
    m = (k >= k_lo) & (k <= k_hi)
    return np.polyfit(np.log(k[m]), np.log(ps[m]), 1)[0]


class TestField:
    def test_sigma_is_exact_and_mean_is_one(self):
        tex = _tex()
        npt.assert_allclose(tex.log_field.std(), 0.5, rtol=1e-12)
        npt.assert_allclose(tex.field.mean(), 1.0, rtol=1e-12)
        assert np.all(tex.field > 0)

    def test_reproducible_and_seed_dependent(self):
        npt.assert_array_equal(_tex().field, _tex().field)
        assert not np.allclose(_tex().field, _tex(rng_seed=8).field)

    def test_zero_sigma_is_uniform(self):
        npt.assert_allclose(_tex(sigma=0).field, 1.0, rtol=1e-12)

    def test_isotropic_1d_slope_is_beta_minus_one(self):
        # A 1D cut through an isotropic 2D field with P ~ k^-beta has slope
        # -(beta - 1) (the conversion used for slit measurements).
        tex = LognormalTexture(sigma=1, beta=3.79, outer_scale=400 * u.arcsec,
                               rng_seed=0, n_grid=1024)
        npt.assert_allclose(_slit_slope(tex, 1 / 200, 1 / 80), -2.79, atol=0.15)

    def test_spectrum_continuous_at_breaks(self):
        tex = _tex()
        for scale in (tex.outer_scale, tex.inner_scale):
            k = 1 / scale.to_value(u.pc)
            lo, hi = tex.power_spectrum([k * (1 - 1e-9), k * (1 + 1e-9)])
            npt.assert_allclose(lo, hi, rtol=1e-6)

    def test_box_and_spacing(self):
        tex = _tex()
        assert tex.spacing == 0.005 * u.pc  # inner_scale / 4
        assert tex.box_size == 256 * 0.005 * u.pc


class TestSampling:
    def test_nodes_are_exact(self):
        tex = _tex()
        dx = tex.spacing
        out = tex.sample(np.arange(5) * dx, np.arange(3) * dx)
        npt.assert_allclose(out, tex.field[:3, :5], rtol=1e-12)

    def test_periodic(self):
        tex = _tex()
        x = np.linspace(0, 0.1, 7) * u.pc
        npt.assert_allclose(tex.sample(x + tex.box_size, x), tex.sample(x, x), rtol=1e-10)

    def test_pixel_average_smooths(self):
        tex = _tex()
        x = np.linspace(0, 1, 50) * u.pc
        fine = np.log(tex.sample(x, x)).std()
        coarse = np.log(tex.sample(x, x, pixel_size=0.08 * u.pc)).std()
        assert coarse < 0.7 * fine

    def test_angular_texture_needs_no_length(self):
        tex = _tex(outer_scale=22 * u.arcsec, inner_scale=8 * u.arcsec)
        assert not tex.is_physical
        assert tex.sample([0, 1] * u.arcmin, [0] * u.arcsec).shape == (1, 2)

    def test_incompatible_sampling_units_raise(self):
        with pytest.raises(u.UnitsError):
            _tex().sample([0, 1] * u.arcsec, [0] * u.arcsec)


class TestValidation:
    @pytest.mark.parametrize(
        "bad",
        [
            {"sigma": -0.1},
            {"beta": 0},
            {"beta_inner": None},
            {"inner_scale": 0.1 * u.pc},
            {"inner_scale": 8 * u.arcsec},
            {"rng_seed": 1.5},
            {"n_grid": 255},
            {"inner_scale": 0.005 * u.pc, "n_grid": 64},  # box < 4 x outer_scale
        ],
        ids=["sigma<0", "beta=0", "inner without beta_inner", "inner>outer",
             "mixed units", "float seed", "odd n_grid", "box too small"],
    )
    def test_bad_parameters_raise(self, bad):
        with pytest.raises(ValueError):
            _tex(**bad)
