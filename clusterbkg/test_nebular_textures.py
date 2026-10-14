# -*- coding: utf-8 -*-
"""Tests for nebular_textures.py (Orion texture statistics, Arthur et al. 2016).

All pure (numpy/scipy/astropy); the calibrations take a few seconds.
"""
import numpy as np
import pytest
from numpy import testing as npt
import astropy.units as u

from . import nebular_textures as nt


def _arcsec_at_scene(angle):
    return angle.to_value(u.rad) * nt.SCENE_DISTANCE.to_value(u.pc)


class TestConstants:
    def test_inner_slope_converts_regime_iii(self):
        npt.assert_allclose(nt.ORION_BETA_INNER, 1 - (-4.26))

    def test_sigma_obs_is_mid_range(self):
        npt.assert_allclose(nt.ORION_SIGMA_OBS, 0.525)

    def test_scales_are_measured_angles_at_scene_distance(self):
        tex = nt.orion_texture(0.3 * u.pc, 0, sigma=0.5, n_grid=256)
        npt.assert_allclose(tex.outer_scale.to_value(u.pc), _arcsec_at_scene(22 * u.arcsec))
        npt.assert_allclose(tex.inner_scale.to_value(u.pc), _arcsec_at_scene(8 * u.arcsec))
        assert tex.is_physical


class TestCalibration:
    def test_beta_constant_matches_calibration(self):
        npt.assert_allclose(nt.calibrate_beta(), nt.ORION_BETA, atol=0.02)

    def test_texture_reproduces_measured_slit_slope(self):
        # Independent seeds, physical units: the 1D slope in the measured
        # window must agree with Table 2 within its stated uncertainty.
        slopes = []
        for seed in range(20, 26):
            tex = nt.orion_texture(0.3 * u.pc, seed, sigma=1.0, n_grid=1024)
            k_lo = 1 / tex.outer_scale.to_value(u.pc)
            k_hi = 1 / tex.inner_scale.to_value(u.pc)
            slopes.append(nt.slit_slope(tex, k_lo, k_hi))
        assert abs(np.mean(slopes) - nt.ORION_GAMMA_T2) < nt.ORION_GAMMA_T2_ERR

    def test_profile_width_decreases_with_core(self):
        widths = [nt.profile_width(r * u.pc) for r in (0.1, 0.2, 0.4)]
        assert widths[0] > widths[1] > widths[2]

    def test_minimum_core(self):
        r_min = nt.minimum_r_core()
        assert 0.15 * u.pc < r_min < 0.25 * u.pc
        npt.assert_allclose(nt.profile_width(r_min), nt.ORION_SIGMA_OBS, atol=1e-3)

    def test_king_core_leaves_no_room_for_texture(self):
        with pytest.warns(UserWarning, match="already more structured"):
            assert nt.calibrate_sigma(0.13 * u.pc) == 0.0

    def test_calibrated_sigma_reproduces_observed_width(self):
        # Forward check with seeds not used in the calibration.
        r_core = 0.3 * u.pc
        sigma = nt.calibrate_sigma(r_core)
        profile, xs, ys = nt._profile_on_field(r_core, nt.DEFAULT_INDEX)
        x = nt._angle_to_length(xs * u.arcsec)
        y = nt._angle_to_length(ys * u.arcsec)
        pix = nt._angle_to_length(nt.ORION_BIN)
        widths = [nt.pdf_core_width(profile * nt.orion_texture(
            r_core, s, sigma=sigma, n_grid=256).sample(x, y, pixel_size=pix)) for s in range(30, 38)]
        assert abs(np.mean(widths) - nt.ORION_SIGMA_OBS) < 0.04

    def test_sigma_does_not_depend_on_seed(self):
        a = nt.orion_texture(0.3 * u.pc, 1, n_grid=256)
        b = nt.orion_texture(0.3 * u.pc, 2, n_grid=256)
        assert a.sigma == b.sigma > 0
