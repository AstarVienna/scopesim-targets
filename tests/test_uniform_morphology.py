# -*- coding: utf-8 -*-
"""Unit tests for UniformRectangularMorphology and its area helper (offline)."""

import pytest
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

from scopesim_targets.stellar.morphology import (
    UniformRectangularMorphology,
    rectangle_area_outside_circle,
)


@pytest.fixture
def center():
    return SkyCoord(0 * u.deg, 0 * u.deg)


class TestRectangleAreaOutsideCircle:
    @pytest.mark.parametrize("radius", (0, 0.5, 1, 2, 3, 3.1))
    def test_matches_monte_carlo(self, radius):
        # 6 x 2 rectangle: the circle stays inside, crosses the long sides,
        # crosses all sides, and nearly covers everything.
        rng = np.random.default_rng(0)
        points = rng.uniform((-3, -1), (3, 1), size=(1_000_000, 2))
        monte_carlo = 12 * np.mean(np.hypot(*points.T) >= radius)
        np.testing.assert_allclose(
            rectangle_area_outside_circle(6, 2, radius), monte_carlo,
            rtol=0.02, atol=5e-3,
        )

    def test_inscribed_circle(self):
        np.testing.assert_allclose(
            rectangle_area_outside_circle(10, 10, 2), 100 - np.pi * 4
        )

    @pytest.mark.parametrize("radius", (np.sqrt(10), 5, 100))
    def test_covered(self, radius):
        assert rectangle_area_outside_circle(6, 2, radius) == 0


class TestUniformRectangularMorphology:
    def test_within_rectangle(self, center):
        x, y = UniformRectangularMorphology(2000, (60, 20), rng_seed=1).sample(center)
        assert len(x) == len(y) == 2000
        assert np.abs(x).max() <= 30 and np.abs(y).max() <= 10
        # roughly uniform: quadrants hold ~ a quarter each
        counts = np.histogram2d(x, y, bins=2, range=((-30, 30), (-10, 10)))[0]
        np.testing.assert_allclose(counts / 2000, 0.25, atol=0.04)

    def test_exclusion_is_exact_count_and_empty_core(self, center):
        x, y = UniformRectangularMorphology(
            500, (60, 20), exclude_radius=8, rng_seed=1
        ).sample(center)
        assert len(x) == 500
        assert np.hypot(x, y).min() >= 8

    def test_square_from_scalar(self, center):
        morph = UniformRectangularMorphology(10, 1 * u.arcmin)
        assert morph.extent_arcsec(center) == (60, 60, 0)

    def test_repeated_sampling_is_one_realization(self, center):
        morph = UniformRectangularMorphology(50, 30, exclude_radius=5, rng_seed=None)
        np.testing.assert_array_equal(morph.sample(center), morph.sample(center))

    def test_seeded(self, center):
        a = UniformRectangularMorphology(50, 30, rng_seed=4).sample(center)
        b = UniformRectangularMorphology(50, 30, rng_seed=4).sample(center)
        c = UniformRectangularMorphology(50, 30, rng_seed=5).sample(center)
        np.testing.assert_array_equal(a, b)
        assert not np.array_equal(a[0], c[0])

    def test_length_size_needs_distance(self, center):
        morph = UniformRectangularMorphology(10, 1 * u.pc)
        with pytest.raises(ValueError, match="needs a field center"):
            morph.sample(center)
        at_1kpc = SkyCoord(0 * u.deg, 0 * u.deg, 1 * u.kpc)
        width, _, _ = morph.extent_arcsec(at_1kpc)
        np.testing.assert_allclose(width, 206.2648, rtol=1e-6)

    def test_area(self, center):
        area = UniformRectangularMorphology(1, 10, exclude_radius=2).area(center)
        np.testing.assert_allclose(area.to_value(u.arcsec**2), 100 - 4 * np.pi)

    def test_zero_stars(self, center):
        x, y = UniformRectangularMorphology(0, 10).sample(center)
        assert x.shape == y.shape == (0,)

    def test_fully_excluded_raises(self, center):
        with pytest.raises(ValueError, match="covers the whole rectangle"):
            UniformRectangularMorphology(5, 10, exclude_radius=8).sample(center)

    @pytest.mark.parametrize(
        ("kwargs", "exc"),
        (
            ({"n_stars": -1, "size": 10}, ValueError),
            ({"n_stars": 1.5, "size": 10}, TypeError),
            ({"n_stars": 1, "size": (1, 2, 3)}, ValueError),
        ),
    )
    def test_invalid_init(self, kwargs, exc):
        with pytest.raises(exc):
            UniformRectangularMorphology(**kwargs)

    @pytest.mark.parametrize(
        "kwargs",
        ({"size": 0}, {"size": (10, -1)}, {"size": 10, "exclude_radius": -1}),
    )
    def test_invalid_extent(self, kwargs, center):
        with pytest.raises(ValueError):
            UniformRectangularMorphology(1, **kwargs).sample(center)
