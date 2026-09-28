# -*- coding: utf-8 -*-
"""Unit tests for stellar/luminosity.py (offline: numpy + scipy only)."""

import pytest
import numpy as np
from scipy import stats

from scopesim_targets.stellar.luminosity import sample_magnitudes


def _exact_cdf(mags, m_min, m_max, slope):
    if slope == 0:
        return (mags - m_min) / (m_max - m_min)
    rate = slope * np.log(10)
    return np.expm1(rate * (mags - m_min)) / np.expm1(rate * (m_max - m_min))


@pytest.mark.parametrize("slope", (0.3, 1.0, 5.0, 0.0, 1e-9, -0.3, -5.0))
def test_matches_exponential_star_counts(slope):
    mags = sample_magnitudes(20_000, (15, 22), slope, np.random.default_rng(2))
    assert mags.min() >= 15 and mags.max() <= 22
    result = stats.kstest(mags, lambda m: _exact_cdf(m, 15, 22, slope))
    assert result.pvalue > 1e-3


def test_positive_slope_favors_faint_end():
    mags = sample_magnitudes(10_000, (15, 22), 0.3, np.random.default_rng(1))
    assert np.median(mags) > 18.5  # uniform would be 18.5


def test_limits_in_either_order():
    a = sample_magnitudes(10, (15, 22), 0.3, np.random.default_rng(1))
    b = sample_magnitudes(10, (22, 15), 0.3, np.random.default_rng(1))
    np.testing.assert_array_equal(a, b)


def test_reproducible_with_same_generator_state():
    a = sample_magnitudes(10, (15, 22), 0.3, np.random.default_rng(7))
    b = sample_magnitudes(10, (15, 22), 0.3, np.random.default_rng(7))
    np.testing.assert_array_equal(a, b)


def test_zero_stars():
    assert sample_magnitudes(0, (15, 22), 0.3, np.random.default_rng()).shape == (0,)


@pytest.mark.parametrize(
    ("args", "exc"),
    (
        ((5, (15, 15), 0.3), ValueError),
        ((5, (15, 22), np.inf), ValueError),
        ((-1, (15, 22), 0.3), ValueError),
        ((2.5, (15, 22), 0.3), TypeError),
        ((True, (15, 22), 0.3), TypeError),
    ),
)
def test_invalid(args, exc):
    with pytest.raises(exc):
        sample_magnitudes(*args, rng=np.random.default_rng())
