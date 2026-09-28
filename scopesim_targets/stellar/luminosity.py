# -*- coding: utf-8 -*-
"""Magnitude distributions (luminosity functions) for star fields."""

from numbers import Integral

import numpy as np
from scipy import stats

__all__ = ["sample_magnitudes"]


def sample_magnitudes(
    n_stars: int,
    mag_range: tuple[float, float],
    slope: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Draw magnitudes from exponential star counts on `mag_range`.

    Star counts follow ``dN/dm ∝ 10**(slope * m)``, i.e.
    ``slope = d log10(N) / dm``. Positive slopes (the usual case: fainter
    stars are more numerous) put most stars near the faint end, ``slope=0``
    is uniform in magnitude.

    Parameters
    ----------
    n_stars : int
        Number of magnitudes to draw (0 gives an empty array).
    mag_range : (float, float)
        Bright and faint limits, in either order.
    slope : float
        Logarithmic count slope ``d log10(N) / dm``. Realistic values depend
        on band and galactic latitude, typically a few tenths.
    rng : numpy.random.Generator
        Source of randomness.

    Returns
    -------
    numpy.ndarray
        `n_stars` magnitudes.

    Notes
    -----
    This is a truncated exponential in magnitude, drawn with scipy's
    ``truncexpon`` (exact analytic inverse CDF). scipy allows no negative
    scale, so a positive slope -- a PDF growing towards the faint end -- is
    drawn as ``m_max - T`` with ``T ~ truncexpon``; a negative slope is
    ``m_min + T``. ``slope=0`` is degenerate for ``truncexpon`` and uses
    ``uniform``.
    """
    if not isinstance(n_stars, Integral) or isinstance(n_stars, bool):
        raise TypeError(f"n_stars must be an int, got {n_stars!r}")
    if n_stars < 0:
        raise ValueError(f"n_stars must be non-negative, got {n_stars}")
    m_min, m_max = sorted(float(m) for m in mag_range)
    if not m_min < m_max:
        raise ValueError(f"mag_range must have two distinct limits, got {mag_range!r}")
    if not np.isfinite(slope):
        raise ValueError(f"slope must be finite, got {slope!r}")

    width = m_max - m_min
    if slope == 0:
        return stats.uniform(loc=m_min, scale=width).rvs(n_stars, random_state=rng)

    rate = abs(slope) * np.log(10)  # 10**(slope * m) == exp(slope ln10 * m)
    offsets = stats.truncexpon(b=rate * width, scale=1 / rate).rvs(
        n_stars, random_state=rng
    )
    return m_max - offsets if slope > 0 else m_min + offsets
