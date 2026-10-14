# -*- coding: utf-8 -*-
"""Multiplicative random textures for extended targets.

A texture multiplies the smooth (mean) profile of a
:class:`~scopesim_targets.extended_source.BrightnessProfile` with a
reproducible, strictly positive random field of unit mean. The smooth profile
keeps its meaning: its surface brightness is the *expected* level, and its
analytic total stays the expected total.
"""
from dataclasses import dataclass
from functools import cached_property

import numpy as np
import astropy.units as u
from scipy.ndimage import map_coordinates


@dataclass(frozen=True)
class LognormalTexture:
    """Lognormal random field with a broken power-law power spectrum.

    ``f = exp(delta) / <exp(delta)>``, where ``delta`` is a Gaussian random
    field with standard deviation ``sigma`` and a 2D power spectrum

    * flat for ``k < 1/outer_scale`` (uncorrelated on large scales),
    * ``k**-beta`` between ``1/outer_scale`` and ``1/inner_scale``,
    * ``k**-beta_inner`` above ``1/inner_scale`` (continuous at the break).

    Wavenumbers are in cycles per unit length (``k = 1/lambda``). Without
    ``inner_scale``/``beta_inner`` the power law extends to the grid
    resolution. ``sigma`` is the standard deviation of ``ln f``; ``f`` has unit
    mean over its (periodic) generation box, so it redistributes flux but does
    not change the expected total.

    Scales are lengths (e.g. pc; sampling then needs the target distance) or
    angles. The field is generated once, in these units, on a fixed periodic
    grid of ``n_grid**2`` cells spaced ``inner_scale/4`` (or
    ``outer_scale/16``); the same seed therefore gives the same structure at
    every distance and pixel scale. ``n_grid`` is a numerical choice: the box
    must comfortably exceed the region of interest, beyond which the field
    repeats.
    """

    sigma: float
    beta: float
    outer_scale: u.Quantity
    rng_seed: int
    inner_scale: u.Quantity | None = None
    beta_inner: float | None = None
    n_grid: int = 2048

    def __post_init__(self):
        if not self.sigma >= 0:
            raise ValueError("sigma must be >= 0")
        if not self.beta > 0:
            raise ValueError("beta must be > 0")
        if (self.inner_scale is None) != (self.beta_inner is None):
            raise ValueError("inner_scale and beta_inner must be given together")
        if not (isinstance(self.rng_seed, (int, np.integer)) and self.rng_seed >= 0):
            raise ValueError("rng_seed must be a non-negative integer")
        if self.n_grid < 64 or self.n_grid % 2:
            raise ValueError("n_grid must be an even number >= 64")
        kind = self.outer_scale.unit.physical_type
        if kind not in ("length", "angle"):
            raise ValueError(f"outer_scale must be a length or an angle, got {self.outer_scale!r}")
        if self.inner_scale is not None:
            if self.inner_scale.unit.physical_type != kind:
                raise ValueError("inner_scale and outer_scale must both be lengths or both angles")
            if not 0 < self.inner_scale < self.outer_scale:
                raise ValueError("need 0 < inner_scale < outer_scale")
        if self.n_grid * self.spacing < 4 * self.outer_scale:
            raise ValueError("generation box smaller than 4 x outer_scale; increase n_grid")

    # ---- geometry ---------------------------------------------------------
    @property
    def unit(self) -> u.Unit:
        return self.outer_scale.unit

    @property
    def is_physical(self) -> bool:
        """True if the scales are lengths (sampling needs a distance)."""
        return self.outer_scale.unit.physical_type == "length"

    @property
    def spacing(self) -> u.Quantity:
        scale = self.inner_scale if self.inner_scale is not None else self.outer_scale / 4
        return (scale / 4).to(self.unit)

    @property
    def box_size(self) -> u.Quantity:
        """Period of the field."""
        return self.n_grid * self.spacing

    # ---- field ------------------------------------------------------------
    def power_spectrum(self, k):
        """Unnormalized 2D power spectrum at wavenumber k [cycles/unit]."""
        k_out = 1 / self.outer_scale.to_value(self.unit)
        k = np.asarray(k, float)
        p = np.where(k < k_out, 1.0, (np.maximum(k, k_out) / k_out) ** -self.beta)
        if self.inner_scale is not None:
            k_in = 1 / self.inner_scale.to_value(self.unit)
            tail = (k_in / k_out) ** -self.beta * (np.maximum(k, k_in) / k_in) ** -self.beta_inner
            p = np.where(k >= k_in, tail, p)
        return p

    @cached_property
    def log_field(self) -> np.ndarray:
        """``ln f`` on the periodic generation grid, std exactly ``sigma``."""
        n, dx = self.n_grid, self.spacing.to_value(self.unit)
        freq = np.fft.fftfreq(n, d=dx)
        k = np.hypot(*np.meshgrid(freq, freq, indexing="ij"))
        amp = np.sqrt(self.power_spectrum(k))
        amp[0, 0] = 0.0
        noise = np.random.default_rng(self.rng_seed).standard_normal((n, n))
        delta = np.real(np.fft.ifft2(np.fft.fft2(noise) * amp))
        delta -= delta.mean()
        std = delta.std()
        delta = delta * (self.sigma / std) if std > 0 else delta
        return delta - np.log(np.mean(np.exp(delta)))

    @cached_property
    def field(self) -> np.ndarray:
        """``f`` on the periodic generation grid (unit mean, strictly positive)."""
        return np.exp(self.log_field)

    # ---- sampling ---------------------------------------------------------
    def sample(self, x, y, pixel_size=None) -> np.ndarray:
        """Texture on the grid spanned by the 1D axes ``x`` (columns) and
        ``y`` (rows), as a ``(len(y), len(x))`` array.

        ``x``, ``y`` and ``pixel_size`` are Quantities in units compatible with
        the texture scales. With ``pixel_size``, each value is the mean over
        the pixel (sub-sampled at the field resolution, at most 16 x 16).
        """
        dx = self.spacing.to_value(self.unit)
        xs = np.asarray(x.to_value(self.unit), float)
        ys = np.asarray(y.to_value(self.unit), float)
        k = 1
        if pixel_size is not None:
            k = int(np.clip(np.ceil(pixel_size.to_value(self.unit) / dx), 1, 16))
            pix = pixel_size.to_value(self.unit)
        out = np.zeros((ys.size, xs.size))
        offs = ((np.arange(k) + 0.5) / k - 0.5) * (pix if k > 1 else 0.0)
        for oy in offs:
            for ox in offs:
                jj, ii = np.meshgrid((ys + oy) / dx, (xs + ox) / dx, indexing="ij")
                out += map_coordinates(self.field, [jj, ii], order=1, mode="grid-wrap")
        return out / k**2
