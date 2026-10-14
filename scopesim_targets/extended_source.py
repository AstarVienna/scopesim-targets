# -*- coding: utf-8 -*-
"""Parametrized and discrete 2+1D and 3D target models."""

from typing import ClassVar
from dataclasses import replace
from collections.abc import Mapping
from numbers import Number  # matches int, float and all the numpy scalars

import numpy as np
from scipy.special import gamma, gammaincinv
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import Distance
from astropy.modeling import Fittable2DModel, Parameter
from astropy.modeling.functional_models import (
    GeneralSersic2D,
    Gaussian2D,
    Ring2D,
    Box2D,
    Disk2D,
    Const2D,
)

from scopesim import Source
from scopesim.source.source_fields import ImageSourceField

from .typing_utils import POSITION_TYPE, SPECTRUM_TYPE, BRIGHTNESS_TYPE
from .target import SpectrumTarget, length_angle_context
from .texture import LognormalTexture
from .brightness import BrightnessError, AmountKind, Brightness


def _as_arcsec(param) -> u.Quantity:
    """Interpret an astropy model parameter as an angle in arcsec.

    Rendering happens on an arcsec grid, so a bare (unitless) parameter is taken
    to be in arcsec by that convention; a parameter carrying a unit is converted.
    """
    unit = param.unit if param.unit is not None else u.arcsec
    return (param.value * unit).to(u.arcsec)


class CoredPowerLaw2D(Fittable2DModel):
    """Circular cored power law, normalized to ``amplitude`` at ``r_ref``.

    .. math::

        I(r) = A \\left[\\frac{1 + (r/r_c)^2}{1 + (r_{ref}/r_c)^2}\\right]^{-p/2}

    Flat (``I ~ const``) inside the core radius ``r_c``, falling as
    ``r**-p`` outside it. ``amplitude`` is the surface brightness *at*
    ``r_ref`` -- not at the peak -- so a brightness measured at a known
    distance from the illuminating source can be used directly. ``p = 0``
    degenerates to a constant (cf. :class:`Flat`). Same functional family as
    ``Moffat2D`` (``p = 2 alpha``, ``r_c = gamma``), reparametrized so the
    power-law index and the reference point are explicit.

    ``r_out`` (default ``inf``) truncates the profile: ``I = 0`` for
    ``r > r_out``. The normalization at ``r_ref`` is that of the untruncated
    law, so truncation never changes the surface brightness inside ``r_out``.
    """

    amplitude = Parameter(default=1.0)
    x_0 = Parameter(default=0.0)
    y_0 = Parameter(default=0.0)
    r_core = Parameter(default=1.0)
    index = Parameter(default=2.0)
    r_ref = Parameter(default=0.0)
    r_out = Parameter(default=np.inf)

    @staticmethod
    def evaluate(x, y, amplitude, x_0, y_0, r_core, index, r_ref, r_out):
        d2 = (x - x_0) ** 2 + (y - y_0) ** 2
        ref = 1.0 + (r_ref / r_core) ** 2
        prof = amplitude * ((1.0 + d2 / r_core**2) / ref) ** (-index / 2.0)
        return np.where(d2 <= r_out**2, prof, 0.0)


class ExtendedSourceTarget(SpectrumTarget):
    """Base class for Extended Source Targets."""


class ParametrizedTarget(ExtendedSourceTarget):
    """Base class for Targets defined via an astropy model class.

    Holds the astropy model (``_model``) and the rendering/WCS helpers. It does
    not define ``to_source`` -- the flux-normalization contract lives on
    :class:`BrightnessProfile`.
    """

    _model_cls: ClassVar[type]

    # Sub-pixel sampling factor for :meth:`_render_image`. 1 = plain pixel-
    # center sampling (exact enough for smooth, well-resolved profiles); k > 1
    # averages k*k midpoint samples per pixel, approximating the *pixel
    # integral*. Profiles with a sharp edge (Disk, Ring) or a cusp (Sersic)
    # set this: pixel-center sampling of a discontinuity mis-weights every
    # edge pixel by up to a full pixel, so the weight-map total quantizes with
    # the pixel phase (up to tens of percent for coarse scales) instead of
    # matching the analytic normalization. See T-QUANT in the tests.
    _oversample: ClassVar[int] = 1

    @staticmethod
    def _grid_centers(npix: int, scale: float) -> np.ndarray:
        """Pixel-center coordinates [arcsec] of an npix-wide centered axis."""
        return (np.arange(npix) - (npix / 2 - 0.5)) * scale

    @staticmethod
    def _scale_arcsec(optical_train) -> float:
        """Pixel scale as a bare arcsec/pixel number (the render space)."""
        return (1 * u.pix).to_value(
            u.arcsec, u.pixel_scale(optical_train["pixel_scale"])
        )

    def _render_image(self, optical_train):
        """Unit-amplitude profile as rendered: the smooth shape
        (:meth:`_render_smooth`) times the texture, if one is set.

        Everything downstream (weight map, effective area) uses this, so a
        texture enters the flux bookkeeping in one place only.
        """
        img = self._render_smooth(optical_train)
        if getattr(self, "texture", None) is not None:
            img = img * self._render_texture(optical_train)
        return img

    def _render_texture(self, optical_train):
        """Pixel-averaged texture on the render grid, centered like the model.

        Pixel coordinates are converted from the arcsec render space to the
        texture's own units -- with the target distance for physical textures
        -- so the same texture gives the same structure at every distance.
        """
        scale = self._scale_arcsec(optical_train)
        x = self._grid_centers(int(optical_train["width"]), scale) * u.arcsec
        y = self._grid_centers(int(optical_train["height"]), scale) * u.arcsec
        pix = scale * u.arcsec
        if self.texture.is_physical:
            with length_angle_context(self._distance):
                x, y, pix = (q.to(self.texture.unit) for q in (x, y, pix))
        return self.texture.sample(x, y, pixel_size=pix)

    def _render_smooth(self, optical_train):
        """Render the unit-amplitude model, pixel-averaged, FITS-oriented.

        Rendering happens in a bare arcsec numeric space: model params live in
        arcsec by convention (see ``_coerce_param`` / ``_as_arcsec``), and
        mixing Quantity coords with bare params breaks model evaluation.

        Orientation: ``Model.render(coords=...)`` expects ``np.indices``-style
        ordering (slow axis first) and evaluates ``model(*coords[::-1])``, so
        the y grid goes in slot 0 and the returned array is ``(height, width)``
        -- numpy's last axis is FITS NAXIS1, matching the WCS header built in
        :meth:`_create_wcs`. (The previous ``meshgrid(x, y)`` ordering handed
        the model swapped axes; on the square grids used so far the swap
        cancelled against the numpy->FITS axis order, but with
        ``width != height`` it put ``x_width`` along the height axis.)

        Each pixel's value is the *mean* of ``_oversample**2`` midpoint samples
        -- the midpoint-rule approximation of the pixel integral of the
        profile. With ``_oversample = 1`` this reduces to the previous
        pixel-center sampling. The subpixel offsets are accumulated one shifted
        full-resolution render at a time, so peak memory stays at one
        ``(height, width)`` array regardless of the factor.
        """
        scale = self._scale_arcsec(optical_train)
        width = int(optical_train["width"])
        height = int(optical_train["height"])
        x = self._grid_centers(width, scale)
        y = self._grid_centers(height, scale)

        k = int(self._oversample)
        offsets = ((np.arange(k) + 0.5) / k - 0.5) * scale
        img = np.zeros((height, width))
        for dy in offsets:
            for dx in offsets:
                img += self._model.render(
                    coords=np.meshgrid(y + dy, x + dx, indexing="ij")
                )
        return img / k**2

    @staticmethod
    def _pixel_area(optical_train) -> u.Quantity:
        """Solid angle subtended by one pixel [arcsec**2]."""
        scale = (1 * u.pix).to_value(
            u.arcsec, u.pixel_scale(optical_train["pixel_scale"])
        )
        return (scale * u.arcsec) ** 2

    def _create_wcs(self, optical_train) -> WCS:
        naxis = np.array([optical_train["width"], optical_train["height"]])
        crpix = (naxis + 1) / 2
        crval = np.array([0, 0])  # TODO: Add support for position here
        # TODO: Use proper u.pixel_scale equivalency
        pxsc = u.pixel_scale(optical_train["pixel_scale"])
        cdelt = np.array(2 * [(1 * u.pix).to_value(u.arcsec, pxsc)])

        wcs = WCS(naxis=2)
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        wcs.wcs.cunit = 2 * [u.arcsec]
        wcs.wcs.crpix = crpix
        wcs.wcs.crval = crval
        wcs.wcs.cdelt = cdelt
        return wcs


class BrightnessProfile(ParametrizedTarget):
    """A parametrized extended target whose model is a *brightness profile*.

    Every concrete profile declares three traits:

    * ``has_finite_total`` -- whether a closed-form spatial integral exists;
    * ``sb_reference`` -- documented point at which a surface brightness sets
      the amplitude ("uniform", "at r_eff", ...);
    * ``total_flux_factor()`` -- that analytic integral for unit amplitude.

    ``brightness`` is the sole flux authority. The profile amplitude is managed
    internally (unit amplitude for the integrated weight map; set from the
    surface-brightness value in the SB path), never supplied by the user -- the
    reserved ``params`` keys enforce this (E8).
    """

    has_finite_total: ClassVar[bool]
    sb_reference: ClassVar[str]

    _RESERVED: ClassVar[dict[str, str]] = {
        "amplitude": "flux is owned by 'brightness' (see defining_brightness.md)",
        "x_0": "position is owned by 'position'/'offset'",
        "y_0": "position is owned by 'position'/'offset'",
    }

    @staticmethod
    def _coerce_param(key: str, value, distance: Distance | None = None):
        """Normalize a param to the bare arcsec render space.

        Lengths/on-sky angles -> arcsec value; the rotation angle ``theta`` ->
        radians value; dimensionless quantities -> their value; bare numbers are
        left untouched. Keeps rendering (bare coords) and ``total_flux_factor``
        (via ``_as_arcsec``) on one consistent footing regardless of whether the
        user or the YAML resolver supplied units.
        """
        if not isinstance(value, u.Quantity):
            return value
        if key == "theta":
            return value.to_value(u.rad)
        if value.unit.physical_type == "angle":
            return value.to_value(u.arcsec)
        if value.unit.physical_type == "length":
            if distance is None:
                raise ValueError(f"Physical length for {key} needs distance.")
            with length_angle_context(distance):
                return value.to_value(u.arcsec)
        return value.value

    def __init__(
        self,
        position: POSITION_TYPE | None = None,
        spectrum: SPECTRUM_TYPE | None = None,
        brightness: BRIGHTNESS_TYPE | None = None,
        params: Mapping[str, u.Quantity | Number] | None = None,
        anchor: str | None = None,
        extinction=None,
        texture: LognormalTexture | None = None,
    ) -> None:
        if position is not None:
            self.position = position
        if spectrum is not None:
            self.spectrum = spectrum
        if brightness is not None:
            self.brightness = brightness
        if anchor is not None:
            self.anchor = anchor
        if extinction is not None:
            self.extinction = extinction

        distance = None
        try:
            if self.position.distance.unit.physical_type == "length":
                distance = self.position.distance
        except AttributeError:
            pass
        self._distance = distance
        # Multiplicative random structure on top of the smooth profile. It has
        # unit mean, so ``brightness`` keeps referring to the smooth (expected)
        # profile and its analytic total stays the expected total.
        self.texture = texture
        if texture is not None and texture.is_physical and distance is None:
            raise ValueError("Physical length for texture needs distance.")
        params = {
            k: self._coerce_param(k, v, distance)
            for k, v in dict(params or {}).items()
        }
        self._validate_params(params)
        # Amplitude is internal, never user-supplied: unit amplitude here, so the
        # rendered image is a pure shape that the weight-map step normalizes.
        self._model = self._model_cls(amplitude=1.0, **params)
        # Required, not optional: render(coords=...) treats a bounding_box as a
        # 0-based pixel-index window, but we render on physical (arcsec) coords
        # centered on zero. With the default bbox that Box2D/Disk2D/Ring2D carry,
        # the carve lands off-grid and render() returns zeros (or raises). None
        # also prevents truncating a profile's wings inside the FOV -- flux
        # beyond the window is accounted for analytically by total_flux_factor.
        self._model.bounding_box = None

    @classmethod
    def _validate_params(cls, params: Mapping) -> None:
        """Reject reserved (E8) and unknown ``params`` keys (plan section 3.2)."""
        accepted = set(cls._model_cls.param_names) - set(cls._RESERVED)
        for key in params:
            if key in cls._RESERVED:
                raise BrightnessError(
                    "E8", f"'{key}' in params: {cls._RESERVED[key]}"
                )
            if key not in accepted:
                raise ValueError(
                    f"unknown parameter '{key}' for {cls.__name__}; "
                    f"accepted: {sorted(accepted)}"
                )

    def total_flux_factor(self) -> u.Quantity:
        """Analytic integral of the unit-amplitude profile over all space.

        Returns a solid angle [arcsec**2]. Raises for non-integrable profiles.
        Single source of truth for the profile's total -- the normalizer calls
        this instead of summing the rendered grid (which would silently
        redistribute flux clipped by the window).
        """
        raise NotImplementedError  # provided by each concrete profile

    def _effective_area(self, optical_train) -> u.Quantity:
        """Solid angle that multiplies a surface brightness to give the total
        flux carried by the spectrum, and that normalizes the weight map.

        For a finite-total profile this is the intrinsic analytic total
        ``P`` (:meth:`total_flux_factor`, field-of-view independent).
        Non-integrable profiles (:class:`Flat`, power laws with index <= 2)
        have no intrinsic total; the rendered field of view closes the open
        integral, ``A_eff = int_FOV p dOmega = Omega_pixel * sum(img)``, giving
        a legitimately FOV-dependent total. :class:`Flat` (``img == 1``) is the
        special case ``Omega_pixel * N_pixels`` and overrides this exactly.
        Everything downstream (weight map ``w_ij = Omega * p_ij / A_eff``,
        surface-brightness reduction ``F = S * A_eff``) is identical; only this
        area differs. Either way the flux landing in a pixel is
        ``S * Omega * p_ij``. With a texture ``f``, ``p`` is the textured
        profile ``p * f``: finite profiles keep the smooth analytic ``P`` (the
        expected total, since ``f`` has unit mean), non-integrable ones close
        the integral over ``p * f``.
        """
        if self.has_finite_total:
            return self.total_flux_factor()
        img = self._render_image(optical_train)
        return self._pixel_area(optical_train) * img.sum()

    def _weightmap(self, optical_train):
        """Render the analytic model as a unit-normalized weight map + header.

        The image is normalized by :meth:`_effective_area`: for finite
        profiles that is the analytic total (never a grid sum), so flux clipped
        by the window is not redistributed; for non-integrable profiles it is
        the profile's integral over the rendered field of view -- a grid sum by
        definition, since nothing outside the FOV is counted -- so the weight
        map sums to 1. Separated from
        :meth:`to_source` so the stored-artifact image is testable without a
        ScopeSim ``Source`` or a resolved spectrum. Carries no flux itself --
        the spectrum does. (E6, integrated-on-non-integrable, is a spectrum-side
        error and is raised in :meth:`to_source`, not here.)
        """
        img = self._render_image(optical_train)
        weightmap = u.Quantity(
            img
            * self._pixel_area(optical_train)
            / self._effective_area(optical_train)
        ).to_value(u.dimensionless_unscaled)

        wcs = self._create_wcs(optical_train)
        hdr = wcs.to_header()
        # Dimensionless weight map: the flux lives in the spectrum. A per-solid-
        # angle BUNIT on a stored weight map is a category error (flux spec 3.4).
        hdr["BUNIT"] = ""
        return weightmap, hdr

    def to_source(self, optical_train) -> Source:
        """Convert to a ScopeSim Source.

        Integrated brightness (Case I) and surface brightness (Case II) both
        produce the canonical stored artifact -- a dimensionless weight map plus
        a flux-calibrated spectrum -- for every profile, finite or not. A
        surface brightness is reduced to the equivalent integrated flux over the
        effective area (:meth:`_effective_integrated_brightness`); the only
        difference between profiles is whether that area is the analytic total
        or the field of view. An *integrated* brightness on a non-integrable
        profile has no finite total to carry, and is the one error (E6).
        """
        if (
            not self.brightness.is_surface_brightness
            and not self.has_finite_total
        ):
            raise BrightnessError(
                "E6",
                f"{type(self).__name__} has no finite analytic total, so an "
                "integrated brightness is undefined -- specify a surface "
                "brightness instead.",
            )

        weightmap, hdr = self._weightmap(optical_train)
        hdu = fits.ImageHDU(header=hdr, data=weightmap)

        # _scale_spectrum resolves, redshifts, and flux-scales in the normative
        # order (redshift before the anchor scale).
        spectrum = self._scale_spectrum(optical_train)

        return Source(field=ImageSourceField(hdu, spectra={0: spectrum}))

    def _effective_integrated_brightness(self, optical_train) -> Brightness:
        """The integrated brightness the spectrum is scaled to.

        Integrated brightness passes through unchanged. A surface brightness is
        reduced to the equivalent integrated amount over the effective area,
        with ``A_eff`` the analytic total for finite profiles or the field of
        view for non-integrable ones (see :meth:`_effective_area`):

        * magnitude: ``SB_mag - 2.5 log10(A_eff / solid_angle)`` -- the log
          keeps the conversion in linear space (never per-pixel magnitude
          area-scaling);
        * linear (flux density / energy flux per solid angle):
          ``total = SB * A_eff``, which strips the per-solid-angle divisor.

        Both express the one rule ``total = SB * A_eff`` at the profile's
        reference point (the amplitude *is* the surface brightness there); they
        differ only because magnitudes are logarithmic. The result is an
        integrated amount the shared scaler handles uniformly.
        """
        b = self.brightness
        if not b.is_surface_brightness:
            return b
        area = self._effective_area(optical_train)
        if b.amount_kind is AmountKind.MAG:
            ratio = (area / b.solid_angle).to_value(u.dimensionless_unscaled)
            implied = b.value - 2.5 * np.log10(ratio) * u.mag
        else:
            # e.g. (MJy / sr) * arcsec2 -> MJy; (W / (m2 arcsec2)) * arcsec2 -> W/m2
            implied = (b.value * area).to(b.value.unit * b.solid_angle)
        return replace(b, value=implied, solid_angle=None)

    def _scale_spectrum(self, optical_train):
        """Resolve, redshift, and flux-scale the spectrum to the brightness.

        Order matches the point-source path and the normative pipeline
        (rest SED -> redshift -> [screens] -> anchor): the flux-authority scale
        is applied last, on the redshifted spectrum. A surface brightness has
        already been reduced to an integrated amount by
        :meth:`_effective_integrated_brightness`, so the shared
        :meth:`~.target.SpectrumTarget._get_spectrum_scale` sees an integrated
        brightness and never raises E7.

        Blackbody spectra are the one exception: spextra bakes the flux scale
        into construction (see
        :meth:`~.target.SpectrumTarget._blackbody_amplitude`), so they are
        already normalized and only need redshifting.
        """
        brightness = self._effective_integrated_brightness(optical_train)
        spectrum = self.resolve_spectrum(self.spectrum, brightness)
        if getattr(self, "_position", None) is not None:
            spectrum = self.redshift_spectrum(spectrum, self.position)
        if isinstance(self.spectrum, str) and self.spectrum.startswith(
            "blackbody:"
        ):
            scaled = spectrum  # spextra baked the scale in at construction
        else:
            scaled = spectrum * self._anchored_spectrum_scale(spectrum, brightness)
        # Extended sources have a single SED (no dedup concern), so the screen is
        # applied in place with the true curve. The scale is a scalar, so this
        # commutes: w * (S * T) == (w * S) * T, and the anchor choice already
        # lives in w (observed scales the reddened SED, intrinsic the bare one).
        return self._redden(scaled)


class Box(BrightnessProfile):
    """Uniform rectangular box profile."""

    _model_cls = Box2D
    has_finite_total = True
    sb_reference = "uniform"

    def total_flux_factor(self) -> u.Quantity:
        return (
            (self._model.x_width << u.arcsec) *
            (self._model.y_width << u.arcsec)
        )

    def _render_smooth(self, optical_train):
        """Exact per-pixel coverage fractions (analytic pixel integral).

        ``Box2D`` is axis-aligned, so the pixel integral separates into two 1D
        box/pixel overlaps -- computed in closed form instead of sampled
        (``_oversample`` is irrelevant here). This makes the weight-map total
        *exactly* the visible fraction of the box at any pixel scale: 1 for a
        fully contained box regardless of whether the widths are commensurate
        with the pixel grid (no edge quantization), and exactly
        ``area(box & FOV) / area(box)`` when the box overflows the window --
        the spectrum keeps carrying the whole-box flux either way.
        """
        scale = self._scale_arcsec(optical_train)

        def coverage(npix: int, full_width: float) -> np.ndarray:
            """Overlap fraction of each pixel with ``[-w/2, +w/2]``."""
            centers = self._grid_centers(npix, scale)
            half = full_width / 2
            lo = np.maximum(centers - scale / 2, -half)
            hi = np.minimum(centers + scale / 2, half)
            return np.clip(hi - lo, 0.0, None) / scale

        cov_x = coverage(
            int(optical_train["width"]),
            _as_arcsec(self._model.x_width).to_value(u.arcsec),
        )
        cov_y = coverage(
            int(optical_train["height"]),
            _as_arcsec(self._model.y_width).to_value(u.arcsec),
        )
        # (height, width): FITS-oriented, same as the base render.
        return np.outer(cov_y, cov_x)


class Disk(BrightnessProfile):
    """Uniform filled disk profile (radius ``R_0``)."""

    _model_cls = Disk2D
    has_finite_total = True
    sb_reference = "uniform"
    # Sharp circular edge: no separable closed-form pixel integral (unlike
    # Box), so approximate it by oversampling; the residual rim error is
    # O(scale / (k * R_0)) instead of O(scale / R_0).
    _oversample = 10

    def total_flux_factor(self) -> u.Quantity:
        return np.pi * (self._model.R_0 << u.arcsec)**2


class Ring(BrightnessProfile):
    """Uniform annulus profile (inner radius ``r_in``, ``width``).

    This is the model formerly exposed as ``Disk`` (an ``astropy`` ``Ring2D``);
    it is a ring, not a filled disk, so it has been renamed. Use :class:`Disk`
    for a uniform filled disk.
    """

    _model_cls = Ring2D
    has_finite_total = True
    sb_reference = "uniform"
    _oversample = 10  # two sharp circular edges; see Disk

    def total_flux_factor(self) -> u.Quantity:
        r_in = self._model.r_in << u.arcsec
        r_out = r_in + (self._model.width << u.arcsec)
        return np.pi * (r_out**2 - r_in**2)


class Sersic(BrightnessProfile):
    """Single-component Sersic profile (``GeneralSersic2D``)."""

    _model_cls = GeneralSersic2D
    has_finite_total = True
    sb_reference = "at r_eff"
    # Smooth except for the central cusp, which pixel-center sampling
    # systematically under-weights for large n (~0.2-0.6% of the total at
    # n = 4 and 0.2-0.4 arcsec/px, drifting with scale); modest oversampling
    # makes the carried fraction sampling-stable at the 1e-4 level.
    _oversample = 4

    def total_flux_factor(self) -> u.Quantity:
        model = self._model
        if float(model.c.value) != 0.0:
            # Boxy isophotes (c != 0) rescale the enclosed area by a Gamma-
            # function factor; that correction is a later refinement.
            raise NotImplementedError(
                "analytic total for a boxy Sersic (c != 0) is not implemented"
            )
        n = float(model.n.value)
        ellip = float(model.ellip.value)
        b_n = gammaincinv(2 * n, 0.5)
        factor = 2 * np.pi * n * np.exp(b_n) * gamma(2 * n) * b_n ** (-2 * n)
        return factor * (1 - ellip) * _as_arcsec(model.r_eff) ** 2


class Gaussian(BrightnessProfile):
    """Elliptical 2D Gaussian profile.

    The clean cross-check among the profiles: smooth everywhere (no sharp edge,
    so pixel-center sampling is already pixel-integral-accurate and there is no
    rasterization quantization to speak of), with a trivial closed-form total
    ``2 pi A sigma_x sigma_y`` (rotation-invariant), and -- for ``theta = 0`` --
    a closed-form *window* fraction too (a product of error functions), so both
    the fully-contained and the clipped weight-map sums have independent
    analytic oracles. See T-QUANT-GAUSS in the tests.

    ``Gaussian2D`` names its center ``x_mean`` / ``y_mean`` (not ``x_0`` /
    ``y_0``), hence the class-specific reserved-parameter table.
    """

    _model_cls = Gaussian2D
    has_finite_total = True
    sb_reference = "at the peak"

    _RESERVED: ClassVar[dict[str, str]] = {
        "amplitude": "flux is owned by 'brightness' (see defining_brightness.md)",
        "x_mean": "position is owned by 'position'/'offset'",
        "y_mean": "position is owned by 'position'/'offset'",
    }

    def total_flux_factor(self) -> u.Quantity:
        return (
            2*np.pi
            * (self._model.x_stddev << u.arcsec)
            * (self._model.y_stddev << u.arcsec)
        )


class Flat(BrightnessProfile):
    """Infinite constant surface brightness (no finite analytic total).

    Only a surface brightness is valid (an integrated brightness is E6). The
    total flux is realized from the rendered field of view -- the FOV closes the
    otherwise-open integral -- so the stored artifact is the usual weight map
    (uniform, summing to 1) plus a flux-calibrated spectrum, and the total is
    legitimately field-of-view dependent.
    """

    _model_cls = Const2D
    has_finite_total = False
    sb_reference = "the constant"

    def total_flux_factor(self) -> u.Quantity:
        raise BrightnessError(
            "E6", "Flat has no finite analytic total (the integral diverges)"
        )

    def _effective_area(self, optical_train) -> u.Quantity:
        # The FOV closes the open integral: A_eff = Omega_pixel * N_pixels.
        # That shortcut needs img == 1; a texture makes the image non-uniform.
        if self.texture is not None:
            return super()._effective_area(optical_train)
        n_pixels = optical_train["width"] * optical_train["height"]
        return self._pixel_area(optical_train) * n_pixels


class CoredPowerLaw(BrightnessProfile):
    """Cored power-law surface brightness, ``I ~ [1 + (r/r_core)^2]^(-index/2)``.

    Meant for diffuse nebular backgrounds that fall off with distance from an
    illuminating star, e.g. the Orion Nebula, whose surface brightness falls
    roughly as ``D**-2`` near theta1 Ori C and closer to ``D**-3`` further out
    (O'Dell & Harris 2010). ``index = 0`` reproduces :class:`Flat`.

    The surface brightness in ``brightness`` applies *at* ``r_ref`` (default 0,
    i.e. the peak), which must lie inside ``r_out``. Finiteness depends on the
    instance: with an outer radius ``r_out`` (e.g. matched to a cluster's tidal
    radius) the total is always finite and FOV-independent. Untruncated
    (``r_out = inf``, the default), the integral converges only for
    ``index > 2``; for ``index <= 2`` the total is closed by the field of view,
    exactly as for :class:`Flat`, and an integrated brightness is E6.

    The profile center is meant to be the illuminating star, which is
    generally *not* the field center. Until extended targets honor
    ``position``/``offset`` in rendering, it is placed at the field center.
    """

    _model_cls = CoredPowerLaw2D
    sb_reference = "at r_ref"
    # Smooth for r_core >> pixel; modest oversampling keeps an unresolved
    # core (r_core ~ pixel) pixel-integral-accurate.
    _oversample = 4

    @property
    def has_finite_total(self) -> bool:  # instance-level: index and r_out
        return bool(
            np.isfinite(float(self._model.r_out.value))
            or float(self._model.index.value) > 2.0
        )

    @classmethod
    def _validate_params(cls, params: Mapping) -> None:
        super()._validate_params(params)
        if "r_core" in params and not float(params["r_core"]) > 0:
            raise ValueError("r_core must be > 0")
        if "index" in params and not float(params["index"]) >= 0:
            raise ValueError("index must be >= 0 (a falling or flat profile)")
        if "r_ref" in params and not float(params["r_ref"]) >= 0:
            raise ValueError("r_ref must be >= 0")
        if "r_out" in params and not float(params["r_out"]) > 0:
            raise ValueError("r_out must be > 0")
        if float(params.get("r_ref", 0.0)) > float(params.get("r_out", np.inf)):
            raise ValueError(
                "r_ref must lie inside r_out (the surface brightness would be "
                "anchored where the profile is zero)"
            )

    def total_flux_factor(self) -> u.Quantity:
        p = float(self._model.index.value)
        r_out = float(self._model.r_out.value)
        if not self.has_finite_total:
            raise BrightnessError(
                "E6", f"untruncated CoredPowerLaw with index={p} <= 2 has no "
                "finite total (the integral diverges); set r_out to truncate"
            )
        r_c = _as_arcsec(self._model.r_core)
        r_ref = _as_arcsec(self._model.r_ref)
        # With u = (r/r_c)^2:  int_0^R 2 pi r [1+(r/r_c)^2]^(-p/2) dr
        #   = pi r_c^2 G(X),  G(X) = int_0^X (1+u)^(-p/2) du,  X = (R/r_c)^2.
        # R = inf (p > 2): G = 2/(p-2).  Finite R: G = expm1(a L)/a with
        # a = 1 - p/2, L = log1p(X), and G = L in the p = 2 limit (a -> 0).
        # Dividing by I(r_ref)/I(0) moves the reference: factor ref**(p/2).
        ref = 1.0 + (r_ref / r_c).to_value(u.dimensionless_unscaled) ** 2
        if np.isfinite(r_out):
            x = (r_out * u.arcsec / r_c).to_value(u.dimensionless_unscaled) ** 2
            a, ln1px = 1.0 - p / 2.0, np.log1p(x)
            g = ln1px if abs(a * ln1px) < 1e-12 else np.expm1(a * ln1px) / a
        else:
            g = 2.0 / (p - 2.0)
        return np.pi * r_c**2 * g * ref ** (p / 2.0)
