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
from astropy.modeling.functional_models import (
    GeneralSersic2D,
    Ring2D,
    Box2D,
    Disk2D,
    Const2D,
)

from scopesim import Source
from scopesim.source.source_fields import ImageSourceField

from .typing_utils import POSITION_TYPE, SPECTRUM_TYPE, BRIGHTNESS_TYPE
from .target import SpectrumTarget
from .brightness import BrightnessError, AmountKind, Brightness


def _as_arcsec(param) -> u.Quantity:
    """Interpret an astropy model parameter as an angle in arcsec.

    Rendering happens on an arcsec grid, so a bare (unitless) parameter is taken
    to be in arcsec by that convention; a parameter carrying a unit is converted.
    """
    unit = param.unit if param.unit is not None else u.arcsec
    return (param.value * unit).to(u.arcsec)


class ExtendedSourceTarget(SpectrumTarget):
    """Base class for Extended Source Targets."""


class ParametrizedTarget(ExtendedSourceTarget):
    """Base class for Targets defined via an astropy model class.

    Holds the astropy model (``_model``) and the rendering/WCS helpers. It does
    not define ``to_source`` -- the flux-normalization contract lives on
    :class:`BrightnessProfile`.
    """

    _model_cls: ClassVar[type]

    def _render_image(self, optical_train):
        pixel_scale = u.pixel_scale(optical_train["pixel_scale"])
        with u.set_enabled_equivalencies(pixel_scale):
            width = optical_train["width"] * u.pixel
            height = optical_train["height"] * u.pixel
            # Render in a bare arcsec numeric space: model params live in arcsec
            # by convention (see _coerce_param / _as_arcsec), and mixing Quantity
            # coords with bare params breaks model evaluation.
            x = (
                np.arange(width.value) * u.pixel - (width / 2 - 0.5 * u.pixel)
            ).to_value(u.arcsec)
            y = (
                np.arange(height.value) * u.pixel
                - (height / 2 - 0.5 * u.pixel)
            ).to_value(u.arcsec)
        coords = np.meshgrid(x, y)
        return self._model.render(coords=coords)

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
        cdelt = np.array(
            2
            * [
                (1 * u.pix).to_value(
                    u.arcsec, u.pixel_scale(optical_train["pixel_scale"])
                )
            ]
        )

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
    def _coerce_param(key: str, value):
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
        return value.value

    def __init__(
        self,
        position: POSITION_TYPE | None = None,
        spectrum: SPECTRUM_TYPE | None = None,
        brightness: BRIGHTNESS_TYPE | None = None,
        params: Mapping[str, u.Quantity | Number] | None = None,
        anchor: str | None = None,
    ) -> None:
        if position is not None:
            self.position = position
        if spectrum is not None:
            self.spectrum = spectrum
        if brightness is not None:
            self.brightness = brightness
        if anchor is not None:
            self.anchor = anchor

        params = {
            k: self._coerce_param(k, v) for k, v in dict(params or {}).items()
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
        (:meth:`total_flux_factor`, field-of-view independent). Non-integrable
        profiles (:class:`Flat`, untruncated power laws) have no intrinsic
        total, so they override this with the *rendered* field of view
        (``Omega_pixel * N_pixels``) -- the FOV closes the open integral, giving
        a legitimately FOV-dependent total. Everything downstream (weight map,
        surface-brightness reduction) is identical; only this area differs.
        """
        return self.total_flux_factor()

    def _weightmap(self, optical_train):
        """Render the analytic model as a unit-normalized weight map + header.

        The image is normalized by :meth:`_effective_area` (never a grid sum):
        for finite profiles that is the analytic total, so flux clipped by the
        window is not redistributed; for non-integrable profiles it is the field
        of view, giving a uniform weight map that sums to 1. Separated from
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
            return spectrum
        return spectrum * self._anchored_spectrum_scale(spectrum, brightness)


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


class Disk(BrightnessProfile):
    """Uniform filled disk profile (radius ``R_0``)."""

    _model_cls = Disk2D
    has_finite_total = True
    sb_reference = "uniform"

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

    def total_flux_factor(self) -> u.Quantity:
        r_in = self._model.r_in << u.arcsec
        r_out = r_in + (self._model.width << u.arcsec)
        return np.pi * (r_out**2 - r_in**2)


class Sersic(BrightnessProfile):
    """Single-component Sersic profile (``GeneralSersic2D``)."""

    _model_cls = GeneralSersic2D
    has_finite_total = True
    sb_reference = "at r_eff"

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
        n_pixels = optical_train["width"] * optical_train["height"]
        return self._pixel_area(optical_train) * n_pixels
