# -*- coding: utf-8 -*-
"""Parametrized and discrete 2+1D and 3D target models."""

from collections.abc import Mapping
from numbers import Number  # matches int, float and all the numpy scalars

import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.modeling import Model
from astropy.modeling.functional_models import GeneralSersic2D

from scopesim import Source
from scopesim.source.source_fields import ImageSourceField

from .typing_utils import POSITION_TYPE, SPECTRUM_TYPE, BRIGHTNESS_TYPE
from .yaml_constructors import register_target_constructor
from .target import SpectrumTarget


class ExtendedSourceTarget(SpectrumTarget):
    """Base class for Extended Source Targets."""


class ParametrizedTarget(ExtendedSourceTarget):
    """Base class for Targets defined via an astropy model class."""

    # _model_cls: Model | None = None

    def to_source(self, optical_train) -> Source:
        """Convert to ScopeSim Source object."""
        img = self._render_image(optical_train)
        img /= img.sum()  # normalize to weightmap

        wcs = self._create_wcs(optical_train)
        hdr = wcs.to_header()
        # TODO: Revisit this!!
        # hdr["BUNIT"] = u.Unit("photlam arcsec-2").to_string("fits")
        hdr["BUNIT"] = u.Unit("photlam").to_string("fits")

        hdu = fits.ImageHDU(header=hdr, data=img)

        spectrum = self.resolve_spectrum(
            self.spectrum).scale_to_magnitude(
                self.brightness.mag, self.brightness.band)

        source = Source(field=ImageSourceField(
            hdu, spectra={0: spectrum}
        ))
        return source

    def _render_image(self, optical_train):
        pixel_scale = u.pixel_scale(optical_train["pixel_scale"])
        with u.set_enabled_equivalencies(pixel_scale):
            width = optical_train["width"] * u.pixel
            height = optical_train["height"] * u.pixel
            # TODO: Find a more elegant implementation here
            coords = np.meshgrid(
                np.arange(width.value) * u.pixel - (width / 2 - .5 * u.pixel),
                np.arange(height.value) * u.pixel - (height / 2 - .5 * u.pixel),
            ) << u.arcsec
        return self._model.render(coords=coords)

    def _create_wcs(self, optical_train) -> WCS:
        naxis = np.array([optical_train["width"], optical_train["height"]])
        crpix = (naxis + 1) / 2
        crval = np.array([0, 0])  # TODO: Add support for position here
        # TODO: Use proper u.pixel_scale equivalency
        cdelt = np.array(2 * [(1*u.pix).to_value(
            u.arcsec, u.pixel_scale(optical_train["pixel_scale"])
        )])

        wcs = WCS(naxis=2)
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        wcs.wcs.cunit = 2 * [u.arcsec]
        wcs.wcs.crpix = crpix
        wcs.wcs.crval = crval
        wcs.wcs.cdelt = cdelt
        return wcs
