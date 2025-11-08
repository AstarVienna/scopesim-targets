# -*- coding: utf-8 -*-
"""Contains main ``Target`` class."""

from abc import ABCMeta, abstractmethod
from collections.abc import Mapping

from astropy import units as u
from astropy.coordinates import SkyCoord, Angle
from synphot import SourceSpectrum

from astar_utils import SpectralType
from spextra import Spextrum, SpecLibrary


DEFAULT_LIBRARY = SpecLibrary("bosz/lr")


class Target(metaclass=ABCMeta):
    """Main class in scopesim-targets."""

    @abstractmethod
    def to_source(self):
        """Convert to ScopeSim Source object."""
        raise NotImplementedError()

    @property
    def position(self):
        """Target position (center) as SkyCoord."""
        # TODO: Consider adding default (with logging) here if
        #       self._position is None and self._offset is None
        #       But consider also how that might interact with parent position
        #       and offset frame from that.
        return self._position

    # TODO: add typing
    @position.setter
    def position(self, position):
        match position:
            case SkyCoord():
                self._position = position
            case (x_arcsec, y_arcsec) | {"x": x_arcsec, "y": y_arcsec}:
                x_arcsec <<= u.arcsec
                y_arcsec <<= u.arcsec
                self._position = SkyCoord(x_arcsec, y_arcsec)
            case _:
                raise TypeError("Unkown postition format.")

    @property
    def offset(self) -> dict:
        """Target offset from parent."""
        return self._offset

    @offset.setter
    def offset(self, offset: Mapping):
        if not isinstance(offset, Mapping):
            raise TypeError("Unkown offset format")

        # TODO: Consider adding warning when self._position is not None, because
        #       that would take precedence over any offset.

        self._offset = {
            "separation": offset["separation"],
            "position_angle": Angle(offset.get("position_angle", 0*u.deg)),
        }


class SpectrumTarget(Target):
    """Base class for Targets with separate spectrum (non-cube)."""

    @property
    def spectrum(self):
        """Target spectral information."""
        return self._spectrum

    # TODO: add typing
    @spectrum.setter
    def spectrum(self, spectrum):
        match spectrum:
            case SourceSpectrum():
                self._spectrum = spectrum
            case str() | SpectralType():
                self._spectrum = SpectralType(spectrum)
            case _:
                raise TypeError("Unkown spectrum format.")
