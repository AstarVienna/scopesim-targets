# -*- coding: utf-8 -*-
"""Contains main ``Target`` class."""

from abc import ABCMeta, abstractmethod

from astropy import units as u
from astropy.coordinates import SkyCoord


class Target(metaclass=ABCMeta):
    """Main class in scopesim-targets."""

    @abstractmethod
    def to_source(self):
        """Convert to ScopeSim Source object."""
        raise NotImplementedError()

    @property
    def position(self):
        """Target position (center) as SkyCoord."""
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
