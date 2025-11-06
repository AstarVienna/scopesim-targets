# -*- coding: utf-8 -*-
"""Contains main ``Target`` class."""

from abc import ABCMeta, abstractmethod


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

    @position.setter
    def position(self, position):
        self._position = position
