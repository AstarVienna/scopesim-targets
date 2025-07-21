# -*- coding: utf-8 -*-
"""Contains main ``Target`` class."""


class Target:
    """Main class in scopesim-targets."""

    def to_source(self):
        """Convert to ScopeSim Source object."""
        raise NotImplementedError()
