# -*- coding: utf-8 -*-
"""Importing this package has the side-effect of adding custom YAML tags."""

from .target import Target
from . import point_source
from . import extended_source

from .yaml_constructors import (
    register_qty,
    register_coord,
    register_target_constructor,
)

# Run YAML registrations
register_qty()
register_coord()

register_target_constructor(point_source.Star)
register_target_constructor(point_source.Binary)
register_target_constructor(point_source.Exoplanet)
register_target_constructor(point_source.PlanetarySystem)
register_target_constructor(point_source.StarField)

register_target_constructor(extended_source.Sersic)
register_target_constructor(extended_source.Disk)
