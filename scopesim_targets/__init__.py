# -*- coding: utf-8 -*-
"""Importing this package has the side-effect of adding custom YAML tags."""

from pathlib import Path

# Needs to be before the other imports to avoid circular issues.
PKG_DIR = Path(__file__).parent
DATA_DIR = PKG_DIR.parent / "data"


from .target import Target
from . import point_source
from . import extended_source
from . import cluster

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

register_target_constructor(cluster.ZeroAgeCluster)
