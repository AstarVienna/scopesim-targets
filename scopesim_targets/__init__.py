# -*- coding: utf-8 -*-
"""Importing this package has the side-effect of adding custom YAML tags."""

from .target import Target

from .yaml_constructors import register_qty, register_coord


# Run YAML registrations
register_qty()
register_coord()
