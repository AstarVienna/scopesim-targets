# -*- coding: utf-8 -*-
"""Custom YAML class constructors (loading) and representations (dumping).

Currently contains representers and constructors for `astropy.units.Quantity`
and `astropy.coordinates.SkyCoord`. Both of them have YAML representers and
constructors defined by astropy, but those are a bit too complex and not easily
readable, so we're defining simplified versions here which fit the needs of this
package. The original YAML tags can still be used if desired, by using the
astropy YAML loader instead of the default one.

The way this works is that the `register_foo` functions are called once in the
top-level `__init__.py`, which makes the representers, constructors and implicit
resolvers available everywhere for the standard YAML loaders and dumpers.
"""

import re

import yaml
import astropy.units as u
from astropy.coordinates import SkyCoord


def register_qty() -> None:
    """Register representer, constructor and implicit resolver for Quantity."""
    # Regex for implicit resolver
    quantity_pattern = re.compile(
        r"""^                    # start of string
            \s*                  # optional leading whitespace
            [-+]?                # optional leading sign
            (\d+(\.\d*)?|\.\d+)  # number (int or float)
            (e[-+]\d+)?          # optional exponent with mandatory +/-
            \s                   # exactly one space between number and unit
            ([a-zA-Z]+(-?\d+)?   # first unit, possibly with negative exponent
            (\s[a-zA-Z]+(-?\d+)?)*)  # optional additional units
            \s*                  # optional trailing whitespace
        $""",
        re.VERBOSE,
    )

    def qty_representer(dumper, data):
        """Convert ``Quantity`` to scalar string with custom tag."""
        return dumper.represent_scalar("!qty", data.to_string())

    def qty_constructor(loader, node):
        """Convert scalar string to ``Quantity``."""
        return u.Quantity(loader.construct_scalar(node))

    yaml.add_representer(u.Quantity, qty_representer)
    yaml.add_constructor("!qty", qty_constructor)
    yaml.add_implicit_resolver("!qty", quantity_pattern)


def register_coord() -> None:
    """Register simplified representer and constructor for SkyCoord."""
    def coord_representer(dumper, data):
        """Convert ``SkyCoord`` to mapping with custom tag."""
        # TODO: Consider adding distance to the representation
        representation = {"ra": f"{data.ra!s}", "dec": f"{data.dec!s}"}
        return dumper.represent_mapping("!Coord", representation)

    def coord_constructor(loader, node):
        """Convert mapping node to ``SkyCoord``."""
        return SkyCoord(**loader.construct_mapping(node))

    yaml.add_representer(SkyCoord, coord_representer)
    yaml.add_constructor("!Coord", coord_constructor)


def register_target_constructor(target_cls) -> None:
    """Register mapping constructor for `target_cls`."""
    def target_constructor(loader, node):
        return target_cls(**loader.construct_mapping(node, deep=True))
    yaml.add_constructor(f"!{target_cls.__name__}", target_constructor)
