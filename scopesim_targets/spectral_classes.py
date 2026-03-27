# -*- coding: utf-8 -*-
"""Stellar parameters lookup table and related functionality.

.. todo:: Expand docstring.

"""

from typing import Any, NamedTuple
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field

from more_itertools import always_iterable
import numpy as np
from astropy import units as u
from astropy.table import Table, QTable, Row, join
from astropy.utils.masked import combine_masks

from astar_utils import SpectralType

from .data_utils import fetch_data_file


_temperature_type = u.get_physical_type("temperature")
TeffRange = NamedTuple("TeffRange", [
    ("min", u.Quantity[_temperature_type] | float),
    ("max", u.Quantity[_temperature_type] | float),
])


@dataclass(order=True, frozen=True, slots=True)
class SpectralClass:
    """Container for main (OBAFGKM) spectral class Teff range and color.

    .. todo:: Add source for colors (Wiki).

    .. todo:: Expand docstring, document creator methods.

    """

    name: str = field(compare=False)
    teff_range: TeffRange
    color: str = field(compare=False)

    colors = {
        "O": "#92B5FF",
        "B": "#A2C0FF",
        "A": "#D5E0FF",
        "F": "#F9F5FF",
        "G": "#FFEDE3",
        "K": "#FFDAB5",
        "M": "#FFB56C",
    }
    _fallback_color = "#EAEAF2"  # default seaborn plot background

    @property
    def midpoint(self):
        """Return mean of `teff_range`."""
        return sum(self.teff_range) / 2

    def is_in_range(self, teff_range):
        """Return True if overlap, False otherwise."""
        try:
            teff_range_overlap(self.teff_range, teff_range)
        except ValueError:
            return False
        return True

    @classmethod
    def from_parameters_table(cls, params_table, name: str):
        teffs = params_table.loc[f"{name}0.0":f"{name}9.9"]["teff"]
        teff_range = TeffRange(teffs.min().value, teffs.max().value)
        return cls(name, teff_range, cls.colors.get(name, cls._fallback_color))

    @classmethod
    def from_table_row(cls, row: Row):
        name = str(row["spectral_class"])
        teff_range = TeffRange(row["teff_min"].value, row["teff_max"].value)
        return cls(name, teff_range, cls.colors.get(name, cls._fallback_color))


def teff_range_overlap(*teff_ranges) -> TeffRange:
    """
    Return overlap between `teff_ranges`.

    Parameters
    ----------
    *teff_ranges : TeffRange
        Two or more TeffRange tuples.

    Raises
    ------
    ValueError
        Raised if no overlap is found, i.e. the ranges are disjoint.

    Returns
    -------
    overlap_range : TeffRange
        Overlapping TeffRange.

    """
    overlap_range = TeffRange(
        max(teff_range.min for teff_range in teff_ranges),
        min(teff_range.max for teff_range in teff_ranges),
    )
    if overlap_range.min >= overlap_range.max:
        raise ValueError("no overlap between Teff ranges")
    return overlap_range


class StellarParameters:
    """Wrapper for Mamajek table of spectral types and physical parameters [1]_.

    Acts as a static lookup table for the parameters listed below.

    Currently only uses a subset of columns from the original table, omitting
    the ones currently not needed in the ScopeSim ecosystem. Some columns have
    been renamed to match nomenclature elsewhere. The following columns are
    currently included:

    - spectral_type (primary key, ``astar_utils.SpectralType`` objects)
    - teff (effective temperature in K)
    - mass (stellar mass in solar masses)
    - radius (stellar radius in solar radii)
    - M_V (absolute V band magnitude)
    - M_J (absolute J band magnitude)
    - M_Ks (absolute Ks band magnitude)
    - U-B color
    - B-V color
    - V-Ks color
    - J-H color
    - H-Ks color

    .. todo:: Consider including the remaining columns if useful.

    .. todo:: Add interpolation functions for main columns.

    .. todo:: Add information about ``.loc[]``-style indexing with spectral
       types, including examples.

    Parameters
    ----------
    required_columns : str | Iterable[str] | None
        Column(s) that may not be NaN. Any rows containing a missing value in
        those columns are excluded from the table. Useful for e.g. clipping the
        table to contain only spectral classes with available IR photometry.

    Attributes
    ----------
    table : astropy.table.QTable
        Processed parameter lookup table.

    Notes
    -----
    Originally taken from:
    https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt,
    initially published in [1]_, partially also in [2]_.

    References
    ----------
    .. [1] Pecaut, Mark J. and Mamajek, Eric E., "Intrinsic Colors,
       Temperatures, and Bolometric Corrections of Pre-main-sequence Stars",
       The Astrophysical Journal Supplement, Volume 208, Issue 1,
       article id. 9, 22 pp., 2013; `2013ApJS..208....9P
       <https://ui.adsabs.harvard.edu/abs/2013ApJS..208....9P>`_

    .. [2] Pecaut, Mark J., Mamajek, Eric E., and Bubar, Eric J., "A Revised Age
       for Upper Scorpius and the Star Formation History among the F-type
       Members of the Scorpius-Centaurus OB Association", The Astrophysical
       Journal, Volume 746, Issue 2, article id. 154, 22 pp., 2012;
       `2012ApJ...746..154P
       <https://ui.adsabs.harvard.edu/abs/2012ApJ...746..154P>`_

    """

    _filename = "stellar/stellar_parameters.ecsv"

    def __init__(self, required_columns: str | Iterable[str] | None = None):
        tbl = self._load_stellar_parameters_table()
        if required_columns is not None:
            combined_mask = combine_masks(
                tbl[col].mask for col in always_iterable(required_columns)
            )
            tbl = tbl[~combined_mask]
        self.table = tbl

        self._closest_indices: dict[str, Any] = {}

    def _load_stellar_parameters_table(self) -> QTable:
        """
        Load and prepare Mamajek table.

        Can't directly load as QTable because of logarithmic units.

        Returns
        -------
        mamajek_redux : QTable
            Reduced Mamajek table.

        """
        mamajek_full = Table.read(fetch_data_file(self._filename))

        mamajek_redux = QTable(mamajek_full[[
            "spectral_type",
            "teff",
            "mass",
            "radius",

            "M_V",
            "M_J",
            "M_Ks",

            "U-B",
            "B-V",
            "V-Ks",
            "J-H",
            "H-Ks",
        ]])

        # Convert primary key column to SpectralType objects
        mamajek_redux["spectral_type"] = [
            SpectralType(spectype)
            for spectype in mamajek_redux["spectral_type"]
        ]
        mamajek_redux.add_index("spectral_type", unique=True)

        return mamajek_redux

    def group_spectral_classes(self) -> Iterator[SpectralClass]:
        """
        Generate SpectralClass objects from grouped parameters table.

        Min and max of Teff are calculated for each spectral class (OBAFGKMLTY)
        and then any "gaps" between classes are closed by the average between
        the previous min and next max.

        Yields
        ------
        spectral_class : SpectralClass
            Instances of SpectralClass created from each row.

        """
        tbl = self.table[["spectral_type", "teff"]]
        tbl["spectral_class"] = [
            SpectralType(spt.spectral_class) for spt in tbl["spectral_type"]
        ]
        by_spec_cls = tbl.group_by("spectral_class")
        minmax = join(
            by_spec_cls["spectral_class", "teff"].groups.aggregate(np.min),
            by_spec_cls["spectral_class", "teff"].groups.aggregate(np.max),
            keys="spectral_class",
            table_names=["min", "max"],
        )
        midpoints = (minmax["teff_max"][1:] + minmax["teff_min"][:-1]) / 2.0
        minmax["teff_max"][1:] = midpoints
        minmax["teff_min"][:-1] = midpoints
        for row in minmax:
            yield SpectralClass.from_table_row(row)
