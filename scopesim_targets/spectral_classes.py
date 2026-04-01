# -*- coding: utf-8 -*-
"""Stellar parameters lookup table and related functionality.

.. todo:: Expand docstring.

"""

from typing import Any, NamedTuple
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field

from more_itertools import always_iterable
import numpy as np
from numpy.lib.recfunctions import structured_to_unstructured
from scipy.spatial import KDTree
from scipy.interpolate import PchipInterpolator, CubicSpline
from astropy import units as u
from astropy.table import Table, QTable, Row, join
from astropy.utils.masked import Masked, combine_masks

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

    .. todo:: Add information about ``.loc[]``-style indexing with spectral
       types, including examples.

    .. todo:: Document interpolation and lookup methods with examples.

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
        self._col_units = {
            colname: column.unit
            for colname, column in self.table.columns.items()
        }

        self._closest_indices: dict[str, Any] = {}
        self._lookup_tree: KDTree | None = None

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

    @staticmethod
    def _calc_halfways_points(array: np.array) -> np.ndarray:
        return array[:-1] / 2.0 + array[1:] / 2.0

    def _sort_col_idx(self, colname: str) -> tuple[np.ndarray, np.ndarray]:
        # Implementation adapted from (deprecated) scipy.interpolate.interp1d
        sorted_indices = self.table[colname].argsort(kind="mergesort")
        sorted_column = self.table[colname][sorted_indices]
        return sorted_indices, sorted_column

    def _setup_closest(self, colname: str) -> tuple[np.ndarray, np.ndarray]:
        sorted_indices, sorted_column = self._sort_col_idx(colname)
        return sorted_indices, self._calc_halfways_points(sorted_column)

    def _get_closest_indices(self, colname: str):
        if (closest_indices := self._closest_indices.get(colname)) is not None:
            return closest_indices
        # Only actually run setup method if not already set
        closest_indices = self._closest_indices.setdefault(
            colname, self._setup_closest(colname)
        )
        return closest_indices

    def _search_halfways_points(
        self,
        search_value: np.ndarray,
        colname: str,
    ) -> np.ndarray:
        _, halfway_points = self._get_closest_indices(colname)
        # Implementation adapted from (deprecated) scipy.interpolate.interp1d
        closest_indices = halfway_points.searchsorted(search_value, side="left")
        # Clip to avoid matching masked values
        if hasattr(halfway_points, "mask"):  # e.g. mass
            max_good_index = len(halfway_points) - halfway_points.mask.sum()
        else:  # e.g. teff
            max_good_index = len(halfway_points)
        return closest_indices.clip(0, max_good_index).astype(np.intp)

    @u.quantity_input
    def closest_mass(self, mass: u.Quantity[u.solMass]) -> QTable | Row:
        """
        Lookup the row(s) closest to a given stellar mass(es).

        Parameters
        ----------
        mass : u.Quantity[u.solMass]
            Stellar mass(es) to look up. Can be scalar or array, but must be
            quantity convertible to mass.

        Returns
        -------
        closest_rows : QTable | Row
            Resulting row(s). If `mass` was provided as an array, the result
            will be a QTable, subset of the original table. If `mass` is scalar,
            the result will be a single Row object.

        """
        sorted_indices, _ = self._get_closest_indices("mass")
        closest_sorted_indices = self._search_halfways_points(mass, "mass")
        actual_indices = sorted_indices[closest_sorted_indices]
        return self.table[actual_indices]

    @u.quantity_input
    def closest_teff(self, teff: u.Quantity[u.K]) -> QTable | Row:
        """
        Lookup the row(s) closest to a given effective temperature(s).

        Parameters
        ----------
        teff : u.Quantity[u.K]
            Stellar effective temperature(s) to look up. Can be scalar or array,
            but must be quantity convertible to temperature.

        Returns
        -------
        closest_rows : QTable | Row
            Resulting row(s). If `teff` was provided as an array, the result
            will be a QTable, subset of the original table. If `teff` is scalar,
            the result will be a single Row object.

        """
        sorted_indices, _ = self._get_closest_indices("teff")
        closest_sorted_indices = self._search_halfways_points(teff, "teff")
        actual_indices = sorted_indices[closest_sorted_indices]
        return self.table[actual_indices]

    def _make_lookup_tree(self) -> KDTree:
        # HACK: Using the existing index should ensure things are already
        #       sorted, but this is untested because we initially sort this
        #       column anyway. Whatever the case, there must be a better way to
        #       get the values out than this double indexing...
        spectral_type_array = [
            spectype.to_array() for spectype in
            self.table.indices["spectral_type"].data.data["spectral_type"]
        ]
        return KDTree(spectral_type_array)

    def closest_spectral_type(self, spectral_type: SpectralType) -> QTable | Row:
        # TODO: docstring
        if self._lookup_tree is None:
            self._lookup_tree = self._make_lookup_tree()

        # Deal with both scalar and array inputs
        spectral_type_array = [
            SpectralType(spectype).to_array() for spectype in
            always_iterable(spectral_type)
        ]
        closest_spectypes = [
            SpectralType.from_array(spectype) for spectype in
            np.atleast_2d(self._lookup_tree.data[
                self._lookup_tree.query(spectral_type_array)[1]
            ])
        ]
        return self.table.loc[closest_spectypes]

    def _get_remaining_colnames(self, colname: str) -> list[str]:
        """Return all column names other than `colname` and "spectral_type"."""
        exclude = {"spectral_type", colname}
        return [col for col in self.table.colnames if col not in exclude]

    def interpolate(self, colname: str, values: u.Quantity, extrapolate_phot: bool = False) -> QTable:
        # TODO: store interpolator
        # TODO: optionally supply output colnames (set intersection stuff)
        if self.table[colname].unit != values.unit:
            raise ValueError("units in values must match column")

        columns = self._get_remaining_colnames(colname)
        sorted_indices, sorted_column = self._sort_col_idx(colname)

        if hasattr(sorted_column, "mask"):
            sorted_indices = sorted_indices[~sorted_column.mask]
            sorted_column = sorted_column[~sorted_column.mask]

        # TODO: find a more efficient way to do this
        mask = structured_to_unstructured(
            self.table.mask[columns][sorted_indices].as_array()
        )
        expanded = Masked(np.repeat(sorted_column[:, None], 10, 1), mask=mask)
        mins = expanded.min(axis=0)
        maxs = expanded.max(axis=0)
        output_mask = (
            (np.repeat(values[:, None], 10, 1) < mins) |
            (np.repeat(values[:, None], 10, 1) > maxs)
        )

        splines = []
        for colname in columns:
            col = self.table[colname][sorted_indices]
            colmask = col.mask if hasattr(col, "mask") else np.zeros_like(col.value).astype(bool)
            if colname in ("mass", "teff", "radius"):
                spline = PchipInterpolator(sorted_column[~colmask], col[~colmask], extrapolate=False)
            else:
                spline = CubicSpline(sorted_column[~colmask], col[~colmask], extrapolate=True)
                # see https://docs.scipy.org/doc/scipy/tutorial/interpolate/extrapolation_examples.html#cubicspline-extend-the-boundary-conditions
                # TODO: Find a better way to extrapolate that works with the (otherwise much better) PchipInterpolator, get rid of that function.
                _add_boundary_knots(spline)
            splines.append(spline)

        # output_data = spline(mass).round(3)  # HACK: use table format instead?
        # Output has to be Quantity BEFORE masking, otherwise QTable doesn't
        # return quantities upon [colname]...
        # TODO: See if this can be done more efficiently and/or robust, as the
        #       current implementation relies on correct column order...
        masked_output = [
            Masked(
                splines[i](values).round(3)*self._col_units.get(col),
                mask=(
                    output_mask[:, i]
                    if col in ("mass", "teff", "radius")
                    or not extrapolate_phot
                    else None
                ),
            ) for i, col in enumerate(columns)
        ]
        result = QTable(
            names=columns,
            data=masked_output,
        )
        return result


def _add_boundary_knots(spline):
    """
    Add knots infinitesimally to the left and right.

    Additional intervals are added to have zero 2nd and 3rd derivatives,
    and to maintain the first derivative from whatever boundary condition
    was selected. The spline is modified in place.

    Taken directly from https://docs.scipy.org/doc/scipy/tutorial/interpolate/extrapolation_examples.html#cubicspline-extend-the-boundary-conditions
    """
    # determine the slope at the left edge
    leftx = spline.x[0]
    lefty = spline(leftx)
    leftslope = spline(leftx, nu=1)

    # add a new breakpoint just to the left and use the
    # known slope to construct the PPoly coefficients.
    leftxnext = np.nextafter(leftx, leftx - 1)
    leftynext = lefty + leftslope*(leftxnext - leftx)
    leftcoeffs = np.array([0, 0, leftslope, leftynext])
    spline.extend(leftcoeffs[..., None], np.r_[leftxnext])

    # repeat with additional knots to the right
    rightx = spline.x[-1]
    righty = spline(rightx)
    rightslope = spline(rightx, nu=1)
    rightxnext = np.nextafter(rightx, rightx + 1)
    rightynext = righty + rightslope * (rightxnext - rightx)
    rightcoeffs = np.array([0, 0, rightslope, rightynext])
    spline.extend(rightcoeffs[..., None], np.r_[rightxnext])
