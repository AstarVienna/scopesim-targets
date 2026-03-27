# -*- coding: utf-8 -*-
"""TBA."""

from collections.abc import Iterable

from more_itertools import always_iterable
from astropy.table import Table, QTable
from astropy.utils.masked import combine_masks

from astar_utils import SpectralType

from .data_utils import fetch_data_file


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

        self._closest_indices = {}

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
