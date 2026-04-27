# -*- coding: utf-8 -*-
"""Stellar populations."""

import numpy as np
from scipy.stats import rv_continuous
from scipy.stats.sampling import NumericalInversePolynomial
from astropy import units as u
from matplotlib import axes

from astar_utils import SpectralType
from spextra import SpecLibrary, Spextrum

from ..spectral_classes import StellarParameters
from ..plot_utils import figure_factory
from .imf import DEFAULT_IMFS

# Split at 1.07 Msol, F/G border
DEFAULT_LIBRARY_LOW_MASS = SpecLibrary("irtf")
DEFAULT_LIBRARY_HIGH_MASS = SpecLibrary("kurucz")
HIGH_LOW_MASS_LIMIT = 1.07*u.solMass


class Population:
    """Base class for stellar populations."""

    def __init__(self, n_stars: int):
        self._n_stars = n_stars
        # TODO: Consider using a singelton-ish thing here
        self._stellar_params = StellarParameters()  # Default lookup table


class ZeroAgePopulation(Population):
    """Stellar population with no consideration about stellar evolution."""


class IMFPopulation(ZeroAgePopulation):
    """Zero-age stellar population sampled from an IMF interpreted as a PDF."""

    imf: rv_continuous = DEFAULT_IMFS["kroupa02"]

    def __init__(self, n_stars: int, imf: rv_continuous | None = None):
        super().__init__(n_stars)
        if imf is not None:
            self.imf = imf

    @classmethod
    @u.quantity_input
    def from_total_mass(cls, total_mass: u.Quantity[u.solMass], imf: rv_continuous | None = None):
        """Generate population for total (cluster) mass.

        For non-continous distributions (e.g. broken powerlaw) this can deviate
        by about 3 %, otherwise within 1 % (all empirical).
        """
        imf = imf or cls.imf  # default if None
        n_stars = int(total_mass.to_value(u.solMass) / imf.expect())
        return cls(n_stars, imf)

    def sample_imf(self) -> u.Quantity[u.solMass]:
        urng = np.random.default_rng()
        rng = NumericalInversePolynomial(self.imf, center=0.1, random_state=urng)
        return rng.rvs(self._n_stars).round(3) * u.solMass

    def _masses_to_brightness(self, masses, absmag_col: str):
        absmags = (
            self._stellar_params
            .interpolate("mass", masses, extrapolate_phot=True)[absmag_col]
            .round(2)
        )
        if hasattr(absmags, "mask") and not absmags.mask.any():
            # ditch any remaining empty masks
            return absmags.unmasked
        return absmags

    def _masses_to_spectra(self, masses):
        # HACK: This is to crop the stellar parameters table to only include
        #       spectral types found in a given spextra library. There are
        #       multiple competing ideas on how to properly implement this, but
        #       I'm not yet sure which is best so I didn't want to commit to
        #       implementing any of those for now. This hack is easiest to
        #       remove again once a proper solution exists...
        # Note: This hack requires that stp is freshly created, otherwise there
        #       might be interpolators and lookup tables already setup which are
        #       not limited to the overwritten table...
        # Filter to those found in both the library and the lookup table,
        # sorted() turns it into the required list.
        # I tried to do some fancy set intersection here, which did work, but
        # not for the Brown Dwarfs, because those are listed as e.g. L5V in the
        # Mamajek table, but as e.g. L5 in the IRTF library, and while e.g. L5
        # is considered equal to L5V and works in table indexing, it does not
        # work in the comparison in a set. But this here is also fine.
        stp_low_mass = StellarParameters("M_J")  # need to override default
        common_spectypes = set()
        for spectype in DEFAULT_LIBRARY_LOW_MASS:
            try:
                spectype = SpectralType(spectype)
                if spectype in stp_low_mass.table["spectral_type"]:
                    common_spectypes.add(spectype)
            except ValueError:
                # Catch and ignore non-standard names in library
                continue
        stp_low_mass.table = stp_low_mass.table.loc[sorted(common_spectypes)]

        stp_high_mass = StellarParameters("M_J")  # need to override default
        common_spectypes = set()
        for spectype in DEFAULT_LIBRARY_HIGH_MASS:
            try:
                spectype = SpectralType(spectype)
                if spectype in stp_high_mass.table["spectral_type"]:
                    common_spectypes.add(spectype)
            except ValueError:
                # Catch and ignore non-standard names in library
                continue
        stp_high_mass.table = stp_high_mass.table.loc[sorted(common_spectypes)]

        stp_low_mass.table = stp_low_mass.table.loc["G0":]
        stp_high_mass.table = stp_high_mass.table.loc[:"F9.9"]

        spectypes = np.where(
            masses < HIGH_LOW_MASS_LIMIT,
            stp_low_mass.closest_mass(masses)["spectral_type"],
            stp_high_mass.closest_mass(masses)["spectral_type"],
        )

        # TODO: forcing M_J now cuts us off at B0V on the high end, which isn't
        #       brilliant, although fine for now
        spectra = {}
        for row in stp_high_mass.table:
            spectype = row["spectral_type"]
            libname = DEFAULT_LIBRARY_HIGH_MASS.name
            spec = Spextrum(f"{libname}/{str(spectype).lower()}")
            absmag = row["M_J"]
            if absmag.mask:
                continue
            spectra[spectype] = spec.scale_to_magnitude(absmag, "J")
        for row in stp_low_mass.table:
            spectype = row["spectral_type"]
            libname = DEFAULT_LIBRARY_LOW_MASS.name
            specname = str(spectype)
            # LTY have no "V" in that library -.-
            if specname.startswith(("L", "T", "Y")):
                specname = specname.removesuffix("V")
            spec = Spextrum(f"{libname}/{specname}")
            absmag = row["M_J"]
            if absmag.mask:
                continue
            spectra[spectype] = spec.scale_to_magnitude(absmag, "J")

        # HACK: While specref still only works with ints in ScopeSim, this seems
        #       to be the quickest way to turn spectypes into IDs...
        specref_lookup = {spec: i for i, spec in enumerate(spectra)}
        specref = np.array([specref_lookup[spec] for spec in spectypes])

        spectra = dict(enumerate(spectra.values()))
        return specref, spectra

    def to_source_columns(self, parent_position, absmag_col: str = "M_J"):
        masses = self.sample_imf()
        specref, spectra = self._masses_to_spectra(masses)
        absmags = self._masses_to_brightness(masses, absmag_col)

        specmags_lookup = {
            ref: self._stellar_params.table.loc[spec.basename][absmag_col].value
            for ref, spec in spectra.items()
        }
        specmags = np.array([specmags_lookup[ref] for ref in specref]) << u.mag

        distmod = parent_position.distance.distmod
        delta_mag = distmod + absmags - specmags
        weights = (delta_mag.value * delta_mag.unit()).physical

        return {"ref": specref, "weight": weights}, spectra

    def plot(
        self,
        samples: np.ndarray | None = None,  # array-like?
        mass_range: tuple[float, float] = (0., np.inf),
        ax: axes.Axes | None = None,
        label: str = "IMF",
    ) -> axes.Axes:
        min_mass, max_mass = sorted(mass_range)
        min_mass = max(min_mass, self.imf.a)
        max_mass = min(max_mass, self.imf.b)

        x = np.geomspace(min_mass + 1e-5, max_mass)
        bins = np.geomspace(min_mass, max_mass, 30)

        if ax is None:
            _, ax = figure_factory()

        if samples is None:
            samples = self.sample_imf()

        ax.loglog(x, self.imf.pdf(x), label=label)
        ax.hist(samples.value, bins=bins, density=True, alpha=.8)
        ax.axvline(samples.mean().value, ls="--")

        if (breakpoints := [
            value for key, value in self.imf.kwds.items() if key.startswith("m")
        ]):
            for value in breakpoints:
                ax.axvline(value, ls=":", c="k")
        ax.legend()

        return ax
