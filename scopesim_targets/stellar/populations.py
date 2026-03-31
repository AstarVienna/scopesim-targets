# -*- coding: utf-8 -*-
"""Stellar populations."""

import numpy as np
from scipy.stats import rv_continuous
from scipy.stats.sampling import NumericalInversePolynomial
from astropy import units as u

from ..spectral_classes import StellarParameters
from .imf import DEFAULT_IMFS


class Population:
    """Base class for stellar populations."""

    def __init__(self, n_stars: int):
        self._n_stars = n_stars


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

    def sample_imf(self) -> np.ndarray:
        urng = np.random.default_rng()
        rng = NumericalInversePolynomial(self.imf, center=0.1, random_state=urng)
        return rng.rvs(self._n_stars)
