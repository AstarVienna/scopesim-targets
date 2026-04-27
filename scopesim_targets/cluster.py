# -*- coding: utf-8 -*-
"""TBA."""

from typing import Any
from collections.abc import Mapping

from astropy import units as u
from astropy.table import Table

from scopesim import Source
from scopesim.source.source_fields import TableSourceField

from .typing_utils import POSITION_TYPE
from .stellar.populations import Population, ZeroAgePopulation
from .stellar.morphology import Morphology


class Cluster:
    def __init__(
        self,
        position: POSITION_TYPE,
        population: Population,
        morphology: Morphology,
    ) -> None:
        self.position = position
        self.population = population
        self.morphology = morphology


class ZeroAgeCluster(Cluster):
    def __init__(
        self,
        position: POSITION_TYPE,
        pop_class: ZeroAgePopulation,
        pop_params: Mapping[str, Any],
        morph_class: Morphology,
        morph_params: Mapping[str, Any],
    ) -> None:
        super().__init__(
            position,
            pop_class(**pop_params),
            morph_class(**morph_params),
        )

    def to_source(self):
        src_coldict, spectra = self.population.to_source_columns(self.position)
        src_coldict.update(self.morphology.to_source_columns(self.position))

        tbl = Table(
            data=src_coldict,
            names=["x", "y", "ref", "weight"],
            units={"x": u.arcsec, "y": u.arcsec},
        )
        src = Source(field=TableSourceField(tbl, spectra=spectra))
        return src

# StarField.from_??(Population, Morphology) ???

# but where to put parameters such that variable??
