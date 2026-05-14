# -*- coding: utf-8 -*-
"""Star cluster targets.

Current issues
--------------

The current `ZeroAgeCluster` class is limited to point sources (or rather stars
really), but that isn't reflected in its name. The parent class `Cluster`
technically isn't limited to point sources.

There is now some overlap in concept between `Cluster` and `StarField`, as both
deal with an internal list of point sources. This should be resolved in the
future. In the context of this, there's also the not-yet-implemented
`TargetCollection` mentioned in the example YAMLs.
"""

from typing import Any
from collections.abc import Mapping

from astropy import units as u
from astropy.table import Table

from scopesim import Source
from scopesim.source.source_fields import TableSourceField

from .typing_utils import POSITION_TYPE
from .target import Target
from .stellar.populations import Population, ZeroAgePopulation
from .stellar.morphology import Morphology


# This look very much like a dataclass now...
class Cluster(Target):
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
        return Source(field=TableSourceField(tbl, spectra=spectra))
