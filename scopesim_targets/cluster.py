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
from .stellar import populations, morphology


# This look very much like a dataclass now...
class Cluster(Target):
    def __init__(
        self,
        position: POSITION_TYPE,
        population: populations.Population,
        morphology: morphology.Morphology,
    ) -> None:
        self.position = position
        self.population = population
        self.morphology = morphology


class ZeroAgeCluster(Cluster):
    def __init__(
        self,
        position: POSITION_TYPE,
        pop_class: populations.ZeroAgePopulation | str,
        pop_params: Mapping[str, Any],
        morph_class: morphology.Morphology | str,
        morph_params: Mapping[str, Any],
    ) -> None:
        # Required for YAML definitions, which provide only strings...
        if isinstance(pop_class, str):
            pop_class = getattr(populations, pop_class)
        if isinstance(morph_class, str):
            morph_class = getattr(morphology, morph_class)

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
            names=["x", "y", "ref", "weight", "absmag", "mass"],
            units={"x": u.arcsec, "y": u.arcsec},
        )
        return Source(field=TableSourceField(tbl, spectra=spectra))
