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

import numpy as np
from astropy import units as u
from astropy.table import Table

from scopesim import Source
from scopesim.source.source_fields import TableSourceField

from .typing_utils import POSITION_TYPE
from .target import Target
from .extinction import (
    parse_extinction,
    parse_extinction_distribution,
    resolve_extinction,
    ExtinctionDistribution,
)
from .stellar import populations, morphology


# Index of the extinction stream in the cluster's spawned seed tree
# (0 = population, 1 = morphology, 2 = extinction).
_N_SEED_STREAMS = 3
_EXT_STREAM = 2


# This look very much like a dataclass now...
class Cluster(Target):
    def __init__(
        self,
        position: POSITION_TYPE,
        population: populations.Population,
        morphology: morphology.Morphology,
        extinction=None,
        rng_seed: int | None = None,
    ) -> None:
        self.position = position
        self.population = population
        self.morphology = morphology
        self.rng_seed = rng_seed
        # A generated cluster's fluxes come from isochrone/absolute magnitudes
        # with the distance modulus already applied by the population, i.e. they
        # are unreddened apparent fluxes that extinction should dim. That is the
        # 'intrinsic' contract; 'observed' has no meaning here (no single anchor
        # band to normalise against), so the frame is fixed, not exposed.
        self.anchor = "intrinsic"
        if extinction is not None:
            self.extinction = extinction

    @property
    def extinction(self):
        """Line-of-sight dust screen(s), or ``None`` for no extinction.

        Accepts the schema forms on assignment (bare ``"2.3 mag"`` A_V sugar, a
        canonical mapping, or a list composed multiplicatively; ``[]`` is an
        explicit opt-out). Stored normalized as a tuple of
        :class:`~.extinction.ExtinctionScreen` (or a deferred
        :class:`~.extinction.FromMap` marker). This is the fourth cascading
        attribute in the spec, but is treated as a plain local attribute for now
        (no inheritance until the scene/role machinery exists).
        """
        return getattr(self, "_extinction", None)

    @extinction.setter
    def extinction(self, value):
        if value is None:
            self._extinction = None
        elif isinstance(value, Mapping) and "distribution" in value:
            # generative per-star form (valid only on sampling targets, e.g.
            # ZeroAgeCluster); non-generative targets reject it at resolve time.
            self._extinction = parse_extinction_distribution(value)
        else:
            self._extinction = parse_extinction(value)

    def _extinction_rng(self) -> "np.random.Generator":
        """RNG for the per-star A_V draw: the extinction child of the seed tree.

        Re-derived from ``rng_seed`` so it matches the child that
        ``ZeroAgeCluster`` handed to the pop/morph samplers -- reproducible and
        independent of them.
        """
        seed = getattr(self, "_ext_seed", None)
        if seed is None:
            seed = np.random.SeedSequence(self.rng_seed).spawn(_N_SEED_STREAMS)[_EXT_STREAM]
        return np.random.default_rng(seed)

    def _apply_extinction(self, tbl: Table) -> Table:
        """Add the per-row ``Av`` column (+ law/rv/anchor meta), if any.

        A distribution draws one A_V per star; a fixed screen broadcasts a
        single value (the N=1, zero-variance case). Absent extinction, the table
        is returned untouched -- byte-identical, no new column.
        """
        ext = self.extinction
        if ext is None:
            return tbl

        if isinstance(ext, ExtinctionDistribution):
            if self.rng_seed is None:
                raise ValueError(
                    "a per-star extinction distribution needs rng_seed for a "
                    "reproducible realization"
                )
            av = ext.sample(len(tbl), self._extinction_rng())
            law, rv = ext.law.value, ext.rv
        else:
            resolved = resolve_extinction(ext)  # V / E(B-V) screens; no bandpass
            if resolved is None:
                return tbl
            av_scalar, law_enum, rv = resolved
            av = av_scalar  # broadcast to every row
            law = law_enum.value

        tbl["Av"] = av
        tbl.meta["extinction_law"] = law
        tbl.meta["extinction_rv"] = rv
        tbl.meta["anchor"] = self.anchor.value
        return tbl


class ZeroAgeCluster(Cluster):
    def __init__(
        self,
        position: POSITION_TYPE,
        pop_class: populations.ZeroAgePopulation | str,
        pop_params: Mapping[str, Any],
        morph_class: morphology.Morphology | str,
        morph_params: Mapping[str, Any],
        extinction=None,
        rng_seed: int | None = None,
    ) -> None:
        # Required for YAML definitions, which provide only strings...
        if isinstance(pop_class, str):
            pop_class = getattr(populations, pop_class)
        if isinstance(morph_class, str):
            morph_class = getattr(morphology, morph_class)

        # One master seed -> independent, reproducible streams for population,
        # morphology and extinction (schema: rng_seed seeds all three).
        pop_seed, morph_seed, ext_seed = np.random.SeedSequence(
            rng_seed
        ).spawn(_N_SEED_STREAMS)

        super().__init__(
            position,
            pop_class(**pop_params, rng=np.random.default_rng(pop_seed)),
            morph_class(**morph_params, rng=np.random.default_rng(morph_seed)),
            extinction=extinction,
            rng_seed=rng_seed,
        )
        self._ext_seed = ext_seed

    def to_source(self):
        src_coldict, spectra = self.population.to_source_columns(self.position)
        src_coldict.update(self.morphology.to_source_columns(self.position))

        tbl = Table(
            data=src_coldict,
            names=["x", "y", "ref", "weight", "absmag", "mass"],
            units={"x": u.arcsec, "y": u.arcsec},
        )
        tbl = self._apply_extinction(tbl)
        return Source(field=TableSourceField(tbl, spectra=spectra))
