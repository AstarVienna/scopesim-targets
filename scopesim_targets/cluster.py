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

Dust
----

A cluster carries the two dust attributes of the spec (``defining_extinction.md``):

* ``extinction`` -- the *sightline*: one foreground screen, identical for every
  member. It is resolved here and applied to the shared spectral templates
  **in place**, so the deduplicated spectrum set stays collapsed and no table
  column is needed. The cluster frame is fixed to ``intrinsic`` (isochrone
  magnitudes with the distance modulus already applied), so the screen dims as
  well as reddens.
* ``attenuation`` -- the *source*: the natal cloud, one column per star (a fixed
  column broadcasts; a ``distribution`` is sampled under the cluster's seed).
  It cannot be folded into the shared templates and is therefore *emitted* as
  the per-row ``attenuation`` column + model meta for the simulation engine
  (see :func:`~.attenuation.apply_attenuation_to_table` for the contract).
"""

from typing import Any
from collections.abc import Mapping

import numpy as np
from astropy import units as u
from astropy.table import Table

from scopesim import Source
from scopesim.source.source_fields import TableSourceField

from .typing_utils import POSITION_TYPE
from .seeding import resolve_seed, stream, new_rng
from .target import Target, _passband
from .extinction import parse_extinction, resolve_extinction, redden
from .attenuation import parse_attenuation, apply_attenuation_to_table
from .stellar import populations, morphology


# This look very much like a dataclass now...
class Cluster(Target):
    def __init__(
        self,
        position: POSITION_TYPE,
        population: populations.Population,
        morphology: morphology.Morphology,
        extinction=None,
        attenuation=None,
        rng_seed: int | None = None,
    ) -> None:
        self.position = position
        self.population = population
        self.morphology = morphology
        # Master seed actually used (fresh entropy for None), recorded so every
        # realization -- including the sampled attenuation -- is reproducible
        # from the target. Idempotent for an int (ZeroAgeCluster passes its
        # already-resolved seed).
        self.rng_seed = resolve_seed(rng_seed)
        # A generated cluster's fluxes come from isochrone/absolute magnitudes
        # with the distance modulus already applied by the population, i.e. they
        # are unreddened apparent fluxes that dust should dim. That is the
        # 'intrinsic' contract; 'observed' has no meaning here (no single anchor
        # band to normalise against), so the frame is fixed, not exposed.
        self.anchor = "intrinsic"
        if extinction is not None:
            self.extinction = extinction
        if attenuation is not None:
            self.attenuation = attenuation

    # -- extinction: the foreground screen ---------------------------------------

    @property
    def extinction(self):
        """Foreground dust screen, or ``None`` for none.

        Accepts the schema forms on assignment (bare ``"2.3 mag"`` A_V sugar or
        the canonical mapping; ``"0 mag"`` is the explicit opt-out). Stored as
        an :class:`~.extinction.ExtinctionScreen` (or a deferred
        :class:`~.extinction.FromMap` marker). Cascades from the scene in the
        spec; treated as a plain local attribute until the scene/role machinery
        exists.
        """
        return getattr(self, "_extinction", None)

    @extinction.setter
    def extinction(self, value):
        self._extinction = None if value is None else parse_extinction(value)

    # -- attenuation: the natal cloud --------------------------------------------

    @property
    def attenuation(self):
        """Source-side dust column per star, or ``None`` for none.

        Accepts bare ``"8.5 mag"`` sugar (screen A_V) or the mapping form incl.
        the sampled ``distribution`` variant; stored as an
        :class:`~.attenuation.Attenuation`. Never inherited.
        """
        return getattr(self, "_attenuation", None)

    @attenuation.setter
    def attenuation(self, value):
        self._attenuation = None if value is None else parse_attenuation(value)

    # -- helpers -----------------------------------------------------------------

    @staticmethod
    def _band_wavelength_lookup(band: str) -> u.Quantity:
        """Effective wavelength of a named band, for non-V dust columns.

        Same lookup as ``SpectrumTarget._band_effective_wavelength`` (no shared
        base class yet; unify once ``Cluster`` joins that hierarchy).
        """
        return _passband(band).avgwave()

    def _attenuation_rng(self) -> np.random.Generator:
        """Fresh generator for the per-star column draw.

        The ``attenuation`` stream of the master seed: reproducible, the same
        on every call, and independent of the population/morphology draws.
        """
        return new_rng(stream(self.rng_seed, "attenuation"))

    def _redden_spectra(self, spectra):
        """Apply the foreground screen to the shared templates, in place.

        ``spectra`` is whatever the population hands to ``TableSourceField``
        (a mapping ``ref -> SourceSpectrum`` or a sequence); the structure is
        preserved so the table's ``ref`` column stays valid. No screen (or an
        explicit ``0 mag``) returns the input untouched.
        """
        resolved = resolve_extinction(
            self.extinction, band_wavelength_lookup=self._band_wavelength_lookup
        )
        if resolved is None or resolved[0] == 0:
            return spectra
        av, law, rv = resolved
        if isinstance(spectra, Mapping):
            return {key: redden(spec, av, law, rv) for key, spec in spectra.items()}
        return type(spectra)(redden(spec, av, law, rv) for spec in spectra)

    def _apply_attenuation(self, tbl: Table) -> Table:
        """Add the per-row ``attenuation`` column (+ model meta), if any.

        A distribution draws one column per star; a fixed column broadcasts (the
        N=1, zero-variance case). Absent attenuation, the table is returned
        untouched -- byte-identical, no new column. The ``anchor`` meta is
        always the cluster's fixed ``intrinsic``.
        """
        att = self.attenuation
        if att is None:
            return tbl
        # rng_seed is always resolved (fresh entropy is recorded), so a sampled
        # column is reproducible even without an explicit seed.
        rng = self._attenuation_rng() if att.is_sampled else None
        tbl = apply_attenuation_to_table(
            tbl, att, rng=rng, band_wavelength_lookup=self._band_wavelength_lookup
        )
        # Target.anchor coerces to the AnchorFrame enum; the table wants the name.
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
        attenuation=None,
        rng_seed: int | None = None,
    ) -> None:
        # Required for YAML definitions, which provide only strings...
        if isinstance(pop_class, str):
            pop_class = getattr(populations, pop_class)
        if isinstance(morph_class, str):
            morph_class = getattr(morphology, morph_class)

        # Resolve the master seed once (fresh entropy for None); population,
        # morphology and attenuation each draw from their own named stream of
        # it, see `seeding.STREAMS`. Passed on resolved, so Cluster records
        # the same value.
        master_seed = resolve_seed(rng_seed)

        super().__init__(
            position,
            pop_class(**pop_params, rng_seed=stream(master_seed, "population")),
            morph_class(
                **morph_params, rng_seed=stream(master_seed, "morphology")
            ),
            extinction=extinction,
            attenuation=attenuation,
            rng_seed=master_seed,
        )

    def to_source(self, optical_train=None):
        src_coldict, spectra = self.population.to_source_columns(self.position)
        src_coldict.update(self.morphology.to_source_columns(self.position))

        tbl = Table(
            data=src_coldict,
            names=["x", "y", "ref", "weight", "absmag", "mass"],
            units={"x": u.arcsec, "y": u.arcsec},
        )
        spectra = self._redden_spectra(spectra)   # foreground: shared, in place
        tbl = self._apply_attenuation(tbl)        # natal cloud: per row, for the engine
        return Source(field=TableSourceField(tbl, spectra=spectra))
