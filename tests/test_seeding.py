# -*- coding: utf-8 -*-
"""Unit tests for seeding.py and the seeded samplers that use it.

All offline: seeding is numpy only, King sampling is scipy/astropy only and
IMF sampling needs no spectra. King-based tests are marked slow: building the
King sampler takes a few seconds per instance.
"""

import pytest
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

from scopesim_targets.seeding import (
    STREAMS,
    resolve_seed,
    stream,
    as_seed_sequence,
    new_rng,
)
from scopesim_targets.stellar.morphology import KingProfileMorphology
from scopesim_targets.stellar.populations import IMFPopulation
from scopesim_targets.cluster import ZeroAgeCluster


def _state(seed_sequence):
    return seed_sequence.generate_state(4)


class TestResolveSeed:
    @pytest.mark.parametrize("seed", (0, 42, 2**70))
    def test_int_passes_through(self, seed):
        assert resolve_seed(seed) == seed

    def test_none_draws_recordable_entropy(self):
        first, second = resolve_seed(None), resolve_seed(None)
        assert isinstance(first, int)
        assert first != second
        # the recorded value reproduces the draw
        np.testing.assert_array_equal(
            new_rng(as_seed_sequence(first)).random(3),
            new_rng(as_seed_sequence(first)).random(3),
        )

    @pytest.mark.parametrize(
        ("seed", "exc"),
        ((-1, ValueError), (True, TypeError), (1.5, TypeError), ("42", TypeError)),
    )
    def test_invalid(self, seed, exc):
        with pytest.raises(exc):
            resolve_seed(seed)


class TestStreams:
    def test_registry_is_pinned(self):
        # APPEND-ONLY contract: existing indices must never change.
        assert STREAMS["population"] == 0
        assert STREAMS["morphology"] == 1
        assert STREAMS["extinction"] == 2
        assert STREAMS["count"] == 3
        assert STREAMS["magnitudes"] == 4
        assert STREAMS["spectra"] == 5
        assert len(set(STREAMS.values())) == len(STREAMS)

    def test_equals_spawned_children(self):
        # Same streams as SeedSequence(seed).spawn(n) -- independent of n.
        children = np.random.SeedSequence(42).spawn(len(STREAMS) + 3)
        for name, index in STREAMS.items():
            np.testing.assert_array_equal(
                _state(stream(42, name)), _state(children[index])
            )

    def test_streams_are_distinct(self):
        states = [_state(stream(42, name)) for name in STREAMS]
        assert len({tuple(state) for state in states}) == len(STREAMS)

    def test_unknown_stream(self):
        with pytest.raises(ValueError, match="unknown stream"):
            stream(42, "bogus")

    def test_as_seed_sequence_keeps_given_sequence(self):
        child = stream(42, "morphology")
        assert as_seed_sequence(child) is child

    def test_new_rng_does_not_advance(self):
        seed_sequence = stream(42, "population")
        np.testing.assert_array_equal(
            new_rng(seed_sequence).random(3), new_rng(seed_sequence).random(3)
        )


@pytest.fixture
def king_kwargs():
    return {"n_stars": 20, "r_core": 1 * u.pc, "r_tide": 10 * u.pc}


@pytest.fixture
def parent_position():
    return SkyCoord(0 * u.deg, 0 * u.deg, 1 * u.kpc)


@pytest.mark.slow
class TestSeededMorphology:
    def test_same_seed_same_draw(self, king_kwargs, parent_position):
        a = KingProfileMorphology(rng_seed=9, **king_kwargs).sample(parent_position)
        b = KingProfileMorphology(rng_seed=9, **king_kwargs).sample(parent_position)
        np.testing.assert_array_equal(a, b)

    def test_different_seed_differs(self, king_kwargs, parent_position):
        a = KingProfileMorphology(rng_seed=9, **king_kwargs).sample(parent_position)
        c = KingProfileMorphology(rng_seed=10, **king_kwargs).sample(parent_position)
        assert not np.array_equal(a[0], c[0])

    @pytest.mark.parametrize("seed", (9, None))
    def test_repeated_sampling_is_one_realization(
        self, seed, king_kwargs, parent_position
    ):
        # e.g. plot() before to_source() must not change the result
        morph = KingProfileMorphology(rng_seed=seed, **king_kwargs)
        np.testing.assert_array_equal(
            morph.sample(parent_position), morph.sample(parent_position)
        )

    def test_one_generator_for_angle_and_radius(
        self, king_kwargs, parent_position
    ):
        # Pins the draw order: angle, then radius, from ONE generator (two
        # fresh generators would replay the same stream -> correlated draws).
        morph = KingProfileMorphology(rng_seed=9, **king_kwargs)
        rng = new_rng(as_seed_sequence(9))
        phi = morph._sample_phi(rng)
        radius = morph._sample_radius(rng)
        rng_replay = new_rng(as_seed_sequence(9))
        rng_replay.uniform(size=king_kwargs["n_stars"])  # skip the angle draw
        replayed = morph._sampler.rvs(king_kwargs["n_stars"], random_state=rng_replay)
        np.testing.assert_array_equal(radius.value, replayed)
        assert phi.shape == radius.shape


class TestSeededPopulation:
    def test_same_seed_same_masses(self):
        a = IMFPopulation(50, rng_seed=3).sample_imf()
        b = IMFPopulation(50, rng_seed=3).sample_imf()
        np.testing.assert_array_equal(a, b)

    def test_different_seed_differs(self):
        a = IMFPopulation(50, rng_seed=3).sample_imf()
        c = IMFPopulation(50, rng_seed=4).sample_imf()
        assert not np.array_equal(a, c)

    @pytest.mark.parametrize("seed", (3, None))
    def test_repeated_sampling_is_one_realization(self, seed):
        pop = IMFPopulation(50, rng_seed=seed)
        np.testing.assert_array_equal(pop.sample_imf(), pop.sample_imf())

    def test_from_total_mass_passes_seed(self):
        a = IMFPopulation.from_total_mass(20 * u.solMass, rng_seed=3)
        b = IMFPopulation.from_total_mass(20 * u.solMass, rng_seed=3)
        np.testing.assert_array_equal(a.sample_imf(), b.sample_imf())


@pytest.mark.slow
class TestSeededZeroAgeCluster:
    @staticmethod
    def _cluster(rng_seed):
        return ZeroAgeCluster(
            position=SkyCoord(0 * u.deg, 0 * u.deg, 1 * u.kpc),
            pop_class="IMFPopulation",
            pop_params={"n_stars": 20},
            morph_class="KingProfileMorphology",
            morph_params={"n_stars": 20, "r_core": 1 * u.pc, "r_tide": 10 * u.pc},
            rng_seed=rng_seed,
        )

    def test_records_given_seed(self):
        assert self._cluster(42).rng_seed == 42

    def test_records_resolved_entropy_for_none(self):
        # An unseeded cluster is reproducible from its recorded seed.
        original = self._cluster(None)
        assert isinstance(original.rng_seed, int)
        replay = self._cluster(original.rng_seed)
        position = original.position
        np.testing.assert_array_equal(
            original.morphology.sample(position), replay.morphology.sample(position)
        )
        np.testing.assert_array_equal(
            original.population.sample_imf(), replay.population.sample_imf()
        )

    def test_components_use_their_named_streams(self):
        cluster = self._cluster(42)
        np.testing.assert_array_equal(
            _state(cluster.population._seed), _state(stream(42, "population"))
        )
        np.testing.assert_array_equal(
            _state(cluster.morphology._seed), _state(stream(42, "morphology"))
        )
