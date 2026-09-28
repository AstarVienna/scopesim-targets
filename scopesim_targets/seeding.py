# -*- coding: utf-8 -*-
"""Reproducible random streams derived from one master seed.

Generative targets (clusters, random star fields) take a single ``rng_seed``.
It is resolved once into a master seed -- the given int, or fresh OS entropy
for ``None`` -- which the target records, so every realization, seeded or not,
can be reproduced from the target alone.

Each independent sampling step (population, morphology, extinction, ...) draws
from its own named stream, a child of the master seed. Streams are addressed by
name via :data:`STREAMS`, so a target can use any subset of them and adding a
new stream never changes the draws of existing ones.

Samplers keep a :class:`~numpy.random.SeedSequence`, not a
:class:`~numpy.random.Generator`, and create a fresh generator per draw
(:func:`new_rng`). A generator is stateful, so a stored one would make results
depend on call history (e.g. ``plot()`` before ``to_source()``); a seed
sequence is not, so every draw of the same object gives the same realization.
"""

from numbers import Integral

import numpy as np

__all__ = [
    "SeedLike",
    "STREAMS",
    "resolve_seed",
    "stream",
    "as_seed_sequence",
    "new_rng",
]


SeedLike = int | np.random.SeedSequence | None
"""What samplers accept as ``rng_seed``."""

# APPEND-ONLY: these indices are part of the reproducibility contract. Changing
# or reusing one silently changes every recorded realization that uses it.
STREAMS: dict[str, int] = {
    "population": 0,
    "morphology": 1,
    "extinction": 2,
    "count": 3,  # Poisson star counts (random star fields)
    "magnitudes": 4,  # per-star magnitudes (random star fields)
    "spectra": 5,  # per-star spectrum choice (random star fields)
}
"""Named child streams of a master seed (index = SeedSequence spawn key)."""


def resolve_seed(rng_seed: int | None) -> int:
    """Master seed to use and record: `rng_seed`, or fresh entropy if None.

    Parameters
    ----------
    rng_seed : int or None
        Non-negative int, or None for a fresh, unpredictable seed.

    Returns
    -------
    int
        The master seed. For ``None`` this is the OS entropy actually drawn
        (a 128-bit int); passing it back as `rng_seed` reproduces the draw.
    """
    if rng_seed is not None and (
        not isinstance(rng_seed, Integral) or isinstance(rng_seed, bool)
    ):
        raise TypeError(f"rng_seed must be an int or None, got {rng_seed!r}")
    if rng_seed is not None and rng_seed < 0:
        raise ValueError(f"rng_seed must be non-negative, got {rng_seed}")
    return int(np.random.SeedSequence(rng_seed).entropy)


def stream(master_seed: int, name: str) -> np.random.SeedSequence:
    """Named, independent child stream of `master_seed`.

    Identical to child ``STREAMS[name]`` of
    ``SeedSequence(master_seed).spawn(...)``, but stateless: it does not depend
    on how many children were spawned before.
    """
    try:
        index = STREAMS[name]
    except KeyError:
        raise ValueError(
            f"unknown stream {name!r}; known streams: {sorted(STREAMS)}"
        ) from None
    return np.random.SeedSequence(master_seed, spawn_key=(index,))


def as_seed_sequence(rng_seed: SeedLike) -> np.random.SeedSequence:
    """Normalize a sampler's `rng_seed` argument into a SeedSequence.

    A SeedSequence (typically a :func:`stream` handed down by a target) is
    used as is; an int or None is resolved via :func:`resolve_seed`, which
    fixes None to one draw of entropy at this point.
    """
    if isinstance(rng_seed, np.random.SeedSequence):
        return rng_seed
    return np.random.SeedSequence(resolve_seed(rng_seed))


def new_rng(seed_sequence: np.random.SeedSequence) -> np.random.Generator:
    """Fresh generator at the start of `seed_sequence`'s stream.

    ``default_rng(SeedSequence)`` does not advance the sequence, so repeated
    calls return generators that produce identical draws.
    """
    return np.random.default_rng(seed_sequence)
