# -*- coding: utf-8 -*-
"""Star cluster and star field morphologies."""

from numbers import Integral

import numpy as np
from scipy.stats.sampling import NumericalInversePolynomial
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle
from astropy.modeling.functional_models import KingProjectedAnalytic1D
from matplotlib import axes

from ..target import length_angle_context
from ..seeding import SeedLike, as_seed_sequence, new_rng
from ..plot_utils import figure_factory, draw_circle


class Morphology:
    """Base class for stellar cluster morphologies.

    `rng_seed` fixes the realization: every :meth:`sample` of the same object
    returns the same positions (see :mod:`~scopesim_targets.seeding`). The
    default ``None`` draws fresh entropy once, at construction.
    """

    def __init__(self, n_stars: int, rng_seed: SeedLike = None):
        self._n_stars = n_stars
        self._seed = as_seed_sequence(rng_seed)

    def _new_rng(self) -> np.random.Generator:
        """Fresh generator at the start of this object's stream."""
        return new_rng(self._seed)


class SphericallySymmetricalMorphology(Morphology):
    def _sample_phi(self, rng: np.random.Generator) -> Angle:
        return Angle(rng.uniform(-np.pi, np.pi, self._n_stars) * u.rad)


def _to_arcsec(value, parent_position: SkyCoord) -> float:
    """Angle [arcsec] of an angular or (with a distance) length Quantity.

    Plain numbers are arcsec. Lengths need `parent_position` to have a
    distance, and are converted at that distance.
    """
    value = value if isinstance(value, u.Quantity) else value * u.arcsec
    if value.unit.physical_type == "angle":
        return value.to_value(u.arcsec)
    distance = parent_position.distance
    if distance.unit.physical_type != "length":
        raise ValueError(
            f"{value} is a length; converting it to an angle needs a field "
            "center position with a distance"
        )
    with length_angle_context(distance):
        return value.to_value(u.arcsec)


def rectangle_area_outside_circle(
    width: float, height: float, radius: float
) -> float:
    """Area of a rectangle minus a circle, both centered on the origin.

    Exact (analytic), for any radius: also when the circle reaches past the
    edges or covers the whole rectangle (area 0). All inputs in the same
    unit; the result is in that unit squared.
    """
    half_w, half_h, radius = width / 2, height / 2, max(radius, 0.0)

    def _quarter_disc_primitive(x):
        # integral of sqrt(r^2 - x^2) dx
        return 0.5 * (
            x * np.sqrt(max(radius**2 - x**2, 0.0))
            + radius**2 * np.arcsin(min(x / radius, 1.0))
        )

    if radius == 0:
        return width * height
    # One quadrant of the overlap: integral of min(half_h, sqrt(r^2 - x^2))
    # over [0, min(half_w, r)]; the circle is above half_h for x < x_cap.
    x_cap = min(np.sqrt(max(radius**2 - half_h**2, 0.0)), half_w)
    x_end = min(radius, half_w)
    quadrant = half_h * x_cap + (
        _quarter_disc_primitive(x_end) - _quarter_disc_primitive(x_cap)
    )
    return max(width * height - 4 * quadrant, 0.0)


class UniformRectangularMorphology(Morphology):
    """Stars spread uniformly over a rectangle, optionally with a hole.

    The rectangle is centered on the parent position, with sides parallel to
    x (East) and y (North). An optional exclusion radius keeps the central
    circle free, e.g. around a science target. `n_stars` is exact: positions
    in the excluded circle are redrawn.

    Parameters
    ----------
    n_stars : int
        Number of stars (0 is allowed).
    size : float, Quantity or (width, height)
        Rectangle extent; a scalar gives a square. Plain numbers are arcsec;
        lengths are converted at the parent's distance when sampling.
    exclude_radius : float or Quantity, optional
        Radius of the central circle kept free of stars.
    rng_seed : optional
        See :class:`Morphology`.
    """

    def __init__(
        self,
        n_stars: int,
        size,
        exclude_radius=None,
        rng_seed: SeedLike = None,
    ):
        if not isinstance(n_stars, Integral) or isinstance(n_stars, bool):
            raise TypeError(f"n_stars must be an int, got {n_stars!r}")
        if n_stars < 0:
            raise ValueError(f"n_stars must be non-negative, got {n_stars}")
        super().__init__(n_stars=int(n_stars), rng_seed=rng_seed)

        if np.ndim(size) == 0:
            size = (size, size)
        if len(size) != 2:
            raise ValueError(f"size must be a scalar or (width, height), got {size!r}")
        self.size = tuple(size)
        self.exclude_radius = exclude_radius

    def extent_arcsec(self, parent_position: SkyCoord) -> tuple[float, float, float]:
        """Width, height and exclusion radius [arcsec] at the parent."""
        width, height = (_to_arcsec(side, parent_position) for side in self.size)
        if not (width > 0 and height > 0):
            raise ValueError(f"size must be positive, got {self.size!r}")
        radius = (
            0.0 if self.exclude_radius is None
            else _to_arcsec(self.exclude_radius, parent_position)
        )
        if radius < 0:
            raise ValueError(
                f"exclude_radius must be non-negative, got {self.exclude_radius!r}"
            )
        return width, height, radius

    def area(self, parent_position: SkyCoord) -> u.Quantity:
        """Area stars can occupy (rectangle minus exclusion) [arcsec2]."""
        width, height, radius = self.extent_arcsec(parent_position)
        return rectangle_area_outside_circle(width, height, radius) * u.arcsec**2

    def sample(self, parent_position: SkyCoord) -> tuple[np.ndarray, np.ndarray]:
        """Offsets from the parent position ``(x_arcsec, y_arcsec)``."""
        width, height, radius = self.extent_arcsec(parent_position)
        if (
            self._n_stars
            and rectangle_area_outside_circle(width, height, radius) <= 0
        ):
            raise ValueError("exclude_radius covers the whole rectangle")

        rng = self._new_rng()
        low, high = (-width / 2, -height / 2), (width / 2, height / 2)
        if radius == 0:
            xy = rng.uniform(low, high, size=(self._n_stars, 2))
        else:
            # Rejection sampling in deterministic batches (same rng stream).
            acceptance = (
                rectangle_area_outside_circle(width, height, radius)
                / (width * height)
            )
            kept = []
            n_kept = 0
            while n_kept < self._n_stars:
                n_draw = max(
                    int(np.ceil(1.1 * (self._n_stars - n_kept) / acceptance)), 16
                )
                batch = rng.uniform(low, high, size=(n_draw, 2))
                batch = batch[np.hypot(batch[:, 0], batch[:, 1]) >= radius]
                kept.append(batch)
                n_kept += len(batch)
            xy = np.concatenate(kept)[: self._n_stars]

        return xy[:, 0].round(6), xy[:, 1].round(6)

    def to_source_columns(self, parent_position):
        x_arcsec, y_arcsec = self.sample(parent_position)
        return {"x": x_arcsec, "y": y_arcsec}


class KingProfileMorphology(SphericallySymmetricalMorphology):
    def __init__(
        self,
        n_stars: int,
        r_core: float,
        r_tide: float,
        rng_seed: SeedLike = None,
    ):
        super().__init__(n_stars=n_stars, rng_seed=rng_seed)

        class KingRadialProfile(KingProjectedAnalytic1D):
            def pdf(self, x):
                return self(x << self.input_units["x"])

        self.radial_profile = KingRadialProfile(
            amplitude=1,  # PDF sampler doesn't need scaling
            r_core=r_core,
            r_tide=r_tide,
        )
        if self.radial_profile.concentration > 2:
            raise ValueError(
                "concentration should be < 2, but is "
                f"{self.radial_profile.concentration}"
            )

        # Built once (setup is the expensive part); the generator is passed
        # per draw, see `_sample_radius`.
        self._sampler = NumericalInversePolynomial(self.radial_profile)

    @property
    def r_unit(self) -> u.Unit:
        return self.radial_profile.input_units["x"]

    def _sample_radius(self, rng: np.random.Generator):
        return self._sampler.rvs(self._n_stars, random_state=rng) << self.r_unit

    def sample(self, parent_position: SkyCoord) -> SkyCoord:
        # HACK: This is WET with PointSourceTarget.....
        local_frame = parent_position.skyoffset_frame()
        # One generator for both draws: separate fresh generators would replay
        # the same stream and correlate angle with radius.
        rng = self._new_rng()
        phi = self._sample_phi(rng)
        radius = self._sample_radius(rng)
        with length_angle_context(parent_position.distance):
            local_positions = parent_position.directional_offset_by(
                phi, radius
            ).transform_to(local_frame)
            x_arcsec = local_positions.lon.to_value(u.arcsec).round(6)
            y_arcsec = local_positions.lat.to_value(u.arcsec).round(6)
        return x_arcsec, y_arcsec

    def to_source_columns(self, parent_position):
        x_arcsec, y_arcsec = self.sample(parent_position)
        return {"x": x_arcsec, "y": y_arcsec}

    def plot(
        self,
        parent_position: SkyCoord,
        samples: np.ndarray | None = None,  # array-like?
        ax: axes.Axes | None = None,
    ) -> axes.Axes:
        if ax is None:
            _, ax = figure_factory()

        if samples is None:
            samples = self.sample(parent_position)

        center_coord = parent_position.transform_to(
            parent_position.skyoffset_frame()
        )
        center = (
            center_coord.lon.to_value(u.arcsec).round(6),
            center_coord.lat.to_value(u.arcsec).round(6),
        )
        x_arcsec, y_arcsec = samples
        ax.scatter(x_arcsec, y_arcsec, s=3, alpha=.8)

        with length_angle_context(parent_position.distance):
            r_core = self.radial_profile.r_core.quantity.to_value(u.arcsec)
            r_tide = self.radial_profile.r_tide.quantity.to_value(u.arcsec)

        draw_circle(
            ax,
            center,
            r_core,
            label=r"$r_\mathrm{core}$",
            linewidth=2,
            linestyle="dotted",
        )
        draw_circle(
            ax,
            center,
            r_tide,
            label=r"$r_\mathrm{tide}$",
            linewidth=2,
            linestyle="dashdot",
        )
        ax.set_xlim(center[0] - 1.05 * r_tide, center[0] + 1.05 * r_tide)
        ax.set_ylim(center[1] - 1.05 * r_tide, center[1] + 1.05 * r_tide)
        ax.set_aspect("equal")
        ax.set_xlabel("x [arcsec]")
        ax.set_ylabel("y [arcsec]")
        ax.set_title("Cluster Morphology")
        ax.legend()

        return ax
