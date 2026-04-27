# -*- coding: utf-8 -*-
"""Star cluster morphologies."""

import numpy as np
from scipy.stats.sampling import NumericalInversePolynomial
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle
from astropy.modeling.functional_models import KingProjectedAnalytic1D
from matplotlib import axes

from ..target import length_angle_context
from ..plot_utils import figure_factory, draw_circle


class Morphology:
    """Base class for stellar cluster morphologies."""

    def __init__(self, n_stars: int):
        self._n_stars = n_stars
        self._rng = np.random.default_rng()


class SphericallySymmetricalMorphology(Morphology):
    def _sample_phi(self) -> Angle:
        return Angle(self._rng.uniform(-np.pi, np.pi, self._n_stars) * u.rad)


class KingProfileMorphology(SphericallySymmetricalMorphology):
    def __init__(self, n_stars: int, r_core: float, r_tide: float):
        super().__init__(n_stars=n_stars)

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

        self._sampler = NumericalInversePolynomial(
            self.radial_profile,
            random_state=self._rng,
        )

    @property
    def r_unit(self) -> u.Unit:
        return self.radial_profile.input_units["x"]

    def _sample_radius(self):
        return self._sampler.rvs(self._n_stars) << self.r_unit

    def sample(self, parent_position: SkyCoord) -> SkyCoord:
        # HACK: This is WET with PointSourceTarget.....
        local_frame = parent_position.skyoffset_frame()
        with length_angle_context(parent_position.distance):
            local_positions = parent_position.directional_offset_by(
                self._sample_phi(),
                self._sample_radius(),
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
