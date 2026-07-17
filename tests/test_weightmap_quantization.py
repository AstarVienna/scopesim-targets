# -*- coding: utf-8 -*-
"""Rasterization-quantization and orientation tests for profile weight maps.

Closes the gap the flux suite left open: the construction and engine tests all
sample pixel scales where the ``Box`` edges happen to land exactly between
pixel centers (6x4 arcsec at 0.05/0.1/0.2/0.4 arcsec/px), so the pixel-center
rasterization error was invisible -- at incommensurate scales it reached tens
of percent (e.g. +10.7% at 0.25, -8.7% at 0.37, +28% at 0.8 arcsec/px). The
fix renders *pixel-averaged* values: exact separable coverage for ``Box``,
midpoint oversampling for the curved sharp edges (``Disk``, ``Ring``) and the
Sersic cusp. These tests pin that behaviour at deliberately awkward scales.

Naming: T-QUANT-* (quantization), T-CLIP-* (profile larger than the FOV, where
the weight sum must be the *visible fraction* while the spectrum keeps the
whole-profile flux -- the reason a grid-sum normalization is never allowed),
T-ORIENT (FITS axis order on non-square grids).

``Gaussian`` is the smooth cross-check: its total (``2 pi sx sy``) and, for
``theta = 0``, its in-window fraction (an erf product) both have closed forms,
so contained *and* clipped weight sums are checked against independent
analytic oracles. Its only discretization error is the midpoint rule on a
smooth integrand, O(scale**2) -- verified to shrink 4x per scale halving --
which is why its tolerances are tighter than the sharp-edged profiles'.

Tiering: everything here is pure (no network, no ScopeSim ``Source``): the
profiles are constructed without ``spectrum``/``brightness`` and only
``_weightmap`` is exercised. The Gaussian ``to_source`` flux recovery lives at
the bottom as ``webtest``, mirroring ``test_source_flux_construction.py``.
"""

import numpy as np
from numpy import testing as npt
from scipy.special import erf
import pytest
import astropy.units as u

from scopesim_targets.extended_source import Box, Disk, Ring, Sersic, Gaussian
from scopesim_targets.brightness import BrightnessError


# Scales chosen to be incommensurate with the 6x4 box (edges land mid-pixel);
# 0.2 is kept as the one commensurate control.
AWKWARD_SCALES = (0.15, 0.2, 0.25, 0.37, 0.8)


def _grid(scale, fov=51.2, fov_y=None):
    npix_w = int(round(fov / scale))
    npix_h = int(round((fov_y if fov_y is not None else fov) / scale))
    return {
        "pixel_scale": scale * u.arcsec / u.pixel,
        "width": npix_w,
        "height": npix_h,
    }


def _weight_sum(target, grid):
    weightmap, _ = target._weightmap(grid)
    return float(weightmap.sum())


class TestBoxQuantization:
    """T-QUANT-BOX: exact coverage kills the pixel-phase quantization."""

    @pytest.mark.parametrize("scale", AWKWARD_SCALES)
    def test_contained_box_sums_to_one(self, scale):
        # A fully contained box must carry exactly its stated total at ANY
        # pixel scale -- including scales where the 6x4 edges land mid-pixel.
        # Pixel-center sampling was off by up to tens of percent here; exact
        # separable coverage is off by nothing.
        box = Box(params={"x_width": 6, "y_width": 4})
        npt.assert_allclose(_weight_sum(box, _grid(scale)), 1.0, rtol=1e-12)

    @pytest.mark.parametrize("scale", (0.2, 0.37))
    def test_overflowing_box_carries_visible_fraction(self, scale):
        # T-CLIP-BOX: box larger than the FOV. The weight sum must be exactly
        # area(box & FOV) / area(box) -- NOT 1 (a grid-sum normalization would
        # force 1 and silently redistribute the unseen flux into the window).
        # The spectrum keeps carrying the whole-box flux.
        grid = _grid(scale)
        fov_x = grid["width"] * scale
        fov_y = grid["height"] * scale
        x_w, y_w = 12.0, 60.0  # overflows in y only (fov ~51.2)
        box = Box(params={"x_width": x_w, "y_width": y_w})
        visible = (min(x_w, fov_x) * min(y_w, fov_y)) / (x_w * y_w)
        npt.assert_allclose(_weight_sum(box, grid), visible, rtol=1e-12)

    def test_subpixel_box(self):
        # Degenerate limit: a box smaller than one pixel still carries exactly
        # its total (one partially-covered pixel), where pixel-center sampling
        # returns either 0 or a full-pixel overestimate.
        box = Box(params={"x_width": 0.3, "y_width": 0.3})
        npt.assert_allclose(_weight_sum(box, _grid(0.8, fov=25.6)), 1.0,
                            rtol=1e-12)


class TestCurvedEdges:
    """T-QUANT-DISK/RING: oversampled sharp curved edges."""

    @pytest.mark.parametrize("scale", AWKWARD_SCALES)
    def test_disk(self, scale):
        disk = Disk(params={"R_0": 3})
        npt.assert_allclose(_weight_sum(disk, _grid(scale)), 1.0, rtol=5e-3)

    @pytest.mark.parametrize("scale", (0.2, 0.37))
    def test_ring(self, scale):
        ring = Ring(params={"r_in": 2, "width": 1.5})
        npt.assert_allclose(_weight_sum(ring, _grid(scale)), 1.0, rtol=5e-3)


class TestSersicCusp:
    """T-QUANT-SERSIC: the cusp no longer drifts with the input scale."""

    def test_carried_fraction_sampling_stable(self):
        # n = 4 is the worst cusp; pixel-center sampling drifted ~4e-3 between
        # these two scales (a systematic under-weighting of the core that
        # shrinks as the grid refines). Oversampled, the carried fraction is a
        # geometric property again.
        fracs = [
            _weight_sum(Sersic(params={"r_eff": 6, "n": 4}), _grid(scale))
            for scale in (0.2, 0.4)
        ]
        npt.assert_allclose(fracs[0], fracs[1], rtol=1e-3)


class TestGaussianOracles:
    """T-QUANT-GAUSS: the smooth profile against its closed forms."""

    def test_total_flux_factor(self):
        gauss = Gaussian(params={"x_stddev": 2, "y_stddev": 3})
        assert u.isclose(
            gauss.total_flux_factor(), 2 * np.pi * 6 * u.arcsec**2
        )

    def test_total_flux_factor_rotation_invariant(self):
        gauss = Gaussian(params={"x_stddev": 2, "y_stddev": 3, "theta": 0.7})
        assert u.isclose(
            gauss.total_flux_factor(), 2 * np.pi * 6 * u.arcsec**2
        )

    @pytest.mark.parametrize("scale", AWKWARD_SCALES)
    def test_contained_matches_erf_oracle(self, scale):
        # Independent oracle: for theta = 0 the in-window fraction of a
        # centered Gaussian in a centered square window of half-width h is
        # erf(h / (sx sqrt(2))) * erf(h / (sy sqrt(2))). Well contained
        # (h ~ 25.6 as >> sigma), the fraction is 1 to double precision, and a
        # smooth profile has no edge quantization -- so this holds tightly at
        # EVERY scale, unlike the sharp-edged profiles above.
        grid = _grid(scale)
        half = grid["width"] * scale / 2
        gauss = Gaussian(params={"x_stddev": 2, "y_stddev": 3})
        oracle = erf(half / (2 * np.sqrt(2))) * erf(half / (3 * np.sqrt(2)))
        npt.assert_allclose(_weight_sum(gauss, grid), oracle, rtol=1e-7)

    def test_clipped_matches_erf_oracle(self):
        # T-CLIP-GAUSS: sigma comparable to the window, so a sizable fraction
        # of the total lies outside the FOV -- the no-sharp-cutoff analogue of
        # the overflowing box. The weight sum must equal the erf window
        # fraction (~0.79), with only the O(scale**2) midpoint-rule residual
        # of a smooth integrand (measured +6.6e-5 at 0.2 arcsec/px, shrinking
        # 4x per scale halving), hence rtol 2e-4.
        grid = _grid(0.2, fov=12.8)
        half = grid["width"] * 0.2 / 2
        gauss = Gaussian(params={"x_stddev": 4, "y_stddev": 4})
        oracle = erf(half / (4 * np.sqrt(2))) ** 2
        assert _weight_sum(gauss, grid) < 0.999  # genuinely clipped
        npt.assert_allclose(_weight_sum(gauss, grid), oracle, rtol=2e-4)

    def test_rotated_total_invariant(self):
        # The carried total of a contained Gaussian is rotation-invariant
        # (2 pi sx sy has no theta); rotation must not leak flux.
        w0 = _weight_sum(
            Gaussian(params={"x_stddev": 2, "y_stddev": 3}), _grid(0.2)
        )
        w1 = _weight_sum(
            Gaussian(params={"x_stddev": 2, "y_stddev": 3, "theta": 0.7}),
            _grid(0.2),
        )
        npt.assert_allclose(w1, w0, rtol=1e-9)

    def test_quantity_params_coerced(self):
        # Units on stddevs/theta normalize to the arcsec/radian render space.
        w_bare = _weight_sum(
            Gaussian(params={"x_stddev": 2, "y_stddev": 3, "theta": 0.7}),
            _grid(0.2),
        )
        w_units = _weight_sum(
            Gaussian(params={
                "x_stddev": 2000 * u.mas,
                "y_stddev": 3 * u.arcsec,
                "theta": (0.7 * u.rad).to(u.deg),
            }),
            _grid(0.2),
        )
        npt.assert_allclose(w_units, w_bare, rtol=1e-9)

    @pytest.mark.parametrize("key", ("amplitude", "x_mean", "y_mean"))
    def test_reserved_params_raise_e8(self, key):
        # Gaussian2D centers are x_mean/y_mean (not x_0/y_0); the class-level
        # _RESERVED override must catch them, or position would be settable
        # through params.
        with pytest.raises(BrightnessError) as excinfo:
            Gaussian(params={"x_stddev": 2, "y_stddev": 3, key: 1.0})
        assert "E8" in str(excinfo.value)


class TestOrientation:
    """T-ORIENT: numpy (height, width) arrays, FITS NAXIS1 = x = width.

    ``Model.render(coords=...)`` expects np.indices-style ordering and
    evaluates ``model(*coords[::-1])``; on the square grids used elsewhere a
    coords-order mistake cancels against the numpy->FITS axis order, so only a
    non-square grid with an asymmetric profile can catch it.
    """

    def test_box_axes_on_nonsquare_grid(self):
        scale = 0.2
        grid = _grid(scale, fov=8.0, fov_y=4.0)  # 40 x 20 pixels
        box = Box(params={"x_width": 6, "y_width": 1})
        weightmap, _ = box._weightmap(grid)

        assert weightmap.shape == (grid["height"], grid["width"])
        # Coverage-weighted extents: sum of per-column (per-row) peak coverage
        # times the pixel scale recovers the box widths on the right axes.
        coverage = weightmap * (6 * 1) / scale**2  # undo weight normalization
        x_extent = coverage.max(axis=0).sum() * scale  # along NAXIS1
        y_extent = coverage.max(axis=1).sum() * scale  # along NAXIS2
        npt.assert_allclose(x_extent, 6.0, rtol=1e-12)
        npt.assert_allclose(y_extent, 1.0, rtol=1e-12)


@pytest.mark.webtest
class TestGaussianConstruction:
    """Gaussian through the full ``to_source`` path (mirrors the Box case in
    ``test_source_flux_construction.py``): stored artifact = weight map summing
    to the erf window fraction (~1 here) + spectrum carrying the stated flux.
    """

    FLUX = 3.5 * u.mJy

    @pytest.fixture
    def flat_spec(self):
        from spextra import Spextrum
        return Spextrum.flat_spectrum(
            5 * u.ABmag, waves=np.linspace(3000, 35000, 6000) * u.AA
        )

    def test_integrated_flux_recovered(self, flat_spec):
        from synphot import Observation
        from spextra import Passband
        from scopesim_targets.target import FILTER_SYSTEM

        grid = _grid(0.2)
        source = Gaussian(
            spectrum=flat_spec,
            brightness=("V", self.FLUX),
            params={"x_stddev": 2, "y_stddev": 3},
        ).to_source(grid)

        w_sum = float(np.asarray(source.fields[0].field.data).sum())
        band = Passband(f"{FILTER_SYSTEM.name}/V")
        spec_flux = Observation(source.fields[0].spectrum, band).effstim(u.mJy)
        npt.assert_allclose(
            w_sum * spec_flux.to_value(u.mJy),
            self.FLUX.to_value(u.mJy),
            rtol=1e-6,
        )
