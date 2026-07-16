# -*- coding: utf-8 -*-
"""Engine-layer flux-conservation tests (flux spec 13.7), PSF disabled.

Where the construction tests (``test_source_flux_construction.py``) check the
*stored artifact* a target produces, these run the same physical source all the
way through a real ScopeSim ``OpticalTrain`` and assert the *detected* flux is
conserved and equal across representations and across pixel scales. The PSF is
switched off so this isolates flux bookkeeping and detector binning from
convolution spreading; flux conservation across PSF convolution is already
covered by ScopeSim's own tests and is not duplicated here.

These integration tests arguably belong in ScopeSim itself (they exercise the
simulation engine, not the target-description framework). But scopesim-targets
depends on ScopeSim and not vice versa, and the tests need the target classes to
build the sources -- so they live here. The point is that the consistency is
tested *somewhere*.

Two independent pixel scales, always labelled (flux spec 9):
  * ``theta_src``  -- WCS pixel scale of the input ImageHDU (a target's grid);
                      varying it exercises the engine's input resampling.
  * ``theta_inst`` -- detector pixel scale of the simulated instrument;
                      varying it exercises detector binning of a fixed source.

Tiering
-------
``webtest``: the sources resolve real spextra templates / CALSPEC. (The example
optical train itself is bundled with ScopeSim and needs no network.) The
detector is read out with PSF and all stochastic/additive effects disabled, so
detected counts are a deterministic function of the input flux.

.. note::
   Tolerances below are first-run starting points per the spec ("start near
   source-level tightness, loosen only if binning/interpolation demands it").
   The upsample case (T10b) is the one most likely to leak and is loosest.
   The cube / spectroscopy round-trip (R4 through SPEC mode) is a heavier,
   separate piece and is deferred -- this file covers the imaging path (R1-R3).
"""

import numpy as np
from numpy import testing as npt
import pytest
import astropy.units as u
from spextra import Spextrum

from scopesim import load_example_optical_train

from scopesim_targets.point_source import Star
from scopesim_targets.extended_source import Box


pytestmark = [
    pytest.mark.webtest,
    # The suite runs with filterwarnings=error. ScopeSim's instrument/detector
    # WCS assembly emits a benign astropy FITSFixedWarning ("END keyrecords...")
    # when the detector pixel scale is overridden -- astropy auto-fixing a header
    # inside ScopeSim, unrelated to our source WCS. Downgrade just that category
    # here (mirrors how pyproject downgrades other third-party warnings) so it
    # doesn't fail the run; the numeric flux assertions still guard correctness.
    pytest.mark.filterwarnings("ignore::astropy.wcs.FITSFixedWarning"),
]

DET_FOV = 51.2 * u.arcsec  # detector field of view, held constant across sweeps
FLUX = 3.5 * u.mJy         # total flux shared by every representation
BOX_PARAMS = {"x_width": 6, "y_width": 4}   # arcsec; A_eff = 24 arcsec**2


@pytest.fixture
def flat_spec():
    return Spextrum.flat_spectrum(
        5 * u.ABmag, waves=np.linspace(3000, 35000, 6000) * u.AA
    )


def _train(theta_inst):
    """A clean example OpticalTrain at detector scale ``theta_inst``."""
    scale = theta_inst.to_value(u.arcsec/u.pixel)
    npix = int(round((DET_FOV / theta_inst).to_value(u.pixel)))
    opt = load_example_optical_train(properties={
        "!INST.pixel_scale": scale,
        "!DET.width": npix,
        "!DET.height": npix,
    })
    opt["psf"].include = False
    opt["atmospheric_radiometry"].include = False
    return opt


def _detected(opt, source):
    """Total detected counts = the engine's F for this source."""
    opt.observe(source, update=True)
    return float(np.nansum(opt.image_planes[0].data))

def _carried(source):
    """Sum of the source's stored image = the input total the engine must
    conserve. A sharp-edged box quantises differently per theta_src (its edges
    land off-grid), so this shifts with the input scale; dividing detected flux
    by it isolates engine resampling fidelity from that template effect (the
    split the flux spec 10 tolerance table draws)."""
    return float(np.nansum(source.fields[0].field.data))


def _src_grid(theta_src):
    """Input ImageHDU grid (theta_src) sized to the constant physical FOV."""
    npix = int(round((DET_FOV / theta_src).to_value(u.pixel)))
    return {"pixel_scale": theta_src, "width": npix, "height": npix}


def _point(flat_spec):
    # R1: point source carries the whole flux; no input grid (table entry).
    return Star(position=(0, 0), spectrum=flat_spec, brightness=("V", FLUX)).to_source()


def _box_integrated(flat_spec, theta_src):
    # R2: 2D image, flux stated as an integrated band flux density.
    return Box(
        spectrum=flat_spec, brightness=("V", FLUX), params=BOX_PARAMS,
    ).to_source(_src_grid(theta_src))


def _box_sb(flat_spec, theta_src):
    # R3: 2D image, the SAME source stated as the equivalent surface brightness.
    A_eff = BOX_PARAMS["x_width"] * u.arcsec * BOX_PARAMS["y_width"] * u.arcsec
    sb = (FLUX / A_eff).to(u.Jy / u.sr)
    return Box(
        spectrum=flat_spec, brightness=("V", sb), params=BOX_PARAMS,
    ).to_source(_src_grid(theta_src))


class TestEngineFluxConservation:
    """T9 / T10 / T11 (flux spec 13.7), PSF off, imaging path."""

    def test_representations_agree(self, flat_spec):
        # T9: point (R1), 2D integrated (R2) and 2D surface brightness (R3) of
        # the same physical flux give equal detected F. theta_src == theta_inst
        # (no input resampling), so a difference points at the flux bookkeeping,
        # not interpolation. R1 is the independent anchor (its flux never touches
        # the image path), so a uniform image-path error shows up as R1 != R2/R3.
        theta = 0.2 * u.arcsec / u.pixel
        opt = _train(theta)
        f_point = _detected(opt, _point(flat_spec))
        f_box = _detected(opt, _box_integrated(flat_spec, theta))
        f_sb = _detected(opt, _box_sb(flat_spec, theta))

        npt.assert_allclose(f_box, f_point, rtol=1e-10)
        npt.assert_allclose(f_sb, f_point, rtol=1e-10)

    def test_downsample_conserves(self, flat_spec):
        # T10a: theta_src FINER than theta_inst -> the engine downsamples/rebins
        # (essentially summation). Per flux spec 10 the engine must conserve the
        # input's *own* total, so compare detected flux per unit carried flux:
        # detected / sum(image) must match the theta_src == theta_inst reference.
        theta_inst = 0.2 * u.arcsec / u.pixel
        opt = _train(theta_inst)
        ref = _box_integrated(flat_spec, theta_inst)
        fine = _box_integrated(flat_spec, 0.1 * u.arcsec / u.pixel)
        r_ref = _detected(opt, ref) / _carried(ref)
        r_fine = _detected(opt, fine) / _carried(fine)
        npt.assert_allclose(r_fine, r_ref, rtol=1e-10)

    def test_upsample_conserves(self, flat_spec):
        # T10b: theta_src COARSER than theta_inst -> the engine interpolates/
        # upsamples, the case most likely to leak. Again the invariant is engine
        # conservation of the input's own total (flux spec 10): detected /
        # sum(image) must be scale-invariant. (Comparing raw detected flux here
        # would instead fold in the box's edge-quantisation at the coarse grid,
        # which is a template effect, not an engine leak.)
        theta_inst = 0.2 * u.arcsec / u.pixel
        opt = _train(theta_inst)
        ref = _box_integrated(flat_spec, theta_inst)
        coarse = _box_integrated(flat_spec, 0.4 * u.arcsec / u.pixel)
        r_ref = _detected(opt, ref) / _carried(ref)
        r_coarse = _detected(opt, coarse) / _carried(coarse)
        npt.assert_allclose(r_coarse, r_ref, rtol=1e-10)

    @pytest.mark.parametrize("theta_inst_mas", (100, 200, 400))
    def test_theta_inst_sweep(self, flat_spec, theta_inst_mas):
        # T11: source fixed, sweep the detector scale -> detected F invariant.
        # The point source is the cleanest fixed source (no input grid, so no
        # theta_src resampling confound); its detected counts should be
        # essentially scale-independent.
        theta_inst = (theta_inst_mas * u.mas / u.pixel).to(u.arcsec / u.pixel)
        opt = _train(theta_inst)
        f = _detected(opt, _point(flat_spec))
        # reference at 200 mas (rebuild to keep each parametrisation independent)
        ref = _detected(_train(0.2 * u.arcsec / u.pixel), _point(flat_spec))
        npt.assert_allclose(f, ref, rtol=1e-10)

    def test_unresolved_extended_equals_point(self, flat_spec):
        # An extended source much smaller than the detector pixel lands
        # (essentially) in one pixel, so it must be indistinguishable from a
        # point source of the same brightness at the ImagePlane. Large detector
        # pixel (2") vs a 0.4" box rendered on a fine source grid (theta_src =
        # 0.05", so the box itself is well sampled, sum(w) = 1). This is the
        # UNRESOLVED limit of the T9 point == extended agreement, stressing
        # sub-pixel placement / detector binning rather than flux bookkeeping.
        #
        # CAUTION: this configuration passes because ScopeSim's decimating
        # downsampler (see TestDownsampleGridPhase below) happens to HIT the
        # centered 8-src-px block with its sparse sample lattice (zoom = 1/40,
        # 128 px grid). Do not "simplify" the grid to a different npix or
        # shrink the box below ~1/zoom source pixels without checking the
        # phase class, or the flux may silently vanish for engine reasons
        # unrelated to what this test asserts.
        opt = _train(2.0 * u.arcsec / u.pixel)
        tiny = Box(
            spectrum=flat_spec, brightness=("V", FLUX),
            params={"x_width": 0.4, "y_width": 0.4},
        ).to_source({
            "pixel_scale": 0.05 * u.arcsec / u.pixel, "width": 128, "height": 128,
        })
        f_point = _detected(opt, _point(flat_spec))
        f_tiny = _detected(opt, tiny)
        npt.assert_allclose(f_tiny, f_point, rtol=1e-10)


class TestDownsampleGridPhase:
    """Deep-sub-detector-pixel sources vs the engine's decimating downsampler.

    ScopeSim resamples an input image onto the FOV grid via
    ``rescale_imagehdu -> scipy.ndimage.zoom(zoom=theta_src/theta_inst,
    order=1)``. For zoom < 1 that is *decimation*, not integration: the output
    lattice point-samples the input, so a source narrower than ~1/zoom input
    pixels is either hit or missed depending purely on the grid phase. For
    zoom = 1/4 the survivor classes are ``npix mod 8 in {3, 4, 5}`` (verified
    against ``ndi.zoom`` directly); every other class returns an all-zero
    array, the ``conserve_flux`` repair is skipped (``if sum_new != 0``), and
    the source silently vanishes -- detected flux exactly 0, no warning.

    This is also why every broad/smooth case in this file conserves at 1e-10:
    any *partial* hit is renormalized to the exact input total by
    ``conserve_flux``, so the aliasing is invisible until the source is narrow
    enough for a clean miss. The invariant these tests state is the T10 one --
    the detected total is a property of the physical source, not of the input
    gridding -- pushed to the unresolved limit where the engine currently
    breaks it.

    The tiny box is 0.5 x 0.5 arcsec on a 0.5 arcsec/px source grid: exactly
    one lit source pixel (odd npix) or a 2x2 quarter-weight block (even npix),
    against theta_inst = 2 arcsec -> zoom = 1/4. The source image must stay
    smaller than the FOV window, or ``extract_from`` crops it and shifts the
    phase class.
    """

    THETA_INST = 2.0 * u.arcsec / u.pixel
    THETA_SRC = 0.5 * u.arcsec / u.pixel  # zoom = 1/4 inside the engine
    TINY_PARAMS = {"x_width": 0.5, "y_width": 0.5}

    def _tiny_box(self, flat_spec, npix):
        return Box(
            spectrum=flat_spec, brightness=("V", FLUX), params=self.TINY_PARAMS,
        ).to_source({
            "pixel_scale": self.THETA_SRC, "width": npix, "height": npix,
        })

    def test_survivor_phase_equals_point(self, flat_spec):
        # npix = 60 (mod 8 = 4): the decimation lattice hits the source, and
        # conserve_flux then pins the total exactly -- the tiny box equals the
        # point source at the ImagePlane. Positive control for the xfail below
        # (same physics, different grid phase), and a regression guard on the
        # unresolved limit itself.
        opt = _train(self.THETA_INST)
        f_point = _detected(opt, _point(flat_spec))
        f_tiny = _detected(opt, self._tiny_box(flat_spec, 60))
        npt.assert_allclose(f_tiny, f_point, rtol=1e-10)

    @pytest.mark.xfail(
        strict=True,
        reason="ScopeSim rescale_imagehdu downsamples by decimation "
               "(ndi.zoom, order=1): on a vanishing-phase source grid "
               "(npix mod 8 not in {3, 4, 5} at zoom = 1/4) the sample "
               "lattice misses a deep-sub-pixel source entirely, the "
               "sum_new == 0 guard skips the conserve_flux repair, and the "
               "detected flux is exactly 0. Remove this marker once the "
               "engine integrates (rebins) instead of decimating.",
    )
    def test_vanishing_phase_equals_point(self, flat_spec):
        # npix = 15 (mod 8 = 7; the documentation example that exposed this):
        # identical physical source and identical brightness as above -- only
        # the input grid size differs -- yet the detected flux is exactly zero.
        # The assertion is the correct physical statement (same as the
        # survivor test); strict xfail turns an upstream fix into a visible
        # pass so the marker gets deleted.
        opt = _train(self.THETA_INST)
        f_point = _detected(opt, _point(flat_spec))
        f_tiny = _detected(opt, self._tiny_box(flat_spec, 15))
        npt.assert_allclose(f_tiny, f_point, rtol=1e-10)


class TestExtendedPlacement:
    """Extended sources land +1 px in x AND y relative to a point source at
    the same position, at matched pixel scales (zoom = 1, so the rescale path
    is bypassed -- this is independent of TestDownsampleGridPhase).

    Mechanism: the WCS handoff computes the canvas pixel coordinate of the
    image centre exactly up to ~1e-8 px of float dust (created when
    ``FieldOfView.extract_from`` recomputes the field header in degrees:
    CRVAL 0.0 -> 3.3e-12 deg). ``overlay_image`` cleans dust only below 5e-11
    (``coords.round(10)``) before an intentional ``np.ceil``, which amplifies
    the surviving +2e-8 px into a full-pixel paste offset. Point sources use
    the table path (``val2pix`` + int truncation) and are unaffected, which is
    why the error shows as a point-vs-extended disagreement.

    Odd grids everywhere so the expected geometry is unambiguous: 15 px
    detector -> centre pixel (7, 7); the 5 x 3 arcsec box at 0.5 arcsec/px has
    exact expected extents x 2..12, y 4..10 (half-covered edge pixels).
    """

    THETA = 0.5 * u.arcsec / u.pixel
    NPIX = 15

    def _train_matched(self):
        opt = load_example_optical_train(properties={
            "!INST.pixel_scale": 0.5,
            "!DET.width": self.NPIX,
            "!DET.height": self.NPIX,
        })
        opt["psf"].include = False
        opt["atmospheric_radiometry"].include = False
        return opt

    def _image(self, opt, source):
        opt.observe(source, update=True)
        return np.asarray(opt.image_planes[0].data)

    @pytest.mark.xfail(
        strict=True,
        reason="ScopeSim overlay_image: np.ceil(coords.round(10)) amplifies "
               "~1e-8 px WCS float dust (from extract_from's deg-space header "
               "recomputation) into a +1 px shift of image fields in x and y; "
               "table sources are placed by a different path and don't shift. "
               "Remove this marker once the dust tolerance upstream is "
               "loosened (round(4) before the ceil).",
    )
    def test_box_centroid_matches_point(self, flat_spec):
        from scipy.ndimage import center_of_mass

        opt = self._train_matched()
        img_pt = self._image(opt, _point(flat_spec))
        py, px = np.unravel_index(np.nanargmax(img_pt), img_pt.shape)
        assert (px, py) == (7, 7), "point-source sanity anchor moved"

        box = Box(
            spectrum=flat_spec, brightness=("V", FLUX),
            params={"x_width": 5, "y_width": 3},
        ).to_source({
            "pixel_scale": self.THETA, "width": self.NPIX, "height": self.NPIX,
        })
        cy, cx = center_of_mass(np.nan_to_num(self._image(opt, box)))
        npt.assert_allclose((cx, cy), (px, py), atol=1e-6)
