# -*- coding: utf-8 -*-
"""Engine-layer flux conservation for datacube sources (flux spec 13.7).

The construction cube tests (``test_cube_flux_construction.py``) check the
stored differential cube; these push it through a real ScopeSim spectroscopic
engine and assert the flux survives. The mode is ``simple_ifu`` -- the
"virtual 3D detector" (``FluxBinning3D``) that re-images the cube onto the
instrument's own (spatial, spectral) grid WITHOUT a slit, so there are no slit
losses to confound a flux-conservation check (the long-slit ``spectroscopy``
mode deliberately clips the source and is not a conservation test).

The invariant, the cube analogue of the imaging T10 conservation:
the detected total is a property of the *physical* source, invariant to how the
input cube is gridded -- both spatially (theta_src) and spectrally (lambda
binning). The engine only rebins; it must not create or lose flux doing so.

Reading the flux
----------------
As in ``test_source_flux_engine.py`` the flux is read at the ImagePlane -- the
layer where it is final, before (downstream, flux-preserving) detector noise. In
``simple_ifu`` the ImagePlane is itself 3D (a cube), which doubles as the
structural check that the source round-trips *as a cube* (R4 has no
image/spectrum split). ``readout()`` is deliberately not called: the 3D-detector
readout path has a noise-frame broadcasting bug in this ScopeSim version, and it
is downstream of the flux anyway.

BUNIT
-----
The cube must carry ``BUNIT = 'ph Angstrom-1 s-1 cm-2 arcsec-2'`` -- the explicit
FITS spelling of PHOTLAM/arcsec2 (per solid angle, per wavelength). ScopeSim's
``prepare_source`` parses BUNIT with a bare ``u.Unit()``, which rejects the
astropy round-tripped ``'PHOTLAM'`` spelling; the explicit form parses and
round-trips cleanly (the same spelling the construction cube mock uses).

"""

import numpy as np
from numpy import testing as npt
import pytest
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS

from scopesim import Source, load_example_optical_train
from scopesim.source.source_fields import CubeSourceField
from synphot import SourceSpectrum
from synphot.models import ConstFlux1D
import synphot.units as su


# simple_ifu's WCS/detector assembly emits the same benign astropy
# FITSFixedWarning as the imaging engine; the suite runs filterwarnings=error.
# (No module-level webtest marker: the pure-cube tests need no network -- no
# spextra template, and basic_instrument is bundled. The spextra-based
# cross-geometry class below carries its own webtest marker.)
pytestmark = pytest.mark.filterwarnings("ignore::astropy.wcs.FITSFixedWarning")

# PHOTLAM/arcsec2 in explicit FITS spelling (parses AND round-trips in astropy).
BUNIT = "ph Angstrom-1 s-1 cm-2 arcsec-2"
# simple_ifu spans 1.7-2.5 um; keep the cube well inside so nothing is clipped.
WAVE_LO, WAVE_HI = 1.9, 2.3  # um


def _sed(waves_um):
    """A simple LINEAR SED, so any spectral sampling reconstructs it exactly --
    input spectral-grid invariance is then a clean equality, not an accuracy
    tolerance."""
    return 10.0 - 6.0 * (waves_um - WAVE_LO) / (WAVE_HI - WAVE_LO)


def _cube_source(theta_src=0.2, nlam=5, fov=6.4, spec_samples=None):
    """Differential datacube of a fixed physical source (a small on-sky box with
    unit-normalised weights) on a theta_src grid with nlam wavelength bins.
    ``spec_samples`` overrides the default linear SED (used to match an
    independent point source's spectrum in the cross-path check below)."""
    waves = np.linspace(WAVE_LO, WAVE_HI, nlam)
    spec = _sed(waves) if spec_samples is None else np.asarray(spec_samples)
    npix = int(round(fov / theta_src))
    w = np.zeros((npix, npix))
    c = npix // 2
    hx, hy = int(round(0.8 / theta_src)), int(round(1.0 / theta_src))
    w[c - hy:c + hy, c - hx:c + hx] = 1.0
    w /= w.sum()  # sum(w) = 1 at any theta_src -> carried total is grid-invariant
    omega = theta_src**2  # arcsec**2
    data = (w[None, :, :] * spec[:, None, None]) / omega  # PHOTLAM/arcsec2

    wcs = WCS(naxis=3)
    wcs.wcs.ctype = "RA---TAN", "DEC--TAN", "WAVE"
    wcs.wcs.cunit = "arcsec", "arcsec", "um"
    wcs.wcs.crpix = (npix + 1) / 2, (npix + 1) / 2, 1.0
    wcs.wcs.crval = 0.0, 0.0, waves[0]
    wcs.wcs.cdelt = theta_src, theta_src, waves[1] - waves[0]
    hdr = wcs.to_header()
    hdr["NAXIS1"], hdr["NAXIS2"], hdr["NAXIS3"] = npix, npix, nlam
    hdr["BUNIT"] = BUNIT
    return Source(field=CubeSourceField(field=fits.ImageHDU(header=hdr, data=data)))


def _observe_cube_total(source):
    """Total flux at the (3D) ImagePlane after the simple_ifu engine."""
    opt = load_example_optical_train(set_modes=["simple_ifu"])
    opt["psf"].include = False
    opt["atmospheric_radiometry"].include = False
    opt.observe(source, update=True)
    return opt.image_planes[0]


class TestCubeEngine:
    """T-CUBE-EQ at the engine layer (flux spec 8, 13.7), simple_ifu mode."""

    def test_round_trips_as_cube(self):
        # Structural: the source stays a cube through the engine -- the ImagePlane
        # is 3D with a spectral axis (R4 has no image/spectrum split), not a
        # collapsed 2D frame.
        ip = _observe_cube_total(_cube_source())
        data = np.asarray(ip.data)
        assert data.ndim == 3, f"expected a 3D ImagePlane, got shape {data.shape}"
        assert data.shape[0] >= 2, "spectral axis should have >1 bin"
        assert ip.hdu.header.get("CTYPE3", "").startswith("WAVE")

    def test_spatial_sampling_conserves(self):
        # Cube analogue of T10: the SAME physical source (same on-sky footprint,
        # same SED) built at two spatial samplings gives the SAME detected total.
        # sum(w) = 1 by construction, so the carried total is theta_src-invariant;
        # the engine's spatial rebinning must not change it.
        f_coarse = np.nansum(_observe_cube_total(_cube_source(theta_src=0.2)).data)
        f_fine = np.nansum(_observe_cube_total(_cube_source(theta_src=0.1)).data)
        npt.assert_allclose(f_fine, f_coarse, rtol=1e-10)

    def test_spectral_sampling_invariant(self):
        # The detected total depends on the physical SED, not on the input
        # wavelength binning: the engine rebins the input lambda grid onto its
        # own without creating or losing flux. Linear SED -> exact equality.
        f5 = np.nansum(_observe_cube_total(_cube_source(nlam=5)).data)
        f9 = np.nansum(_observe_cube_total(_cube_source(nlam=9)).data)
        npt.assert_allclose(f9, f5, rtol=1e-10)


class TestCubeVsPoint:
    """Absolute cross-path normalisation (the axis the self-consistency tests
    above cannot see): a datacube and a POINT source of the SAME SED should give
    the same IFU total. They do not -- ScopeSim under-counts the FIRST spectral
    bin of the cube path by half (per-lambda ratio [0.5, 1, 1, ...]), so
    cube_total = point_total * (N - 0.5) / N  (0.90 for the N=5 output grid).

    This is an engine bug, not a mock artefact: it survives over-covering the
    input spectral range far beyond the instrument band, so the cube does supply
    the first bin's wavelengths -- ScopeSim's FluxBinning3D drops half of it. It
    cancels in the cube-vs-cube tests above (both sides lose the same half-bin),
    which is exactly why those pass at 1e-10 and why nothing caught it until a
    cube was compared against a non-cube reference. A native ScopeSim point
    source is used (no spextra) so the check stays offline and isolates the
    engine paths.

    Marked xfail(strict): when the first-bin handling is fixed upstream this
    starts passing, the strict marker turns that into a failure, and the fix is
    to delete the marker -- the assertion is already the correct one.
    """

    @staticmethod
    def _ifu_total(source):
        opt = load_example_optical_train(set_modes=["simple_ifu"])
        opt["psf"].include = False
        opt["atmospheric_radiometry"].include = False
        opt.observe(source, update=True)
        return float(np.nansum(opt.image_planes[0].data))

    @pytest.mark.xfail(
        strict=True,
        reason="ScopeSim under-counts the first cube spectral bin by half "
               "(FluxBinning3D edge): cube total = point total * (N-0.5)/N. "
               "Remove this marker once fixed upstream.",
    )
    def test_cube_matches_point(self):
        spec = SourceSpectrum(ConstFlux1D, amplitude=15*u.ABmag)
        point = Source(
            x=[0]*u.arcsec, y=[0]*u.arcsec, ref=[0], weight=[1], spectra=[spec],
        )
        samples = spec(
            np.linspace(WAVE_LO, WAVE_HI, 5)*u.um, flux_unit=su.PHOTLAM
        ).value
        cube = _cube_source(spec_samples=samples)
        npt.assert_allclose(
            self._ifu_total(cube), self._ifu_total(point), rtol=1e-3
        )


@pytest.mark.webtest
@pytest.mark.filterwarnings("ignore::astropy.wcs.FITSFixedWarning")
class TestCrossGeometryIFU:
    """Representation-independence in the 3D FOV (flux spec 13.7): a point source
    and a 2D extended source of the same brightness give the same total at the
    (3D) IFU ImagePlane -- '2D blown up to 3D' agreeing with 'point in 3D'.

    Uses the target classes (Star, Box), so unlike the pure-cube tests above this
    needs spextra -> webtest.

    NOTE: the datacube geometry is compared against the point/2D path separately,
    in TestCubeVsPoint -- it currently disagrees by a factor (N-0.5)/N from a
    first-spectral-bin engine bug, so it lives there as a strict xfail rather
    than here. The cube's flux *conservation* (invariance to input sampling) is
    locked by TestCubeEngine.
    """

    FLUX = 3.5*u.mJy

    @staticmethod
    def _ifu_total(source):
        opt = load_example_optical_train(set_modes=["simple_ifu"])
        opt["psf"].include = False
        opt["atmospheric_radiometry"].include = False
        opt.observe(source, update=True)
        return float(np.nansum(opt.image_planes[0].data))

    def test_point_equals_2d(self, flat_spec):
        from scopesim_targets.point_source import Star
        from scopesim_targets.extended_source import Box

        point = Star(
            position=(0, 0), spectrum=flat_spec, brightness=("V", self.FLUX),
        ).to_source()
        box = Box(
            spectrum=flat_spec, brightness=("V", self.FLUX),
            params={"x_width": 1.6, "y_width": 2.0},
        ).to_source({
            "pixel_scale": 0.2 * u.arcsec / u.pixel, "width": 32, "height": 32,
        })
        npt.assert_allclose(
            self._ifu_total(box), self._ifu_total(point), rtol=1e-10
        )


@pytest.fixture
def flat_spec():
    from spextra import Spextrum
    return Spextrum.flat_spectrum(
        5 * u.ABmag, waves=np.linspace(3000, 35000, 6000) * u.AA
    )
