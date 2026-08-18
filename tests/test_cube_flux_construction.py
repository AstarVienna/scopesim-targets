# -*- coding: utf-8 -*-
"""Construction-layer tests for spectroscopic datacube sources (flux spec 8).

A datacube is the convergent internal representation for spectroscopy: unlike
the 2+1D case there is no image/spectrum split and no provenance axis -- a cube
is *always* differential on every axis (per solid angle spatially, per bin
spectrally), i.e. ``PHOTLAM/arcsec2``, and that unit lives in ``BUNIT`` (flux
spec 8). The Omega round-trip is the crux (flux spec 8.1): forming the cube from
a 2+1D source divides each column by Omega to make it per-solid-angle, and the
collapse multiplies Omega back:

    column(i,j,lam) = w(i,j) * spectrum(lam) / Omega        # PHOTLAM/arcsec2
    F = sum_lam sum_ij column * Omega * dlam                # -> original F

With uniform tangent-plane Omega the factor cancels numerically, so a total-flux
check alone cannot catch a mislabeled cube -- the values themselves must be
verified per-solid-angle (T-CUBE-OMEGA), not just that F collapses correctly.

Tiering
-------
``webtest``: builds a real ScopeSim ``CubeSourceField`` (needs the WCS and BUNIT
machinery). The flux carried is treated as a conservation *invariant* (equality
across representations); absolute photon-flux unit correctness -- reconciling
PHOTLAM's per-Angstrom against a per-micron wavelength grid -- is an engine-layer
concern, so here dlam is taken from ``waveset`` and used consistently on both
sides (the equality is independent of the wavelength unit).
"""

import numpy as np
from numpy import testing as npt
import pytest
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from synphot import units as su

from scopesim import Source
from scopesim.source.source_fields import CubeSourceField


pytestmark = pytest.mark.webtest


@pytest.fixture
def optical_train():
    return {"pixel_scale": 0.2 * u.arcsec / u.pixel, "width": 8, "height": 6}


class MockCubeTarget:
    """Stand-in for a future cube Target (none exists in the package yet).

    Holds a physical source as a unit-normalised weight map ``w`` (sum = 1) plus
    a spectrum sampled on a wavelength grid, and renders the canonical
    differential datacube (flux spec 8): every voxel is ``w * S / Omega`` in
    ``PHOTLAM/arcsec2`` with a full 3D WCS. Omega is divided out here at build
    time and must be multiplied back at collapse -- the round-trip the tests
    assert. ``ctype3`` selects a linear (``WAVE``) or logarithmic (``WAVE-LOG``,
    non-uniform dlam) dispersion axis.
    """

    def __init__(self, weights, spectrum_photlam, *,
                 ctype3="WAVE", crval3=1.0, cdelt3=0.1, cunit3="um"):
        w = np.asarray(weights, dtype=float)
        self.weights = w / w.sum()  # enforce sum(w) = 1
        self.spectrum = np.asarray(spectrum_photlam, dtype=float)
        self.ctype3 = ctype3
        self.crval3 = crval3
        self.cdelt3 = cdelt3
        self.cunit3 = cunit3

    def _pixel_scale(self, optical_train):
        return optical_train["pixel_scale"].to_value(u.arcsec / u.pixel)

    def omega(self, optical_train):
        """Pixel solid angle [arcsec**2] used to differentiate the cube."""
        scale = self._pixel_scale(optical_train)
        return (scale * u.arcsec) ** 2

    def _header(self, optical_train):
        ny, nx = self.weights.shape
        nlam = self.spectrum.size
        scale = self._pixel_scale(optical_train)
        wcs = WCS(naxis=3)
        wcs.wcs.ctype = "RA---TAN", "DEC--TAN", self.ctype3
        wcs.wcs.cunit = "arcsec", "arcsec", self.cunit3
        wcs.wcs.crpix = (nx + 1) / 2, (ny + 1) / 2, 1.0
        wcs.wcs.crval = 0.0, 0.0, self.crval3
        wcs.wcs.cdelt = scale, scale, self.cdelt3
        hdr = wcs.to_header()
        hdr["NAXIS1"], hdr["NAXIS2"], hdr["NAXIS3"] = nx, ny, nlam
        hdr["BUNIT"] = "ph Angstrom-1 s-1 cm-2 arcsec-2"
        return hdr

    def to_source(self, optical_train):
        """Build the differential cube and wrap it in a ScopeSim Source."""
        hdr = self._header(optical_train)
        omega = self.omega(optical_train).to_value(u.arcsec**2)
        # voxel = w(i,j) * S(lam) / Omega -> PHOTLAM/arcsec2, axis order [lam,y,x]
        data = (self.weights[None, :, :] * self.spectrum[:, None, None]) / omega
        hdu = fits.ImageHDU(header=hdr, data=data)
        return Source(field=CubeSourceField(field=hdu))


def _cube(source):
    return source.fields[0]


def _dlam(field):
    """Per-bin wavelength widths [um] from the cube's spectral WCS."""
    return np.gradient(field.waveset.to_value(u.um)) * u.um


def _collapse(field):
    """F = sum_lam sum_ij cube * Omega * dlam, with Omega and dlam from the WCS.

    Omega comes from the spatial WCS (``pixel_area``) and dlam from the spectral
    WCS (``waveset``) -- each read back from the stored artifact, so a wrong WCS
    encoding of either breaks the round-trip.
    """
    data = np.asarray(field.data)  # [lam, y, x], PHOTLAM/arcsec2
    omega = field.pixel_area.to_value(u.arcsec**2)
    dlam = _dlam(field).to_value(u.um)
    per_lam = data.sum(axis=(1, 2))  # collapse space
    return float((per_lam * omega * dlam).sum())


def _analytic_total(target, field):
    """Independent F for the same source: sum_lam S(lam) * dlam (sum(w) = 1)."""
    dlam = _dlam(field).to_value(u.um)
    return float((target.spectrum * dlam).sum())


# Non-trivial, non-symmetric physical source reused across the cube tests.
WEIGHTS = np.zeros((6, 8))
WEIGHTS[2:4, 3:6] = np.array([[1.0, 2.0, 1.0], [2.0, 3.0, 1.0]])
SPECTRUM = np.array([10.0, 8.0, 6.0, 5.0, 4.0, 3.5])  # non-flat PHOTLAM samples


class TestCube:
    """T-CUBE-OMEGA / T-CUBE-EQ / T-CUBE-DISP (flux spec 8, 8.1, 13.4)."""

    def test_cube_is_per_solid_angle(self, optical_train):
        # T-CUBE-OMEGA (flux spec 8.1): the stored values must be genuinely
        # per-solid-angle, not merely collapse to the right F. Assert (a) BUNIT
        # is differential spatially AND spectrally (PHOTLAM/arcsec2), and (b) a
        # sample voxel equals w * S / Omega with Omega read from the WCS.
        target = MockCubeTarget(WEIGHTS, SPECTRUM)
        field = _cube(target.to_source(optical_train))

        assert field.is_bunit_spatially_differential
        assert field.bunit.is_equivalent(su.PHOTLAM / u.arcsec**2)

        omega = field.pixel_area.to_value(u.arcsec**2)
        w = target.weights
        data = np.asarray(field.data)
        for (k, j, i) in [(0, 2, 3), (3, 3, 4), (5, 3, 5)]:
            npt.assert_allclose(
                data[k, j, i], w[j, i] * target.spectrum[k] / omega, rtol=1e-12
            )
        # spatial collapse of one plane returns S(lam)/Omega (since sum(w) = 1)
        npt.assert_allclose(
            data[0].sum(), target.spectrum[0] / omega, rtol=1e-12
        )

    def test_collapse_conserves_flux(self, optical_train):
        # T-CUBE-EQ (flux spec 8): the Omega/dlam round-trip recovers the source
        # total. Omega divides out (in at build, out at collapse) and the spatial
        # weights sum to 1, so sum_lam sum_ij cube * Omega * dlam == sum_lam S*dlam.
        target = MockCubeTarget(WEIGHTS, SPECTRUM)
        field = _cube(target.to_source(optical_train))
        npt.assert_allclose(
            _collapse(field), _analytic_total(target, field), rtol=1e-12
        )

    def test_nonuniform_dispersion(self, optical_train):
        # T-CUBE-DISP (flux spec 11): with a logarithmic (non-uniform) dispersion
        # axis, dlam varies per bin. The collapse must read the per-bin width from
        # waveset -- a single-CDELT3 (uniform) assumption would mis-integrate. F
        # is still conserved.
        target = MockCubeTarget(WEIGHTS, SPECTRUM, ctype3="WAVE-LOG")
        field = _cube(target.to_source(optical_train))

        dlam = _dlam(field).to_value(u.um)
        assert not np.allclose(dlam, dlam[0]), "dispersion axis is not non-uniform"

        npt.assert_allclose(
            _collapse(field), _analytic_total(target, field), rtol=1e-12
        )
