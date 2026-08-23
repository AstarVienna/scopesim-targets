# -*- coding: utf-8 -*-
"""Targets-level realization of the extended-source flux test spec.

This is the ``test_source_flux_construction.py`` the flux-authority increment
unlocks: it checks that the *stored artifact* a target produces (a dimensionless
weight map + a flux-calibrated spectrum, per the flux spec's canonical boundary
form) carries the flux the ``brightness`` says it should, and that every
representation of the same physical source agrees.

Tiering
-------
These are ``webtest``: ``to_source`` resolves a real spextra template and builds
a ScopeSim ``Source``, and the Vega zero point needs CALSPEC. The *scaling*
invariants underneath them (four-form equivalence T-B1, mag/Jy SB equivalence
T-B6, the magnitude closed form T8) are proven offline against synphot in
``test_flux_scaling.py``; this module adds the Source-construction layer.

.. note::
   ``_source_total_flux`` reconstructs F from the stored artifact = ``sum(weight
   map) * band-integrated spectrum flux`` (flux spec 3, the "spectrum owns the
   flux, image is weights" boundary form). The two attribute accesses marked
   ``TODO(scopesim-api)`` -- pulling the weight-map array and the spectrum out
   of the built ``Source`` -- follow the shapes used elsewhere in the point- and
   extended-source paths but are not exercised in the offline sandbox; confirm
   them on first run and adjust if the ``ImageSourceField`` accessors differ.
"""

import numpy as np
from numpy import testing as npt
from scipy.special import gammaincinv, gamma
import pytest
import astropy.units as u
from astropy.modeling.models import Sersic2D
from synphot import Observation
from synphot.units import FLAM
from spextra import Passband, Spextrum

from scopesim_targets.extended_source import Box, Sersic
from scopesim_targets.target import FILTER_SYSTEM
from scopesim_targets.brightness import BrightnessError


pytestmark = pytest.mark.webtest

BOX_PARAMS = {"x_width": 6, "y_width": 4}


@pytest.fixture
def optical_train():
    """Padded field so profile wings are not clipped."""
    return {
        "pixel_scale": 0.2 * u.arcsec / u.pixel,
        "width": 256,
        "height": 256,
    }


@pytest.fixture
def flat_spec():
    """Simple flat spectrum for testing, produces ~10 PHOTLAM @ 550 nm."""
    return Spextrum.flat_spectrum(
        5*u.ABmag,
        waves=np.linspace(3000, 35000, 6000)*u.AA,
    )

def _band(name):
    return Passband(f"{FILTER_SYSTEM.name}/{name}")


def _weight_sum(source):
    """Sum of the stored dimensionless weight map (the in-window fraction)."""
    return float(np.asarray(source.fields[0].field.data).sum())


def _spectrum_band_flux(source, band, flux_unit=u.mJy):
    """Band-integrated flux of the stored spectrum ALONE (no weight factor).

    This is the absolute flux the spectrum owns -- for a finite profile it is
    the analytic *total*, independent of how much of it the window contains.
    """
    return Observation(source.fields[0].spectrum, band).effstim(flux_unit)


def _source_total_flux(source, band, flux_unit=u.Jy):
    """Reconstruct total flux from a built Source's stored artifact.

    F = (sum of the dimensionless weight map) * (band-integrated flux of the
    spectrum). This is the boundary form the flux spec fixes: the spectrum owns
    the absolute flux, the image is weights.
    """
    return _weight_sum(source) * _spectrum_band_flux(source, band, flux_unit)


def _sersic_window_fraction(r_eff, n, scale, npix, ellip=0.0):
    """Fraction of a Sersic's analytic total that falls inside an npix**2 window.

    Independent oracle for T-B3: astropy's stock ``Sersic2D`` (amplitude=I_e)
    rendered on the same centred arcsec grid, quadratured and divided by the
    closed-form total. It does not touch the target code, so a weight map that
    was (wrongly) grid-sum-normalised -- forcing the fraction to 1 -- disagrees
    with it.
    """
    b_n = gammaincinv(2 * n, 0.5)
    total = (
        2 * np.pi * n * np.exp(b_n) * gamma(2 * n) * b_n ** (-2 * n)
        * (1 - ellip) * r_eff**2
    )
    xs = (np.arange(npix) - (npix - 1) / 2) * scale
    xx, yy = np.meshgrid(xs, xs)
    img = Sersic2D(amplitude=1.0, r_eff=r_eff, n=n, ellip=ellip, x_0=0, y_0=0)(xx, yy)
    return img.sum() * scale**2 / total


class TestIntegratedBrightness:
    """R3-style: a finite profile with an integrated brightness recovers F."""

    def test_box_integrated_flux_density(self, optical_train, flat_spec):
        # A box whose band flux density is stated directly; the stored artifact
        # must reproduce it (weights sum to ~1 for a contained box, spectrum
        # carries the flux).
        box = Box(
            spectrum=flat_spec,
            brightness=("V", 3.5*u.mJy),
            params=BOX_PARAMS,
        )
        f = _source_total_flux(box.to_source(optical_train), _band("V"))
        npt.assert_allclose(f, 3.5*u.mJy, rtol=1e-10)  # accurate

    def test_sersic_analytic_not_grid(self, subtests, optical_train, flat_spec):
        # T-B3 (regression for the grid-sum-normalisation bug). A Sersic whose
        # wings leave the window is normalised by its ANALYTIC total, so:
        #   * the stored spectrum carries the full stated total (3.5 mJy), and
        #   * the weight map sums to the in-window FRACTION (< 1) -- the clipped
        #     wings are NOT redistributed back into the window.
        # Therefore the in-window flux (sum(w) * spectrum) is LESS than the
        # stated brightness. The bug this guards against grid-sum-normalises the
        # map, forcing sum(w) == 1 and the in-window flux == 3.5 mJy always.
        scale = optical_train["pixel_scale"].to_value(u.arcsec / u.pixel)
        npix = optical_train["width"]
        r_eff, n = 8.0, 4.0  # r_eff in arcsec (bare params are arcsec by convention)

        ser = Sersic(
            spectrum=flat_spec,
            brightness=("V", 3.5*u.mJy),
            params={"r_eff": r_eff, "n": n},
        )
        band = _band("V")
        src = ser.to_source(optical_train)

        # (1) the spectrum owns the full analytic total, exactly
        with subtests.test(msg="flux fully in spec"):
            npt.assert_allclose(
                _spectrum_band_flux(src, band), 3.5*u.mJy, rtol=1e-10
            )

        # (2) the weight map is analytically normalised -> clipped fraction < 1,
        #     matching an independent Sersic quadrature (~0.823 here). A grid-sum
        #     bug would give exactly 1.0.
        with subtests.test(msg="weight sum"):
            w = _weight_sum(src)
            frac = _sersic_window_fraction(r_eff, n, scale, npix)
            assert w < 0.95, f"weight sum {w:.4f} looks grid-sum-normalised (should be < 1)"
            npt.assert_allclose(w, frac, rtol=2e-2)

        # (3) consistency: in-window flux == fraction * total, and < the stated total
        with subtests.test(msg="in-window flux"):
            f = _source_total_flux(src, band)
            npt.assert_allclose(f, 3.5*u.mJy * w, rtol=1e-10)
            assert f < 3.5*u.mJy


class TestFourFormEquivalence:
    """T-B1: the same physical box built via all four 5.3.6 pairs -> same F.

    The four amounts are derived from one reference build so they describe the
    identical source; equivalence then tests only the input-form dispatch.
    """

    def test_four_forms(self, subtests, optical_train, flat_spec):
        ref = Box(
            spectrum=flat_spec,
            brightness=("V", 10*u.ABmag),
            params=BOX_PARAMS,
        )
        band = _band("V")
        src = ref.to_source(optical_train)
        f_ref = _source_total_flux(src, band)  # band flux density, in Jy

        # Derive the wavelength- and frequency-locator amounts from the SAME
        # built SED so all four forms describe the identical physical source;
        # equivalence then tests only the locator/amount dispatch in
        # flux_scaling (band-mag, band-density, monochromatic-per-lambda,
        # monochromatic-per-nu), not four different sources.
        spectrum = src.fields[0].spectrum
        lam0 = band.pivot()  # a representative in-band wavelength
        nu0 = lam0.to(u.GHz, equivalencies=u.spectral())
        f_lam0 = spectrum(lam0, flux_unit=FLAM)  # F_lambda at lam0
        f_nu0 = spectrum(lam0, flux_unit=u.Jy)   # F_nu at lam0

        forms = [
            ("V", f_ref.to(u.Jy)),  # band + spectral flux density (per nu)
            ("V", 10*u.ABmag),      # band + magnitude
            (lam0, f_lam0),         # wavelength + flux density (per lambda)
            (nu0, f_nu0),           # frequency  + flux density (per nu)
        ]
        for brightness in forms:
            with subtests.test(msg="ref flux", brightness=brightness):
                box = Box(
                    spectrum=flat_spec,
                    brightness=brightness,
                    params=BOX_PARAMS,
                )
                f = _source_total_flux(box.to_source(optical_train), band)
                npt.assert_allclose(f, f_ref, rtol=1e-3)


class TestSurfaceBrightness:
    """T7/T8/T-B6: SB units and mag/arcsec2 reduce to the right total flux."""

    @pytest.mark.parametrize(
        "amount", (5*u.MJy/u.sr, 3.5*u.Jy/u.arcsec**2),
    )
    def test_linear_sb_units_equal(self, optical_train, flat_spec, amount):
        # T7: linear SB reduces to total = SB * A_eff over the analytic area.
        box = Box(
            spectrum=flat_spec,
            brightness=("V", amount),
            params=BOX_PARAMS,
        )
        A_eff = box.total_flux_factor()
        b = box._effective_integrated_brightness(optical_train)
        # The reduced integrated amount equals SB * A_eff independent of units.
        expected = (u.Quantity(amount) * A_eff).to(b.value.unit)
        npt.assert_allclose(b.value, expected, rtol=1e-9)

    @pytest.mark.parametrize("unit", ("mag", "mag(Vega)", "mag(AB)", "mag(ST)"))
    def test_mag_arcsec2_closed_form(self, optical_train, flat_spec, unit):
        # T8 (spec 4 / 13.2): the logarithmic analogue of the linear T7 check.
        # A uniform box at surface brightness mu [mag/arcsec2] reduces to the
        # integrated magnitude by the closed form
        #     m_int = mu - 2.5 * log10(A_arcsec2),
        # i.e. delog -> area-weight -> relog, done in LINEAR space. A naive
        # per-pixel magnitude area-scaling (summing/scaling in mag space) gets
        # the mu/area split wrong; this pins it to 1e-10. The closed form is a
        # magnitude difference, so it is independent of the photometric system
        # (the zero point cancels) -- AB here to be explicit; T8b will sweep
        # {AB, ST, Vega}.
        mu = 21.5  # mag(AB) / arcsec2
        box = Box(
            spectrum=flat_spec,
            brightness=("V", f"{mu} {unit} / arcsec2"),
            params=BOX_PARAMS,
        )
        b = box._effective_integrated_brightness(optical_train)
        A_arcsec2 = box.total_flux_factor().to_value(u.arcsec**2)
        expected = (mu - 2.5 * np.log10(A_arcsec2)) * u.mag
        npt.assert_allclose(
            b.value.to_value(u.mag), expected.to_value(u.mag), rtol=1e-10
        )

    def test_mag_arcsec2_matches_jy_sr(self, optical_train, flat_spec):
        # T-B6 at the Source level: mag/arcsec2 and the equivalent Jy/sr give the
        # same total flux (conversion to linear precedes area scaling).
        band = _band("V")
        box_mag = Box(
            spectrum=flat_spec,
            brightness=("V", "21.5 mag / arcsec2"),
            params=BOX_PARAMS,
        )
        f_mag = _source_total_flux(box_mag.to_source(optical_train), band)
        # the Jy/sr that carries the same total over the same area
        A_eff = box_mag.total_flux_factor()
        box_jy = Box(
            spectrum=flat_spec,
            brightness=("V", (f_mag / A_eff).to(u.Jy/u.sr)),
            params=BOX_PARAMS,
        )
        f_jy = _source_total_flux(box_jy.to_source(optical_train), band)
        npt.assert_allclose(f_mag, f_jy, rtol=1e-3)


class TestDispatchAndDoubleCounting:
    """Flux authority lives in exactly one place (flux spec 3.3-3.4).

    T-D1  -- the input flux authority (total vs surface brightness) is a
             dispatch, not a different answer: both build the same F.
    T-DBL -- the input SED's own absolute level is divided out, never multiplied
             back in.
    E8    -- an independent second flux authority (``amplitude``) is over-
             determination and must raise, not silently win.
    """

    def test_integrated_and_sb_agree(self, optical_train, flat_spec):
        # T-D1 (spec 3.4): the same physical box, built from a total flux and
        # from the equivalent surface brightness, gives identical F. The SB path
        # folds Omega in (F = SB * A_eff) while the integrated path passes the
        # flux straight through; the two must converge. Box is fully contained
        # (sum(w) = 1), so F also equals the independent point-source value 3.5
        # mJy -- the R1 anchor that would expose a uniform image-path error.
        band = _band("V")

        integrated = Box(
            spectrum=flat_spec,
            brightness=("V", 3.5*u.mJy),
            params=BOX_PARAMS,
        )
        A_eff = integrated.total_flux_factor()  # analytic area [arcsec**2]
        sb_amount = (3.5*u.mJy / A_eff).to(u.Jy/u.sr)  # same object, as SB
        surface = Box(
            spectrum=flat_spec,
            brightness=("V", sb_amount),
            params=BOX_PARAMS,
        )

        f_int = _source_total_flux(integrated.to_source(optical_train), band)
        f_sb = _source_total_flux(surface.to_source(optical_train), band)
        npt.assert_allclose(f_sb, f_int, rtol=1e-9)
        npt.assert_allclose(f_int, 3.5*u.mJy, rtol=1e-9)

    def test_double_counting_level_divided_out(self, optical_train, flat_spec):
        # T-DBL (spec 3.3): the flux scale is set by the brightness, not by the
        # input SED's normalisation. Building the same SB source from a spectrum
        # and from that spectrum * k (k != 1) must give identical F -- if the
        # level leaked through instead of being divided out, F would be off by k.
        band = _band("V")
        # Box analytic area
        A_eff = BOX_PARAMS["x_width"]*u.arcsec * BOX_PARAMS["y_width"]*u.arcsec
        sb = ("V", (3.5*u.mJy / A_eff).to(u.Jy/u.sr))

        base = Box(
            spectrum=flat_spec,
            brightness=sb,
            params=BOX_PARAMS,
        )
        rescaled = Box(
            spectrum=flat_spec * 1000.0,
            brightness=sb,
            params=BOX_PARAMS,
        )

        f_base = _source_total_flux(base.to_source(optical_train), band)
        f_scaled = _source_total_flux(rescaled.to_source(optical_train), band)
        npt.assert_allclose(f_scaled, f_base, rtol=1e-9)

    def test_amplitude_param_is_over_determined(self, flat_spec):
        # E8 (spec 3.4): 'brightness' is the sole flux authority. A model
        # 'amplitude' independently sets the flux, so supplying it alongside a
        # brightness is over-determined -- the construction must raise, never
        # silently honour one and discard the other. (With a single-brightness
        # API, 'amplitude' is the only route to over-determination, so E8 is the
        # realisation of the spec's "supplying both must not pass silently".)
        with pytest.raises(BrightnessError) as excinfo:
            Box(
                spectrum=flat_spec,
                brightness=("V", 3.5*u.mJy),
                params=BOX_PARAMS | {"amplitude": 2.0},
            )
        assert "E8" in str(excinfo.value)


class TestSampling:
    """Carried flux of a smooth profile is a geometric property, robust to the
    input sampling once the profile is resolved. This is the counterpart to T-B3
    (which pins the carried fraction at one scale): here the fraction must
    (a) match the independent analytic window fraction at each scale, and
    (b) be stable across scales. It is exactly the property a sharp-edged box
    lacks -- its edges land off-grid, so its carried fraction drifts with
    theta_src -- which is why the engine layer conserves the input's *own* total
    rather than comparing raw detected flux across scales (flux spec 10
    tolerance table: template quadrature is sampling-dependent, engine
    conservation is tight).
    """

    def test_smooth_profile_sampling_invariant(self, subtests, flat_spec):
        # Same Sersic (n=3, smooth) well inside a physically FIXED window, sampled
        # at two well-resolved scales (r_eff = 6 arcsec = 60 px @ 0.1, 30 px @ 0.2).
        # Only the input pixel scale changes; the window fraction is geometric,
        # so the carried weight-map sum should barely move.
        r_eff, n = 6.0, 3.0
        fov = 51.2  # arcsec; constant physical window (window ~= 4.3 r_eff)
        carried = {}
        for scale in (0.1, 0.2):
            npix = int(round(fov / scale))
            grid = {
                "pixel_scale": scale * u.arcsec / u.pixel,
                "width": npix,
                "height": npix,
            }
            ser = Sersic(
                spectrum=flat_spec,
                brightness=("V", 3.5 * u.mJy),
                params={"r_eff": r_eff, "n": n},
            )
            w = _weight_sum(ser.to_source(grid))
            # (a) equals the independent analytic window fraction at this scale
            with subtests.test(msg="matches analytic fraction", scale=scale):
                frac = _sersic_window_fraction(r_eff, n, scale, npix)
                npt.assert_allclose(w, frac, rtol=1e-2)
            carried[scale] = w
        # (b) sampling-invariance: the two well-sampled scales agree tightly
        # (the analytic oracle drifts < 1e-3 here, so this is the real teeth).
        with subtests.test(msg="sampling-invariant across scales"):
            npt.assert_allclose(carried[0.1], carried[0.2], rtol=2e-3)
