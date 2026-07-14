# -*- coding: utf-8 -*-
"""Unit tests for extended_source.py.

Two tiers:

* Geometry / normalization (``total_flux_factor``, ``_weightmap``, param
  validation, traits) -- pure: needs only the module import and astropy, no
  ScopeSim ``Source`` and no spectrum/filter data.
* Full ``to_source`` -- marked ``webtest``: constructs a real ScopeSim
  ``Source`` and resolves/scales a spectrum via spextra.

``optical_train`` is duck-typed: ``to_source`` only reads ``pixel_scale``,
``width`` and ``height``, so a plain dict suffices for tests.
"""

import pytest
import numpy as np
from numpy import testing as npt
import astropy.units as u

from scopesim_targets.extended_source import Box, Disk, Ring, Sersic, Flat
from scopesim_targets.brightness import AnchorFrame
from scopesim_targets.brightness import (
    BrightnessError,
    parse_brightness,
    AmountKind,
)


@pytest.fixture
def optical_train():
    return {"pixel_scale": 1 * u.arcsec / u.pixel, "width": 64, "height": 64}


class TestTraits:
    @pytest.mark.parametrize(
        "cls, finite, ref",
        [
            (Box, True, "uniform"),
            (Disk, True, "uniform"),
            (Ring, True, "uniform"),
            (Sersic, True, "at r_eff"),
            (Flat, False, "the constant"),
        ],
    )
    def test_traits(self, cls, finite, ref):
        assert cls.has_finite_total is finite
        assert cls.sb_reference == ref


class TestTotalFluxFactor:
    def test_box(self):
        assert (
            Box(params={"x_width": 6, "y_width": 4}).total_flux_factor()
            == 24 * u.arcsec**2
        )

    def test_disk(self):
        tff = Disk(params={"R_0": 3}).total_flux_factor()
        npt.assert_allclose(tff.to_value(u.arcsec**2), np.pi * 9)

    def test_ring(self):
        tff = Ring(params={"r_in": 2, "width": 1}).total_flux_factor()
        npt.assert_allclose(tff.to_value(u.arcsec**2), np.pi * (9 - 4))

    def test_sersic(self):
        # cross-checked against a 1-D radial quadrature
        tff = Sersic(
            params={"r_eff": 3, "n": 4, "ellip": 0.3}
        ).total_flux_factor()
        npt.assert_allclose(tff.to_value(u.arcsec**2), 142.7910, rtol=1e-5)

    def test_unit_and_bare_params_agree(self):
        bare = Sersic(params={"r_eff": 3, "n": 4}).total_flux_factor()
        unit = Sersic(
            params={"r_eff": 3 * u.arcsec, "n": 4}
        ).total_flux_factor()
        assert bare == unit

    def test_boxy_sersic_not_implemented(self):
        with pytest.raises(NotImplementedError):
            Sersic(params={"r_eff": 3, "n": 4, "c": 0.5}).total_flux_factor()

    def test_flat_raises_E6(self):
        with pytest.raises(BrightnessError) as exc:
            Flat().total_flux_factor()
        assert exc.value.code == "E6"


class TestParamValidation:
    @pytest.mark.parametrize("reserved", ["amplitude", "x_0", "y_0"])
    def test_reserved_params_raise_E8(self, reserved):
        with pytest.raises(BrightnessError) as exc:
            Box(params={"x_width": 6, "y_width": 4, reserved: 1.0})
        assert exc.value.code == "E8"

    def test_unknown_param_raises_valueerror(self):
        with pytest.raises(ValueError):
            Box(params={"x_width": 6, "bogus": 1})

    def test_valid_params_ok(self):
        assert isinstance(Box(params={"x_width": 6, "y_width": 4}), Box)

    def test_flat_takes_no_params(self):
        assert isinstance(Flat(), Flat)
        with pytest.raises(ValueError):
            Flat(params={"R_0": 3})


class TestWeightMap:
    """Rendering + analytic normalization: no Source, no spectrum, no network."""

    def test_box_weightmap_sums_to_one(self, optical_train):
        # 6x4 box on a 1 arcsec grid, fully contained and pixel-aligned: the
        # midpoint rule is exact, so the weight map sums to exactly 1 (T-BOX).
        wm, _ = Box(params={"x_width": 6, "y_width": 4})._weightmap(
            optical_train
        )
        npt.assert_allclose(wm.sum(), 1.0, rtol=1e-12)

    def test_bunit_is_dimensionless(self, optical_train):
        _, hdr = Box(params={"x_width": 6, "y_width": 4})._weightmap(
            optical_train
        )
        assert hdr["BUNIT"] == ""

    def test_clipped_profile_sums_below_one(self):
        # Sersic wings clipped by a small FOV: analytic normalization does NOT
        # redistribute the missing flux (regression against grid-sum norm).
        small = {
            "pixel_scale": 1 * u.arcsec / u.pixel,
            "width": 16,
            "height": 16,
        }
        wm, _ = Sersic(params={"r_eff": 5, "n": 4})._weightmap(small)
        assert wm.sum() < 0.95

    def test_contained_profile_pixel_scale_invariant(self):
        # A fully-contained box has in-window fraction 1 at any pixel scale.
        box = Box(params={"x_width": 6, "y_width": 4})
        coarse = {
            "pixel_scale": 1.0 * u.arcsec / u.pixel,
            "width": 64,
            "height": 64,
        }
        fine = {
            "pixel_scale": 0.5 * u.arcsec / u.pixel,
            "width": 128,
            "height": 128,
        }
        npt.assert_allclose(
            box._weightmap(coarse)[0].sum(),
            box._weightmap(fine)[0].sum(),
            rtol=1e-12,
        )

    def test_flat_weightmap_is_uniform(self, optical_train):
        # Non-integrable profile: A_eff is the field of view, so the weight map
        # is uniform and sums to 1 (the FOV closes the open integral).
        wm, _ = Flat()._weightmap(optical_train)
        npt.assert_allclose(wm.sum(), 1.0, rtol=1e-12)
        assert np.allclose(wm, wm.flat[0])


class TestSurfaceBrightness:
    """Case II reduction (pure): SB -> implied integrated mag over the effective
    area. Brightness is set directly (bypassing the FilterSystem band check) so
    these stay offline; the reduction itself is synphot-free.
    """

    def test_finite_sb_reduces_to_integrated_mag(self, optical_train):
        box = Box(
            params={"x_width": 6, "y_width": 4}
        )  # total_flux_factor = 24 arcsec2
        box._brightness = parse_brightness(["V", "21.5 mag / arcsec2"])
        eff = box._effective_integrated_brightness(optical_train)
        assert eff.solid_angle is None  # now integrated
        npt.assert_allclose(
            eff.value.to_value(u.mag), 21.5 - 2.5 * np.log10(24), rtol=1e-9
        )

    def test_integrated_passes_through_unchanged(self, optical_train):
        box = Box(params={"x_width": 6, "y_width": 4})
        box._brightness = parse_brightness(["V", "15 mag"])
        assert (
            box._effective_integrated_brightness(optical_train)
            is box._brightness
        )

    def test_sb_equivalent_to_matching_integrated(self, optical_train):
        # A Sersic given SB must reduce to the integrated mag that carries the
        # same total flux (T-B2 at the brightness level).
        ser = Sersic(params={"r_eff": 3, "n": 4})
        tff = ser.total_flux_factor().to_value(u.arcsec**2)
        ser._brightness = parse_brightness(["V", "20 mag / arcsec2"])
        eff = ser._effective_integrated_brightness(optical_train)
        npt.assert_allclose(
            eff.value.to_value(u.mag), 20 - 2.5 * np.log10(tff), rtol=1e-9
        )

    def test_nonintegrable_sb_uses_fov_area(self, optical_train):
        # Flat has no analytic total, so A_eff is the FOV: 64x64 @ 1 arcsec/pix
        # -> A_FOV = 4096 arcsec2.
        flat = Flat()
        flat._brightness = parse_brightness(["V", "21.5 mag / arcsec2"])
        eff = flat._effective_integrated_brightness(optical_train)
        npt.assert_allclose(
            eff.value.to_value(u.mag),
            21.5 - 2.5 * np.log10(64 * 64),
            rtol=1e-9,
        )

    def test_nonintegrable_sb_is_fov_dependent(self):
        # Same SB in a bigger field -> more flux -> brighter total (smaller mag).
        flat = Flat()
        flat._brightness = parse_brightness(["V", "21.5 mag / arcsec2"])
        small = {
            "pixel_scale": 1 * u.arcsec / u.pixel,
            "width": 64,
            "height": 64,
        }
        big = {
            "pixel_scale": 1 * u.arcsec / u.pixel,
            "width": 128,
            "height": 128,
        }
        m_small = flat._effective_integrated_brightness(small).value
        m_big = flat._effective_integrated_brightness(big).value
        assert m_big < m_small

    def test_nonmagnitude_sb_reduces_to_integrated(self, optical_train):
        # Linear surface brightness now reduces too: total = SB * A_eff, with
        # the per-solid-angle divisor stripped (Box total_flux_factor = 24
        # arcsec2). Retires the former Jy/sr deferral.
        box = Box(params={"x_width": 6, "y_width": 4})
        box._brightness = parse_brightness(["V", "5 MJy / sr"])
        eff = box._effective_integrated_brightness(optical_train)
        assert eff.solid_angle is None
        assert eff.amount_kind is AmountKind.FLUX_DENSITY_NU
        npt.assert_allclose(
            eff.value.to_value(u.Jy),
            (5 * u.MJy / u.sr * 24 * u.arcsec**2).to_value(u.Jy),
            rtol=1e-9,
        )

    def test_radiance_sb_reduces_to_energy_flux(self, optical_train):
        # Per-solid-angle energy flux (radiance) reduces to a band-integrated
        # energy flux over the effective area.
        box = Box(params={"x_width": 6, "y_width": 4})
        box._brightness = parse_brightness(["V", "3e-8 W / (m2 arcsec2)"])
        eff = box._effective_integrated_brightness(optical_train)
        assert eff.solid_angle is None
        assert eff.amount_kind is AmountKind.ENERGY_FLUX
        npt.assert_allclose(
            eff.value.to_value(u.W / u.m**2),
            (
                3e-8 * u.W / (u.m**2 * u.arcsec**2) * 24 * u.arcsec**2
            ).to_value(u.W / u.m**2),
            rtol=1e-9,
        )

    def test_nonfinite_integrated_raises_E6(self, optical_train):
        # An *integrated* brightness on a non-integrable profile has no finite
        # total to carry: E6 (raised before any rendering/scaling).
        flat = Flat()
        flat._brightness = parse_brightness(["V", "15 mag"])
        with pytest.raises(BrightnessError) as exc:
            flat.to_source(optical_train)
        assert exc.value.code == "E6"


class TestToSource:
    """Full conversion: builds a real ScopeSim Source and scales a spectrum."""

    @pytest.mark.webtest
    def test_returns_source(self, optical_train):
        box = Box(
            spectrum="G2V",
            brightness=["V", 15],
            params={"x_width": 6, "y_width": 4},
        )
        assert box.to_source(optical_train) is not None

    @pytest.mark.webtest
    def test_position_guard_without_position(self, optical_train):
        # Must not raise AttributeError when position was never set.
        box = Box(
            spectrum="G2V",
            brightness=["V", 15],
            params={"x_width": 6, "y_width": 4},
        )
        box.to_source(optical_train)

    @pytest.mark.webtest
    def test_finite_surface_brightness_ok(self, optical_train):
        # Finite profile + surface brightness renders via the Case II reduction.
        box = Box(
            spectrum="G2V",
            brightness=["V", "21.5 mag / arcsec2"],
            params={"x_width": 6, "y_width": 4},
        )
        assert box.to_source(optical_train) is not None

    @pytest.mark.webtest
    def test_nonintegrable_surface_brightness_ok(self, optical_train):
        # Flat + surface brightness now renders: total flux from the FOV, uniform
        # weight map, flux in the spectrum (no per-pixel image-owns-flux path).
        flat = Flat(spectrum="G2V", brightness=["V", "21.5 mag / arcsec2"])
        assert flat.to_source(optical_train) is not None


class TestAnchorOnProfiles:
    """Anchor-frame behaviour specific to extended (surface-brightness) targets."""

    def test_absolute_surface_brightness_raises_E11(self):
        # A surface brightness is distance-invariant, so anchor: absolute on an
        # SB amount is a category error -- even for a finite profile, and even
        # with a distance present. The guard keys on the *original* SB
        # brightness, before the SB->integrated reduction.
        prof = Sersic(
            spectrum="G2V",
            brightness=("V", "21 mag(AB) / arcsec2"),
            params={"r_eff": 2, "n": 1},
            position={"distance": 25 * u.pc},
            anchor="absolute",
        )
        with pytest.raises(BrightnessError) as exc:
            prof._anchored_spectrum_scale(None, prof.brightness)
        assert exc.value.code == "E11"

    @pytest.mark.webtest
    def test_observed_vs_intrinsic_surface_brightness_with_screen(
        self, optical_train
    ):
        # T-B8. With a 2-mag (A_V) extinction screen:
        #   * observed: the reddened SED is scaled so band photometry matches the
        #     SB value -> the anchored scale is set against the reddened SED.
        #   * intrinsic: the unextincted SED is scaled to the SB value first,
        #     then the screen dims it -> realized band flux is fainter by the
        #     band transmission.
        # So intrinsic / observed == T_band at the brightness band. This is the
        # *law-computed* transmission (CCM89, R_V=3.1), ~0.1583 for A_V=2, which
        # differs at the third digit from the idealized 10**(-0.8)=0.15849
        # because A(V_eff)/A(V) != 1 over a real V bandpass -- assert the law.
        common = dict(
            spectrum="G2V",
            brightness=("V", "18 mag / arcsec2"),
            params={"r_eff": 3, "n": 4},
            extinction="2 mag",
        )
        observed = Sersic(anchor="observed", **common)
        intrinsic = Sersic(anchor="intrinsic", **common)
        s_obs = observed._scale_spectrum(optical_train)
        s_int = intrinsic._scale_spectrum(optical_train)
        npt.assert_allclose(
            _band_flux(s_int, "V") / _band_flux(s_obs, "V"),
            0.1583,
            rtol=1.5e-2,
        )


def _band_flux(spectrum, band):
    """Realized flux of ``spectrum`` through ``band`` (webtest helper)."""
    from synphot import Observation
    from spextra import Passband
    from scopesim_targets.target import FILTER_SYSTEM

    passband = Passband(f"{FILTER_SYSTEM.name}/{band}")
    return Observation(spectrum, passband).effstim(
        u.Unit("erg / (s cm2 AA)")
    ).value
