# -*- coding: utf-8 -*-
"""Unit tests for extended_source.py.

Two tiers:

* Geometry / normalization (``total_flux_factor``, ``_weightmap``, param
  validation, traits) -- pure: needs only the module import and astropy, no
  ScopeSim ``Source`` and no spectrum/filter data.
* Full ``to_source`` -- marked ``webtest``: constructs a real ScopeSim
  ``Source`` and resolves/scales a spectrum via spextra.

"""

import pytest
import numpy as np
from numpy import testing as npt
import astropy.units as u
from astropy.coordinates import SkyCoord
from scipy.integrate import quad

from scopesim_targets.extended_source import (
    Box,
    Disk,
    Ring,
    Sersic,
    Flat,
    CoredPowerLaw,
    CoredPowerLaw2D,
)
from scopesim_targets.brightness import (
    BrightnessError,
    parse_brightness,
    AmountKind,
)


@pytest.fixture
def grid():
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


    @pytest.mark.parametrize(
        "index, finite",
        [(3.0, True), (2.5, True), (2.0, False), (1.0, False), (0.0, False)],
    )
    def test_cored_power_law_traits(self, index, finite):
        # Finiteness depends on the instance (converges only for index > 2), so
        # it is a property: check instances, not the class (cf. test_traits).
        prof = CoredPowerLaw(params={"r_core": 2, "index": index})
        assert prof.has_finite_total is finite
        assert CoredPowerLaw.sb_reference == "at r_ref"

    @pytest.mark.parametrize("index", [0.0, 1.0, 2.0, 3.0])
    def test_truncated_cored_power_law_is_finite(self, index):
        prof = CoredPowerLaw(params={"r_core": 2, "index": index, "r_out": 20})
        assert prof.has_finite_total is True

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


    @pytest.mark.parametrize("r_ref", [0, 10])
    def test_cored_power_law(self, r_ref):
        # Closed form vs 1-D radial quadrature of the r_ref-normalized profile.
        r_c, p = 2.0, 3.0
        tff = CoredPowerLaw(
            params={"r_core": r_c, "index": p, "r_ref": r_ref}
        ).total_flux_factor()
        ref = 1 + (r_ref / r_c) ** 2

        def integrand(r):
            return 2 * np.pi * r * ((1 + (r / r_c) ** 2) / ref) ** (-p / 2)

        expected = quad(integrand, 0, 50)[0] + quad(integrand, 50, np.inf)[0]
        npt.assert_allclose(tff.to_value(u.arcsec**2), expected, rtol=1e-8)

    def test_cored_power_law_closed_form(self):
        # r_ref = 0 (peak): 2 pi r_core^2 / (index - 2) = 2 pi * 4 / 2 = 4 pi.
        tff = CoredPowerLaw(params={"r_core": 2, "index": 4}).total_flux_factor()
        npt.assert_allclose(tff.to_value(u.arcsec**2), 4 * np.pi)

    def test_cored_power_law_unit_and_bare_params_agree(self):
        bare = CoredPowerLaw(
            params={"r_core": 2, "index": 3, "r_ref": 10}
        ).total_flux_factor()
        unit = CoredPowerLaw(
            params={"r_core": 2 * u.arcsec, "index": 3, "r_ref": 1 / 6 * u.arcmin}
        ).total_flux_factor()
        npt.assert_allclose(
            bare.to_value(u.arcsec**2), unit.to_value(u.arcsec**2), rtol=1e-12
        )

    @pytest.mark.parametrize("index", [2.0, 1.0, 0.0])
    def test_cored_power_law_nonfinite_raises_E6(self, index):
        with pytest.raises(BrightnessError) as exc:
            CoredPowerLaw(params={"r_core": 2, "index": index}).total_flux_factor()
        assert exc.value.code == "E6"

    @pytest.mark.parametrize("index", [0.0, 1.0, 2.0, 2.0 + 1e-9, 3.0, 4.5])
    def test_truncated_cored_power_law(self, index):
        # Closed form (incl. the log branch at index 2) vs radial quadrature
        # out to r_out, with the reference point off the peak.
        r_c, r_ref, r_out = 2.0, 5.0, 15.0
        tff = CoredPowerLaw(
            params={"r_core": r_c, "index": index, "r_ref": r_ref, "r_out": r_out}
        ).total_flux_factor()
        ref = 1 + (r_ref / r_c) ** 2

        def integrand(r):
            return 2 * np.pi * r * ((1 + (r / r_c) ** 2) / ref) ** (-index / 2)

        npt.assert_allclose(
            tff.to_value(u.arcsec**2), quad(integrand, 0, r_out)[0], rtol=1e-9
        )

    def test_truncated_index_zero_is_disk(self):
        # Index 0 truncated at r_out is a uniform disk: pi r_out^2.
        tff = CoredPowerLaw(
            params={"r_core": 2, "index": 0, "r_out": 7}
        ).total_flux_factor()
        npt.assert_allclose(tff.to_value(u.arcsec**2), np.pi * 49)

    def test_r_out_beyond_convergence_matches_untruncated(self):
        # index > 2: a very large r_out recovers the untruncated total.
        base = {"r_core": 2, "index": 3, "r_ref": 10}
        inf = CoredPowerLaw(params=base).total_flux_factor()
        big = CoredPowerLaw(params=base | {"r_out": 1e9}).total_flux_factor()
        npt.assert_allclose(big.to_value(u.arcsec**2), inf.to_value(u.arcsec**2), rtol=1e-6)

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


    @pytest.mark.parametrize("reserved", ["amplitude", "x_0", "y_0"])
    def test_cored_power_law_reserved_params_raise_E8(self, reserved):
        with pytest.raises(BrightnessError) as exc:
            CoredPowerLaw(params={"r_core": 2, reserved: 1.0})
        assert exc.value.code == "E8"

    @pytest.mark.parametrize(
        "bad",
        [{"r_core": 0}, {"r_core": -1}, {"index": -0.5}, {"r_ref": -1}],
    )
    def test_cored_power_law_out_of_range_raises(self, bad):
        with pytest.raises(ValueError):
            CoredPowerLaw(params=bad)

    @pytest.mark.parametrize(
        "bad",
        [{"r_out": 0}, {"r_out": -5}, {"r_ref": 20, "r_out": 10}],
        ids=["r_out zero", "r_out negative", "r_ref outside r_out"],
    )
    def test_cored_power_law_bad_truncation_raises(self, bad):
        with pytest.raises(ValueError):
            CoredPowerLaw(params={"r_core": 2} | bad)

    @pytest.mark.parametrize(
        "cls, params",
        [
            (CoredPowerLaw, {"r_core": 0.1 * u.pc}),
            (Sersic, {"r_eff": 0.2 * u.pc, "n": 1}),
        ],
    )
    def test_physical_length_param_raises(self, cls, params):
        # A physical length needs a distance position to resolve
        with pytest.raises(ValueError, match="needs distance"):
            cls(params=params)

class TestWeightMap:
    """Rendering + analytic normalization: no Source, no spectrum, no network."""

    def test_box_weightmap_sums_to_one(self, grid):
        # 6x4 box on a 1 arcsec grid, fully contained and pixel-aligned: the
        # midpoint rule is exact, so the weight map sums to exactly 1 (T-BOX).
        wm, _ = Box(params={"x_width": 6, "y_width": 4})._weightmap(grid)
        npt.assert_allclose(wm.sum(), 1.0, rtol=1e-12)

    def test_bunit_is_dimensionless(self, grid):
        _, hdr = Box(params={"x_width": 6, "y_width": 4})._weightmap(grid)
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

    def test_flat_weightmap_is_uniform(self, grid):
        # Non-integrable profile: A_eff is the field of view, so the weight map
        # is uniform and sums to 1 (the FOV closes the open integral).
        wm, _ = Flat()._weightmap(grid)
        npt.assert_allclose(wm.sum(), 1.0, rtol=1e-12)
        assert np.allclose(wm, wm.flat[0])


    def test_nonuniform_nonintegrable_sums_to_one(self, grid):
        # The FOV closes the integral over the *profile*, A_eff = Omega*sum(p)
        # (flux_scaling.md), not over a constant (Omega*N, Flat only), which
        # would make the map sum to mean(p). With r_ref outside the FOV, p > 1
        # almost everywhere, so the old rule would fail by a large factor.
        prof = CoredPowerLaw(params={"r_core": 3, "index": 2, "r_ref": 40})
        wm, _ = prof._weightmap(grid)
        npt.assert_allclose(wm.sum(), 1.0, rtol=1e-12)
        assert prof._render_image(grid).mean() > 1.5

    def test_cored_power_law_index_zero_is_flat(self, grid):
        prof = CoredPowerLaw(params={"r_core": 1, "index": 0})
        npt.assert_array_equal(prof._weightmap(grid)[0], Flat()._weightmap(grid)[0])
        npt.assert_allclose(
            prof._effective_area(grid).to_value(u.arcsec**2),
            Flat()._effective_area(grid).to_value(u.arcsec**2),
            rtol=1e-15,
        )

    def test_surface_brightness_applies_at_r_ref(self):
        # sb_reference "at r_ref": w_ij * A_eff / Omega is the unit-amplitude
        # profile, which is 1 at r = r_ref (pixel-averaged, hence rtol).
        fine = {"pixel_scale": 0.25 * u.arcsec / u.pixel, "width": 161, "height": 161}
        prof = CoredPowerLaw(params={"r_core": 5, "index": 2, "r_ref": 10})
        wm, _ = prof._weightmap(fine)
        rel = wm * (
            prof._effective_area(fine) / prof._pixel_area(fine)
        ).to_value(u.dimensionless_unscaled)
        c = 80  # odd grid: the central pixel is centered on 0
        npt.assert_allclose(rel[c, c + 40], 1.0, rtol=1e-3)  # 40 px = 10 arcsec
        npt.assert_allclose(rel[c + 40, c], 1.0, rtol=1e-3)

    def test_pixel_flux_independent_of_branch(self, grid, monkeypatch):
        # flux_scaling.md: w_ij = Omega p_ij / A_eff and F = S A_eff, so the
        # flux per pixel (S Omega p_ij) does not depend on A_eff. Forcing a
        # finite profile through the FOV-closed branch changes only the split
        # between weight-map sum and spectrum, never w_ij * A_eff.
        prof = CoredPowerLaw(params={"r_core": 3, "index": 3, "r_ref": 10})
        w_fin, a_fin = prof._weightmap(grid)[0], prof._effective_area(grid)
        monkeypatch.setattr(
            CoredPowerLaw, "has_finite_total", property(lambda self: False)
        )
        w_fov, a_fov = prof._weightmap(grid)[0], prof._effective_area(grid)
        assert w_fin.sum() < 1
        npt.assert_allclose(w_fov.sum(), 1.0, rtol=1e-12)
        npt.assert_allclose(
            (w_fin * a_fin).to_value(u.arcsec**2),
            (w_fov * a_fov).to_value(u.arcsec**2),
            rtol=1e-12,
        )

    def test_contained_truncated_profile_sums_to_one(self, grid):
        # Truncated inside the FOV: the analytic total is all in the window, so
        # the map sums to 1 up to the pixelization of the circular edge (a
        # fully contained Disk behaves the same way). Index 2 would be
        # non-integrable without r_out.
        prof = CoredPowerLaw(params={"r_core": 3, "index": 2, "r_out": 25})
        wm, _ = prof._weightmap(grid)
        npt.assert_allclose(wm.sum(), 1.0, rtol=5e-3)

    def test_truncated_profile_zero_outside_r_out(self, grid):
        prof = CoredPowerLaw(params={"r_core": 3, "index": 2, "r_out": 20})
        img = prof._render_image(grid)
        yy, xx = np.indices(img.shape)
        r = np.hypot(xx - 31.5, yy - 31.5)  # even 64 px grid, 1 arcsec/px
        assert np.all(img[r > 21] == 0)
        assert np.all(img[r < 19] > 0)

class TestSurfaceBrightness:
    """Case II reduction (pure): SB -> implied integrated mag over the effective
    area. Brightness is set directly (bypassing the FilterSystem band check) so
    these stay offline; the reduction itself is synphot-free.
    """

    def test_finite_sb_reduces_to_integrated_mag(self, grid):
        box = Box(
            params={"x_width": 6, "y_width": 4}
        )  # total_flux_factor = 24 arcsec2
        box._brightness = parse_brightness(["V", "21.5 mag / arcsec2"])
        eff = box._effective_integrated_brightness(grid)
        assert eff.solid_angle is None  # now integrated
        npt.assert_allclose(
            eff.value.to_value(u.mag), 21.5 - 2.5 * np.log10(24), rtol=1e-9
        )

    def test_integrated_passes_through_unchanged(self, grid):
        box = Box(params={"x_width": 6, "y_width": 4})
        box._brightness = parse_brightness(["V", "15 mag"])
        assert box._effective_integrated_brightness(grid) is box._brightness

    def test_sb_equivalent_to_matching_integrated(self, grid):
        # A Sersic given SB must reduce to the integrated mag that carries the
        # same total flux (T-B2 at the brightness level).
        ser = Sersic(params={"r_eff": 3, "n": 4})
        tff = ser.total_flux_factor().to_value(u.arcsec**2)
        ser._brightness = parse_brightness(["V", "20 mag / arcsec2"])
        eff = ser._effective_integrated_brightness(grid)
        npt.assert_allclose(
            eff.value.to_value(u.mag), 20 - 2.5 * np.log10(tff), rtol=1e-9
        )

    def test_nonintegrable_sb_uses_fov_area(self, grid):
        # Flat has no analytic total, so A_eff is the FOV: 64x64 @ 1 arcsec/pix
        # -> A_FOV = 4096 arcsec2.
        flat = Flat()
        flat._brightness = parse_brightness(["V", "21.5 mag / arcsec2"])
        eff = flat._effective_integrated_brightness(grid)
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

    def test_nonmagnitude_sb_reduces_to_integrated(self, grid):
        # Linear surface brightness now reduces too: total = SB * A_eff, with
        # the per-solid-angle divisor stripped (Box total_flux_factor = 24
        # arcsec2). Retires the former Jy/sr deferral.
        box = Box(params={"x_width": 6, "y_width": 4})
        box._brightness = parse_brightness(["V", "5 MJy / sr"])
        eff = box._effective_integrated_brightness(grid)
        assert eff.solid_angle is None
        assert eff.amount_kind is AmountKind.FLUX_DENSITY_NU
        npt.assert_allclose(
            eff.value.to_value(u.Jy),
            (5 * u.MJy / u.sr * 24 * u.arcsec**2).to_value(u.Jy),
            rtol=1e-9,
        )

    def test_radiance_sb_reduces_to_energy_flux(self, grid):
        # Per-solid-angle energy flux (radiance) reduces to a band-integrated
        # energy flux over the effective area.
        box = Box(params={"x_width": 6, "y_width": 4})
        box._brightness = parse_brightness(["V", "3e-8 W / (m2 arcsec2)"])
        eff = box._effective_integrated_brightness(grid)
        assert eff.solid_angle is None
        assert eff.amount_kind is AmountKind.ENERGY_FLUX
        npt.assert_allclose(
            eff.value.to_value(u.W / u.m**2),
            (
                3e-8 * u.W / (u.m**2 * u.arcsec**2) * 24 * u.arcsec**2
            ).to_value(u.W / u.m**2),
            rtol=1e-9,
        )

    def test_nonfinite_integrated_raises_E6(self, grid):
        # An *integrated* brightness on a non-integrable profile has no finite
        # total to carry: E6 (raised before any rendering/scaling).
        flat = Flat()
        flat._brightness = parse_brightness(["V", "15 mag"])
        with pytest.raises(BrightnessError) as exc:
            flat.to_source(grid)
        assert exc.value.code == "E6"


    def test_nonuniform_nonintegrable_sb_uses_rendered_integral(self, grid):
        # A_eff = Omega * sum(p) (1 arcsec2 pixels here), not Flat's Omega * N.
        prof = CoredPowerLaw(params={"r_core": 3, "index": 2, "r_ref": 40})
        prof._brightness = parse_brightness(["V", "21.5 mag / arcsec2"])
        eff = prof._effective_integrated_brightness(grid)
        a_eff = prof._render_image(grid).sum()
        npt.assert_allclose(
            eff.value.to_value(u.mag), 21.5 - 2.5 * np.log10(a_eff), rtol=1e-9
        )
        assert not np.isclose(a_eff, 64 * 64)

    def test_cored_power_law_finite_sb_uses_analytic_total(self, grid):
        prof = CoredPowerLaw(params={"r_core": 2, "index": 3, "r_ref": 10})
        tff = prof.total_flux_factor().to_value(u.arcsec**2)
        prof._brightness = parse_brightness(["V", "20 mag / arcsec2"])
        eff = prof._effective_integrated_brightness(grid)
        npt.assert_allclose(
            eff.value.to_value(u.mag), 20 - 2.5 * np.log10(tff), rtol=1e-9
        )

    def test_cored_power_law_nonfinite_integrated_raises_E6(self, grid):
        prof = CoredPowerLaw(params={"r_core": 2, "index": 2})
        prof._brightness = parse_brightness(["V", "15 mag"])
        with pytest.raises(BrightnessError) as exc:
            prof.to_source(grid)
        assert exc.value.code == "E6"

class TestToSource:
    """Full conversion: builds a real ScopeSim Source and scales a spectrum."""

    @pytest.mark.webtest
    def test_returns_source(self, grid):
        box = Box(
            spectrum="G2V",
            brightness=["V", 15],
            params={"x_width": 6, "y_width": 4},
        )
        assert box.to_source(grid) is not None

    @pytest.mark.webtest
    def test_position_guard_without_position(self, grid):
        # Must not raise AttributeError when position was never set.
        box = Box(
            spectrum="G2V",
            brightness=["V", 15],
            params={"x_width": 6, "y_width": 4},
        )
        box.to_source(grid)

    @pytest.mark.webtest
    def test_finite_surface_brightness_ok(self, grid):
        # Finite profile + surface brightness renders via the Case II reduction.
        box = Box(
            spectrum="G2V",
            brightness=["V", "21.5 mag / arcsec2"],
            params={"x_width": 6, "y_width": 4},
        )
        assert box.to_source(grid) is not None

    @pytest.mark.webtest
    def test_nonintegrable_surface_brightness_ok(self, grid):
        # Flat + surface brightness now renders: total flux from the FOV, uniform
        # weight map, flux in the spectrum (no per-pixel image-owns-flux path).
        flat = Flat(spectrum="G2V", brightness=["V", "21.5 mag / arcsec2"])
        assert flat.to_source(grid) is not None


    @pytest.mark.webtest
    def test_cored_power_law_surface_brightness_ok(self, grid):
        prof = CoredPowerLaw(
            spectrum="G2V",
            brightness=["V", "21.5 mag / arcsec2"],
            params={"r_core": 3, "index": 2, "r_ref": 20},
        )
        assert prof.to_source(grid) is not None


class TestCoredPowerLaw2D:
    """The bare astropy model, without the target machinery."""

    def test_amplitude_applies_at_r_ref(self):
        model = CoredPowerLaw2D(amplitude=3.0, r_core=2, index=2.5, r_ref=7)
        npt.assert_allclose(model(7, 0), 3.0)
        npt.assert_allclose(model(0, -7), 3.0)
        npt.assert_allclose(model(7 / np.sqrt(2), 7 / np.sqrt(2)), 3.0)

    def test_power_law_tail(self):
        model = CoredPowerLaw2D(r_core=1, index=3)
        r = np.array([1e3, 1e4])
        slope = np.diff(np.log(model(r, 0))) / np.diff(np.log(r))
        npt.assert_allclose(slope, -3, rtol=1e-5)

    def test_flat_core(self):
        model = CoredPowerLaw2D(r_core=10, index=2)
        npt.assert_allclose(model(0.1, 0), 1.0, rtol=2e-4)

    def test_index_zero_is_constant(self):
        model = CoredPowerLaw2D(amplitude=2.0, r_core=1, index=0, r_ref=5)
        npt.assert_array_equal(model(np.array([0.0, 3.0, 100.0]), 0), 2.0)

    def test_truncation(self):
        # Zero beyond r_out; inside, identical to the untruncated law (the
        # r_ref normalization is not affected by truncation).
        cut = CoredPowerLaw2D(r_core=2, index=2, r_ref=5, r_out=10)
        full = CoredPowerLaw2D(r_core=2, index=2, r_ref=5)
        r = np.array([0.0, 5.0, 9.99])
        npt.assert_array_equal(cut(r, 0), full(r, 0))
        npt.assert_array_equal(cut(np.array([10.01, 50.0]), 0), 0.0)


def _sky(distance=None):
    """Arbitrary sky position, optionally carrying a distance."""
    kwargs = {} if distance is None else {"distance": distance}
    return SkyCoord(10 * u.deg, 20 * u.deg, **kwargs)


def _subtended(length, distance):
    """Angle subtended by a projected length (small-angle)."""
    return (length / distance * u.rad).to(u.arcsec)


class TestPhysicalLengthParams:
    """Radii given as physical lengths, converted with ``position.distance``."""

    @pytest.mark.parametrize("distance", [414 * u.pc, 0.414 * u.kpc])
    def test_pc_radii_match_angular_twin(self, grid, distance):
        # A profile in pc must render exactly like its angular equivalent,
        # independent of the unit the distance is given in.
        phys = CoredPowerLaw(
            params={"r_core": 0.1 * u.pc, "index": 2, "r_ref": 0.05 * u.pc},
            position=_sky(distance),
        )
        twin = CoredPowerLaw(
            params={
                "r_core": _subtended(0.1 * u.pc, 414 * u.pc),
                "index": 2,
                "r_ref": _subtended(0.05 * u.pc, 414 * u.pc),
            }
        )
        npt.assert_allclose(
            phys._weightmap(grid)[0], twin._weightmap(grid)[0], rtol=1e-6
        )
        npt.assert_allclose(
            phys._effective_area(grid).to_value(u.arcsec**2),
            twin._effective_area(grid).to_value(u.arcsec**2),
            rtol=1e-6,
        )

    def test_mixed_length_and_angle_params(self, grid):
        # Each radius is converted on its own; angles pass through unchanged.
        phys = CoredPowerLaw(
            params={"r_core": 0.1 * u.pc, "index": 2, "r_ref": 20 * u.arcsec},
            position=_sky(414 * u.pc),
        )
        twin = CoredPowerLaw(
            params={
                "r_core": _subtended(0.1 * u.pc, 414 * u.pc),
                "index": 2,
                "r_ref": 20 * u.arcsec,
            }
        )
        npt.assert_allclose(
            phys._weightmap(grid)[0], twin._weightmap(grid)[0], rtol=1e-6
        )

    def test_finite_total_scales_as_inverse_distance_squared(self):
        # Pushing a target of fixed physical size away: its solid angle, and
        # hence the total implied by a given surface brightness, falls as d^-2
        # (the surface brightness itself is distance-invariant).
        params = {"r_core": 0.05 * u.pc, "index": 3, "r_ref": 0.2 * u.pc}
        near = CoredPowerLaw(params=params, position=_sky(400 * u.pc))
        far = CoredPowerLaw(params=params, position=_sky(4 * u.kpc))
        ratio = near.total_flux_factor() / far.total_flux_factor()
        npt.assert_allclose(ratio.decompose().value, 100.0, rtol=1e-6)

    def test_king_matched_background_scales_as_inverse_distance_squared(self):
        # All radii in pc (core/tidal radius of the reference cluster), index
        # 2: finite only because of r_out, and still exactly d^-2 when pushed.
        params = {"r_core": 0.13 * u.pc, "index": 2, "r_ref": 0.224 * u.pc,
                  "r_out": 2.0 * u.pc}
        near = CoredPowerLaw(params=params, position=_sky(414 * u.pc))
        far = CoredPowerLaw(params=params, position=_sky(4140 * u.pc))
        ratio = near.total_flux_factor() / far.total_flux_factor()
        npt.assert_allclose(ratio.decompose().value, 100.0, rtol=1e-9)

    @pytest.mark.parametrize(
        "position",
        [None, _sky()],
        ids=["no position", "SkyCoord without distance"],
    )
    def test_pc_radius_without_distance_raises(self, position):
        # Not xfail: holds now and must keep holding. A SkyCoord without a
        # distance reports a *dimensionless* distance of 1, not None, so the
        # implementation must check the unit and raise a clear error naming the
        # distance, rather than fail with an opaque unit-conversion error.
        with pytest.raises(ValueError, match="distance"):
            CoredPowerLaw(params={"r_core": 0.1 * u.pc}, position=position)

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

    @pytest.mark.skip(
        reason="T-B8: needs the extinction attribute/screens, which are 'to be "
        "implemented' (defining_extinction.md). Written in final form; enable "
        "once screens are applied in the normative pipeline."
    )
    def test_observed_vs_intrinsic_surface_brightness_with_screen(self, grid):
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
        s_obs = observed._scale_spectrum(grid)
        s_int = intrinsic._scale_spectrum(grid)

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
