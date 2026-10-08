# -*- coding: utf-8 -*-
"""Unit tests for spectral_classes.py."""

import pytest
import numpy as np
from astropy import units as u
from astropy.table import QTable, Row

from scopesim_targets.spectral_classes import (
    BROWN_DWARF_TEFF_MASS,
    StellarParameters,
    SpectralClass,
    TeffRange,
    teff_range_overlap,
    _end_secants,
    _extrapolation_bounds,
)


class TestStellarParameters:
    @pytest.mark.parametrize(
        ("req", "tbl_len"), ((None, 118), ("J-H", 101), (["J-H", "B-V"], 76))
    )
    def test_required_cols(self, req, tbl_len):
        stp = StellarParameters(req)
        assert len(stp.table) == tbl_len

    @pytest.mark.parametrize(
        ("mass", "desired"),
        (
            (1 * u.solMass, "G2V"),
            ([1, 1] * u.solMass, ("G2V", "G2V")),
            # brown-dwarf masses are filled down to Y4V (~0.024 solMass)
            ([1e4, 1e-4] * u.solMass, ("O3V", "Y4V")),
            (0.05 * u.solMass, "T0V"),
        ),
    )
    def test_closest_mass(self, mass, desired):
        stp = StellarParameters()
        closest = stp.closest_mass(mass)
        np.testing.assert_array_equal(closest["spectral_type"], desired)
        assert isinstance(closest, (Row, QTable))

    @pytest.mark.parametrize(
        ("teff", "desired"),
        (
            (5778 * u.K, "G2V"),
            ([5778, 5778] * u.K, ("G2V", "G2V")),
            ([45000, 2000] * u.K, ("O3V", "L2V")),
        ),
    )
    def test_closest_teff(self, teff, desired):
        stp = StellarParameters()
        closest = stp.closest_teff(teff)
        np.testing.assert_array_equal(closest["spectral_type"], desired)
        assert isinstance(closest, (Row, QTable))

# TODO: Add tests to check if wrong or missing input units throw


@pytest.fixture(scope="module")
def stellar_params():
    return StellarParameters()


class TestBrownDwarfMasses:
    def test_no_masses_missing(self, stellar_params):
        assert not stellar_params.table["mass"].mask.any()

    @pytest.mark.parametrize(
        ("spectral_type", "mass"),
        (("G2V", 1.0), ("M9V", 0.079), ("L2V", 0.075)),
    )
    def test_tabulated_masses_untouched(self, stellar_params, spectral_type, mass):
        assert stellar_params.table.loc[spectral_type]["mass"] == mass * u.solMass

    @pytest.mark.parametrize("spectral_type", ("L3V", "T0V", "Y0V", "Y4V"))
    def test_filled_masses_follow_relation(self, stellar_params, spectral_type):
        row = stellar_params.table.loc[spectral_type]
        expected = BROWN_DWARF_TEFF_MASS(row["teff"]).to(u.solMass)
        assert u.isclose(row["mass"], expected, atol=1e-4 * u.solMass)

    def test_filled_masses_decrease(self, stellar_params):
        # L2V (last tabulated mass) through Y4V
        masses = stellar_params.table["mass"][-29:].to_value(u.solMass)
        assert (np.diff(masses) < 0).all()

    def test_provenance_recorded(self, stellar_params):
        fill = stellar_params.table.meta["mass_fill"]
        assert fill["model"] == BROWN_DWARF_TEFF_MASS.name
        assert len(fill["rows"]) == 28
        assert fill["rows"][0] == "L3V" and fill["rows"][-1] == "Y4V"


class TestInterpolate:
    @staticmethod
    def _masses(stellar_params, *spectral_types):
        return u.Quantity(
            [stellar_params.table.loc[spt]["mass"] for spt in spectral_types]
        )

    def test_reproduces_table_at_knots(self, stellar_params):
        tbl = stellar_params.table
        result = stellar_params.interpolate("teff", tbl["teff"])
        np.testing.assert_allclose(
            result["M_V"].filled(np.nan).value,
            tbl["M_V"].filled(np.nan).value,
            atol=1e-3,
        )

    def test_keeps_column_order(self, stellar_params):
        result = stellar_params.interpolate("mass", [1.0] * u.solMass)
        expected = [
            col for col in stellar_params.table.colnames
            if col not in ("spectral_type", "mass")
        ]
        assert result.colnames == expected

    def test_scalar_input(self, stellar_params):
        result = stellar_params.interpolate("mass", 1 * u.solMass)
        assert len(result) == 1
        assert u.isclose(result["teff"][0], 5770 * u.K, atol=10 * u.K)

    def test_unit_mismatch_raises(self, stellar_params):
        with pytest.raises(ValueError):
            stellar_params.interpolate("mass", [2e30] * u.kg)

    def test_non_positive_input_masked(self, stellar_params):
        result = stellar_params.interpolate("mass", [0, -1, 1] * u.solMass)
        np.testing.assert_array_equal(result["teff"].mask, [True, True, False])

    @pytest.mark.parametrize("colname", ("mass", "teff"))
    @pytest.mark.parametrize(
        "output", ("M_V", "M_J", "M_Ks", "B-V", "V-Ks", "J-H", "H-Ks")
    )
    def test_no_overshoot(self, stellar_params, colname, output):
        tbl = stellar_params.table
        valid = ~tbl[output].mask
        knots = np.sort(tbl[colname][valid].to_value(tbl[colname].unit))
        values = tbl[output][valid][np.argsort(tbl[colname][valid])].value
        # dense sampling of every interval between neighbouring rows
        frac = np.linspace(0, 1, 21)[1:-1]
        search = (knots[:-1, None] + np.diff(knots)[:, None] * frac).ravel()
        result = stellar_params.interpolate(
            colname, search * tbl[colname].unit
        )[output].value.unmasked.reshape(-1, len(frac))
        lower = np.minimum(values[:-1], values[1:])[:, None]
        upper = np.maximum(values[:-1], values[1:])[:, None]
        tol = 1e-3  # output rounding
        assert (result >= lower - tol).all() and (result <= upper + tol).all()

    def test_y_dwarf_mj_regression(self, stellar_params):
        # Former CubicSpline gave M_J ~ 29.9 here, beyond Y4V's 28.2 mag.
        result = stellar_params.interpolate("mass", [0.025] * u.solMass)
        assert 23.6 * u.mag <= result["M_J"][0] <= 28.2 * u.mag

    def test_physical_columns_never_extrapolated(self, stellar_params):
        result = stellar_params.interpolate(
            "mass", [0.01, 100] * u.solMass, extrapolate_phot=True
        )
        assert result["teff"].mask.all()
        assert result["radius"].mask.all()

    def test_masked_outside_data_by_default(self, stellar_params):
        # Y4V has no M_Ks (data ends at Y2V)
        masses = self._masses(stellar_params, "Y2V", "Y4V")
        result = stellar_params.interpolate("mass", masses)
        np.testing.assert_array_equal(result["M_Ks"].mask, [False, True])

    def test_extrapolation_continues_trend(self, stellar_params):
        masses = self._masses(stellar_params, "Y2V", "Y4V")
        m_ks = stellar_params.interpolate(
            "mass", masses, extrapolate_phot=True
        )["M_Ks"]
        assert not m_ks.mask.any()
        # fainter at lower mass, but no runaway (old code: ~31 mag)
        assert m_ks[0] < m_ks[1] < m_ks[0] + 3 * u.mag

    @pytest.mark.parametrize(
        ("steps", "masked"),
        (
            # J-H data ends at T9V; T9.5V, Y0V, Y0.5V are 1, 2, 3 steps beyond
            (None, [False, False, False]),
            (1, [False, True, True]),
            (2, [False, False, True]),
        ),
    )
    def test_max_extrapolation_steps(self, stellar_params, steps, masked):
        masses = self._masses(stellar_params, "T9.5V", "Y0V", "Y0.5V")
        result = stellar_params.interpolate(
            "mass", masses, extrapolate_phot=True,
            max_extrapolation_steps=steps,
        )
        np.testing.assert_array_equal(result["J-H"].mask, masked)

    def test_max_extrapolation_steps_ignored_without_extrapolation(
        self, stellar_params
    ):
        masses = self._masses(stellar_params, "T9.5V")
        result = stellar_params.interpolate(
            "mass", masses, max_extrapolation_steps=3
        )
        assert result["J-H"].mask.all()


class TestExtrapolationHelpers:
    @pytest.mark.parametrize(
        ("steps", "bounds"),
        ((None, (-np.inf, np.inf)), (1, (2, 7)), (2, (1, 8)), (5, (-2, 11))),
    )
    def test_bounds(self, steps, bounds):
        knots = np.arange(10.0)
        valid = np.zeros(10, dtype=bool)
        valid[3:7] = True
        assert _extrapolation_bounds(knots, valid, steps) == bounds

    def test_bounds_beyond_table_use_outer_interval_width(self):
        knots = np.array([0.0, 1.0, 3.0, 7.0])
        valid = np.ones(4, dtype=bool)
        assert _extrapolation_bounds(knots, valid, 2) == (-2.0, 15.0)

    def test_end_secants(self):
        x = np.array([0.0, 0.1, 1.0, 2.0, 2.1])
        y = np.array([0.0, 5.0, 1.0, 2.0, 9.0])
        # spans two intervals, so a short, steep end interval is damped
        assert _end_secants(x, y) == pytest.approx((1.0, 8.0 / 1.1))

    def test_end_secants_two_points(self):
        assert _end_secants(np.array([0.0, 2.0]), np.array([1.0, 5.0])) == (2.0, 2.0)


class TestTeffRange:
    def test_is_still_tuple(self):
        tr = TeffRange(4000, 5000)
        assert isinstance(tr, tuple)

    def test_can_unpack_like_tuple(self):
        teff_min, teff_max = TeffRange(4000, 5000)
        assert teff_min == 4000
        assert teff_max == 5000

    def test_has_named_attrs(self):
        tr = TeffRange(4000, 5000)
        assert tr.min == 4000
        assert tr.max == 5000


class TestTeffRangeOverlap:
    @pytest.mark.parametrize(
        ("teff_range_a", "teff_range_b", "overlap"),
        (
            (TeffRange(400, 600), TeffRange(500, 700), TeffRange(500, 600)),
            (TeffRange(500, 700), TeffRange(400, 600), TeffRange(500, 600)),
            (TeffRange(400, 600), TeffRange(300, 700), TeffRange(400, 600)),
            (TeffRange(400, 800), TeffRange(500, 700), TeffRange(500, 700)),
        ),
    )
    def test_finds_overlap_of_two(self, teff_range_a, teff_range_b, overlap):
        assert teff_range_overlap(teff_range_a, teff_range_b) == overlap

    def test_finds_overlap_of_three(self):
        teff_range_a = TeffRange(400, 900)
        teff_range_b = TeffRange(300, 700)
        teff_range_c = TeffRange(600, 800)
        overlap = teff_range_overlap(teff_range_a, teff_range_b, teff_range_c)
        assert overlap == TeffRange(600, 700)

    def test_throws_for_disjoint(self):
        with pytest.raises(ValueError):
            teff_range_overlap(TeffRange(400, 600), TeffRange(700, 800))


@pytest.fixture(scope="class")
def spectral_class_a():
    return SpectralClass("A", TeffRange(7400, 9700), "#D5E0FF")


class TestSpectralClass:
    def test_instantiate_manually(self, spectral_class_a):
        assert isinstance(spectral_class_a, SpectralClass)
        assert spectral_class_a.name == "A"
        assert spectral_class_a.color == "#D5E0FF"
        assert spectral_class_a.teff_range.min == 7400
        assert spectral_class_a.teff_range.max == 9700

    def test_from_params_table(self, spectral_class_a):
        tbl = StellarParameters().table
        sc = SpectralClass.from_parameters_table(tbl, "A")
        assert sc == SpectralClass("A", TeffRange(7400, 9700), "#D5E0FF")

    def test_from_params_table_grouped(self):
        stp = StellarParameters()
        scs = stp.group_spectral_classes()  # generator
        assert "".join(sc.name for sc in scs) == "OBAFGKMLTY"

    @pytest.mark.parametrize(
        ("teff_range", "result"),
        ((TeffRange(9000, 10000), True), (TeffRange(5000, 6000), False)),
    )
    def test_is_in_range(self, spectral_class_a, teff_range, result):
        assert spectral_class_a.is_in_range(teff_range) is result

    def test_midpoint(self, spectral_class_a):
        assert spectral_class_a.midpoint == 8550.0
