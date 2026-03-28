# -*- coding: utf-8 -*-
"""Unit tests for spectral_classes.py."""

import pytest

from scopesim_targets.spectral_classes import (
    StellarParameters,
    SpectralClass,
    TeffRange,
    teff_range_overlap,
)


class TestStellarParameters:
    @pytest.mark.parametrize(
        ("req", "tbl_len"), ((None, 118), ("J-H", 101), (["J-H", "B-V"], 76))
    )
    def test_required_cols(self, req, tbl_len):
        stp = StellarParameters(req)
        assert len(stp.table) == tbl_len


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
