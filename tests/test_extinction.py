# -*- coding: utf-8 -*-
"""Pure tests for :mod:`scopesim_targets.extinction` (no scopesim/spextra)."""

import numpy as np
from numpy import testing as npt
import pytest
import astropy.units as u
from synphot import SourceSpectrum, SpectralElement, Observation
from synphot.models import Empirical1D, Box1D

from scopesim_targets.extinction import (
    parse_extinction,
    resolve_extinction,
    redden,
    transmission_element,
    ExtinctionScreen,
    ExtinctionLaw,
    FromMap,
    ExtinctionError,
    ExtinctionDistribution,
    parse_extinction_distribution,
)

FLAM = u.Unit("erg / (s cm2 AA)")


@pytest.fixture
def sed():
    wave = np.geomspace(1000, 100000, 6000) * u.AA
    return SourceSpectrum(
        Empirical1D, points=wave,
        lookup_table=(wave.value / 5500.0) ** -1.0 * 1e-14 * FLAM,
    )


class TestParseGrammar:
    def test_bare_av_sugar(self):
        (screen,) = parse_extinction("2.3 mag")
        assert screen == ExtinctionScreen(av=2.3, band="V")

    def test_bare_number_is_av(self):
        (screen,) = parse_extinction(1.5)
        assert screen.av == 1.5 and screen.band == "V"

    def test_mapping_full(self):
        (screen,) = parse_extinction(
            {"value": "8.5 mag", "band": "Ks", "law": "f99", "rv": 4.0}
        )
        assert screen == ExtinctionScreen(
            av=8.5, band="Ks", law=ExtinctionLaw.F99, rv=4.0
        )

    def test_ebv_currency(self):
        (screen,) = parse_extinction({"ebv": 0.3})
        assert screen.ebv == 0.3 and screen.av is None

    def test_list_is_ordered_tuple(self):
        screens = parse_extinction([{"value": "1 mag"}, {"value": "0.5 mag"}])
        assert [s.av for s in screens] == [1.0, 0.5]

    def test_empty_list_is_opt_out(self):
        assert parse_extinction([]) == ()

    def test_from_map_marker(self):
        marker = parse_extinction({"from_map": "bayestar2019"})
        assert isinstance(marker, FromMap) and marker.map_id == "bayestar2019"


class TestParseErrors:
    @pytest.mark.parametrize(
        "spec, code",
        [
            ({"value": "1 mag", "ebv": 0.2}, "E13"),
            ({}, "E14"),
            ("not a mag", "E14"),
            ({"value": "1 mag", "law": "bogus"}, "E15"),
            ({"ebv": -0.1}, "E16"),
            ({"value": "-1 mag"}, "E16"),
            ({"value": "1 mag", "rv": 0}, "E16"),
            ([{"from_map": "x"}, {"value": "1 mag"}], "E13"),
        ],
    )
    def test_error_codes(self, spec, code):
        with pytest.raises(ExtinctionError) as exc:
            parse_extinction(spec)
        assert exc.value.code == code


class TestResolve:
    def test_av_sugar(self):
        assert resolve_extinction(parse_extinction("2.3 mag")) == (
            2.3, ExtinctionLaw.CCM89, 3.1,
        )

    def test_ebv_to_av(self):
        av, law, rv = resolve_extinction(parse_extinction({"ebv": 0.3}))
        assert av == pytest.approx(0.93) and rv == 3.1

    def test_composition_sums_av(self):
        av, _, _ = resolve_extinction(
            parse_extinction([{"value": "1 mag"}, {"value": "0.5 mag"}])
        )
        assert av == pytest.approx(1.5)

    def test_opt_out_and_none(self):
        assert resolve_extinction(parse_extinction([])) is None
        assert resolve_extinction(None) is None

    def test_mixed_law_raises_E17(self):
        screens = parse_extinction(
            [{"value": "1 mag", "law": "ccm89"}, {"value": "1 mag", "law": "f99"}]
        )
        with pytest.raises(ExtinctionError) as exc:
            resolve_extinction(screens)
        assert exc.value.code == "E17"

    def test_from_map_resolve_raises_E18(self):
        with pytest.raises(ExtinctionError) as exc:
            resolve_extinction(parse_extinction({"from_map": "x"}))
        assert exc.value.code == "E18"

    def test_non_v_band_needs_lookup(self):
        screens = parse_extinction({"value": "8.5 mag", "band": "Ks"})
        with pytest.raises(ExtinctionError) as exc:
            resolve_extinction(screens)  # no lookup supplied
        assert exc.value.code == "E14"

    def test_non_v_band_with_lookup(self):
        screens = parse_extinction({"value": "8.5 mag", "band": "Ks"})
        av, _, _ = resolve_extinction(
            screens, band_wavelength_lookup=lambda band: 21600 * u.AA
        )
        # A_Ks = 8.5 -> A_V = 8.5 / (A(Ks)/A(V)) > 8.5 (Ks extinction is small)
        assert av > 8.5


class TestRedden:
    def test_av_zero_is_identity(self, sed):
        assert redden(sed, 0.0, ExtinctionLaw.CCM89, 3.1) is sed

    def test_dims_v_band_by_about_av(self, sed):
        V = SpectralElement(Box1D, amplitude=1, x_0=5500 * u.AA, width=300 * u.AA)
        red = redden(sed, 2.0, ExtinctionLaw.CCM89, 3.1)
        dim = -2.5 * np.log10(
            Observation(red, V).effstim(FLAM).value
            / Observation(sed, V).effstim(FLAM).value
        )
        assert dim == pytest.approx(2.0, abs=0.02)  # ~A_V, not exactly (band eff wav)

    def test_out_of_range_wavelengths_do_not_raise(self):
        # SED spanning far outside CCM89 validity (0.1-3.33 micron): clamp, no raise
        wave = np.geomspace(300, 1_000_000, 5000) * u.AA
        wide = SourceSpectrum(Empirical1D, points=wave,
                              lookup_table=np.ones_like(wave.value) * 1e-14 * FLAM)
        assert redden(wide, 1.0, ExtinctionLaw.CCM89, 3.1) is not None

    @pytest.mark.parametrize("law", list(ExtinctionLaw))
    def test_all_laws_build_a_transmission(self, law):
        wave = np.geomspace(1500, 30000, 1000) * u.AA
        el = transmission_element(wave, 1.0, law, 3.1)
        assert isinstance(el, SpectralElement)


class TestDistributionParse:
    def test_parse_full(self):
        d = parse_extinction_distribution({
            "distribution": "column_density_pdf",
            "params": {"av_median": 1.0, "sigma": 0.4, "av_break": 3.0, "slope": 2.5},
            "law": "f99", "rv": 4.5,
        })
        assert d.kind == "column_density_pdf"
        assert d.law is ExtinctionLaw.F99 and d.rv == 4.5

    @pytest.mark.parametrize(
        "spec, code",
        [
            ({"distribution": "bogus", "params": {"low": 0, "high": 1}}, "E15"),
            ({"distribution": "uniform", "params": {}}, "E14"),
            ({"params": {"low": 0, "high": 1}}, "E14"),          # no distribution key
            ({"distribution": "uniform", "params": {"low": 0, "high": 1}, "rv": 0}, "E16"),
        ],
    )
    def test_parse_errors(self, spec, code):
        with pytest.raises(ExtinctionError) as exc:
            parse_extinction_distribution(spec)
        assert exc.value.code == code


class TestDistributionSample:
    def _rng(self, seed=0):
        return np.random.default_rng(seed)

    def test_uniform_bounds(self):
        d = ExtinctionDistribution("uniform", {"low": 0.5, "high": 1.5})
        av = d.sample(50000, self._rng())
        assert av.min() >= 0.5 and av.max() <= 1.5
        assert av.mean() == pytest.approx(1.0, abs=0.02)

    def test_lognormal_median(self):
        d = ExtinctionDistribution("lognormal", {"av_median": 2.0, "sigma": 0.3})
        av = d.sample(200000, self._rng())
        assert (av > 0).all()
        assert np.median(av) == pytest.approx(2.0, rel=0.02)

    def test_column_density_pdf_core_and_tail(self):
        d = ExtinctionDistribution(
            "column_density_pdf",
            {"av_median": 1.0, "sigma": 0.4, "av_break": 3.0, "slope": 2.5},
        )
        av = d.sample(200000, self._rng())
        assert (av >= 0).all()
        assert np.median(av) == pytest.approx(1.0, rel=0.05)   # core-dominated
        assert (av > 3.0).any() and av.max() > 10               # power-law tail

    def test_reproducible_under_seed(self):
        d = ExtinctionDistribution(
            "column_density_pdf",
            {"av_median": 1.0, "sigma": 0.4, "av_break": 3.0, "slope": 2.5},
        )
        npt.assert_array_equal(d.sample(1000, self._rng(7)), d.sample(1000, self._rng(7)))

    def test_quantity_string_params(self):
        d = ExtinctionDistribution("uniform", {"low": "0 mag", "high": "2 mag"})
        assert d.sample(20000, self._rng()).mean() == pytest.approx(1.0, abs=0.05)

    def test_bad_kind_raises_E15(self):
        with pytest.raises(ExtinctionError) as exc:
            ExtinctionDistribution("nope", {"low": 0, "high": 1})
        assert exc.value.code == "E15"

    def test_tail_slope_must_exceed_one(self):
        d = ExtinctionDistribution(
            "column_density_pdf",
            {"av_median": 1.0, "sigma": 0.4, "av_break": 3.0, "slope": 0.5},
        )
        with pytest.raises(ExtinctionError) as exc:
            d.sample(10, self._rng())
        assert exc.value.code == "E16"
