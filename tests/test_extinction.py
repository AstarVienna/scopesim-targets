# -*- coding: utf-8 -*-
"""Pure tests for :mod:`scopesim_targets.extinction` (no scopesim/spextra).

``extinction`` is the single foreground screen; dust at the source lives in
:mod:`scopesim_targets.attenuation` (see ``test_attenuation.py``).
"""

import numpy as np
import pytest
import astropy.units as u
from synphot import SourceSpectrum, SpectralElement, Observation
from synphot.models import Empirical1D, Box1D

from scopesim_targets.extinction import (
    parse_extinction,
    resolve_extinction,
    redden,
    transmission_element,
    as_plain_mag,
    ExtinctionScreen,
    ExtinctionLaw,
    FromMap,
    ExtinctionError,
    UnphysExtError,
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
        assert parse_extinction("2.3 mag") == ExtinctionScreen(av=2.3, band="V")

    def test_bare_number_is_av(self):
        screen = parse_extinction(1.5)
        assert screen.av == 1.5 and screen.band == "V"

    def test_zero_mag_is_explicit_opt_out_screen(self):
        assert parse_extinction("0 mag") == ExtinctionScreen(av=0.0)

    def test_mapping_full(self):
        screen = parse_extinction(
            {"value": "8.5 mag", "band": "Ks", "law": "f99", "rv": 4.0}
        )
        assert screen == ExtinctionScreen(
            av=8.5, band="Ks", law=ExtinctionLaw.F99, rv=4.0
        )

    def test_ebv_currency(self):
        screen = parse_extinction({"ebv": 0.3})
        assert screen.ebv == 0.3 and screen.av is None

    def test_from_map_marker(self):
        marker = parse_extinction({"from_map": "bayestar2019"})
        assert isinstance(marker, FromMap) and marker.map_id == "bayestar2019"

    @pytest.mark.parametrize("spec", [[], ["1 mag"], [{"value": "1 mag"}, {"value": "0.5 mag"}], ()])
    def test_list_form_is_rejected_E17(self, spec):
        with pytest.raises(ExtinctionError):
            parse_extinction(spec)


class TestParseErrors:
    @pytest.mark.parametrize(
        "spec",
        [
            ({"value": "1 mag", "ebv": 0.2}),
            ({"ebv": 0.2, "from_map": "x"}),
            ({}),
            ("not a mag"),
            ({"value": "1 mag", "law": "bogus"}),
            ({"ebv": -0.1}),
            ({"value": "-1 mag"}),
            ({"value": "1 mag", "rv": 0}),
        ],
    )
    def test_error_codes(self, spec):
        with pytest.raises(ExtinctionError):
            parse_extinction(spec)


class TestResolve:
    def test_av_sugar(self):
        assert resolve_extinction(parse_extinction("2.3 mag")) == (
            2.3, ExtinctionLaw.CCM89, 3.1,
        )

    def test_ebv_to_av(self):
        av, law, rv = resolve_extinction(parse_extinction({"ebv": 0.3}))
        assert av == pytest.approx(0.93) and rv == 3.1

    def test_opt_out_is_zero_triple_and_none_is_none(self):
        # '0 mag' is an explicit zero screen (matters for the cascade); None is
        # "nothing set". Both leave the SED alone in redden().
        assert resolve_extinction(parse_extinction("0 mag")) == (0.0, ExtinctionLaw.CCM89, 3.1)
        assert resolve_extinction(None) is None

    def test_from_map_resolve_raises_E18(self):
        with pytest.raises(ExtinctionError):
            resolve_extinction(parse_extinction({"from_map": "x"}))

    def test_non_v_band_needs_lookup(self):
        screen = parse_extinction({"value": "8.5 mag", "band": "Ks"})
        with pytest.raises(ExtinctionError):
            resolve_extinction(screen)  # no lookup supplied

    def test_non_v_band_with_lookup(self):
        screen = parse_extinction({"value": "8.5 mag", "band": "Ks"})
        av, _, _ = resolve_extinction(
            screen, band_wavelength_lookup=lambda band: 21600 * u.AA
        )
        # A_Ks = 8.5 -> A_V = 8.5 / (A(Ks)/A(V)) > 8.5 (Ks extinction is small)
        assert av > 8.5

    def test_two_screens_compose_into_one(self):
        # The rationale for dropping the list form: for one law/R_V, A(lambda)
        # is linear in A_V, so screens 1 + 0.5 == one 1.5 mag screen exactly.
        wave = np.geomspace(1500, 30000, 500) * u.AA
        t1 = transmission_element(wave, 1.0, ExtinctionLaw.CCM89, 3.1)(wave)
        t2 = transmission_element(wave, 0.5, ExtinctionLaw.CCM89, 3.1)(wave)
        t12 = transmission_element(wave, 1.5, ExtinctionLaw.CCM89, 3.1)(wave)
        np.testing.assert_allclose(t1 * t2, t12, rtol=1e-12)


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
        # SED spanning far outside CCM89 validity (0.1-3.33 um): clamp, no raise
        wave = np.geomspace(300, 1_000_000, 5000) * u.AA
        wide = SourceSpectrum(Empirical1D, points=wave,
                              lookup_table=np.ones_like(wave.value) * 1e-14 * FLAM)
        assert redden(wide, 1.0, ExtinctionLaw.CCM89, 3.1) is not None

    @pytest.mark.parametrize("law", list(ExtinctionLaw))
    def test_all_laws_build_a_transmission(self, law):
        wave = np.geomspace(1500, 30000, 1000) * u.AA
        el = transmission_element(wave, 1.0, law, 3.1)
        assert isinstance(el, SpectralElement)


class TestAsPlainMag:
    def test_accepts_number_and_mag_string(self):
        assert as_plain_mag(1.5, "x") == 1.5
        assert as_plain_mag("2 mag", "x") == 2.0

    @pytest.mark.parametrize("bad", ["2 Jy", "abc", True, None])
    def test_rejects_others_E16(self, bad):
        with pytest.raises(UnphysExtError):
            as_plain_mag(bad, "x")
