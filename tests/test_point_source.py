# -*- coding: utf-8 -*-
"""Unit tests for point_source.py."""

from unittest.mock import patch

import pytest
import yaml
import numpy as np
from astropy import units as u

from scopesim_targets.brightness import parse_brightness
from scopesim_targets.point_source import (
    PointSourceTarget,
    Star,
    Binary,
    Exoplanet,
    PlanetarySystem,
    StarField,
)


class TestStar:
    def test_basic(self):
        tgt = Star()
        assert isinstance(tgt, PointSourceTarget)

    # Webtest??
    @pytest.mark.parametrize("position", ((0, 0), (10, -4.2)))
    def test_to_source(self, position):
        # Note: Without any additional info, single point source will be placed
        #       in the center of the FOV for ScopeSim.
        src = Star(
            position=position, spectrum="A0V", brightness=("V", 10)
        ).to_source()
        assert src.fields[0].field["x"] == 0
        assert src.fields[0].field["y"] == 0
        assert src.fields[0].field["ref"][0] in src.fields[0].spectra

    # Webtest??
    def test_loads_yaml(self):
        tgt = yaml.full_load(
            """
            !Star
            position: [2, 3]
            spectrum: A0V
            brightness: ["R", 15 mag]
        """
        )
        assert isinstance(tgt, Star)


@pytest.fixture
def basic_binary():
    return Binary(brightness=("R", 10))


class TestBinary:
    def test_two_brightnesses(self):
        tgt = Binary(brightness=(("R", 10), ("V", 15 * u.mag)))
        assert tgt.brightness == parse_brightness(["R", 10])
        assert tgt.brightness_secondary == parse_brightness(["V", 15])
        with pytest.raises(AttributeError):
            tgt.contrast

    def test_brightness_and_contrast(self):
        tgt = Binary(brightness=("R", 10), contrast=100.0)
        assert tgt.brightness == parse_brightness(["R", 10])
        assert tgt.contrast == 100
        with pytest.raises(AttributeError):
            tgt.brightness_secondary

    def test_invalid_contrast_throws(self):
        with pytest.raises(TypeError):
            Binary(brightness=("R", 10), contrast="bogus")

    def test_two_brightnesses_and_contrast_throws(self):
        with pytest.raises(TypeError):
            Binary(brightness=(("R", 10), ("V", 15 * u.mag)), contrast=100.0)

    def test_invalid_brightness_and_contrast_throws(self):
        with pytest.raises(TypeError):
            Binary(brightness="bogus", contrast="bogus")

    @pytest.mark.parametrize(
        "single",
        (
            ("R", "15 mag"),  # string amount
            ("R", 15 * u.mag),  # Quantity amount
            ("R", 15),  # bare number amount
            ("230 GHz", "5 mJy"),  # frequency string locator
            (230 * u.GHz, "5 mJy"),  # Quantity frequency locator
            {"band": "V", "value": 15},  # canonical mapping form
        ),
    )
    def test_single_spec_not_misrouted_as_pair(self, single):
        # The old (str(), Quantity()|Number()) match silently treated every one
        # of these as two brightnesses (or raised). They must route to a single
        # primary brightness with no secondary.
        tgt = Binary(brightness=single)
        assert tgt.brightness == parse_brightness(single)
        with pytest.raises(AttributeError):
            tgt.brightness_secondary

    @pytest.mark.parametrize(
        "pair",
        (
            (("R", "10 mag"), ("V", "15 mag")),  # string amounts
            (("656.3 nm", "5 mJy"), ("230 GHz", "3 mJy")),  # non-band locators
            (
                {"band": "R", "value": 10},
                {"band": "V", "value": 15},
            ),  # mappings
        ),
    )
    def test_pair_spec_routes_to_two_brightnesses(self, pair):
        tgt = Binary(brightness=pair)
        assert tgt.brightness == parse_brightness(pair[0])
        assert tgt.brightness_secondary == parse_brightness(pair[1])
        with pytest.raises(AttributeError):
            tgt.contrast

    @pytest.mark.parametrize(
        ("brightness", "contrast"),
        (
            (None, None),
            (("R", 10), None),
            (None, 10.0),
        ),
    )
    def test_other_cases(self, brightness, contrast):
        # TODO: Replace this with more meaningful tests!
        tgt = Binary(brightness=brightness, contrast=contrast)
        assert isinstance(tgt, Binary)

    def test_to_source(self):
        # TODO: cover more possible cases
        tgt = Binary(
            brightness=("R", 10),
            contrast=100.0,
            spectra=["A0V", "M2V"],
            offset={"separation": 5 * u.arcsec},
        )
        src = tgt.to_source()
        np.testing.assert_array_equal(src.fields[0].field["x"], [0, 0])
        np.testing.assert_array_equal(src.fields[0].field["y"], [0, 5])

    def test_distance_separtation_resolves(self):
        # 1 AU at 1 pc should produce 1 arcsec separation
        tgt = Binary(
            brightness=("R", 10),
            contrast=100.0,
            spectra=["A0V", "M2V"],
            position={"distance": 1 * u.pc},
            offset={"separation": 1 * u.AU},
        )
        src = tgt.to_source()
        np.testing.assert_array_equal(src.fields[0].field["x"], [0, 0])
        np.testing.assert_array_equal(src.fields[0].field["y"], [0, 1])

    def test_throws_if_no_contrast_or_secondary_brightness(self):
        tgt = Binary(
            spectra=("A0V", "M2V"),
            brightness=("R", 10),
            offset={"separation": 2 * u.arcsec},
        )
        with pytest.raises(ValueError):
            tgt.to_source()

    @pytest.mark.parametrize(
        ("spectra", "refs", "called"),
        (
            (None, None, True),
            (None, 42, True),
            ({0: "foo", 1: "bar"}, None, False),
            ({0: "foo", 1: "bar"}, [0, 1], False),
        ),
    )
    def test_resolve_spectra_refs(self, basic_binary, spectra, refs, called):
        start = refs if isinstance(refs, int) else 0
        expected_spec = {start: "foo", start + 1: "bar"}
        expected_refs = (start, start + 1)

        with patch.object(basic_binary, "source_spectra") as mock_spectra:
            mock_spectra.return_value = expected_spec

            result = basic_binary._resolve_spectra_refs(spectra, refs)

            if called:
                if refs is not None:
                    mock_spectra.assert_called_once_with(refs)
                else:
                    mock_spectra.assert_called_once()
            else:
                mock_spectra.assert_not_called()

        assert result == (expected_spec, expected_refs)

    @pytest.mark.parametrize(
        ("spectra", "refs", "exc", "msg"),
        (
            (
                None,
                [0, 1],
                ValueError,
                "refs sequence must have matching spectra",
            ),
            ({0: "foo"}, [0, 1], ValueError, "not all refs found in spectra"),
            ("bogus", 42, TypeError, "refs and spectra not understood"),
        ),
    )
    def test_resolve_spectra_refs_throws(
        self, basic_binary, spectra, refs, exc, msg
    ):
        with pytest.raises(exc, match=msg):
            basic_binary._resolve_spectra_refs(spectra, refs)


class TestExoplanet:
    def test_basic(self):
        tgt = Exoplanet()
        assert isinstance(tgt, Exoplanet)

    @pytest.mark.webtest
    def test_all_attributes(self):
        # TODO: Replace this with a more meaningful test
        Exoplanet(
            position=(0, 1),
            offset={"separation": 2 * u.arcsec},
            brightness=("K", 23),
            spectrum="spex:irtf/Jupiter",
            contrast=1e3,
        )

    @pytest.mark.webtest
    def test_default_spectrum(self):
        tgt = Exoplanet(offset={"separation": 2 * u.arcsec}, contrast=1e3)
        assert str(tgt.spectrum) == "Spextrum(irtf/Neptune)"


class TestPlanetarySystem:
    def test_to_source(self):
        src = PlanetarySystem(
            position=(0, 0),
            primary=Star(
                spectrum="A0V",
                brightness=("V", 10),
            ),
            components=[
                Exoplanet(contrast=1e5),
            ],
        ).to_source()
        assert len(src.fields[0]) == 2  # primary and one planet


class TestStarField:
    def test_to_source(self):
        src = StarField(
            positions=[(0, 0), (0, 1), (1, 0)],
            spectra=["A0V", "G2V", "A0V"],
            brightnesses=[5 * u.mag, 8 * u.mag, 6 * u.mag],
            band="R",
        ).to_source()
        assert len(src.fields[0]) == 3
        assert len(src.fields[0].spectra) == 2  # two A0V stars share spectrum
        np.testing.assert_array_equal(src.fields[0].field["x"], [0, 0, 1])
        np.testing.assert_array_equal(src.fields[0].field["y"], [0, 1, 0])

    def test_len_mismatch_throws(self):
        tgt = StarField(
            positions=[(0, 0), (0, 1), (1, 0)],
            spectra=["A0V", "G2V", "A0V"],
            brightnesses=[5 * u.mag, 8 * u.mag, 6 * u.mag],
            band="R",
        )
        with pytest.raises(ValueError):
            tgt.positions = [(0, 0), (1, 1)]
        with pytest.raises(ValueError):
            tgt.spectra = ["A0V", "G2V"]
        with pytest.raises(ValueError):
            tgt.brightnesses = [5 * u.mag, 6 * u.mag]

    def test_extinction_kwarg_resolves(self):
        # E(B-V) -> A_V = R_V * E(B-V) needs no bandpass, so this stays offline.
        tgt = StarField(
            positions=[(0, 0)],
            spectra=["A0V"],
            brightnesses=[5 * u.mag],
            band="R",
            extinction={"ebv": 0.3},
        )
        av, law, rv = tgt._extinction_av_meta()
        np.testing.assert_allclose(av, 0.93)
        assert (law, rv) == ("ccm89", 3.1)

    def test_extinction_column_broadcast_and_dedup(self):
        # A single shared screen -> one broadcast Av column + law/rv/anchor meta,
        # and the spectrum dedup is untouched (two A0V stars still share one SED).
        src = StarField(
            positions=[(0, 0), (0, 1), (1, 0)],
            spectra=["A0V", "G2V", "A0V"],
            brightnesses=[5 * u.mag, 8 * u.mag, 6 * u.mag],
            band="R",
            extinction={"ebv": 0.3},
            anchor="intrinsic",
        ).to_source()
        table = src.fields[0].field
        assert len(src.fields[0].spectra) == 2  # dedup intact despite extinction
        np.testing.assert_allclose(table["Av"], [0.93, 0.93, 0.93])
        assert table.meta["extinction_law"] == "ccm89"
        assert table.meta["extinction_rv"] == 3.1
        assert table.meta["anchor"] == "intrinsic"

    def test_no_extinction_has_no_av_column(self):
        src = StarField(
            positions=[(0, 0), (0, 1)],
            spectra=["A0V", "G2V"],
            brightnesses=[5 * u.mag, 8 * u.mag],
            band="R",
        ).to_source()
        table = src.fields[0].field
        assert "Av" not in table.colnames
        assert "extinction_law" not in table.meta
