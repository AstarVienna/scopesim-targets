# -*- coding: utf-8 -*-
"""Unit tests for point_source.py."""

from unittest.mock import patch

import pytest
import yaml
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from synphot import SpectralElement
from synphot.models import Box1D
from spextra import Spextrum

from scopesim_targets import target as target_module

from scopesim_targets.brightness import parse_brightness, BrightnessError
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


@pytest.fixture
def offline_field(monkeypatch):
    """Make StarField.to_source run offline against real synphot.

    Library spectra are replaced by flat stand-ins (distinct amplitudes, so
    refs are distinguishable), and the network-backed passband lookup by a
    synthetic V-ish boxcar. Returns the name -> stand-in spectrum mapping.
    """
    flat = Spextrum.flat_spectrum(
        5 * u.ABmag, waves=np.linspace(3000, 25000, 600) * u.AA
    )
    stand_ins = {"A0V": flat, "G2V": flat * 2, "M2V": flat * 5}
    monkeypatch.setattr(
        StarField,
        "resolve_spectrum",
        staticmethod(lambda spectrum, brightness=None: stand_ins[str(spectrum)]),
    )
    box = SpectralElement(Box1D, amplitude=1, x_0=5500 * u.AA, width=1000 * u.AA)
    monkeypatch.setattr(target_module, "_passband", lambda band: box)
    return stand_ins


class TestStarFieldBrightnesses:
    """Per-star brightness entries: full specs, or bare amounts + `band`."""

    @pytest.mark.parametrize(
        ("entry", "expected"),
        (
            (15, ("V", 15)),  # bare number
            (np.float64(15), ("V", 15)),  # from iterating a numpy array
            (15 * u.mag, ("V", 15)),  # bare Quantity
            ("15 mag", ("V", "15 mag")),  # bare string (was TypeError)
            ("15 mag(AB)", ("V", "15 mag(AB)")),
            ("3.5 mJy", ("V", "3.5 mJy")),
            (("R", "15 mag"), ("R", "15 mag")),  # full spec overrides band
            ({"band": "K", "value": 12}, ("K", 12)),  # mapping (was misrouted)
            (("230 GHz", "5 mJy"), ("230 GHz", "5 mJy")),
        ),
    )
    def test_entry_forms(self, entry, expected):
        tgt = StarField(
            positions=[(0, 0)], spectra=["A0V"], brightnesses=[entry], band="V"
        )
        assert tgt.brightnesses == [parse_brightness(expected)]

    def test_numpy_array_of_mags(self):
        tgt = StarField(
            positions=[(0, 0), (1, 0), (2, 0)],
            spectra=["A0V"] * 3,
            brightnesses=np.linspace(10, 12, 3),
            band="V",
        )
        assert tgt.brightnesses == [
            parse_brightness(("V", m)) for m in (10, 11, 12)
        ]

    def test_bare_amount_without_band_raises_E2(self):
        with pytest.raises(BrightnessError, match="needs a locator") as exc:
            StarField(positions=[(0, 0)], spectra=["A0V"], brightnesses=[15])
        assert exc.value.code == "E2"

    def test_full_specs_need_no_band(self):
        tgt = StarField(
            positions=[(0, 0)], spectra=["A0V"], brightnesses=[("R", 15)]
        )
        assert tgt.brightnesses == [parse_brightness(("R", 15))]

    def test_from_spectral_type_rejected(self):
        with pytest.raises(TypeError, match="from_spectral_type"):
            StarField(
                positions=[(0, 0)],
                spectra=["A0V"],
                brightnesses=[{"from_spectral_type": "mamajek"}],
            )


class TestStarFieldToSource:
    def test_refs_in_first_seen_order(self, offline_field):
        src = StarField(
            positions=[(0, 0), (1, 0), (2, 0), (3, 0)],
            spectra=["G2V", "A0V", "G2V", "M2V"],
            brightnesses=["15 mag(AB)"] * 4,
            band="V",
        ).to_source()
        field = src.fields[0]
        np.testing.assert_array_equal(field.field["ref"], [0, 1, 0, 2])
        assert field.spectra[0] is offline_field["G2V"]
        assert field.spectra[1] is offline_field["A0V"]
        assert field.spectra[2] is offline_field["M2V"]

    def test_batched_weights_match_per_star(self, offline_field):
        # Mixed groups: two spectra x (AB mag, ST mag, band flux density,
        # monochromatic flux density), several values each.
        spectra = ["A0V", "G2V"] * 6
        brightnesses = [
            "14 mag(AB)", "15 mag(AB)", "16.5 mag(AB)",
            "15 mag(ST)", "12 mag(ST)",
            "3 mJy", "0.5 Jy",
            ("656.3 nm", "5 mJy"), ("656.3 nm", "80 mJy"),
            ("R", "13 mag(AB)"), ("R", "18 mag(AB)"), "13 mag(AB)",
        ]
        tgt = StarField(
            positions=[(i, 0) for i in range(len(spectra))],
            spectra=spectra,
            brightnesses=brightnesses,
            band="V",
        )
        weights = tgt.to_source().fields[0].field["weight"]
        expected = [
            tgt._anchored_spectrum_scale(offline_field[str(spec)], bright)
            for spec, bright in zip(tgt.spectra, tgt.brightnesses)
        ]
        np.testing.assert_allclose(weights, expected, rtol=1e-10)

    def test_photometry_once_per_group(self, offline_field, monkeypatch):
        calls = []

        def spy(spectrum, brightness):
            calls.append(brightness)
            return 1.0

        monkeypatch.setattr(StarField, "_get_spectrum_scale", staticmethod(spy))
        n_stars = 50
        StarField(
            positions=[(i, 0) for i in range(n_stars)],
            spectra=["A0V", "G2V"] * (n_stars // 2),
            brightnesses=np.linspace(10, 20, n_stars),
            band="V",
        ).to_source()
        # Vega mags would need the (network) Vega reference, so the spy stubs
        # the photometry: two spectra x one (band, system) group -> 2 calls.
        assert len(calls) == 2


class TestStarFieldPositions:
    """Vectorized per-star offsets from the field center."""

    EXPECTED = [[0, 0], [0, 1], [60, 0]]  # arcsec

    @pytest.mark.parametrize(
        "positions",
        (
            [(0, 0), (0, 1), (60, 0)],  # plain numbers [arcsec]
            np.array([(0, 0), (0, 1), (60, 0)]),
            np.array([(0, 0), (0, 1), (60, 0)]) * u.arcsec,
            np.array([(0, 0), (0, 1 / 60), (1, 0)]) * u.arcmin,
            [  # rows of angular Quantities: units converted, never stripped
                u.Quantity([0, 0], u.arcsec),
                u.Quantity([0, 1 / 60], u.arcmin),
                u.Quantity([1, 0], u.arcmin),
            ],
            [(0 * u.arcsec, 0 * u.deg), (0, 1), (1 * u.arcmin, 0)],  # mixed
            [{"x": 0, "y": 0}, {"x": 0, "y": 1}, {"x": 60, "y": 0}],
        ),
    )
    def test_input_forms(self, positions):
        tgt = StarField(positions=positions)
        assert tgt.offsets.shape == (3, 2)
        np.testing.assert_allclose(
            tgt.offsets.to_value(u.arcsec), self.EXPECTED, atol=1e-12
        )

    @pytest.mark.parametrize(
        ("positions", "exc"),
        (
            ([(1 * u.AU, 0 * u.AU)], u.UnitConversionError),
            ([(1, 2, 3)], ValueError),
            ([{"x": 1, "y": 2, "distance": 3 * u.pc}], ValueError),
            ("bogus", TypeError),
            ([], ValueError),
        ),
    )
    def test_invalid_forms(self, positions, exc):
        with pytest.raises(exc):
            StarField(positions=positions)

    def test_empty_field_constructs(self):
        tgt = StarField()
        assert tgt.offsets is None
        assert tgt.positions is None
        assert tgt.spectra is None

    def test_absolute_positions_around_default_center(self):
        tgt = StarField(positions=[(3, 4)])
        center = SkyCoord(0 * u.deg, 0 * u.deg)
        np.testing.assert_allclose(center.separation(tgt.positions).arcsec, 5)

    def test_absolute_positions_around_field_center(self):
        center = SkyCoord(150.1 * u.deg, 2.2 * u.deg)
        tgt = StarField(position=center, positions=[(0, 0), (0, 10), (3, 4)])
        positions = tgt.positions
        np.testing.assert_allclose(
            center.separation(positions).arcsec, [0, 10, 5], atol=1e-9
        )
        # +y is North, +x East
        np.testing.assert_allclose(
            center.position_angle(positions[1:]).deg,
            [0, np.degrees(np.arctan2(3, 4))],
            atol=1e-6,
        )

    def test_skycoord_input_round_trips(self):
        center = SkyCoord(150.1 * u.deg, 2.2 * u.deg)
        stars = SkyCoord([150.1, 150.101, 150.0995] * u.deg, [2.2, 2.2, 2.201] * u.deg)
        tgt = StarField(position=center, positions=stars)
        np.testing.assert_allclose(
            tgt.positions.separation(stars).to_value(u.arcsec), 0, atol=1e-9
        )
        # list of scalar SkyCoords is equivalent
        again = StarField(position=center, positions=list(stars))
        np.testing.assert_allclose(again.offsets, tgt.offsets)

    def test_moving_center_moves_field_rigidly(self):
        tgt = StarField(positions=[(0, 0), (0, 10)])
        tgt.position = SkyCoord(150.1 * u.deg, 2.2 * u.deg)
        np.testing.assert_allclose(tgt.offsets.to_value(u.arcsec), [[0, 0], [0, 10]])
        np.testing.assert_allclose(
            tgt.positions[0].separation(tgt.positions[1]).arcsec, 10
        )

    def test_center_must_be_scalar(self):
        with pytest.raises(ValueError, match="single field center"):
            StarField(position=[(0, 0), (1, 1)])

    def test_center_distance_is_field_distance(self):
        tgt = StarField(position={"distance": 50 * u.pc}, positions=[(0, 0)])
        assert tgt._distance_or_none() == 50 * u.pc


class TestStarFieldSharedSpectrum:
    def test_single_spectrum_broadcasts(self):
        tgt = StarField(
            positions=[(0, 0), (1, 0), (2, 0)], spectra="A0V",
            brightnesses=[10, 11, 12], band="V",
        )
        assert len(tgt.spectra) == 3
        assert len(set(tgt.spectra)) == 1
        assert tgt.spectra[0] == "a0v"

    def test_positions_can_change_length_with_shared_spectrum(self):
        tgt = StarField(positions=[(0, 0)], spectra="A0V")
        tgt.positions = [(0, 0), (1, 0), (2, 0)]
        assert len(tgt.spectra) == 3

    def test_brightnesses_still_length_checked(self):
        tgt = StarField(positions=[(0, 0), (1, 0)], spectra="A0V")
        with pytest.raises(ValueError, match="Brightnesses length"):
            tgt.brightnesses = [10, 11, 12]

    def test_switch_between_shared_and_per_star(self):
        tgt = StarField(positions=[(0, 0), (1, 0)], spectra="A0V")
        tgt.spectra = ["A0V", "G2V"]
        assert tgt.spectra == ["a0v", "g2v"]
        with pytest.raises(ValueError, match="Spectra length"):
            tgt.spectra = ["A0V"]
        tgt.spectra = "G2V"
        assert tgt.spectra == ["g2v", "g2v"]

    def test_to_source_single_spectrum(self, offline_field):
        src = StarField(
            positions=[(0, 0), (1, 0), (2, 0)], spectra="G2V",
            brightnesses=["15 mag(AB)"] * 3, band="V",
        ).to_source()
        field = src.fields[0]
        np.testing.assert_array_equal(field.field["ref"], [0, 0, 0])
        assert list(field.spectra) == [0]
        assert field.spectra[0] is offline_field["G2V"]


class TestStarFieldTable:
    def test_xy_are_offsets_independent_of_center(self, offline_field):
        positions = [(0, 0), (-12.5, 3), (7, -0.25)]
        src = StarField(
            position=SkyCoord(150.1 * u.deg, 2.2 * u.deg),
            positions=positions, spectra="A0V",
            brightnesses=["15 mag(AB)"] * 3, band="V",
        ).to_source()
        np.testing.assert_array_equal(
            src.fields[0].field["x"], [p[0] for p in positions]
        )
        np.testing.assert_array_equal(
            src.fields[0].field["y"], [p[1] for p in positions]
        )
        assert src.fields[0].field["x"].unit == u.arcsec


class TestStarFieldFromGrid:
    def test_square_geometry_and_order(self):
        tgt = StarField.from_grid(3, 2, spectra="A0V", brightnesses=15, band="V")
        # x runs fastest, starting bottom-left (smallest x and y)
        np.testing.assert_array_equal(
            tgt.offsets.to_value(u.arcsec),
            [[-2, -2], [0, -2], [2, -2],
             [-2, 0], [0, 0], [2, 0],
             [-2, 2], [0, 2], [2, 2]],
        )

    def test_shape_is_width_height(self):
        tgt = StarField.from_grid((4, 2), 1, spectra="A0V", brightnesses=15, band="V")
        x, y = tgt.offsets.to_value(u.arcsec).T
        np.testing.assert_array_equal(np.unique(x), [-1.5, -0.5, 0.5, 1.5])
        np.testing.assert_array_equal(np.unique(y), [-0.5, 0.5])

    def test_even_grid_has_no_center_star(self):
        tgt = StarField.from_grid(2, 1, spectra="A0V", brightnesses=15, band="V")
        assert not np.any(np.all(tgt.offsets.value == 0, axis=1))
        np.testing.assert_allclose(tgt.offsets.value.mean(axis=0), [0, 0])

    def test_spacing_units(self):
        tgt = StarField.from_grid(2, 1 * u.arcmin, spectra="A0V", brightnesses=15, band="V")
        np.testing.assert_array_equal(np.unique(tgt.offsets.to_value(u.arcsec)), [-30, 30])

    @pytest.mark.parametrize(
        ("shape", "exc"),
        ((0, ValueError), ((3, 0), ValueError), (2.5, TypeError),
         (True, TypeError), ((1, 2, 3), TypeError), ((2, 2.0), TypeError)),
    )
    def test_invalid_shape(self, shape, exc):
        with pytest.raises(exc):
            StarField.from_grid(shape, 1, spectra="A0V", brightnesses=15, band="V")

    @pytest.mark.parametrize(
        ("spacing", "exc"),
        ((0, ValueError), (-1, ValueError), ([1, 2], ValueError),
         (1 * u.pc, u.UnitConversionError)),
    )
    def test_invalid_spacing(self, spacing, exc):
        with pytest.raises(exc):
            StarField.from_grid(2, spacing, spectra="A0V", brightnesses=15, band="V")

    @pytest.mark.parametrize(
        "shared",
        (15, "15 mag", 15 * u.mag, ("V", 15), {"band": "V", "value": 15}),
    )
    def test_shared_brightness(self, shared):
        # ("V", 15) on a 2-star grid: one pair, not two per-star amounts
        tgt = StarField.from_grid((2, 1), 1, spectra="A0V", brightnesses=shared, band="V")
        expected = parse_brightness(("V", 15))
        assert tgt.brightnesses == [expected, expected]

    @pytest.mark.parametrize(
        "per_star",
        ([15, 16], ["15 mag", "16 mag"], np.array([15, 16]), [15, 16] * u.mag,
         [("V", 15), ("V", 16)]),
    )
    def test_per_star_brightness(self, per_star):
        tgt = StarField.from_grid((2, 1), 1, spectra="A0V", brightnesses=per_star, band="V")
        assert tgt.brightnesses == [parse_brightness(("V", m)) for m in (15, 16)]

    def test_brightness_ramp_follows_order(self):
        tgt = StarField.from_grid(
            (3, 2), 1, spectra="A0V", brightnesses=np.linspace(10, 15, 6), band="V"
        )
        # first element bottom-left, last top-right
        assert tgt.brightnesses[0] == parse_brightness(("V", 10))
        np.testing.assert_array_equal(tgt.offsets[0].value, [-1, -0.5])
        np.testing.assert_array_equal(tgt.offsets[-1].value, [1, 0.5])

    @pytest.mark.parametrize(
        ("name", "value"),
        (("brightnesses", np.arange(4).reshape(2, 2)),
         ("brightnesses", [[15, 16], [17, 18]]),
         ("spectra", [["A0V", "G2V"], ["A0V", "G2V"]]),
         ("spectra", np.array([["A0V", "G2V"], ["A0V", "G2V"]]))),
    )
    def test_2d_rejected(self, name, value):
        kwargs = {"spectra": "A0V", "brightnesses": 15, name: value}
        with pytest.raises(ValueError, match="flat sequence"):
            StarField.from_grid(2, 1, band="V", **kwargs)

    @pytest.mark.parametrize("name", ("spectra", "brightnesses"))
    def test_wrong_length(self, name):
        kwargs = {"spectra": "A0V", "brightnesses": 15}
        kwargs[name] = ["A0V"] * 3 if name == "spectra" else [15] * 3
        with pytest.raises(ValueError, match="expected 4 per-star values"):
            StarField.from_grid(2, 1, band="V", **kwargs)

    def test_per_star_spectra(self):
        tgt = StarField.from_grid(
            (2, 1), 1, spectra=np.array(["A0V", "G2V"]), brightnesses=15, band="V"
        )
        assert tgt.spectra == ["a0v", "g2v"]

    def test_field_center(self):
        center = SkyCoord(150.1 * u.deg, 2.2 * u.deg)
        tgt = StarField.from_grid(
            1, 1, spectra="A0V", brightnesses=15, band="V", position=center
        )
        assert tgt.positions[0].separation(center).arcsec < 1e-9

    def test_to_source(self, offline_field):
        src = StarField.from_grid(
            (3, 2), 5, spectra="G2V",
            brightnesses=[f"{m} mag(AB)" for m in (15, 15, 15, 17.5, 17.5, 17.5)],
            band="V",
        ).to_source()
        field = src.fields[0].field
        np.testing.assert_array_equal(field["x"], [-5, 0, 5] * 2)
        np.testing.assert_array_equal(field["y"], [-2.5] * 3 + [2.5] * 3)
        np.testing.assert_allclose(field["weight"][3:] / field["weight"][:3], 0.1)
