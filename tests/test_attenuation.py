# -*- coding: utf-8 -*-
"""Unit tests for attenuation.py."""

import importlib.util

import numpy as np, pytest, astropy.units as u
from astropy.table import Table
from synphot import SourceSpectrum
from synphot.models import ConstFlux1D, Empirical1D

from scopesim_targets.extinction import (
    parse_extinction,
    resolve_extinction,
    ExtinctionScreen,
    ExtinctionError,
    ExtinctionLaw,
    redden,
)
from scopesim_targets.attenuation import (
    parse_attenuation,
    Attenuation,
    AttenuationModel,
    AttenuationError,
    EmptyAttenError,
    UnphysAttenError,
    DustColumn,
    ColumnDistribution,
    apply_attenuation_to_table,
    ATTENUATION_COLUMN,
)

# wg00 evaluates dust_attenuation (optional, GitHub-only); everything else,
# including wg00 parsing and tau_v-currency handling, works without it.
requires_dust_attenuation = pytest.mark.skipif(
    importlib.util.find_spec("dust_attenuation") is None,
    reason="optional dependency: pip install git+https://github.com/karllark/dust_attenuation.git",
)

KS = lambda band: {"Ks": 2.15 * u.um, "V": 0.551 * u.um, "B": 0.445 * u.um}[band]

# ---------------- extinction: single screen ----------------
def test_extinction_list_removed():
    with pytest.raises(ExtinctionError):
        parse_extinction(["0.3 mag", "1 mag"])
    with pytest.raises(ExtinctionError):
        parse_extinction([])

def test_extinction_zero_mag_opt_out_is_explicit_screen():
    s = parse_extinction("0 mag")
    assert s == ExtinctionScreen(av=0.0)
    assert resolve_extinction(s) == (0.0, ExtinctionLaw.CCM89, 3.1)
    assert resolve_extinction(None) is None
    sp = SourceSpectrum(ConstFlux1D, amplitude=1 * u.mJy)
    assert redden(sp, 0.0, ExtinctionLaw.CCM89, 3.1) is sp

def test_extinction_band_screen_resolves():
    av, law, rv = resolve_extinction(
        parse_extinction({"value": "8.5 mag", "band": "Ks", "rv": 4.7}),
        band_wavelength_lookup=KS,
    )
    assert law is ExtinctionLaw.CCM89
    assert rv == 4.7 and 60 < av < 80  # A_Ks/A_V ~ 0.12

# ---------------- attenuation: parsing ----------------
def test_sugar_is_screen_av():
    a = parse_attenuation("8.5 mag")
    assert a.model is AttenuationModel.SCREEN
    assert a.column == DustColumn("mag", 8.5)
    assert not a.is_sampled
    assert a.native_amount() == 8.5
    assert a.native_parameter == "Av"

@pytest.mark.parametrize("spec", [
    ({"value": "2 mag", "tau_v": 2}),
    ({"model": "screen", "law": "ccm89"}),
    ({"model": "wg00", "dust_distribution": "clumpy", "tau_v": 2}),
    ({"model": "wg00", "geometry": "shell", "dust_distribution": "clumpy", "tau_v": 2, "law": "g23"}),
    ({"geometry": "shell", "value": "2 mag"}),
    ({"model": "c00", "value": "1 mag"}),
    ({"value": "1 mag", "params": {"a": 1}}),
    ({"ebv": 0.1, "band": "Ks"}),
    ({"value": "-1 mag"}),
    ({"model": "wg00", "geometry": "shell", "dust_distribution": "clumpy", "tau_v": 60}),
    ({"bogus": 1, "value": "1 mag"}),
    (["1 mag"]),
])
def test_parse_errors(spec):
    with pytest.raises(AttenuationError):
        parse_attenuation(spec).native_amount()

def test_attenuation_error_is_extinction_error():
    assert issubclass(AttenuationError, ExtinctionError)

# ---------------- column -> native ----------------
def test_screen_column_conversions():
    assert parse_attenuation({"ebv": 0.5, "rv": 4.0}).native_amount() == pytest.approx(2.0)
    assert parse_attenuation({"tau_v": 1.0}).native_amount() == pytest.approx(1.0857, abs=1e-3)
    a_ks = parse_attenuation({"value": "1 mag", "band": "Ks"})
    ratio = float(ExtinctionLaw.CCM89.model(3.1)(KS("Ks")))
    assert a_ks.native_amount(band_wavelength_lookup=KS) == pytest.approx(1 / ratio)
    with pytest.raises(EmptyAttenError):
        a_ks.native_amount()

@requires_dust_attenuation
def test_wg00_column_conversions_use_models_own_curve():
    base = {
        "model": "wg00",
        "geometry": "shell",
        "dust_type": "mw",
        "dust_distribution": "homogeneous",
    }
    assert parse_attenuation({**base, "tau_v": 2.0}).native_amount() == 2.0
    a = parse_attenuation({**base, "value": "0.97 mag", "band": "Ks"})
    tau = a.native_amount(band_wavelength_lookup=KS)
    from dust_attenuation.radiative_transfer import WG00
    m = WG00(tau_V=1.0, geometry="shell", dust_type="mw", dust_distribution="homogeneous")
    assert m.get_extinction(KS("Ks"), tau) == pytest.approx(0.97, rel=1e-6)  # round trip
    # V and E(B-V) work without a lookup (built-in fallbacks)
    tau_v = parse_attenuation({**base, "value": "1.109 mag"}).native_amount()
    assert tau_v == pytest.approx(1.0, rel=2e-3)
    ebv = parse_attenuation({**base, "ebv": 0.3}).native_amount()
    assert ebv > 0
    assert a.native_parameter == "tau_V"

# ---------------- curves ----------------
WAVE = np.linspace(0.3, 5.0, 300) * u.um

def test_screen_curve_is_law_times_av_and_clamps():
    a = parse_attenuation({"value": "2 mag", "rv": 4.5})
    c = a.curve(WAVE, 2.0)
    assert np.isclose(c[np.argmin(abs(WAVE.value - 0.551))], 2.0, atol=0.02)
    assert np.all(np.isfinite(c)) and c[-1] == c[-2]  # held constant beyond validity
    assert np.all(a.curve(WAVE, 0.0) == 0)

@requires_dust_attenuation
def test_wg00_curve_geometry_and_thin_limit():
    a = parse_attenuation({
        "model": "wg00",
        "geometry": "shell",
        "dust_distribution": "homogeneous",
        "tau_v": 2.0,
    })
    c = a.curve(WAVE, 2.0)
    assert 1.3 < c[np.argmin(abs(WAVE.value - 0.55))] < 1.5  # scattering: Att_V < 1.086*2
    assert np.all(np.isfinite(c))
    assert c[-1] == c[-2]           # clamp beyond 3 um

    thin = a.curve(WAVE, 0.05); ref = a.curve(WAVE, 0.25)
    assert np.allclose(thin, ref * 0.05 / 0.25)                 # optically-thin scaling

    sp = SourceSpectrum(
        Empirical1D,
        points=WAVE,
        lookup_table=np.ones(WAVE.size) * 1e-15 * u.Unit("erg/(s cm2 AA)"),
    )
    att = a.attenuate(sp, 2.0)
    assert att(2.2 * u.um).value < sp(2.2 * u.um).value

# ---------------- sampling + table ----------------
def test_alves_pdf_in_ak_and_table_contract():
    spec = {
        "model": "screen",
        "law": "ccm89",
        "rv": 4.5,
        "band": "Ks",
        "distribution": "column_density_pdf",
        "params": {
            "median": "0.6 mag",
            "sigma": 0.45,
            "break": "1.2 mag",
            "slope": 2.0,
        },
    }

    a = parse_attenuation(spec)
    assert a.is_sampled

    tbl = Table({"x": np.zeros(500), "y": np.zeros(500)})
    out = apply_attenuation_to_table(
        tbl,
        a,
        rng=np.random.default_rng(1138),
        band_wavelength_lookup=KS,
    )
    assert out is tbl
    assert ATTENUATION_COLUMN in tbl.colnames
    assert tbl[ATTENUATION_COLUMN].unit == u.mag

    ak_over_av = float(ExtinctionLaw.CCM89.model(4.5)(KS("Ks")))
    ak = tbl[ATTENUATION_COLUMN] * ak_over_av          # back to the A_K the PDF was written in
    assert 0.4 < np.median(ak) < 0.9 and np.all(ak > 0)
    assert tbl.meta["attenuation_model"] == "screen"
    assert tbl.meta["attenuation_parameter"] == "Av"
    assert tbl.meta["attenuation_law"] == "ccm89"
    assert tbl.meta["attenuation_rv"] == 4.5

    # reproducible under one seed
    t2 = apply_attenuation_to_table(
        Table({"x": np.zeros(500)}),
        a,
        rng=np.random.default_rng(1138),
        band_wavelength_lookup=KS,
    )
    assert np.array_equal(t2[ATTENUATION_COLUMN], tbl[ATTENUATION_COLUMN])
    # engine side: rebuild from meta and get a curve per row
    eng = Attenuation.from_table_meta(tbl.meta)
    assert eng.model is AttenuationModel.SCREEN
    assert eng.rv == 4.5
    assert eng.curve(WAVE, float(tbl[ATTENUATION_COLUMN][0])).shape == WAVE.shape

def test_fixed_column_broadcasts_and_absent_is_untouched():
    tbl = Table({"x": np.zeros(3)})
    assert apply_attenuation_to_table(tbl, None) is tbl
    assert tbl.colnames == ["x"]
    assert not tbl.meta
    apply_attenuation_to_table(tbl, parse_attenuation("3 mag"))
    assert np.all(tbl[ATTENUATION_COLUMN] == 3.0)

def test_sampled_needs_rng_and_tau_v_currency():
    a = parse_attenuation({
        "distribution": "lognormal",
        "currency": "tau_v",
        "params": {"median": 2.0, "sigma": 0.3},
        "model": "wg00",
        "geometry": "cloudy",
        "dust_distribution": "clumpy",
    })

    with pytest.raises(AttenuationError):
        apply_attenuation_to_table(Table({"x": [0]}), a)

    tbl = apply_attenuation_to_table(
        Table({"x": np.zeros(200)}),
        a,
        rng=np.random.default_rng(0),
    )
    assert tbl[ATTENUATION_COLUMN].unit == u.dimensionless_unscaled
    assert tbl.meta["attenuation_geometry"] == "cloudy"
    assert 1.5 < np.median(tbl[ATTENUATION_COLUMN]) < 2.6
    eng = Attenuation.from_table_meta(tbl.meta)
    assert eng.geometry == "cloudy"
    assert eng.dust_distribution == "clumpy"

def test_uniform_and_param_errors():
    with pytest.raises(UnphysAttenError):
        ColumnDistribution(
            "column_density_pdf",
            {"median": 1, "sigma": 0.3, "break": 2, "slope": 0.5},
        ).sample(3, np.random.default_rng(0))

    with pytest.raises(EmptyAttenError):
        ColumnDistribution(
            "lognormal",
            {"median": 1},
        ).sample(3, np.random.default_rng(0))

    u_ = ColumnDistribution(
        "uniform",
        {"low": "1 mag", "high": "2 mag"},
    ).sample(50, np.random.default_rng(0))
    assert np.all((u_ >= 1) & (u_ <= 2))
