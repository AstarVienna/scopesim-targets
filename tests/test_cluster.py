# -*- coding: utf-8 -*-
"""Unit tests for cluster.py."""

from unittest.mock import Mock

import pytest
import numpy as np
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
from synphot import SourceSpectrum
from synphot.models import Empirical1D

from scopesim_targets.seeding import stream, new_rng
from scopesim_targets.cluster import Cluster, ZeroAgeCluster
from scopesim_targets.stellar.populations import IMFPopulation
from scopesim_targets.stellar.morphology import KingProfileMorphology
from scopesim_targets.extinction import ExtinctionError, EmptyExtError


@pytest.fixture
def basic_cluster():
    tgt = ZeroAgeCluster(
        SkyCoord(0*u.deg, 0*u.deg, 1*u.kpc),
        IMFPopulation,
        {"n_stars": 10},
        KingProfileMorphology,
        {"n_stars": 10, "r_core": 1*u.pc, "r_tide": 10*u.pc},
    )
    return tgt


@pytest.mark.slow
class TestZeroAgeCluster:
    def test_basic(self, basic_cluster):
        assert isinstance(basic_cluster, Cluster)

    @pytest.mark.webtest  # because spextra templates need download
    def test_to_source(self, basic_cluster):
        src = basic_cluster.to_source()
        assert len(src.fields[0].field) == 10


class TestClusterDust:
    """The two dust seams (mocked pop/morph, no isochrone sampling).

    ``extinction`` (the sightline) reddens the shared templates in place;
    ``attenuation`` (the natal cloud) becomes the per-row ``attenuation`` column.
    """

    @staticmethod
    def _cluster(**kwargs):
        class ConcreteCluster(Cluster):
            """`Cluster` is abstract (no `to_source`); these tests only
            exercise its dust helpers, so a minimal concrete subclass will do
            (same pattern as ``MockTargetSubcls`` in test_target.py)."""

            def to_source(self, optical_train=None):
                raise NotImplementedError

        return ConcreteCluster(
            position=SkyCoord(0*u.deg, 0*u.deg, 1*u.kpc),
            population=Mock(),
            morphology=Mock(),
            **kwargs,
        )

    @staticmethod
    def _table(n=5):
        return Table(
            {"x": np.arange(n) * 1.0, "y": np.arange(n) * 1.0,
             "ref": np.zeros(n, int), "weight": np.ones(n),
             "absmag": np.zeros(n), "mass": np.ones(n)}
        )

    @staticmethod
    def _flat_sed():
        wave = np.geomspace(3000, 30000, 300) * u.AA
        return SourceSpectrum(
            Empirical1D, points=wave,
            lookup_table=np.ones(300) * 1e-14 * u.Unit("erg / (s cm2 AA)"),
        )

    # -- attenuation: the per-row column --------------------------------------

    def test_distribution_samples_per_row(self):
        c = self._cluster(
            attenuation={
                "distribution": "column_density_pdf",
                "params": {"median": 1.0, "sigma": 0.4, "break": 3.0, "slope": 2.5},
                "rv": 4.5,
            },
            rng_seed=1,
        )
        out = c._apply_attenuation(self._table(200))
        col = np.asarray(out["attenuation"])
        assert len(col) == 200 and (col >= 0).all()
        assert np.ptp(col) > 0                              # per-row, not broadcast
        assert out["attenuation"].unit == u.mag             # screen model: A_V
        assert out.meta["attenuation_model"] == "screen"
        assert out.meta["attenuation_parameter"] == "Av"
        assert out.meta["attenuation_law"] == "ccm89"
        assert out.meta["attenuation_rv"] == 4.5
        assert out.meta["anchor"] == "intrinsic"            # cluster frame is forced

    def test_distribution_in_a_k_converts_to_av(self):
        # Alves+2017 quote the PDF in A_K: 'band: Ks' makes that native input;
        # the column is emitted in the model's native A_V (A_K/A_V ~ 0.11 for
        # CCM89), so the realized A_V is ~9x the A_K parameters. Offline via
        # an injected band-wavelength lookup (Cluster's real one needs spextra).
        c = self._cluster(
            attenuation={
                "distribution": "lognormal",
                "params": {"median": "1 mag", "sigma": 0.01},
                "band": "Ks",
            },
            rng_seed=1,
        )
        c._band_wavelength_lookup = lambda band: 2.15 * u.micron
        col = np.asarray(c._apply_attenuation(self._table(50))["attenuation"])
        assert 7 < np.median(col) < 11

    def test_fixed_column_broadcasts(self):
        c = self._cluster(attenuation="2 mag", rng_seed=1)
        out = c._apply_attenuation(self._table(4))
        np.testing.assert_array_equal(np.asarray(out["attenuation"]), [2.0, 2.0, 2.0, 2.0])
        assert out.meta["attenuation_law"] == "ccm89"

    def test_no_attenuation_leaves_table_untouched(self):
        # ... even with a foreground screen: that one lives in the spectra.
        for c in (self._cluster(), self._cluster(extinction="2 mag")):
            out = c._apply_attenuation(self._table(3))
            assert "attenuation" not in out.colnames and "Av" not in out.colnames
            assert not out.meta

    def test_distribution_without_seed_is_reproducible_from_record(self):
        # No explicit seed: fresh entropy is drawn once and recorded as
        # rng_seed, so the sampled column can be reproduced from the target.
        attenuation = {"distribution": "uniform", "params": {"low": 0, "high": 1}}
        original = self._cluster(attenuation=attenuation)
        assert isinstance(original.rng_seed, int)
        replay = self._cluster(attenuation=attenuation, rng_seed=original.rng_seed)
        np.testing.assert_array_equal(
            np.asarray(original._apply_attenuation(self._table(20))["attenuation"]),
            np.asarray(replay._apply_attenuation(self._table(20))["attenuation"]),
        )

    def test_rng_reproducible_and_independent(self):
        c = self._cluster(
            attenuation={"distribution": "uniform", "params": {"low": 0, "high": 5}},
            rng_seed=123,
        )
        a = np.asarray(c._apply_attenuation(self._table(50))["attenuation"])
        b = np.asarray(c._apply_attenuation(self._table(50))["attenuation"])
        np.testing.assert_array_equal(a, b)                 # reproducible under seed
        bare = np.random.default_rng(123).uniform(0, 5, 50)
        assert not np.allclose(a, bare)                     # independent sub-stream

    # -- extinction: the shared templates ---------------------------------------

    def test_extinction_distribution_is_rejected(self):
        # A PDF is a column at the source, i.e. attenuation; the screen parser
        # sees a mapping without any currency.
        with pytest.raises(EmptyExtError):
            self._cluster(
                extinction={"distribution": "uniform", "params": {"low": 0, "high": 1}},
            )

    def test_extinction_list_is_rejected(self):
        with pytest.raises(ExtinctionError):
            self._cluster(extinction=["0.3 mag", "1 mag"])

    def test_extinction_reddens_shared_templates_in_place(self):
        sed = self._flat_sed()
        c = self._cluster(extinction="2 mag")
        out = c._redden_spectra({0: sed, 1: sed})
        assert set(out) == {0, 1}                            # refs preserved
        ratio_v = float(out[0](5500 * u.AA) / sed(5500 * u.AA))
        ratio_k = float(out[0](22000 * u.AA) / sed(22000 * u.AA))
        np.testing.assert_allclose(ratio_v, 10 ** (-0.4 * 2.0), rtol=0.03)  # dims (intrinsic frame)
        assert ratio_k > ratio_v                                             # and reddens
        # sequence containers keep their type
        assert isinstance(c._redden_spectra([sed]), list)

    def test_no_or_zero_extinction_returns_input(self):
        spectra = {0: self._flat_sed()}
        assert self._cluster()._redden_spectra(spectra) is spectra
        assert self._cluster(extinction="0 mag")._redden_spectra(spectra) is spectra

    @pytest.mark.slow
    @pytest.mark.webtest  # spextra templates need download
    def test_zeroage_cluster_end_to_end_distribution(self):
        tgt = ZeroAgeCluster(
            SkyCoord(0 * u.deg, 0 * u.deg, 1 * u.kpc),
            IMFPopulation,
            {"n_stars": 10},
            KingProfileMorphology,
            {"n_stars": 10, "r_core": 1 * u.pc, "r_tide": 10 * u.pc},
            extinction="0.3 mag",
            attenuation={
                "distribution": "column_density_pdf",
                "params": {"median": 2.0, "sigma": 0.4, "break": 5.0, "slope": 2.5},
            },
            rng_seed=42,
        )
        table = tgt.to_source().fields[0].field
        assert "attenuation" in table.colnames and "Av" not in table.colnames
        assert (np.asarray(table["attenuation"]) >= 0).all()
        assert np.ptp(np.asarray(table["attenuation"])) > 0
        assert table.meta["attenuation_model"] == "screen"
        assert table.meta["anchor"] == "intrinsic"


@pytest.mark.slow
class TestClusterReproducibility:
    """rng_seed threads through population, morphology AND the attenuation draw."""

    # Morphology/population seeding itself is covered in test_seeding.py.

    def test_attenuation_uses_its_named_stream(self):
        # Offline: the attenuation draw is the "attenuation" stream of the
        # cluster's master seed, independent of population and morphology.
        tgt = ZeroAgeCluster(
            SkyCoord(0 * u.deg, 0 * u.deg, 1 * u.kpc),
            IMFPopulation, {"n_stars": 10},
            KingProfileMorphology,
            {"n_stars": 10, "r_core": 1 * u.pc, "r_tide": 10 * u.pc},
            rng_seed=7,
        )
        np.testing.assert_array_equal(
            tgt._attenuation_rng().random(5),
            new_rng(stream(7, "attenuation")).random(5),
        )

    @pytest.mark.webtest  # spextra templates need download
    def test_same_seed_reproduces_whole_cluster(self):
        def build():
            return ZeroAgeCluster(
                SkyCoord(0 * u.deg, 0 * u.deg, 1 * u.kpc),
                IMFPopulation, {"n_stars": 10},
                KingProfileMorphology,
                {"n_stars": 10, "r_core": 1 * u.pc, "r_tide": 10 * u.pc},
                attenuation={"distribution": "uniform", "params": {"low": 0, "high": 3}},
                rng_seed=7,
            ).to_source().fields[0].field

        t1, t2 = build(), build()
        for col in ("mass", "x", "y", "attenuation"):
            np.testing.assert_array_equal(np.asarray(t1[col]), np.asarray(t2[col]))
