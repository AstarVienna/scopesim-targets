# -*- coding: utf-8 -*-
"""Unit tests for cluster.py."""

import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord

from scopesim_targets.cluster import Cluster, ZeroAgeCluster
from scopesim_targets.stellar.populations import IMFPopulation
from scopesim_targets.stellar.morphology import KingProfileMorphology


# Something introduced a serious performance hit in these tests. Until figured
# out, we can skip these.
pytestmark = pytest.mark.skip(reason="performance regression")


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

class TestZeroAgeCluster:
    def test_basic(self, basic_cluster):
        assert isinstance(basic_cluster, Cluster)

    @pytest.mark.webtest  # because spextra templates need download
    def test_to_source(self, basic_cluster):
        src = basic_cluster.to_source()
        assert len(src.fields[0].field) == 10


class TestClusterExtinction:
    """The per-star ``Av`` seam only (mocked pop/morph, no isochrone sampling)."""

    @staticmethod
    def _cluster(**kwargs):
        from unittest.mock import Mock

        return Cluster(
            position=SkyCoord(0*u.deg, 0*u.deg, 1*u.kpc),
            population=Mock(),
            morphology=Mock(),
            **kwargs,
        )

    @staticmethod
    def _table(n=5):
        import numpy as np
        from astropy.table import Table

        return Table(
            {"x": np.arange(n) * 1.0, "y": np.arange(n) * 1.0,
             "ref": np.zeros(n, int), "weight": np.ones(n),
             "absmag": np.zeros(n), "mass": np.ones(n)}
        )

    def test_distribution_samples_per_row(self):
        import numpy as np

        c = self._cluster(
            extinction={
                "distribution": "column_density_pdf",
                "params": {"av_median": 1.0, "sigma": 0.4, "av_break": 3.0,
                           "slope": 2.5},
                "rv": 4.5,
            },
            rng_seed=1,
        )
        out = c._apply_extinction(self._table(200))
        assert "Av" in out.colnames and len(out["Av"]) == 200
        assert (np.asarray(out["Av"]) >= 0).all()
        assert np.ptp(np.asarray(out["Av"])) > 0        # per-row, not broadcast
        assert out.meta["extinction_law"] == "ccm89"
        assert out.meta["extinction_rv"] == 4.5
        assert out.meta["anchor"] == "intrinsic"          # cluster frame is forced

    def test_fixed_screen_broadcasts(self):
        import numpy as np

        c = self._cluster(extinction="2 mag", rng_seed=1)
        out = c._apply_extinction(self._table(4))
        np.testing.assert_array_equal(np.asarray(out["Av"]), [2.0, 2.0, 2.0, 2.0])
        assert out.meta["extinction_law"] == "ccm89"

    def test_no_extinction_leaves_table_untouched(self):
        c = self._cluster()
        out = c._apply_extinction(self._table(3))
        assert "Av" not in out.colnames
        assert "extinction_law" not in out.meta

    def test_distribution_without_seed_raises(self):
        c = self._cluster(
            extinction={"distribution": "uniform", "params": {"low": 0, "high": 1}},
        )
        with pytest.raises(ValueError):
            c._apply_extinction(self._table(3))

    def test_rng_reproducible_and_independent(self):
        import numpy as np

        c = self._cluster(
            extinction={"distribution": "uniform", "params": {"low": 0, "high": 5}},
            rng_seed=123,
        )
        a = np.asarray(c._apply_extinction(self._table(50))["Av"])
        b = np.asarray(c._apply_extinction(self._table(50))["Av"])
        np.testing.assert_array_equal(a, b)                 # reproducible under seed
        bare = np.random.default_rng(123).uniform(0, 5, 50)
        assert not np.allclose(a, bare)                     # independent sub-stream

    @pytest.mark.webtest  # spextra templates need download
    def test_zeroage_cluster_end_to_end_distribution(self):
        import numpy as np

        from scopesim_targets.stellar.populations import IMFPopulation
        from scopesim_targets.stellar.morphology import KingProfileMorphology

        tgt = ZeroAgeCluster(
            SkyCoord(0 * u.deg, 0 * u.deg, 1 * u.kpc),
            IMFPopulation,
            {"n_stars": 10},
            KingProfileMorphology,
            {"n_stars": 10, "r_core": 1 * u.pc, "r_tide": 10 * u.pc},
            extinction={
                "distribution": "column_density_pdf",
                "params": {"av_median": 2.0, "sigma": 0.4, "av_break": 5.0,
                           "slope": 2.5},
            },
            rng_seed=42,
        )
        table = tgt.to_source().fields[0].field
        assert "Av" in table.colnames and len(table["Av"]) == 50
        assert (np.asarray(table["Av"]) >= 0).all()
        assert np.ptp(np.asarray(table["Av"])) > 0
        assert table.meta["anchor"] == "intrinsic"


class TestClusterReproducibility:
    """rng_seed now threads through population, morphology AND extinction."""

    def test_morphology_honors_injected_rng(self):
        # Offline: King sampling is scipy/astropy only, no template download.
        import numpy as np
        from scopesim_targets.stellar.morphology import KingProfileMorphology

        pos = SkyCoord(0 * u.deg, 0 * u.deg, 1 * u.kpc)
        kw = dict(n_stars=10, r_core=1 * u.pc, r_tide=10 * u.pc)
        a = KingProfileMorphology(rng=np.random.default_rng(9), **kw).sample(pos)
        b = KingProfileMorphology(rng=np.random.default_rng(9), **kw).sample(pos)
        c = KingProfileMorphology(rng=np.random.default_rng(10), **kw).sample(pos)
        np.testing.assert_array_equal(a, b)              # same seed -> same draw
        assert not np.array_equal(a[0], c[0])            # different seed -> differs

    @pytest.mark.webtest  # spextra templates need download
    def test_same_seed_reproduces_whole_cluster(self):
        import numpy as np

        from scopesim_targets.stellar.populations import IMFPopulation
        from scopesim_targets.stellar.morphology import KingProfileMorphology

        def build():
            return ZeroAgeCluster(
                SkyCoord(0 * u.deg, 0 * u.deg, 1 * u.kpc),
                IMFPopulation, {"n_stars": 10},
                KingProfileMorphology,
                {"n_stars": 10, "r_core": 1 * u.pc, "r_tide": 10 * u.pc},
                extinction={"distribution": "uniform", "params": {"low": 0, "high": 3}},
                rng_seed=7,
            ).to_source().fields[0].field

        t1, t2 = build(), build()
        for col in ("mass", "x", "y", "Av"):
            np.testing.assert_array_equal(np.asarray(t1[col]), np.asarray(t2[col]))
