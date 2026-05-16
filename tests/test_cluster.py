# -*- coding: utf-8 -*-
"""Unit tests for cluster.py."""

import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord

from scopesim_targets.cluster import Cluster, ZeroAgeCluster
from scopesim_targets.stellar.populations import IMFPopulation
from scopesim_targets.stellar.morphology import KingProfileMorphology


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
