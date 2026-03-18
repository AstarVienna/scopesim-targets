# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 16:09:31 2026

@author: ghost
"""

import pytest

from scopesim_targets.spectral_classes import (
    StellarParameters,
)


class TestStellarParameters:
    @pytest.mark.parametrize(
        ("req", "tbl_len"), ((None, 118), ("J-H", 101), (["J-H", "B-V"], 76))
    )
    def test_required_cols(self, req, tbl_len):
        stp = StellarParameters(req)
        assert len(stp.table) == tbl_len
