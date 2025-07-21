# -*- coding: utf-8 -*-
"""Unit tests for target.py."""

from scopesim_targets.target import Target


class TestTarget:
    def test_basic(self):
        tgt = Target()
        assert isinstance(tgt, Target)
