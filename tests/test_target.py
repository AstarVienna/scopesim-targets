# -*- coding: utf-8 -*-
"""Unit tests for target.py."""

from scopesim_targets.target import Target


class TestTarget:
    def test_basic(self):
        # Note: This is necessary because Target is an abstract base class which
        #       cannot be instantiated directly. The hack below allows to run a
        #       super basic test on it, without using any of its actual subcls.
        class MockTargetSubcls(Target):
            def to_source(self):
                pass
        tgt = MockTargetSubcls()
        assert isinstance(tgt, Target)
