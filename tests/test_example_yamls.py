# -*- coding: utf-8 -*-
"""Test example YAMLs for RTD."""

from pathlib import Path

import pytest
import yaml

from scopesim_targets.target import Target

# HACK: Temporary to ensue proper import and yaml tag registring
from scopesim_targets.point_source import *
from scopesim_targets.extended_source import *


EXAMPLE_YAML_PATH = Path(__package__).parent / "docs"


@pytest.mark.xfail
def test_examle_yamls_parse(subtests):
    for file in EXAMPLE_YAML_PATH.rglob("*.yaml"):
        with subtests.test(filename=file.name):
            with file.open("r", encoding="utf-8") as stream:
                tgt = yaml.full_load(stream)
            assert isinstance(tgt, Target)


@pytest.mark.xfail
def test_examle_yamls_can_make_source(subtests):
    for file in EXAMPLE_YAML_PATH.rglob("*.yaml"):
        with subtests.test(filename=file.name):
            with file.open("r", encoding="utf-8") as stream:
                tgt = yaml.full_load(stream)
            tgt.to_source()
