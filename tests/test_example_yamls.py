# -*- coding: utf-8 -*-
"""Test example YAMLs for RTD."""

from pathlib import Path

import pytest
import yaml

from scopesim_targets.target import Target


EXAMPLE_YAML_PATH = Path(__package__).parent / "docs"


@pytest.mark.xfail
@pytest.mark.parametrize("file", EXAMPLE_YAML_PATH.rglob("*.yaml"))
def test_examle_yamls_parse(file):
    with file.open("r", encoding="utf-8") as stream:
        tgt = yaml.full_load(stream)
    assert isinstance(tgt, Target)


@pytest.mark.xfail
@pytest.mark.parametrize("file", EXAMPLE_YAML_PATH.rglob("*.yaml"))
def test_examle_yamls_can_make_source(file):
    with file.open("r", encoding="utf-8") as stream:
        tgt = yaml.full_load(stream)
    tgt.to_source()
