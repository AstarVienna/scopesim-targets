# -*- coding: utf-8 -*-
"""Test example YAMLs for RTD."""

from pathlib import Path

import pytest
import yaml
from astropy import units as u

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
    opt_dict = {"pixel_scale": 0.1*u.arcsec/u.pix, "width": 200, "height": 100}
    for file in EXAMPLE_YAML_PATH.rglob("*.yaml"):
        with subtests.test(filename=file.name):
            with file.open("r", encoding="utf-8") as stream:
                tgt = yaml.full_load(stream)
            if isinstance(tgt, ParametrizedTarget):
                tgt.to_source(opt_dict)
            else:
                tgt.to_source()
