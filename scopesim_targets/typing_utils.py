# -*- coding: utf-8 -*-
"""Custom composite types used in this package.

Currently contains `POSITION_TYPE`, `SPECTRUM_TYPE` and `BRIGHTNESS_TYPE`,
neither of which is final. These are kept here so that Target subclasses can
simply import and use them, and when we eventually refine them, the code doesn't
need to be updated everywhere.
"""

import typing
from collections.abc import Mapping

from astropy import units as u
from astropy.coordinates import SkyCoord
from synphot import SourceSpectrum

from astar_utils import SpectralType


# TODO: Properly define POSITION_TYPE
POSITION_TYPE = SkyCoord | tuple[float, float] | Mapping[str, float | u.Quantity]

# TODO: Properly define SPECTRUM_TYPE
SPECTRUM_TYPE = SourceSpectrum | SpectralType | str

# TODO: Properly define BRIGHTNESS_TYPE
BRIGHTNESS_TYPE = typing.Any  # PLACEHOLDER
