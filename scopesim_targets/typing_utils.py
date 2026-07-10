# -*- coding: utf-8 -*-
"""Custom composite types used in this package.

Currently contains `POSITION_TYPE`, `SPECTRUM_TYPE` and `BRIGHTNESS_TYPE`,
neither of which is final. These are kept here so that Target subclasses can
simply import and use them, and when we eventually refine them, the code doesn't
need to be updated everywhere.
"""

from collections.abc import Mapping, Sequence

from astropy import units as u
from astropy.coordinates import SkyCoord
from synphot import SourceSpectrum

from astar_utils import SpectralType


# TODO: Properly define POSITION_TYPE
POSITION_TYPE = SkyCoord | tuple[float, float] | Mapping[str, float | u.Quantity]

# TODO: Properly define SPECTRUM_TYPE
SPECTRUM_TYPE = SourceSpectrum | SpectralType | str

# Accepted `brightness` input: either the ``(locator, amount)`` tuple sugar or
# the canonical mapping form (``{band|wavelength|frequency: ..., value: ...}``
# or ``{from_spectral_type: ...}``). ``brightness.parse_brightness`` normalizes
# either into a ``brightness.Brightness``. Both slots individually accept
# ``str``/``Quantity``/number, so this stays deliberately loose at the top
# level. (Retires the former ``typing.Any`` placeholder.)
BRIGHTNESS_TYPE = Sequence | Mapping
