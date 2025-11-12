# -*- coding: utf-8 -*-
"""Currently only ``Star`` and baseclass."""

from astropy import units as u
from astropy.table import Table

from scopesim import Source
from scopesim.source.source_fields import TableSourceField

from .typing_utils import POSITION_TYPE, SPECTRUM_TYPE, BRIGHTNESS_TYPE
from .yaml_constructors import register_target_constructor
from .target import SpectrumTarget


class PointSourceTarget(SpectrumTarget):
    """Base class for Point Source Targets."""

    def __init__(
        self,
        position: POSITION_TYPE | None = None,
        spectrum: SPECTRUM_TYPE | None = None,
        brightness: BRIGHTNESS_TYPE | None = None,
    ) -> None:
        if position is not None:
            self.position = position
        if spectrum is not None:
            self.spectrum = spectrum
        if brightness is not None:
            self.brightness = brightness

    def to_source(self) -> Source:
        """Convert to ScopeSim Source object."""
        tbl = Table(names=["x", "y", "ref", "weight"],
                    units={"x": u.arcsec, "y": u.arcsec})
        tbl.meta["x_unit"] = "arcsec"  # TODO: are those really needed?
        tbl.meta["y_unit"] = "arcsec"
        tbl.add_row(self._to_table_row())

        source = Source(field=TableSourceField(
            tbl, spectra={0: self.resolve_spectrum()}
        ))
        return source

    def _xy_arcsec_position(self, local_frame) -> tuple[float, float]:
        # Transform to local offset for ScopeSim
        local_position = self.position.transform_to(local_frame)

        # ra, dec turn into lon, lat in offset frame, .round(6) is microarcsec
        x_arcsec = local_position.lon.to(u.arcsec).round(6)
        y_arcsec = local_position.lat.to(u.arcsec).round(6)
        return x_arcsec, y_arcsec

    def _to_table_row(
        self,
        local_frame=None,
        spectrum=None,
        ref: int = 0,
    ) -> dict[str, float]:
        # If not given from parent, position is always (0, 0) locally
        if local_frame is None:
            local_frame = self.position.skyoffset_frame()
        x_arcsec, y_arcsec = self._xy_arcsec_position(local_frame)

        # If not given from parent, resolve now
        if spectrum is None:
            spectrum = self.resolve_spectrum()
        weight = self._get_spectrum_scale(spectrum)

        row = {
            "x": x_arcsec,
            "y": y_arcsec,
            "weight": weight,
            "ref": ref,
        }
        return row


class Star(PointSourceTarget):
    """A single star."""


register_target_constructor(Star)
