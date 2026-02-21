# -*- coding: utf-8 -*-
"""Currently only ``Star`` and baseclass."""

from collections.abc import Sequence, Mapping
from itertools import count
from numbers import Number  # matches int, float and all the numpy scalars

from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
from synphot import SourceSpectrum

from astar_utils.guard_functions import guard_same_len
from spextra import Spextrum
from scopesim import Source
from scopesim.source.source_fields import TableSourceField

from .typing_utils import POSITION_TYPE, SPECTRUM_TYPE, BRIGHTNESS_TYPE
from .yaml_constructors import register_target_constructor
from .target import Brightness, SpectrumTarget


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
        source = Source(field=TableSourceField(
            self.to_table(), spectra=self.source_spectra()
        ))
        return source

    def to_table(self, local_frame=None) -> Table:
        """Convert to table for Source conversion."""
        tbl = self._create_source_table()
        tbl.add_row(self._to_table_row(local_frame))
        return tbl

    def source_spectra(self, start: int = 0) -> dict[int, SourceSpectrum]:
        """Create spectra dict for Source conversion."""
        return {start: self.resolve_spectrum(self.spectrum)}

    def _create_source_table(self) -> Table:
        tbl = Table(names=["x", "y", "ref", "weight"],
                    units={"x": u.arcsec, "y": u.arcsec})

        # TODO: Figure out if those are really needed
        tbl.meta["x_unit"] = "arcsec"
        tbl.meta["y_unit"] = "arcsec"
        return tbl

    @staticmethod
    def _xy_arcsec_position(position, local_frame) -> tuple[float, float]:
        # Transform to local offset for ScopeSim
        local_position = position.transform_to(local_frame)

        # ra, dec turn into lon, lat in offset frame, .round(6) is microarcsec
        x_arcsec = local_position.lon.to_value(u.arcsec).round(6)
        y_arcsec = local_position.lat.to_value(u.arcsec).round(6)
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
        x_arcsec, y_arcsec = self._xy_arcsec_position(self.position, local_frame)

        # If not given from parent, resolve now
        if spectrum is None:
            spectrum = self.resolve_spectrum(self.spectrum)
        weight = self._get_spectrum_scale(spectrum, self.brightness)

        row = {
            "x": x_arcsec,
            "y": y_arcsec,
            "weight": weight,
            "ref": ref,
        }
        return row


class Star(PointSourceTarget):
    """A single star."""


class Binary(PointSourceTarget):
    """Binary star."""

    def __init__(
        self,
        position: POSITION_TYPE | None = None,
        offset: Mapping[str, float | u.Quantity] | None = None,
        spectra: Sequence[SPECTRUM_TYPE] | None = None,
        brightness: BRIGHTNESS_TYPE | Sequence[BRIGHTNESS_TYPE] | None = None,
        contrast: float | None = None,
    ) -> None:
        if position is not None:
            self.position = position
        if offset is not None:
            self.offset = offset

        if spectra is not None:
            self.primary_spectrum, self.secondary_spectrum = spectra

        # or rather put this in contrast setter???
        match brightness, contrast:
            case None, None:
                pass  # TODO: What to do here?
            case (str(), u.Quantity() | Number()) as primary, None:
                self.brightness = primary
                # self.contrast = None  # replace?
            case None, contrast:
                # self.brightness = None  # replace?
                self.contrast = contrast
            case (primary, secondary), None:
                self.brightness = primary
                self.brightness_secondary = secondary
            case (str(), u.Quantity() | Number()) as primary, contrast:
                self.brightness = primary
                self.contrast = contrast
            case (primary, secondary), contrast:
                raise TypeError("Either supply brightness and contrast or two "
                                "brightness tuples, but not both.")
            case _:
                raise TypeError("Unkown brightness format.")

    @property
    def primary_spectrum(self) -> SPECTRUM_TYPE:
        """Spectral information of primary component."""
        return self._primary_spectrum

    @primary_spectrum.setter
    def primary_spectrum(self, spectrum: SPECTRUM_TYPE):
        self._primary_spectrum = self._parse_spectrum(spectrum)

    @property
    def secondary_spectrum(self) -> SPECTRUM_TYPE:
        """Spectral information of secondary component."""
        return self._secondary_spectrum

    @secondary_spectrum.setter
    def secondary_spectrum(self, spectrum: SPECTRUM_TYPE):
        self._secondary_spectrum = self._parse_spectrum(spectrum)

    @property
    def brightness_secondary(self) -> Brightness:
        """Brightness of secondary component, if not set via `contrast`."""
        return self._brightness_secondary

    @brightness_secondary.setter
    def brightness_secondary(self, brightness: BRIGHTNESS_TYPE):
        self._brightness_secondary = self._parse_brightness(brightness)

    @property
    def contrast(self) -> float:
        """Contrast ratio between primary and secondary component.

        The contrast ratio is interpreted as ratio between physical flux of the
        primary and secondary component, not difference in magnitudes.

        Brightness of the secondary can also be specified via the
        `brightness_secondary` attribute instead, which supports magnitudes.
        """
        return self._contrast

    @contrast.setter
    def contrast(self, contrast: float):
        if not isinstance(contrast, float):
            # TODO: Also check for dimensionless Quantity here, and also numpy
            #       float-ish types. Int should also be fine...
            raise TypeError("contrast must be float or dimensionless Quantity")
        self._contrast = contrast

    def _resolve_spectra_refs(
        self,
        spectra: Mapping[int, SourceSpectrum] | None,
        refs: Sequence[int] | int | None,
    ) -> tuple[dict[int, SourceSpectrum], tuple[int, ...]]:
        match spectra, refs:
            case None, None:
                spectra = self.source_spectra()
                ref_pri, ref_sec = 0, 1  # default spectra indices
            case None, (_, _):
                raise ValueError("refs sequence must have matching spectra")
            case dict(spectra), None:
                if len(spectra) == 2:
                    ref_pri, ref_sec = spectra.keys()
            case None, int(start):
                spectra = self.source_spectra(start)
                ref_pri, ref_sec = start, start + 1
            case dict(spectra), (ref_pri, ref_sec):
                if not {ref_pri, ref_sec}.issubset(spectra):
                    raise ValueError("not all refs found in spectra")
            case _:
                raise TypeError("refs and spectra not understood")
        return spectra, (ref_pri, ref_sec)

    def _resolve_secondary_weight(
        self,
        secondary_spectrum: SourceSpectrum,
        primary_weight: float,
    ) -> float:
        if hasattr(self, "_contrast"):
            return primary_weight / self.contrast
        if hasattr(self, "_brightness_secondary"):
            return self._get_spectrum_scale(
                secondary_spectrum, self.brightness_secondary)
        raise ValueError("Either contrast or secondary brightness is needed.")

    def to_table(self, local_frame=None, spectra=None, refs=None) -> Table:
        """Convert to table for Source conversion."""
        tbl = self._create_source_table()

        # TODO: Add support for parent frame as in base class
        try:
            primary_position = self.position
        except AttributeError:
            # Default to (0, 0)
            primary_position = SkyCoord(0*u.deg, 0*u.deg)

        local_frame = primary_position.skyoffset_frame()
        x_arcsec_pri, y_arcsec_pri = self._xy_arcsec_position(primary_position, local_frame)
        x_arcsec_sec, y_arcsec_sec = self._xy_arcsec_position(
            self.resolve_position(primary_position),
            local_frame,
        )

        spectra, (ref_pri, ref_sec) = self._resolve_spectra_refs(spectra, refs)

        primary = {
            "x": x_arcsec_pri,
            "y": y_arcsec_pri,
            "weight": self._get_spectrum_scale(
                spectra[ref_pri], self.brightness),
            "ref": ref_pri,
        }

        secondary = {
            "x": x_arcsec_sec,
            "y": y_arcsec_sec,
            "weight": self._resolve_secondary_weight(
                spectra[ref_sec], primary["weight"]),
            "ref": ref_sec,
        }

        tbl.add_row(primary)
        tbl.add_row(secondary)
        return tbl

    def source_spectra(self, start: int = 0) -> dict[int, SourceSpectrum]:
        """Create spectra dict for Source conversion."""
        spectra = {
            start: self.resolve_spectrum(self.primary_spectrum),
            start + 1: self.resolve_spectrum(self.secondary_spectrum),
        }
        return spectra


class Exoplanet(PointSourceTarget):
    """Exoplanet (point source) with default spectrum of Neptune."""

    def __init__(
        self,
        position: POSITION_TYPE | None = None,
        offset: Mapping[str, float | u.Quantity] | None = None,
        spectrum: SPECTRUM_TYPE | None = None,
        brightness: BRIGHTNESS_TYPE | None = None,
        contrast: float | None = None,
    ) -> None:
        if position is not None:
            self.position = position
        if offset is not None:
            self.offset = offset
        if spectrum is not None:
            self.spectrum = spectrum
        if brightness is not None:
            self.brightness = brightness
        if contrast is not None:
            self.contrast = contrast

    @property
    def spectrum(self):
        # Deal with default here
        try:
            return self._spectrum
        except AttributeError:
            pass
        return Spextrum("irtf/Neptune")

    @spectrum.setter
    def spectrum(self, spectrum: SPECTRUM_TYPE):
        self._spectrum = self._parse_spectrum(spectrum)


# TODO: Common base class for multi-component targets
class PlanetarySystem(PointSourceTarget):
    """Planetary system with primary and components."""

    def __init__(
        self,
        position: POSITION_TYPE | None = None,
        primary: PointSourceTarget | None = None,
        components: Sequence[PointSourceTarget] | None = None,
    ) -> None:
        if position is not None:
            self.position = position
        if primary is not None:
            self.primary = primary
        if components is not None:
            self.components = components

    def to_source(self) -> Source:
        """Convert to ScopeSim Source object."""
        local_frame = self.position.skyoffset_frame()

        # HACK: Should be able to pass this down
        self.primary.position = self.position

        table = self.primary.to_table(local_frame)
        spectra = self.primary.source_spectra()

        for ref, component in enumerate(self.components, start=max(spectra)+1):
            spectrum = component.resolve_spectrum(component.spectrum)

            x_arcsec, y_arcsec = self._xy_arcsec_position(
                component.resolve_position(self.position),
                local_frame,
            )

            row = {
                "x": x_arcsec,
                "y": y_arcsec,
                "weight": table[0]["weight"] / component.contrast,
                "ref": ref,
            }

            table.add_row(row)
            spectra[ref] = spectrum

        source = Source(field=TableSourceField(table, spectra=spectra))
        return source


# TODO: Common base class for multi-component targets
class StarField(PointSourceTarget):
    """Multiple Stars."""

    def __init__(
        self,
        positions: Sequence[POSITION_TYPE] | None = None,
        spectra: Sequence[SPECTRUM_TYPE] | None = None,
        brightnesses: Sequence[BRIGHTNESS_TYPE] | None = None,
        band: str | None = None,  # TODO: Proper typing
    ) -> None:
        guard_same_len(positions, spectra, brightnesses)
        self.band = band
        self.positions = positions
        self.spectra = spectra
        self.brightnesses = brightnesses

    @property
    def positions(self):
        try:
            return self._positions
        except AttributeError:
            pass  # return None

    @positions.setter
    def positions(self, positions: Sequence[POSITION_TYPE]):
        try:
            guard_same_len(positions, self.spectra, self.brightnesses)
        except ValueError as err:
            raise ValueError(
                "Positions length doesn't match other attributes"
            ) from err
        self._positions = [
            self._parse_position(position) for position in positions
        ]

    @property
    def spectra(self):
        try:
            return self._spectra
        except AttributeError:
            pass  # return None

    @spectra.setter
    def spectra(self, spectra: Sequence[SPECTRUM_TYPE]):
        try:
            guard_same_len(self.positions, spectra, self.brightnesses)
        except ValueError as err:
            raise ValueError(
                "Spectra length doesn't match other attributes"
            ) from err
        self._spectra = [self._parse_spectrum(spectrum) for spectrum in spectra]

    @property
    def brightnesses(self):
        try:
            return self._brightnesses
        except AttributeError:
            pass  # return None

    @brightnesses.setter
    def brightnesses(self, brightnesses: Sequence[BRIGHTNESS_TYPE]):
        try:
            guard_same_len(self.positions, self.spectra, brightnesses)
        except ValueError as err:
            raise ValueError(
                "Brightnesses length doesn't match other attributes"
            ) from err
        self._brightnesses = [
            self._parse_brightness(brightness)
            if isinstance(brightness, Sequence) and len(brightness) > 1
            else self._parse_brightness((self.band, brightness))
            for brightness in brightnesses
        ]

    def to_source(self) -> Source:
        """Convert to ScopeSim Source object."""
        # TODO: Consider top-level center coords (somehow...)
        local_frame = SkyCoord(0*u.deg, 0*u.deg).skyoffset_frame()

        xy_positions = []
        for position in self.positions:
            xy_positions.append(self._xy_arcsec_position(position, local_frame))

        from more_itertools import unzip
        x_positions, y_positions = unzip(xy_positions)

        spectra_ids = dict(zip(set(self.spectra), count()))
        resolved_spectra = {
            spectrum_id: self.resolve_spectrum(spectrum)
            for spectrum, spectrum_id in spectra_ids.items()
        }

        spec_refs = [spectra_ids[spectrum] for spectrum in self.spectra]
        weights = [
            self._get_spectrum_scale(resolved_spectra[spectrum_id], brightness)
            for spectrum_id, brightness in zip(spec_refs, self.brightnesses)
        ]

        # TODO: Refactor...
        table = Table(
            names=["x", "y", "ref", "weight"],
            units={"x": u.arcsec, "y": u.arcsec},
            data={
                "x": list(x_positions),
                "y": list(y_positions),
                "ref": spec_refs,
                "weight": weights,
            },
        )

        # TODO: Figure out if those are really needed
        table.meta["x_unit"] = "arcsec"
        table.meta["y_unit"] = "arcsec"

        return Source(field=TableSourceField(table, spectra=resolved_spectra))


# TODO: Move these to __init__.py?
register_target_constructor(Star)
register_target_constructor(Binary)
register_target_constructor(Exoplanet)
register_target_constructor(PlanetarySystem)
register_target_constructor(StarField)
