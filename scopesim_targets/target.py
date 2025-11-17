# -*- coding: utf-8 -*-
"""Contains main ``Target`` class."""

from abc import ABCMeta, abstractmethod
from collections import namedtuple
from collections.abc import Mapping
from numbers import Number  # matches int, float and all the numpy scalars

from astropy import units as u
from astropy.coordinates import SkyCoord, Angle, Distance
from synphot import SourceSpectrum, Observation
from synphot.units import PHOTLAM

from astar_utils import SpectralType
from spextra import Spextrum, SpecLibrary, FilterSystem, Passband
from scopesim import Source

from .typing_utils import POSITION_TYPE, SPECTRUM_TYPE, BRIGHTNESS_TYPE


Brightness = namedtuple("Brightness", ["band", "mag"])

# For now, limit possible bands to ETC filters in SpeXtra
FILTER_SYSTEM = FilterSystem("etc")
DEFAULT_LIBRARY = SpecLibrary("bosz/lr")


class Target(metaclass=ABCMeta):
    """Main class in scopesim-targets."""

    @abstractmethod
    def to_source(self) -> Source:
        """Convert to ScopeSim Source object."""
        raise NotImplementedError()

    @property
    def position(self) -> SkyCoord:
        """Target position (center) as SkyCoord."""
        # TODO: Consider adding default (with logging) here if
        #       self._position is None and self._offset is None
        #       But consider also how that might interact with parent position
        #       and offset frame from that.
        return self._position

    @position.setter
    def position(self, position: POSITION_TYPE):
        match position:
            case SkyCoord():
                self._position = position
            case {"x": x_arcsec, "y": y_arcsec, "distance": distance}:
                x_arcsec <<= u.arcsec
                y_arcsec <<= u.arcsec
                distance = Distance(distance)
                self._position = SkyCoord(x_arcsec, y_arcsec, distance)
            case (x_arcsec, y_arcsec) | {"x": x_arcsec, "y": y_arcsec}:
                x_arcsec <<= u.arcsec
                y_arcsec <<= u.arcsec
                self._position = SkyCoord(x_arcsec, y_arcsec)
            case {"distance": distance}:
                # Assume target in center of field
                self._position = SkyCoord(0*u.deg, 0*u.deg, Distance(distance))
            case _:
                raise TypeError("Unkown postition format.")

    @property
    def offset(self) -> dict:
        """Target offset from parent."""
        return self._offset

    @offset.setter
    def offset(self, offset: Mapping[str, float | u.Quantity]):
        if not isinstance(offset, Mapping):
            raise TypeError("Unkown offset format")

        # TODO: Consider adding warning when self._position is not None, because
        #       that would take precedence over any offset.

        self._offset = {
            "separation": offset["separation"],
            "position_angle": Angle(offset.get("position_angle", 0*u.deg)),
        }

    def resolve_offset(self, parent_position: SkyCoord | None = None):
        if hasattr(self, "_offset") and self.offset is not None:
            if parent_position is None:
                raise ValueError("offset needs parent position to resolve")

            separation = self.offset["separation"]
            if separation.unit.physical_type == "length":
                with u.set_enabled_equivalencies(u.dimensionless_angles()):
                    separation = separation / parent_position.distance
                    separation <<= u.arcsec
            elif separation.unit.physical_type != "angle":
                # TODO: Or move this to offset setter??
                raise ValueError("separation must be length or angle")

            position = parent_position.directional_offset_by(
                self.offset["position_angle"],
                separation,
            )
            return position

        # Default to (0, 0)
        return SkyCoord(0*u.deg, 0*u.deg)


class SpectrumTarget(Target):
    """Base class for Targets with separate spectrum (non-cube)."""

    @property
    def spectrum(self) -> SPECTRUM_TYPE:
        """Target spectral information."""
        return self._spectrum

    @spectrum.setter
    def spectrum(self, spectrum: SPECTRUM_TYPE):
        self._spectrum = self._parse_spectrum(spectrum)

    @staticmethod
    def _parse_spectrum(spectrum: SPECTRUM_TYPE):
        match spectrum:
            case SourceSpectrum():
                return spectrum
            case str(spex) if spex.startswith("spex:"):
                # TODO: Consider adding check at this point if spex exists
                return spex
            case str(file) if file.startswith("file:"):
                # TODO: Consider adding check if file exists already here
                return file
            case str() | SpectralType():
                return SpectralType(spectrum)
            case _:
                raise TypeError("Unkown spectrum format.")

    @staticmethod
    def resolve_spectrum(spectrum: SPECTRUM_TYPE) -> SourceSpectrum:
        """
        Create SpeXtrum instance from `spectrum` identifier.

        Can resolve a ``SpectralType`` instance (next-closest available template
        spectrum) or a string that is a valid entry in the SpeXtrum database.

        .. todo:: Actually implement this "next-closest available template", see
            :issue:`68`.

        Returns
        -------
        Spextrum

        """
        if isinstance(spectrum, str) and spectrum.startswith("spex:"):
            # Explicit SpeXtra identifier
            return Spextrum(spectrum.removeprefix("spex:"))

        if isinstance(spectrum, str) and spectrum.startswith("file:"):
            # Explicit SpeXtra identifier
            # TODO: Use pathlib file URI here
            return SourceSpectrum.from_file(spectrum.removeprefix("file:"))

        # HACK: The current DEFAULT_LIBRARY stores spectral classes in lowercase
        #       letters, while SpectralType converts to uppercase. This needs a
        #       proper fix down the road.
        return Spextrum(f"{DEFAULT_LIBRARY.name}/{str(spectrum).lower()}")

    @property
    def brightness(self) -> Brightness:
        """Target brightness information."""
        return self._brightness

    @brightness.setter
    def brightness(self, brightness: BRIGHTNESS_TYPE):
        self._brightness = self._parse_brightness(brightness)

    @staticmethod
    def _parse_brightness(brightness: BRIGHTNESS_TYPE):
        match brightness:
            case str(band), u.Quantity() | Number() as mag:
                # TODO: Consider adding logging about unit assumptions
                # TODO: Implement support for flux instead of mag
                if band not in FILTER_SYSTEM:
                    raise ValueError(f"Band '{band}' unknown.")
                return Brightness(band, mag << u.mag)
            case _:
                raise TypeError("Unkown brightness format.")

    def _get_spectrum_scale(self, spectrum: SourceSpectrum) -> float:
        band = Passband(f"{FILTER_SYSTEM.name}/{self.brightness.band}")

        # TODO: Carefully check this implementation!
        #       Why does Spextrum.flat_spectrum() not need a band?
        ref_flux = Observation(
            Spextrum.flat_spectrum(amplitude=self.brightness.mag),
            band,
        ).effstim(flux_unit=PHOTLAM)
        real_flux = Observation(spectrum, band).effstim(flux_unit=PHOTLAM)

        return float(ref_flux / real_flux)
