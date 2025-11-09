# -*- coding: utf-8 -*-
"""Contains main ``Target`` class."""

from abc import ABCMeta, abstractmethod
from collections import namedtuple
from collections.abc import Mapping
from numbers import Number  # matches int, float and all the numpy scalars

from astropy import units as u
from astropy.coordinates import SkyCoord, Angle
from synphot import SourceSpectrum, Observation
from synphot.units import PHOTLAM

from astar_utils import SpectralType
from spextra import Spextrum, SpecLibrary, FilterSystem, Passband


Brightness = namedtuple("Brightness", ["band", "mag"])

# For now, limit possible bands to ETC filters in SpeXtra
FILTER_SYSTEM = FilterSystem("etc")
DEFAULT_LIBRARY = SpecLibrary("bosz/lr")


class Target(metaclass=ABCMeta):
    """Main class in scopesim-targets."""

    @abstractmethod
    def to_source(self):
        """Convert to ScopeSim Source object."""
        raise NotImplementedError()

    @property
    def position(self):
        """Target position (center) as SkyCoord."""
        # TODO: Consider adding default (with logging) here if
        #       self._position is None and self._offset is None
        #       But consider also how that might interact with parent position
        #       and offset frame from that.
        return self._position

    # TODO: add typing
    @position.setter
    def position(self, position):
        match position:
            case SkyCoord():
                self._position = position
            case (x_arcsec, y_arcsec) | {"x": x_arcsec, "y": y_arcsec}:
                x_arcsec <<= u.arcsec
                y_arcsec <<= u.arcsec
                self._position = SkyCoord(x_arcsec, y_arcsec)
            case _:
                raise TypeError("Unkown postition format.")

    @property
    def offset(self) -> dict:
        """Target offset from parent."""
        return self._offset

    @offset.setter
    def offset(self, offset: Mapping):
        if not isinstance(offset, Mapping):
            raise TypeError("Unkown offset format")

        # TODO: Consider adding warning when self._position is not None, because
        #       that would take precedence over any offset.

        self._offset = {
            "separation": offset["separation"],
            "position_angle": Angle(offset.get("position_angle", 0*u.deg)),
        }


class SpectrumTarget(Target):
    """Base class for Targets with separate spectrum (non-cube)."""

    @property
    def spectrum(self):
        """Target spectral information."""
        return self._spectrum

    # TODO: add typing
    @spectrum.setter
    def spectrum(self, spectrum):
        match spectrum:
            case SourceSpectrum():
                self._spectrum = spectrum
            case str(spex) if spex.startswith("spex:"):
                # TODO: Consider adding check at this point if spex exists
                self._spectrum = spex
            case str(file) if file.startswith("file:"):
                # TODO: Consider adding check if file exists already here
                self._spectrum = file
            case str() | SpectralType():
                self._spectrum = SpectralType(spectrum)
            case _:
                raise TypeError("Unkown spectrum format.")

    def resolve_spectrum(self) -> Spextrum:
        """
        Create SpeXtrum instance from `self.spectrum` identifier.

        Can resolve a ``SpectralType`` instance (next-closest available template
        spectrum) or a string that is a valid entry in the SpeXtrum database.

        .. todo:: Actually implement this "next-closest available template".

        Returns
        -------
        Spextrum

        """
        if isinstance(self.spectrum, str) and self.spectrum.startswith("spex:"):
            # Explicit SpeXtra identifier
            return Spextrum(self.spectrum.removeprefix("spex:"))

        if isinstance(self.spectrum, str) and self.spectrum.startswith("file:"):
            # Explicit SpeXtra identifier
            # TODO: Use pathlib file URI here
            return SourceSpectrum.from_file(self.spectrum.removeprefix("file:"))

        # HACK: The current DEFAULT_LIBRARY stores spectral classes in lowercase
        #       letters, while SpectralType converts to uppercase. This needs a
        #       proper fix down the road.
        return Spextrum(f"{DEFAULT_LIBRARY.name}/{str(self.spectrum).lower()}")

    @property
    def brightness(self):
        """Target brightness information."""
        return self._brightness

    # TODO: add typing
    @brightness.setter
    def brightness(self, brightness):
        match brightness:
            case str(band), u.Quantity() | Number() as mag:
                # TODO: Consider adding logging about unit assumptions
                # TODO: Implement support for flux instead of mag
                if band not in FILTER_SYSTEM:
                    raise ValueError(f"Band '{band}' unknown.")
                self._brightness = Brightness(band, mag << u.mag)
            case _:
                raise TypeError("Unkown brightness format.")

    def _get_spectrum_scale(self, spectrum: SourceSpectrum) -> float:
        filter_name = f"{FILTER_SYSTEM.name}/{self.brightness.band}"
        band = Passband(filter_name)

        # TODO: Carefully check this implementation!
        #       Why does Spextrum.flat_spectrum() not need a band?
        ref_flux = Observation(
            Spextrum.flat_spectrum(amplitude=self.brightness.mag),
            band,
        ).effstim(flux_unit=PHOTLAM)
        real_flux = Observation(spectrum, band).effstim(flux_unit=PHOTLAM)

        return float(ref_flux / real_flux)
