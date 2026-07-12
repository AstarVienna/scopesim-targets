# -*- coding: utf-8 -*-
"""Contains main ``Target`` class."""

from abc import ABCMeta, abstractmethod
from functools import lru_cache
from collections.abc import Mapping

from astropy import units as u
from astropy.coordinates import SkyCoord, Angle, Distance
from synphot import SourceSpectrum
from synphot.units import VEGAMAG

from astar_utils import SpectralType
from spextra import Spextrum, SpecLibrary, FilterSystem, Passband
from scopesim import Source

from .typing_utils import POSITION_TYPE, SPECTRUM_TYPE, BRIGHTNESS_TYPE
from .flux_scaling import synphot_flux_scale
from .brightness import (
    parse_brightness,
    Brightness,
    LocatorKind,
    AmountKind,
    PhotometricSystem,
    AnchorFrame,
    FromSpectralType,
    BrightnessError,
)


# Magnitude system -> the synphot function unit black_body_spectrum expects.
_BB_MAG_UNIT = {
    PhotometricSystem.VEGA: VEGAMAG,
    PhotometricSystem.AB: u.ABmag,
    PhotometricSystem.ST: u.STmag,
}


@lru_cache(maxsize=1)
def _vega_reference() -> SourceSpectrum:
    """The Vega reference spectrum (CALSPEC), fetched once and cached.

    Used as the zero-point reference for VEGA-system magnitudes. Network-backed
    (hence a webtest in the suite); cached so a StarField of N stars downloads
    it at most once.
    """
    return SourceSpectrum.from_vega()


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
    def anchor(self) -> AnchorFrame:
        """Frame the ``brightness`` value(s) are anchored in.

        One of :class:`~.brightness.AnchorFrame` (``observed`` default,
        ``intrinsic``, ``absolute``). Accepts the enum or the schema strings
        on assignment. The flux-scaling path (:meth:`SpectrumTarget.
        _anchored_spectrum_scale`) reads this to decide which SED is scaled and
        whether the distance modulus enters.
        """
        return getattr(self, "_anchor", AnchorFrame.OBSERVED)

    @anchor.setter
    def anchor(self, value: AnchorFrame | str | None):
        self._anchor = AnchorFrame.coerce(value)

    def _distance_or_none(self) -> u.Quantity | None:
        """Return ``position.distance`` as a length, or ``None`` if unset.

        A :class:`~astropy.coordinates.SkyCoord` with no distance exposes a
        *dimensionless* ``distance`` of ``1``; only a genuine length counts
        (used to raise E10 for ``anchor: absolute`` without a distance).
        """
        position = getattr(self, "_position", None)
        if position is None:
            return None
        distance = getattr(position, "distance", None)
        if distance is None:
            return None
        try:
            return distance.to(u.pc)
        except u.UnitConversionError:
            return None

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
        self._position = self._parse_position(position)

    @staticmethod
    def _parse_position(position: POSITION_TYPE) -> SkyCoord:
        match position:
            case SkyCoord():
                return position
            case {"x": x_arcsec, "y": y_arcsec, "distance": distance}:
                x_arcsec <<= u.arcsec
                y_arcsec <<= u.arcsec
                distance = Distance(distance)
                return SkyCoord(x_arcsec, y_arcsec, distance)
            case (x_arcsec, y_arcsec) | {"x": x_arcsec, "y": y_arcsec}:
                x_arcsec <<= u.arcsec
                y_arcsec <<= u.arcsec
                return SkyCoord(x_arcsec, y_arcsec)
            case {"distance": distance}:
                # Assume target in center of field
                return SkyCoord(0*u.deg, 0*u.deg, Distance(distance))
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

    def _resolve_offset(self, parent_position: SkyCoord | None = None):
        if parent_position is None:
            raise ValueError("If offset is used, parent_position is required.")

        try:
            with length_angle_context(parent_position.distance):
                position = parent_position.directional_offset_by(
                    self.offset["position_angle"],
                    self.offset["separation"] << u.arcsec,
                )
        except u.UnitConversionError as err:
            # TODO: Catch wrong units in offset setter??
            raise ValueError("separation must be length or angle") from err

        return position

    def resolve_position(self, parent_position: SkyCoord | None = None):
        """
        Resolve target position or offset.

        This uses the following lookup order:
        1. `self.position` set? -> use that
        2. parent position present?
          a. `self.offset` set? -> resolve offset to parent position
          b. otherwise use parent position
        3. `self.offset` set, but no parent position present -> Error
        4. default to (0, 0)

        Parameters
        ----------
        parent_position : SkyCoord | None, optional
            Position of any parent target. If None (the default), `self.offset`
            must not be set.

        Raises
        ------
        ValueError
            Raised if `self.offset` set, but `parent_position` is None.

        Returns
        -------
        position : SkyCoord
            Resolved position as SkyCoord object.

        """
        # Offset needs to be checked first in order for binary with separation
        # to be resolved correctly, otherwise would just return primary.
        if hasattr(self, "_offset") and self.offset is not None:
            return self._resolve_offset(parent_position)

        if hasattr(self, "_position") and self.position is not None:
            # parent_position is ignored here...
            return self.position

        if parent_position is not None:
            return parent_position

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
            case str(file) if file.startswith("blackbody:"):
                # TODO: idk
                return file
            case str() | SpectralType():
                return SpectralType(spectrum)
            case _:
                raise TypeError("Unkown spectrum format.")

    @staticmethod
    def resolve_spectrum(
        spectrum: SPECTRUM_TYPE,
        brightness: Brightness | None = None,
    ) -> SourceSpectrum:
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
        if isinstance(spectrum, SourceSpectrum):
            # TODO: Convert to SpeXtrum to get full method access?
            return spectrum

        if isinstance(spectrum, str) and spectrum.startswith("spex:"):
            # Explicit SpeXtra identifier
            return Spextrum(spectrum.removeprefix("spex:"))

        if isinstance(spectrum, str) and spectrum.startswith("file:"):
            # TODO: Convert to SpeXtrum to get full method access?
            # TODO: Use pathlib file URI here
            return SourceSpectrum.from_file(spectrum.removeprefix("file:"))

        if isinstance(spectrum, str) and spectrum.startswith("blackbody:"):
            temp = u.Quantity(spectrum.removeprefix("blackbody:"))
            amplitude, band = SpectrumTarget._blackbody_amplitude(brightness)
            spec = Spextrum.black_body_spectrum(temp, amplitude, band)
            return spec

        # HACK: The current DEFAULT_LIBRARY stores spectral classes in lowercase
        #       letters, while SpectralType converts to uppercase. This needs a
        #       proper fix down the road.
        return Spextrum(f"{DEFAULT_LIBRARY.name}/{str(spectrum).lower()}")

    @staticmethod
    def redshift_spectrum(spectrum: Spextrum, position: SkyCoord) -> Spextrum:
        """Doppler shift spectrum based on position `z` or `v_rad`."""
        # TODO: Add proper unit tests for this!
        if (
            isinstance(position.distance, Distance)  # catch default distance
            and position.distance > 50 * u.Mpc
        ):
            return spectrum.redshift(z=position.distance.z)

        try:
            return spectrum.redshift(vel=position.radial_velocity)
        except ValueError:  # no radial_velocity defined
            pass

        return spectrum  # return as-is

    @property
    def brightness(self) -> Brightness:
        """Target brightness information.

        If a ``{from_spectral_type: ...}`` resolver was supplied, it is resolved
        lazily on first access (it needs the already-set ``spectrum``), memoized
        into a concrete :class:`Brightness`, and its provenance recorded in
        :attr:`brightness_provenance`.
        """
        if getattr(self, "_brightness_resolver", None) is not None:
            self._brightness = self._resolve_from_spectral_type(
                self._brightness_resolver
            )
            self._brightness_resolver = None
        return self._brightness

    @brightness.setter
    def brightness(self, brightness: BRIGHTNESS_TYPE):
        parsed = self._parse_brightness(brightness)
        if isinstance(parsed, FromSpectralType):
            # Deferred: resolved at first `brightness` access, once `spectrum`
            # is known (see the getter).
            self._brightness_resolver = parsed
            self._brightness = None
        else:
            self._brightness_resolver = None
            self._brightness = parsed

    @property
    def brightness_provenance(self) -> dict | None:
        """Resolver provenance (table name/version, resolved value), or ``None``.

        Populated when ``brightness`` was given as ``{from_spectral_type: ...}``;
        recorded so exports remain reproducible when the lookup table is later
        revised (same contract as the extinction ``from_map`` resolver).
        """
        return getattr(self, "_brightness_provenance", None)

    @staticmethod
    def _parse_brightness(
        brightness: BRIGHTNESS_TYPE,
    ) -> Brightness | FromSpectralType:
        """Normalize the ``brightness`` input into a :class:`Brightness`.

        Structural parsing (the two-slot grammar, physical-type dispatch,
        magnitude systems, surface-brightness divisor and the E1-E5 error
        matrix) is delegated to :func:`.brightness.parse_brightness`. This
        wrapper adds the one load-time check the pure parser deliberately omits:
        band membership in the active :data:`FILTER_SYSTEM`. The
        ``{from_spectral_type: ...}`` resolver marker is passed through untouched
        (its band is validated when the resolver fires).
        """
        parsed = parse_brightness(brightness)
        if isinstance(parsed, FromSpectralType):
            return parsed
        if (
            parsed.locator_kind is LocatorKind.BAND
            and parsed.locator not in FILTER_SYSTEM
        ):
            raise ValueError(f"Band '{parsed.locator}' unknown.")
        return parsed

    def _resolve_from_spectral_type(
        self, resolver: FromSpectralType
    ) -> Brightness:
        """Resolve ``{from_spectral_type: <table>}`` to an absolute-magnitude
        :class:`Brightness`, forcing ``anchor: absolute``.

        Looks up the target's spectral type in the named table (currently only
        ``mamajek`` -> :class:`~.spectral_classes.StellarParameters`), reads the
        absolute magnitude in the resolver's band (``M_<band>``), and records
        the table name/version plus the resolved value in
        :attr:`brightness_provenance`. The ordinary machinery then applies the
        distance modulus (missing distance -> E10) and any screens.

        Raises E12 if the spectrum is not a :class:`~astar_utils.SpectralType`
        (a ``file:``/``blackbody:``/``SourceSpectrum`` has no spectral type to
        look up). Never a fallback for missing brightness -- only fires when the
        resolver form was explicitly given.
        """
        spectrum = self.spectrum
        if not isinstance(spectrum, SpectralType):
            raise BrightnessError(
                "E12",
                "from_spectral_type needs a SpectralType spectrum, not "
                f"{spectrum!r}",
            )

        table = self._lookup_table(resolver.table)
        row = table.closest_spectral_type(spectrum)
        column = f"M_{resolver.band}"
        if column not in row.colnames:
            available = [c for c in row.colnames if c.startswith("M_")]
            raise ValueError(
                f"table {resolver.table!r} has no absolute magnitude for band "
                f"{resolver.band!r} (column {column!r}); available: {available}"
            )
        abs_mag = u.Quantity(row[column]).to_value(u.mag)

        # Absolute magnitudes here are on the table's native system (Vega for
        # the Johnson/2MASS columns); anchor: absolute means "unextincted, at
        # 10 pc", so the standard machinery supplies distance modulus + screens.
        self.anchor = AnchorFrame.ABSOLUTE
        self._brightness_provenance = {
            "resolver": "from_spectral_type",
            "table": resolver.table,
            "table_version": self._table_version(table),
            "matched_spectral_type": str(row["spectral_type"]),
            "band": resolver.band,
            "resolved_value": f"{abs_mag} mag",
        }
        return self._parse_brightness((resolver.band, f"{abs_mag} mag"))

    @staticmethod
    def _lookup_table(name: str):
        """Map a resolver table name to its lookup object (E: unknown name).

        The name is validated *before* the backend is imported, so an unknown
        table fails cleanly without pulling the data-file machinery.
        """
        known = {"mamajek"}
        if name not in known:
            raise ValueError(
                f"unknown from_spectral_type table {name!r}; "
                f"available: {sorted(known)}"
            )
        # Imported lazily: StellarParameters pulls a data file, keeping the
        # heavy dependency off the module-import path.
        from .spectral_classes import StellarParameters

        return {"mamajek": StellarParameters}[name]()

    @staticmethod
    def _table_version(table) -> str:
        """Best-effort version string for provenance (table meta, else n/a)."""
        meta = getattr(getattr(table, "table", None), "meta", {}) or {}
        return str(meta.get("version", meta.get("reference", "unversioned")))

    @staticmethod
    def _blackbody_amplitude(
        brightness: Brightness | None,
    ) -> tuple[u.Quantity, str]:
        """Amplitude + reference band for ``Spextrum.black_body_spectrum``.

        Unlike the library-spectrum path (resolve a shape, then scale it via
        :func:`~.flux_scaling.synphot_flux_scale`), spextra bakes the flux scale
        into blackbody *construction*: it wants an ``amplitude`` (a magnitude in
        some system, or a spectral flux density such as ``Jy``) at a reference
        ``band``. So the scaling for a blackbody happens here rather than in the
        shared scaler.

        .. note::
            This asymmetry is a spextra API constraint
            (``black_body_spectrum`` has no shape-only mode). If spextra gains
            one, fold blackbody into the shared resolve-then-scale path and drop
            this helper.
        """
        if brightness is None:
            raise BrightnessError(
                "E9", "a blackbody spectrum requires a brightness to scale it"
            )
        b = brightness
        if b.is_surface_brightness:
            raise BrightnessError(
                "E7", "surface brightness is invalid for a point source"
            )
        if b.locator_kind is not LocatorKind.BAND:
            raise BrightnessError(
                "E3",
                "a blackbody amplitude needs a band locator (spextra's "
                "black_body_spectrum takes a reference band)",
            )
        if b.amount_kind is AmountKind.MAG:
            amplitude = b.value.value * _BB_MAG_UNIT[b.system]
        elif b.amount_kind is AmountKind.ENERGY_FLUX:
            raise BrightnessError(
                "E1",
                "a band-integrated energy flux is not a valid blackbody "
                "amplitude -- give a magnitude or a spectral flux density",
            )
        else:  # spectral flux density (Jy, erg/s/cm2/A): pass through
            amplitude = b.value
        return amplitude, b.locator

    @staticmethod
    def _get_spectrum_scale(
        spectrum: SourceSpectrum,
        brightness: Brightness,
    ) -> float:
        """Factor to scale ``spectrum`` so its photometry matches ``brightness``.

        Resolves the two network-backed inputs -- the bandpass (spextra
        ``Passband``) and, for Vega magnitudes, the Vega reference -- and
        delegates the unit arithmetic to
        :func:`~.flux_scaling.synphot_flux_scale`. Covers every branch of the
        grammar (magnitude in any system; flux density per frequency or
        wavelength; band-integrated energy flux) at band / wavelength /
        frequency locators.
        """
        # A point source has no solid angle, so a surface brightness is invalid.
        # Profiles reduce a surface brightness to an integrated amount *before*
        # calling this, so they never trip E7 here.
        if brightness.is_surface_brightness:
            raise BrightnessError(
                "E7", "surface brightness is invalid for a point source"
            )

        band = None
        if brightness.locator_kind is LocatorKind.BAND:
            band = Passband(f"{FILTER_SYSTEM.name}/{brightness.locator}")

        vegaspec = None
        if (
            brightness.amount_kind is AmountKind.MAG
            and brightness.system is PhotometricSystem.VEGA
        ):
            vegaspec = _vega_reference()

        return synphot_flux_scale(
            spectrum, brightness, band=band, vegaspec=vegaspec
        )

    def _select_anchor_sed(self, spectrum: SourceSpectrum) -> SourceSpectrum:
        """Return the SED the flux scale should be applied to, per ``anchor``.

        ``observed`` scales the reddened SED; ``intrinsic``/``absolute`` scale
        the unextincted SED (screens are then applied afterwards). Extinction
        screens are not yet wired in, so ``spectrum`` here is already the
        (un-screened) SED and every branch returns it unchanged -- this is the
        single seam where the two frames will diverge once the ``extinction``
        attribute lands, so the callers need no further change then.
        """
        # TODO(extinction): observed -> pass the screened SED; intrinsic/
        # absolute -> pass the unextincted SED (screens applied after scaling).
        return spectrum

    def _anchored_spectrum_scale(
        self,
        spectrum: SourceSpectrum,
        brightness: Brightness,
    ) -> float:
        """Flux scale for ``spectrum``, incorporating the ``anchor`` frame.

        Wraps the pure photometric scale (:meth:`_get_spectrum_scale`) with the
        anchor-frame semantics that need target context (``anchor`` and
        ``position.distance``):

        * ``observed`` / ``intrinsic`` -- the photometric scale as-is (they
          differ only in which SED is scaled, handled by
          :meth:`_select_anchor_sed`; identical until extinction is wired in).
        * ``absolute`` -- the value is the *absolute* amount (unextincted, at
          10 pc), so the photometric scale is multiplied by the achromatic
          inverse-square distance factor ``(10 pc / d)**2``. This applies
          uniformly to every amount kind. Guards: a missing distance is E10; a
          surface-brightness amount is E11 (surface brightness is
          distance-invariant).

        ``brightness`` is the (possibly SB-reduced) amount handed to the
        scaler; the E11 surface-brightness guard is keyed on the *original*
        ``self.brightness`` so a profile's SB->integrated reduction does not
        hide it.
        """
        anchor = self.anchor
        if anchor is AnchorFrame.ABSOLUTE:
            if self.brightness.is_surface_brightness:
                raise BrightnessError(
                    "E11",
                    "anchor: absolute with a surface-brightness amount -- "
                    "surface brightness is distance-invariant",
                )
            distance = self._distance_or_none()
            if distance is None:
                raise BrightnessError(
                    "E10",
                    "anchor: absolute requires a position.distance (the "
                    "distance modulus has no default)",
                )

        sed = self._select_anchor_sed(spectrum)
        scale = self._get_spectrum_scale(sed, brightness)

        if anchor is AnchorFrame.ABSOLUTE:
            distance_factor = float(
                ((10 * u.pc) / distance).to_value(u.dimensionless_unscaled) ** 2
            )
            scale *= distance_factor
        return scale


# TODO: docstring
def length_angle_equivalency(distance: Distance | u.Quantity[u.pc]):
    length_unit = u.AU
    angle_unit = u.arcsec

    distance = Distance(distance)

    def length_to_angle(length):
        # Equivalency functions receive and return only values
        angle = (length << length_unit) / distance
        return angle.to_value(angle_unit, u.dimensionless_angles())

    def angle_to_length(angle):
        # Equivalency functions receive and return only values
        length = (angle << angle_unit) * distance
        return length.to_value(length_unit, u.dimensionless_angles())

    return [(length_unit, angle_unit, length_to_angle, angle_to_length)]


# TODO: docstring
# TODO: better name??
def length_angle_context(distance: Distance | u.Quantity[u.pc]):
    return u.set_enabled_equivalencies(length_angle_equivalency(distance))
