# -*- coding: utf-8 -*-
"""Parsing and normalization of the ``brightness`` target attribute.

``brightness`` is the single **flux authority** of a target (see
``defining_brightness.md``). It is a two-slot value -- *where* (a band,
wavelength or frequency locator) and *how much* (a magnitude or flux amount) --
which this module normalizes into one frozen :class:`Brightness` value type.

This module only *classifies and stores*; it invents no flux and discards no
information, so the canonical mapping form can be reconstructed losslessly.
Conversion of magnitudes to linear flux (via synphot) and the actual spectrum
scaling happen later, at ``to_source`` time -- not here.

Dispatch principles (verified against astropy 8.0 / synphot 1.7):

* Locators are structurally disjoint: a band name starts with a letter and
  contains no whitespace, so anything that parses as a quantity is a
  wavelength or frequency, never a band.
* Amounts dispatch on the astropy **physical type** of the unit, never on the
  string shape -- *except* magnitudes, which are parsed here (see
  :func:`_parse_mag_string`) because the photometric *system* is metadata of a
  magnitude, not part of its unit, and because ``mag(Vega ...)`` is not an
  astropy function unit at all.
* A per-solid-angle divisor selects surface brightness; its absence selects
  integrated flux. Detection reuses ScopeSim's ``solid angle`` convention.

Wiring: ``SpectrumTarget._parse_brightness`` delegates to
:func:`parse_brightness`; ``typing_utils.BRIGHTNESS_TYPE`` becomes the accepted
input union below.
"""

import re
from dataclasses import dataclass
from enum import Enum, auto
from numbers import Number
from collections.abc import Mapping, Sequence

import astropy.units as u
from synphot.units import VEGAMAG

from .typing_utils import BRIGHTNESS_TYPE

__all__ = [
    "Brightness",
    "LocatorKind",
    "AmountKind",
    "PhotometricSystem",
    "AnchorFrame",
    "FromSpectralType",
    "BrightnessError",
    "parse_brightness",
    "solid_angle_unit",
]


class LocatorKind(Enum):
    BAND = auto()
    WAVELENGTH = auto()
    FREQUENCY = auto()


class AmountKind(Enum):
    MAG = auto()
    FLUX_DENSITY_NU = auto()  # spectral flux density (per frequency): Jy, ...
    # spectral flux density (per wavelength): erg/s/cm2/A
    FLUX_DENSITY_LAM = auto()
    ENERGY_FLUX = auto()  # band-/line-integrated energy flux: W/m2


class PhotometricSystem(Enum):
    VEGA = "Vega"
    AB = "AB"
    ST = "ST"


class AnchorFrame(Enum):
    """Frame the ``brightness`` value is anchored in (schema ``anchorFrame``).

    Controls which SED the flux-authority scale is applied to and whether a
    distance modulus enters (see ``defining_brightness.md`` /
    ``defining_extinction.md``):

    * ``OBSERVED`` (default) -- scale the *reddened* SED so photometry in the
      band matches the value (catalog/ETC convention; extinction then only
      changes colours).
    * ``INTRINSIC`` -- scale the *unextincted* SED, then apply screens
      (extinction dims *and* reddens).
    * ``ABSOLUTE`` -- like ``INTRINSIC`` but the value refers to 10 pc; the
      achromatic distance modulus from ``position.distance`` is applied before
      the screens.

    Because ``OBSERVED`` and ``INTRINSIC`` differ only in *when* extinction
    screens act, they are behaviourally identical until the extinction
    attribute is wired in; the distinction is kept here so screens slot into
    the documented place without a downstream change.
    """

    OBSERVED = "observed"
    INTRINSIC = "intrinsic"
    ABSOLUTE = "absolute"

    @classmethod
    def coerce(cls, value: "AnchorFrame | str | None") -> "AnchorFrame":
        """Normalize ``None``/str/enum into an :class:`AnchorFrame`."""
        if value is None:
            return cls.OBSERVED
        if isinstance(value, cls):
            return value
        try:
            return cls(str(value))
        except ValueError as exc:
            allowed = ", ".join(repr(m.value) for m in cls)
            raise ValueError(
                f"anchor must be one of {allowed}, not {value!r}"
            ) from exc


@dataclass(frozen=True)
class FromSpectralType:
    """Deferred ``brightness: {from_spectral_type: <table>}`` resolver marker.

    The pure parser cannot resolve this -- it needs the target's *spectrum*
    (a :class:`~astar_utils.SpectralType`) and its ``position.distance``, which
    live on the target. So :func:`parse_brightness` returns this marker and
    :class:`~.target.SpectrumTarget` performs the (network-backed) lookup at
    brightness-access time, recording the table name/version for reproducible
    exports (same resolver pattern as extinction ``from_map``). ``band``
    selects the absolute-magnitude column (default ``V`` -> ``M_V``).
    """

    table: str
    band: str = "V"


class BrightnessError(ValueError):
    """A brightness parse/consistency error, tagged with its matrix code."""

    _ERROR_DOCS = {
        "E6": "integrated amount on a profile without a finite analytic total; see the validity matrix",
        "E7": "surface-brightness amount on a point source; see the validity matrix",
        "E8": "reserved 'params' key (amplitude/x_0/y_0); see defining_brightness.md",
        "E9": "blackbody spectrum without a brightness to scale it; see defining_spectra.md",
        "E10": "'anchor: absolute' without a 'position.distance'; see the absolute-magnitudes section",
        "E11": "'anchor: absolute' with a surface-brightness amount (surface brightness is distance-invariant); see the absolute-magnitudes section",
        "E12": "'from_spectral_type' resolver on a non-SpectralType spectrum; see the spectral-type resolver section",
    }
    code = None

    def __init__(self, code=None, detail="") -> None:
        if self.code is not None:
            super().__init__(code)
            return

        self.code = code
        body = detail or self._ERROR_DOCS[code]
        pointer = "" if detail == "" else f" ({self._ERROR_DOCS[code]})"
        super().__init__(f"[{code}] {body}{pointer}")


class AmountError(BrightnessError):
    """Amount unit of unrecognized physical type; see the dispatch table in
    defining_brightness.md"""

    code = "E1"


class LocatorError(BrightnessError):
    """Locator is neither band-shaped nor a length/frequency quantity; see the
    locator grammar"""

    code = "E2"


class MissingBandError(BrightnessError):
    """Magnitude requires a band locator."""

    code = "E3"


class SystemMismatchError(BrightnessError):
    """'system' is only meaningful for magnitude amounts.

    See the photometric-systems section."""

    code = "E4"


class SystemDuplicateError(BrightnessError):
    """'system' given both as a field and embedded in the amount string.

    See the photometric-systems section."""

    code = "E5"


def solid_angle_unit(unit: u.UnitBase) -> u.UnitBase | None:
    """Return the solid-angle unit (``sr``, ``arcsec2``, ...) of a per-Omega
    unit, or ``None`` if `unit` carries no per-solid-angle divisor."""
    try:
        bases, powers = unit.bases, unit.powers
    except AttributeError:
        return None
    for base, power in zip(bases, powers):
        if 1 / (base**power).physical_type == "solid angle":
            return (base**power) ** -1
    return None


@dataclass(frozen=True)
class Brightness:
    """One normalized brightness specification (single entry in v1).

    Every accepted input form maps into this. ``value`` is always a plain
    linear :class:`~astropy.units.Quantity` for flux amounts, or a plain
    ``... mag`` Quantity for magnitudes -- the photometric ``system`` is carried
    separately (astropy's plain ``u.mag`` has no zero point of its own), and the
    synphot function unit (VEGAMAG/ABmag/STmag) is reconstructed only at
    conversion time.
    """

    locator_kind: LocatorKind
    locator: str | u.Quantity
    amount_kind: AmountKind
    value: u.Quantity
    system: PhotometricSystem = PhotometricSystem.VEGA
    solid_angle: u.UnitBase | None = None

    @property
    def is_surface_brightness(self) -> bool:
        return self.solid_angle is not None


# Band: letter-leading, no whitespace (schema 'band' pattern) -> disjoint from
# any quantity string.
_BAND_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_'*+-]{0,11}$")

# Number token shared by mag/quantity shapes.
_NUM = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"

# Magnitude amount: '<num> mag[(SYS)]' with an optional '/ (arcsec2|sr)' tail.
# The system lives inside mag(...), never combined with the divisor -- we parse
# it ourselves, so mag(Vega)/arcsec2 works even though astropy rejects it.
_MAG_SB_TAIL = re.compile(r"\s*/\s*(arcsec2|sr)\s*$")
_MAG_CORE = re.compile(rf"^\s*({_NUM})\s+mag(?:\((AB|ST|Vega|VEGA)\))?\s*$")

_SYSTEM = {
    "AB": PhotometricSystem.AB,
    "ST": PhotometricSystem.ST,
    "Vega": PhotometricSystem.VEGA,
    "VEGA": PhotometricSystem.VEGA,
}


@dataclass(frozen=True)
class _Amount:
    kind: AmountKind
    value: u.Quantity
    system: PhotometricSystem | None  # None => amount did not state one
    solid_angle: u.UnitBase | None


def _looks_like_magnitude(s: str) -> bool:
    return bool(_MAG_CORE.match(_MAG_SB_TAIL.sub("", s)))


def _parse_mag_string(s: str) -> _Amount:
    """Parse a magnitude amount string into (value, system, solid_angle)."""
    solid = None
    tail = _MAG_SB_TAIL.search(s)
    if tail:
        solid = u.Unit(tail.group(1))
        s = s[: tail.start()]
    core = _MAG_CORE.match(s)
    if core is None:  # pragma: no cover - guarded by caller
        raise AmountError(f"could not parse magnitude amount {s!r}")
    value = float(core.group(1)) * u.mag
    system = _SYSTEM[core.group(2)] if core.group(2) else None
    return _Amount(AmountKind.MAG, value, system, solid)


_FLUX_KINDS: dict[u.physical.PhysicalType, AmountKind] = {
    u.get_physical_type("spectral flux density"): AmountKind.FLUX_DENSITY_NU,
    u.get_physical_type("spectral flux density wav"): AmountKind.FLUX_DENSITY_LAM,
    u.get_physical_type("energy flux"): AmountKind.ENERGY_FLUX,
    u.get_physical_type("surface brightness"): AmountKind.FLUX_DENSITY_NU,  # per-freq SB
    u.get_physical_type("surface brightness wav"): AmountKind.FLUX_DENSITY_LAM,  # per-wav SB
    u.get_physical_type("radiance"): AmountKind.ENERGY_FLUX,  # integrated-flux SB
}


def _flux_kind(unit: u.UnitBase) -> AmountKind:
    kind = _FLUX_KINDS.get(unit.physical_type)
    if kind is None:
        raise AmountError(
            f"unrecognized amount physical type {unit.physical_type!s} for "
            f"unit {unit!s}",
        )
    return kind


def _parse_amount(amount: object) -> _Amount:
    """Normalize a *how-much* value: Number | Quantity | str -> _Amount."""
    # bare number -> Vega magnitude (band-locator cross-check happens later, E3)
    # TODO: use match here again
    if isinstance(amount, Number) and not isinstance(amount, bool):
        return _Amount(AmountKind.MAG, float(amount) * u.mag, None, None)

    if isinstance(amount, str):
        if _looks_like_magnitude(amount):
            return _parse_mag_string(amount)
        try:
            q = u.Quantity(amount)  # rejects function units (by design)
        except (TypeError, ValueError) as exc:
            raise AmountError(f"could not parse amount {amount!r}: {exc}") from exc
        return _parse_amount(q)

    if isinstance(amount, u.Quantity):
        unit = amount.unit
        # a plain-mag Quantity (resolver-fired 'N mag', or explicit mag SB)
        mag_units = {u.mag, VEGAMAG, u.ABmag, u.STmag}
        if u.mag in getattr(unit, "bases", [unit]) or unit in mag_units:
            solid = solid_angle_unit(unit)
            try:
                system = _SYSTEM[unit.physical_unit.to_string()]
            except AttributeError:
                system = PhotometricSystem.VEGA  # plain mag -> assume Vega
            return _Amount(AmountKind.MAG, amount.value * u.mag, system, solid)

        kind = _flux_kind(unit)
        return _Amount(kind, amount, None, solid_angle_unit(unit))

    raise AmountError(f"unsupported amount type {type(amount).__name__}")


def _parse_locator(locator: object) -> tuple[LocatorKind, str | u.Quantity]:
    """Normalize a *where* value: band string | wavelength/frequency."""
    # TODO: use match here again
    if isinstance(locator, str) and _BAND_RE.match(locator):
        return LocatorKind.BAND, locator
    if isinstance(locator, str):
        try:
            locator = u.Quantity(locator)
        except (TypeError, ValueError) as exc:
            raise LocatorError(
                f"locator {locator!r} is not a band or quantity: {exc}"
            ) from exc
    if isinstance(locator, u.Quantity):
        pt = locator.unit.physical_type
        if pt == "length":
            return LocatorKind.WAVELENGTH, locator
        if pt == "frequency":
            return LocatorKind.FREQUENCY, locator
        raise LocatorError(
            f"locator quantity has physical type {str(pt)!r} "
            "(need length or frequency)",
        )
    raise LocatorError(
        f"unsupported locator type {type(locator).__name__}"
    )


def parse_brightness(brightness: BRIGHTNESS_TYPE) -> "Brightness | FromSpectralType":
    """Normalize any accepted `brightness` input into a :class:`Brightness`.

    Accepts the tuple sugar ``(locator, amount)`` and the canonical mapping
    ``{band|wavelength|frequency: ..., value: ..., system: ...}``. Both slots
    accept ``str``, :class:`~astropy.units.Quantity` or (for the amount) a bare
    number, so behaviour is independent of whether the YAML implicit resolver
    fired on a scalar.
    """
    # shape dispatch
    # TODO: use match here again
    if isinstance(brightness, Mapping):
        if (table := brightness.get("from_spectral_type")) is not None:
            if not isinstance(table, str):
                raise LocatorError(
                    "from_spectral_type must name a lookup table (a string), "
                    f"not {type(table).__name__}",
                )
            if extra := brightness.keys() - {"from_spectral_type", "band"}:
                raise LocatorError(
                    "from_spectral_type takes only an optional 'band', "
                    f"got unexpected keys {sorted(extra)}",
                )
            return FromSpectralType(
                table=table, band=brightness.get("band", "V")
            )

        locator_key = {"band", "wavelength", "frequency"} & brightness.keys()
        if len(locator_key) != 1 or "value" not in brightness:
            raise LocatorError(
                "canonical brightness needs exactly one of "
                "band/wavelength/frequency plus a 'value'",
            )
        (key,) = locator_key
        loc_kind, loc = (
            _parse_locator(brightness[key])
            if key != "band"
            else (LocatorKind.BAND, brightness["band"])
        )
        # a band given via the 'band' key is trusted as a band even if 1 char
        amount = _parse_amount(brightness["value"])
        field_system = brightness.get("system")

    elif isinstance(brightness, Sequence) and not isinstance(brightness, str):
        if len(brightness) != 2:
            raise LocatorError(
                "tuple brightness must be exactly (locator, amount)"
            )
        loc_kind, loc = _parse_locator(brightness[0])
        amount = _parse_amount(brightness[1])
        field_system = None

    else:
        raise TypeError(
            "brightness must be a (locator, amount) tuple or a mapping, "
            f"not {type(brightness).__name__}"
        )

    # cross-field consistency (E3/E4/E5)
    system = amount.system
    if field_system is not None:
        if amount.kind is not AmountKind.MAG:
            raise SystemMismatchError()
        if amount.system is not None:
            raise SystemDuplicateError()
        system = (
            field_system
            if isinstance(field_system, PhotometricSystem)
            else _SYSTEM[str(field_system)]
        )

    if amount.kind is AmountKind.MAG:
        if loc_kind is not LocatorKind.BAND:
            raise MissingBandError()
        if system is None:
            system = PhotometricSystem.VEGA
    else:
        system = PhotometricSystem.VEGA  # unused for non-mag; kept as default

    return Brightness(
        locator_kind=loc_kind,
        locator=loc,
        amount_kind=amount.kind,
        value=amount.value,
        system=system,
        solid_angle=amount.solid_angle,
    )
