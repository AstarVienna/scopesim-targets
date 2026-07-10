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


__all__ = [
    "Brightness",
    "LocatorKind",
    "AmountKind",
    "PhotometricSystem",
    "BrightnessError",
    "parse_brightness",
    "unit_includes_per_physical_type",
    "solid_angle_unit",
]


class LocatorKind(Enum):
    BAND = auto()
    WAVELENGTH = auto()
    FREQUENCY = auto()


class AmountKind(Enum):
    MAG = auto()
    FLUX_DENSITY_NU = auto()  # spectral flux density (per frequency): Jy, ...
    FLUX_DENSITY_LAM = (
        auto()
    )  # spectral flux density (per wavelength): erg/s/cm2/A
    ENERGY_FLUX = auto()  # band-/line-integrated energy flux: W/m2


class PhotometricSystem(Enum):
    VEGA = "Vega"
    AB = "AB"
    ST = "ST"


_ERROR_DOCS: dict[str, str] = {
    "E1": "amount unit of unrecognized physical type; see the dispatch table in defining_brightness.md",
    "E2": "locator is neither band-shaped nor a length/frequency quantity; see the locator grammar",
    "E3": "magnitude requires a band locator (monochromatic AB is a deferred discussion item)",
    "E4": "'system' is only meaningful for magnitude amounts; see the photometric-systems section",
    "E5": "'system' given both as a field and embedded in the amount string; see the photometric-systems section",
    "E6": "integrated amount on a profile without a finite analytic total; see the validity matrix",
    "E7": "surface-brightness amount on a point source; see the validity matrix",
    "E8": "reserved 'params' key (amplitude/x_0/y_0); see defining_brightness.md",
}


class BrightnessError(ValueError):
    """A brightness parse/consistency error, tagged with its matrix code."""

    def __init__(self, code: str, detail: str = "") -> None:
        self.code = code
        prefix = f"[{code}] "
        body = detail or _ERROR_DOCS[code]
        pointer = "" if detail == "" else f" ({_ERROR_DOCS[code]})"
        super().__init__(f"{prefix}{body}{pointer}")


def _err(code: str, detail: str = "") -> "BrightnessError":
    return BrightnessError(code, detail)


def unit_includes_per_physical_type(
    unit: u.UnitBase, physical_type: str
) -> bool:
    """True if one of ``unit``'s bases is of ``1 / physical_type``."""
    try:
        bases, powers = unit.bases, unit.powers
    except AttributeError:  # function units (e.g. ABmag) have no .bases
        return False
    return any(
        1 / (base**power).physical_type == physical_type
        for base, power in zip(bases, powers)
    )


def solid_angle_unit(unit: u.UnitBase) -> u.UnitBase | None:
    """Return the solid-angle unit (``sr``, ``arcsec2``, ...) of a per-Omega
    unit, or ``None`` if ``unit`` carries no per-solid-angle divisor."""
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

    # -- temporary compat shims -------------------------------------------
    # The current ParametrizedTarget.to_source reads .band and .mag off the old
    # namedtuple. These keep Sersic/Disk green until the Phase 2 to_source
    # rewrite consumes the full Brightness; remove them then.
    @property
    def band(self) -> str:
        if self.locator_kind is not LocatorKind.BAND:
            raise AttributeError("brightness has no band locator")
        return self.locator

    @property
    def mag(self) -> u.Quantity:
        if self.amount_kind is not AmountKind.MAG:
            raise AttributeError("brightness amount is not a magnitude")
        if self.system is not PhotometricSystem.VEGA:
            raise NotImplementedError(
                f"{self.system.value} magnitude scaling arrives in Phase 2; the "
                "legacy .mag shim represents Vega only (see "
                "brightness_implementation_plan.md)."
            )
        return self.value


# Band: letter-leading, no whitespace (schema 'band' pattern) -> disjoint from
# any quantity string.
_BAND_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_'*+-]{0,11}$")

# Number token shared by mag/quantity shapes.
_NUM = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"

# Magnitude amount: '<num> mag[(SYS)]' with an optional '/ (arcsec2|sr)' tail.
# The system lives inside mag(...), never combined with the divisor -- we parse
# it ourselves, so mag(Vega)/arcsec2 works even though astropy rejects it.
_MAG_SB_TAIL = re.compile(r"\s*/\s*(arcsec2|sr)\s*$")
_MAG_CORE = re.compile(rf"^\s*({_NUM})\s+mag(?:\((AB|ST|Vega)\))?\s*$")

_SYSTEM = {
    "AB": PhotometricSystem.AB,
    "ST": PhotometricSystem.ST,
    "Vega": PhotometricSystem.VEGA,
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
        raise _err("E1", f"could not parse magnitude amount {s!r}")
    value = float(core.group(1)) * u.mag
    system = _SYSTEM[core.group(2)] if core.group(2) else None
    return _Amount(AmountKind.MAG, value, system, solid)


# physical_type -> (AmountKind). Membership via astropy's alias-aware equality.
def _flux_kind(unit: u.UnitBase) -> AmountKind:
    pt = unit.physical_type
    if pt == "spectral flux density":
        return AmountKind.FLUX_DENSITY_NU
    if pt == "spectral flux density wav":
        return AmountKind.FLUX_DENSITY_LAM
    if pt == "energy flux" or pt == "irradiance":
        return AmountKind.ENERGY_FLUX
    if pt == "surface brightness":  # per-frequency SB
        return AmountKind.FLUX_DENSITY_NU
    if pt == "surface brightness wav":  # per-wavelength SB
        return AmountKind.FLUX_DENSITY_LAM
    if pt == "radiance":  # integrated-flux SB
        return AmountKind.ENERGY_FLUX
    raise _err(
        "E1",
        f"unrecognized amount physical type {str(pt)!r} for unit {unit!s}",
    )


def _parse_amount(amount: object) -> _Amount:
    """Normalize a *how-much* value: Number | Quantity | str -> _Amount."""
    # bare number -> Vega magnitude (band-locator cross-check happens later, E3)
    if isinstance(amount, Number) and not isinstance(amount, bool):
        return _Amount(AmountKind.MAG, float(amount) * u.mag, None, None)

    if isinstance(amount, str):
        if _looks_like_magnitude(amount):
            return _parse_mag_string(amount)
        try:
            q = u.Quantity(amount)  # rejects function units (by design)
        except (TypeError, ValueError) as exc:
            raise _err(
                "E1", f"could not parse amount {amount!r}: {exc}"
            ) from exc
        return _parse_amount(q)

    if isinstance(amount, u.Quantity):
        unit = amount.unit
        # a plain-mag Quantity (resolver-fired 'N mag', or explicit mag SB)
        if u.mag in getattr(unit, "bases", [unit]) or unit == u.mag:
            solid = solid_angle_unit(unit)
            return _Amount(AmountKind.MAG, amount.value * u.mag, None, solid)
        kind = _flux_kind(unit)
        return _Amount(kind, amount, None, solid_angle_unit(unit))

    raise _err("E1", f"unsupported amount type {type(amount).__name__}")


def _parse_locator(locator: object) -> tuple[LocatorKind, str | u.Quantity]:
    """Normalize a *where* value: band string | wavelength/frequency."""
    if isinstance(locator, str) and _BAND_RE.match(locator):
        return LocatorKind.BAND, locator
    if isinstance(locator, str):
        try:
            locator = u.Quantity(locator)
        except (TypeError, ValueError) as exc:
            raise _err(
                "E2", f"locator {locator!r} is not a band or quantity: {exc}"
            ) from exc
    if isinstance(locator, u.Quantity):
        pt = locator.unit.physical_type
        if pt == "length":
            return LocatorKind.WAVELENGTH, locator
        if pt == "frequency":
            return LocatorKind.FREQUENCY, locator
        raise _err(
            "E2",
            f"locator quantity has physical type {str(pt)!r} (need length or frequency)",
        )
    raise _err("E2", f"unsupported locator type {type(locator).__name__}")


def parse_brightness(brightness: object) -> Brightness:
    """Normalize any accepted ``brightness`` input into a :class:`Brightness`.

    Accepts the tuple sugar ``(locator, amount)`` and the canonical mapping
    ``{band|wavelength|frequency: ..., value: ..., system: ...}``. Both slots
    accept ``str``, :class:`~astropy.units.Quantity` or (for the amount) a bare
    number, so behaviour is independent of whether the YAML implicit resolver
    fired on a scalar.
    """
    # ---- shape dispatch --------------------------------------------------
    if isinstance(brightness, Mapping):
        if "from_spectral_type" in brightness:
            raise NotImplementedError(
                "brightness: {from_spectral_type: ...} is the Phase 5 resolver "
                "(see defining_brightness.md); not yet implemented."
            )
        locator_key = {"band", "wavelength", "frequency"} & brightness.keys()
        if len(locator_key) != 1 or "value" not in brightness:
            raise _err(
                "E2",
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
            raise _err(
                "E2", "tuple brightness must be exactly (locator, amount)"
            )
        loc_kind, loc = _parse_locator(brightness[0])
        amount = _parse_amount(brightness[1])
        field_system = None
    else:
        raise TypeError(
            "brightness must be a (locator, amount) tuple or a mapping, "
            f"not {type(brightness).__name__}"
        )

    # ---- cross-field consistency (E3/E4/E5) ------------------------------
    system = amount.system
    if field_system is not None:
        if amount.kind is not AmountKind.MAG:
            raise _err("E4")
        if amount.system is not None:
            raise _err("E5")
        system = (
            field_system
            if isinstance(field_system, PhotometricSystem)
            else _SYSTEM[str(field_system)]
        )

    if amount.kind is AmountKind.MAG:
        if loc_kind is not LocatorKind.BAND:
            raise _err("E3")
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
