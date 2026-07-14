# -*- coding: utf-8 -*-
"""Line-of-sight (foreground) dust extinction (see ``docs/defining_extinction.md``).

This module owns the *structural* half of the ``extinction`` attribute: the law
vocabulary, the single-screen grammar (bare ``A_V`` sugar / canonical mapping),
its resolution to one ``(A_V, law, R_V)`` triple, and the ``synphot`` reddening
curve. It deliberately does **not** know about targets, inheritance, or the
source table -- that wiring lives in :mod:`.target` and the point-/extended-
source classes.

``extinction`` is the **sightline**: one foreground screen between the observer
and the scene, uniform over a component, cascading from the scene. There is no
list form -- for a fixed law and R_V two screens compose exactly into one
(``A(lambda)`` is linear in ``A_V``), and dust *at the source* (natal cloud,
circumstellar envelope) is a different physical quantity with its own attribute,
:mod:`.attenuation`. The explicit opt-out is ``"0 mag"``.

Normative pipeline order (``defining_spectra.md``): rest SED -> attenuation
(source frame) -> [distance modulus] -> redshift -> extinction (observer frame)
-> flux anchoring. Because the screen is uniform, it is applied to the spectrum
*in place* for every field type -- point-source spectrum dedup is unaffected
and no table column is emitted for it. The flux *anchor* only decides whether
the scale is computed against the reddened SED (``observed``) or the intrinsic
one (``intrinsic``/``absolute``).
"""

import re
from dataclasses import dataclass
from enum import Enum
from numbers import Real
from collections.abc import Mapping

import numpy as np
import astropy.units as u
from synphot import SourceSpectrum, SpectralElement
from synphot.models import Empirical1D

__all__ = [
    "ExtinctionLaw",
    "ExtinctionScreen",
    "FromMap",
    "ExtinctionError",
    "parse_extinction",
    "resolve_extinction",
    "transmission_element",
    "redden",
    "PLAIN_MAG",
    "as_plain_mag",
]

# Bare plain-magnitude quantity, e.g. "2.3 mag" (schema plainMagString). No
# photometric system: extinction amounts are flux *ratios*, so no zero point.
PLAIN_MAG = re.compile(
    r"^[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?\smag$"
)


class ExtinctionError(ValueError):
    """Structured extinction error with a stable ``code`` (E13-E18).

    :class:`~.attenuation.AttenuationError` extends the table with E19-E22."""


class EmptyExtError(ExtinctionError):
    """E14, empty extinction screen (needs value, ebv, or from_map)"""


class UnphysExtError(ExtinctionError):
    """E16, non-physical extinction parameter (A_V/E(B-V) < 0 or R_V <= 0)"""


class ExtinctionLaw(Enum):
    """Reddening law, pinned to the ``dust_extinction`` package (schema enum)."""

    CCM89 = "ccm89"
    F99 = "f99"
    G23 = "g23"

    @classmethod
    def coerce(cls, value: "ExtinctionLaw | str | None") -> "ExtinctionLaw":
        if value is None:
            return cls.CCM89  # matches ESO-371803 5.3.7 prose
        if isinstance(value, cls):
            return value
        try:
            return cls(str(value).lower())
        except ValueError as exc:
            raise ExtinctionError(
                "law must be one of "
                f"{', '.join(repr(m.value) for m in cls)}, not {value!r}"
            ) from exc

    def model(self, rv: float):
        """Return the ``dust_extinction`` model instance for this law + R_V."""
        # Lazy import: keeps the heavy optional dependency off the module path.
        from dust_extinction.parameter_averages import CCM89, F99, G23

        return {self.CCM89: CCM89, self.F99: F99, self.G23: G23}[self](Rv=rv)


@dataclass(frozen=True)
class FromMap:
    """Deferred ``extinction: {from_map: <id>}`` resolver marker.

    Directional dust maps (e.g. ``bayestar2019``) need the target's *position*
    and external map data, so resolution is deferred exactly like the brightness
    ``from_spectral_type`` marker. Not yet implemented (raises E18 on resolve).
    """

    map_id: str
    law: "ExtinctionLaw" = ExtinctionLaw.CCM89
    rv: float = 3.1


@dataclass(frozen=True)
class ExtinctionScreen:
    """One line-of-sight screen, normalized. Hashable (plain fields only).

    Exactly one currency is set: ``av`` (A in ``band``, magnitudes) or ``ebv``
    (E(B-V)). ``band`` defaults to ``V`` -- with the default, ``av`` is the
    model's V-normalized A_V directly and needs no bandpass. A non-V band needs
    that band's effective wavelength to convert A_band -> A_V, supplied at
    resolve time (see :func:`resolve_extinction`).
    """

    av: float | None = None
    band: str = "V"
    ebv: float | None = None
    law: ExtinctionLaw = ExtinctionLaw.CCM89
    rv: float = 3.1


def _parse_screen(spec: object) -> ExtinctionScreen | FromMap:
    """Normalize one screen (bare A_V string or mapping) to a screen/marker."""
    if isinstance(spec, str):
        if not PLAIN_MAG.match(spec):
            raise EmptyExtError(
                f"bare extinction must be a plain 'N mag' A_V, got {spec!r}",
            )
        return ExtinctionScreen(av=float(spec.split()[0]))  # band V, ccm89, 3.1

    if isinstance(spec, Real):
        # tolerate a bare number as A_V mag (parity with brightness sugar)
        return ExtinctionScreen(av=float(spec))

    if not isinstance(spec, Mapping):
        raise EmptyExtError(f"screen not understood: {spec!r}")

    law = ExtinctionLaw.coerce(spec.get("law"))
    rv = float(spec.get("rv", 3.1))
    if rv <= 0:
        raise UnphysExtError(f"R_V must be > 0, got {rv}")

    ext_keys = {"value", "ebv", "from_map"} & spec.keys()
    try:
        ext_key = ext_keys.pop()
    except KeyError as err:
        raise EmptyExtError(
            "empty extinction screen (needs value, ebv, or from_map)"
        ) from err
    if ext_keys:
        raise ExtinctionError(f"screen sets multiple of {ext_keys}")

    if ext_key == "from_map":
        return FromMap(map_id=str(spec[ext_key]), law=law, rv=rv)

    if ext_key == "ebv":
        ebv = float(spec[ext_key])
        if ebv < 0:
            raise UnphysExtError(f"E(B-V) must be >= 0, got {ebv}")
        return ExtinctionScreen(ebv=ebv, law=law, rv=rv)

    value = spec[ext_key]
    if not (isinstance(value, str) and PLAIN_MAG.match(value)):
        raise EmptyExtError(
            f"screen 'value' must be a plain 'N mag' amount, got {value!r}"
        )
    av = float(value.split()[0])
    if av < 0:
        raise UnphysExtError(f"A_V must be >= 0, got {av}")
    return ExtinctionScreen(av=av, band=str(spec.get("band", "V")), law=law, rv=rv)


# What ``extinction`` normalizes to: one screen, or a single deferred FromMap
# marker. ``ExtinctionScreen(av=0.0)`` is the explicit opt-out ("0 mag").
Extinction = "ExtinctionScreen | FromMap"


def parse_extinction(spec: object) -> "ExtinctionScreen | FromMap":
    """Normalize the ``extinction`` input into one screen (or a marker).

    * bare ``"2.3 mag"`` / number -> A_V screen (band V, ccm89, R_V 3.1);
      ``"0 mag"`` is the explicit no-extinction opt-out
    * mapping -> one screen (``value``+``band`` | ``ebv`` | ``from_map``)
    * list -> E17 (removed; see the module docstring)
    """
    if isinstance(spec, (list, tuple)):
        raise ExtinctionError(f"got a list of {len(spec)} screen(s)")
    return _parse_screen(spec)


def _screen_av(screen: ExtinctionScreen, *, band_wavelength=None) -> float:
    """A_V for one screen. ``ebv`` -> R_V*E(B-V); band V -> av directly; a
    non-V band needs its effective wavelength to divide out A(band)/A(V)."""
    if screen.ebv is not None:
        return screen.rv * screen.ebv
    if screen.band == "V":
        return screen.av
    if band_wavelength is None:
        raise EmptyExtError(
            f"A_{screen.band} needs the band effective wavelength to convert "
            "to A_V; supply band_wavelength_lookup",
        )
    return screen.av / float(screen.law.model(screen.rv)(band_wavelength))


def resolve_extinction(
    screen: "ExtinctionScreen | FromMap | None",
    *,
    band_wavelength_lookup=None,
) -> tuple[float, ExtinctionLaw, float] | None:
    """Resolve to one ``(A_V, law, R_V)`` triple, or ``None`` for *no screen*.

    ``None`` in means "nothing set" (inherit / no extinction); an explicit
    ``"0 mag"`` comes back as ``(0.0, law, R_V)`` -- both leave the SED untouched
    in :func:`redden`, but the distinction matters for the cascade.
    ``from_map`` is not yet resolvable (E18).

    ``band_wavelength_lookup`` is a ``band_name -> Quantity`` callable, needed
    only for a screen given in a non-V band.
    """
    if screen is None:
        return None
    if isinstance(screen, FromMap):
        raise ExtinctionError(
            "from_map extinction resolver is not yet implemented"
        )

    lam = None
    if screen.ebv is None and screen.band != "V":
        if band_wavelength_lookup is None:
            raise EmptyExtError(
                f"A_{screen.band} needs a band_wavelength_lookup to reach A_V"
            )
        lam = band_wavelength_lookup(screen.band)

    return _screen_av(screen, band_wavelength=lam), screen.law, screen.rv


def transmission_element(
    wave: u.Quantity, av: float, law: ExtinctionLaw, rv: float
) -> SpectralElement:
    """Dimensionless transmission T(lambda) as a ``synphot`` throughput.

    Wavelengths outside the law's validity range are clamped to the nearest
    in-range value (transmission held constant beyond validity) so a full-range
    SED can be multiplied without the ``dust_extinction`` model raising.
    """
    model = law.model(rv)
    x = 1.0 / wave.to_value(u.um)  # inverse ums
    wav_clamped = (1.0 / np.clip(x, *model.x_range)) * u.um
    trans = np.asarray(model.extinguish(wav_clamped, Av=av), dtype=float)
    return SpectralElement(Empirical1D, points=wave, lookup_table=trans)


def redden(
    spectrum: SourceSpectrum,
    av: float,
    law: ExtinctionLaw,
    rv: float
) -> SourceSpectrum:
    """Return ``spectrum`` reddened by the true (un-normalized) screen curve."""
    if av == 0:
        return spectrum
    return spectrum * transmission_element(spectrum.waveset, av, law, rv)


def as_plain_mag(value: object, name: str) -> float:
    """Coerce a magnitude-valued parameter to float (number or ``"N mag"``).

    Shared with :mod:`.attenuation` (distribution parameters, columns).
    """
    if isinstance(value, str):
        if not PLAIN_MAG.match(value):
            raise UnphysExtError(
                f"{name} must be a number or 'N mag', got {value!r}"
            )
        return float(value.split()[0])
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return float(value)
    raise UnphysExtError(f"{name} must be numeric, got {value!r}")
