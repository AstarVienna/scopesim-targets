# -*- coding: utf-8 -*-
"""Line-of-sight dust extinction (see ``docs/defining_extinction.md``).

This module owns the *structural* half of the ``extinction`` attribute: the law
vocabulary, the screen grammar (bare ``A_V`` sugar / canonical mapping / list),
resolution of a screen list to a single ``(A_V, law, R_V)`` triple, and the
``synphot`` reddening curve. It deliberately does **not** know about targets,
inheritance, or the source table -- that wiring lives in :mod:`.target` and the
point-/extended-source classes.

Normative pipeline order (``defining_spectra.md``): rest SED -> redshift ->
extinction screens (observer frame) -> flux anchoring. Screens compose
multiplicatively; the reddening curve is applied to the SED, while the flux
*anchor* only decides whether the scale is computed against the reddened SED
(``observed``) or the intrinsic one (``intrinsic``/``absolute``).

Delivery differs by field type (a scale/memory decision, not a physics one):
extended sources redden their single SED in place, while point sources keep the
deduplicated intrinsic spectra and emit a per-row ``Av`` column (+ ``law``/``rv``
in ``table.meta``) for ScopeSim to apply downstream. Either way the numbers here
are the single source of truth.
"""

import re
from dataclasses import dataclass
from enum import Enum
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
    "ExtinctionDistribution",
    "parse_extinction_distribution",
]

# Bare plain-magnitude quantity, e.g. "2.3 mag" (schema plainMagString). No
# photometric system: extinction amounts are flux *ratios*, so no zero point.
_PLAIN_MAG = re.compile(
    r"^[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?\smag$"
)


class ExtinctionError(ValueError):
    """Structured extinction error with a stable ``code`` (E13-E18).

    Mirrors :class:`~.brightness.BrightnessError`: the code is the contract, the
    detail is human-facing, and the doc-pointer suffix is looked up by code.
    """

    _ERROR_DOCS = {
        "E13": "a screen takes exactly one of value / ebv / from_map",
        "E14": "empty extinction screen (needs value, ebv, or from_map)",
        "E15": "unknown extinction law (allowed: ccm89, f99, g23)",
        "E16": "non-physical extinction parameter (A_V/E(B-V) < 0 or R_V <= 0)",
        "E17": "a screen list mixes laws or R_V values (v1: one law/R_V per table)",
        "E18": "from_map extinction resolver is not yet implemented",
    }

    def __init__(self, code: str, detail: str = ""):
        self.code = code
        body = detail or self._ERROR_DOCS[code]
        pointer = "" if detail == "" else f" ({self._ERROR_DOCS[code]})"
        super().__init__(f"[{code}] {body}{pointer}")


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
            allowed = ", ".join(repr(m.value) for m in cls)
            raise ExtinctionError(
                "E15", f"law must be one of {allowed}, not {value!r}"
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


# What ``extinction`` normalizes to: an ordered tuple of screens (empty tuple =
# explicit opt-out), or a single deferred FromMap marker.
Extinction = "tuple[ExtinctionScreen, ...] | FromMap"


def _parse_screen(spec: object) -> ExtinctionScreen | FromMap:
    """Normalize one screen (bare A_V string or mapping) to a screen/marker."""
    if isinstance(spec, str):
        if not _PLAIN_MAG.match(spec):
            raise ExtinctionError(
                "E14",
                f"bare extinction sugar must be a plain 'N mag' A_V, got {spec!r}",
            )
        return ExtinctionScreen(av=float(spec.split()[0]))  # band V, ccm89, 3.1

    if isinstance(spec, (int, float)) and not isinstance(spec, bool):
        # tolerate a bare number as A_V mag (parity with brightness sugar)
        return ExtinctionScreen(av=float(spec))

    if isinstance(spec, dict):
        law = ExtinctionLaw.coerce(spec.get("law"))
        rv = float(spec.get("rv", 3.1))
        if rv <= 0:
            raise ExtinctionError("E16", f"R_V must be > 0, got {rv}")

        currencies = [k for k in ("value", "ebv", "from_map") if k in spec]
        if len(currencies) == 0:
            raise ExtinctionError("E14")
        if len(currencies) > 1:
            raise ExtinctionError(
                "E13", f"screen sets multiple of {currencies}"
            )

        if "from_map" in spec:
            return FromMap(map_id=str(spec["from_map"]), law=law, rv=rv)

        if "ebv" in spec:
            ebv = float(spec["ebv"])
            if ebv < 0:
                raise ExtinctionError("E16", f"E(B-V) must be >= 0, got {ebv}")
            return ExtinctionScreen(ebv=ebv, law=law, rv=rv)

        value = spec["value"]
        if not (isinstance(value, str) and _PLAIN_MAG.match(value)):
            raise ExtinctionError(
                "E14", f"screen 'value' must be a plain 'N mag' amount, got {value!r}"
            )
        av = float(value.split()[0])
        if av < 0:
            raise ExtinctionError("E16", f"A_V must be >= 0, got {av}")
        return ExtinctionScreen(av=av, band=str(spec.get("band", "V")), law=law, rv=rv)

    raise ExtinctionError("E14", f"screen not understood: {spec!r}")


def parse_extinction(spec: object) -> "tuple[ExtinctionScreen, ...] | FromMap":
    """Normalize the ``extinction`` input into a tuple of screens (or a marker).

    * bare ``"2.3 mag"`` / number -> one A_V screen (band V, ccm89, R_V 3.1)
    * mapping -> one screen (``value``+``band`` | ``ebv`` | ``from_map``)
    * list -> the screens in order (``[]`` = explicit opt-out -> empty tuple)

    A ``from_map`` screen is returned as a bare :class:`FromMap` marker (a map
    is a whole-line-of-sight quantity, not one item in a multiplicative stack).
    """
    if isinstance(spec, list):
        screens = tuple(_parse_screen(item) for item in spec)
        maps = [s for s in screens if isinstance(s, FromMap)]
        if maps:
            raise ExtinctionError(
                "E13", "from_map cannot be combined in a screen list"
            )
        return screens
    parsed = _parse_screen(spec)
    return parsed if isinstance(parsed, FromMap) else (parsed,)


def _screen_av(screen: ExtinctionScreen, *, band_wavelength=None) -> float:
    """A_V for one screen. ``ebv`` -> R_V*E(B-V); band V -> av directly; a
    non-V band needs its effective wavelength to divide out A(band)/A(V)."""
    if screen.ebv is not None:
        return screen.rv * screen.ebv
    if screen.band == "V":
        return screen.av
    if band_wavelength is None:
        raise ExtinctionError(
            "E14",
            f"A_{screen.band} needs the band effective wavelength to convert "
            "to A_V; supply band_wavelength_lookup",
        )
    a_over_av = float(screen.law.model(screen.rv)(band_wavelength))
    return screen.av / a_over_av


def resolve_extinction(
    screens: "tuple[ExtinctionScreen, ...] | FromMap | None",
    *,
    band_wavelength_lookup=None,
) -> tuple[float, ExtinctionLaw, float] | None:
    """Collapse screens to one ``(A_V, law, R_V)`` triple, or ``None`` for none.

    v1 restriction (matches the single ``law``/``rv`` stored in ``table.meta``):
    every screen in a list must share one law *and* one R_V; A_V then adds
    (screens compose multiplicatively -> transmissions multiply -> A_V sums).
    Mixed laws/R_V raise E17. ``from_map`` is not yet resolvable (E18).

    ``band_wavelength_lookup`` is a ``band_name -> Quantity`` callable used only
    for screens given in a non-V band (each screen resolves its own band).
    """
    if screens is None:
        return None
    if isinstance(screens, FromMap):
        raise ExtinctionError("E18")
    if len(screens) == 0:
        return None  # explicit opt-out

    laws = {s.law for s in screens}
    rvs = {s.rv for s in screens}
    if len(laws) > 1 or len(rvs) > 1:
        raise ExtinctionError(
            "E17", f"list mixes laws={ {l.value for l in laws} } / rv={rvs}"
        )
    law, rv = screens[0].law, screens[0].rv

    def av_of(screen: ExtinctionScreen) -> float:
        lam = None
        if screen.ebv is None and screen.band != "V":
            if band_wavelength_lookup is None:
                raise ExtinctionError(
                    "E14",
                    f"A_{screen.band} needs a band_wavelength_lookup to reach A_V",
                )
            lam = band_wavelength_lookup(screen.band)
        return _screen_av(screen, band_wavelength=lam)

    return sum(av_of(s) for s in screens), law, rv


def transmission_element(
    wave: u.Quantity, av: float, law: ExtinctionLaw, rv: float
) -> SpectralElement:
    """Dimensionless transmission T(lambda) as a ``synphot`` throughput.

    Wavelengths outside the law's validity range are clamped to the nearest
    in-range value (transmission held constant beyond validity) so a full-range
    SED can be multiplied without the ``dust_extinction`` model raising.
    """
    model = law.model(rv)
    x = 1.0 / wave.to_value(u.micron)  # inverse microns
    lo, hi = model.x_range
    wav_clamped = (1.0 / np.clip(x, lo, hi)) * u.micron
    trans = np.asarray(model.extinguish(wav_clamped, Av=av), dtype=float)
    return SpectralElement(Empirical1D, points=wave, lookup_table=trans)


def redden(
    spectrum: SourceSpectrum, av: float, law: ExtinctionLaw, rv: float
) -> SourceSpectrum:
    """Return ``spectrum`` reddened by the true (un-normalized) screen curve."""
    if av == 0:
        return spectrum
    grid = spectrum.waveset
    return spectrum * transmission_element(grid, av, law, rv)


# --------------------------------------------------------------------------- #
# Per-star extinction distribution (generative targets)                       #
# --------------------------------------------------------------------------- #
# Same draw-from-PDF pattern as the population/morphology sampling: at cluster
# expansion each star gets its own A_V, which becomes an ordinary per-row `Av`
# column (the single-value screen path is the N=1, zero-variance special case).
# Zeroth-order screen model: no differential extinction within the PSF and no
# A_V-vs-position correlation (a future extension shares the same RNG).


def _as_mag(value: object, name: str) -> float:
    """Coerce a distribution parameter to a magnitude float (accepts 'N mag')."""
    if isinstance(value, str):
        if not _PLAIN_MAG.match(value):
            raise ExtinctionError(
                "E16", f"{name} must be a number or 'N mag', got {value!r}"
            )
        return float(value.split()[0])
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return float(value)
    raise ExtinctionError("E16", f"{name} must be numeric, got {value!r}")


@dataclass(frozen=True)
class ExtinctionDistribution:
    """A per-star A_V distribution for generative targets (e.g. clusters).

    ``kind`` is one of ``"column_density_pdf"`` / ``"lognormal"`` / ``"uniform"``;
    ``params`` are interpreted per-kind (all A_V values in magnitudes, plain
    numbers or ``"N mag"`` strings). ``law`` and ``rv`` accompany the realized
    A_V values exactly as for a fixed screen -- one law / R_V per table.

    Parameters by kind
    ------------------
    * ``uniform``: ``low``, ``high`` (A_V bounds).
    * ``lognormal``: ``av_median`` (median A_V), ``sigma`` (width in ln A_V).
    * ``column_density_pdf``: a lognormal core (``av_median``, ``sigma``) that
      hands over to a power-law tail ``p(A_V) ∝ A_V**-slope`` above ``av_break``,
      continuous at the break (Alves, Lombardi & Lada 2017, A&A 606, L2). The
      cloud's gas-to-dust conversion is folded into the choice of these A_V-space
      parameters. Optional ``av_min``/``av_max`` bound the support (defaults:
      ``av_median/100`` and ``20*av_break``); sampling is grid inverse-CDF.
    """

    kind: str
    params: Mapping[str, object]
    law: ExtinctionLaw = ExtinctionLaw.CCM89
    rv: float = 3.1

    _KINDS = ("column_density_pdf", "lognormal", "uniform")

    def __post_init__(self):
        if self.kind not in self._KINDS:
            raise ExtinctionError(
                "E15",
                f"distribution must be one of {self._KINDS}, not {self.kind!r}",
            )
        if self.rv <= 0:
            raise ExtinctionError("E16", f"R_V must be > 0, got {self.rv}")

    def sample(self, n: int, rng: "np.random.Generator") -> "np.ndarray":
        """Draw ``n`` non-negative A_V values (magnitudes) from ``rng``.

        ``rng`` is supplied by the caller so the whole realization stays
        reproducible under the cluster's single seed.
        """
        if n < 0:
            raise ValueError("n must be >= 0")
        p = self.params
        if self.kind == "uniform":
            low = _as_mag(p["low"], "low")
            high = _as_mag(p["high"], "high")
            if not 0 <= low <= high:
                raise ExtinctionError("E16", f"need 0 <= low <= high, got {low}, {high}")
            return rng.uniform(low, high, size=n)

        if self.kind == "lognormal":
            median = _as_mag(p["av_median"], "av_median")
            sigma = _as_mag(p["sigma"], "sigma")
            if median <= 0 or sigma <= 0:
                raise ExtinctionError("E16", "av_median and sigma must be > 0")
            return rng.lognormal(mean=np.log(median), sigma=sigma, size=n)

        return self._sample_column_density_pdf(n, rng)

    def _sample_column_density_pdf(self, n, rng):
        p = self.params
        median = _as_mag(p["av_median"], "av_median")
        sigma = _as_mag(p["sigma"], "sigma")
        av_break = _as_mag(p["av_break"], "av_break")
        slope = _as_mag(p["slope"], "slope")
        if median <= 0 or sigma <= 0 or av_break <= 0:
            raise ExtinctionError("E16", "av_median, sigma, av_break must be > 0")
        if slope <= 1:
            raise ExtinctionError(
                "E16", f"tail slope must be > 1 (integrable), got {slope}"
            )
        av_min = _as_mag(p.get("av_min", median / 100.0), "av_min")
        av_max = _as_mag(p.get("av_max", 20.0 * av_break), "av_max")
        if not 0 < av_min < av_max:
            raise ExtinctionError("E16", "need 0 < av_min < av_max")

        grid = np.linspace(av_min, av_max, 8192)

        def lognormal_pdf(a):
            return (1.0 / (a * sigma * np.sqrt(2 * np.pi))) * np.exp(
                -((np.log(a) - np.log(median)) ** 2) / (2 * sigma**2)
            )

        core = lognormal_pdf(grid)
        # Power-law tail beyond the break, scaled for continuity at av_break.
        tail = lognormal_pdf(av_break) * (grid / av_break) ** (-slope)
        pdf = np.where(grid <= av_break, core, tail)

        cdf = np.cumsum(pdf)
        cdf /= cdf[-1]
        # Deduplicate flat CDF stretches so np.interp stays monotone.
        uniq, idx = np.unique(cdf, return_index=True)
        return np.interp(rng.random(n), uniq, grid[idx])


def parse_extinction_distribution(spec: object) -> ExtinctionDistribution:
    """Normalize the ``extinctionDistribution`` mapping into a sampler.

    Distinguished from a screen by the presence of a ``distribution`` key (the
    schema puts screens and this sampler under one ``extinction`` slot for
    generative targets); the loader/cluster dispatches on that key.
    """
    if not isinstance(spec, Mapping) or "distribution" not in spec:
        raise ExtinctionError(
            "E14", f"extinction distribution needs a 'distribution' key: {spec!r}"
        )
    params = spec.get("params")
    if not isinstance(params, Mapping) or not params:
        raise ExtinctionError("E14", "distribution needs a non-empty 'params'")
    law = ExtinctionLaw.coerce(spec.get("law"))
    rv = float(spec.get("rv", 3.1))
    return ExtinctionDistribution(
        kind=str(spec["distribution"]), params=dict(params), law=law, rv=rv
    )
