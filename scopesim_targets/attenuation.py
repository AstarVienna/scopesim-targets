# -*- coding: utf-8 -*-
"""Source-side dust attenuation (see ``docs/defining_extinction.md``).

``attenuation`` is dust at the source, e.g. the cloud of an embedded cluster,
the circumstellar envelope of a YSO, following the definition of the
``dust_extinction`` package: stars mixed with or shelled by dust, where the
dust column differs from star to star and (for a shell) light scattered
back _into_ the line of sight matters. It is a property of the *component*, not
of the sightline, so it never cascades from a scene.

The grammar has two parts, following column-density literature conventions:

* a dust **column**, in extinction "currency": ``value`` + ``band`` (e.g. A_K),
  ``ebv``, or ``tau_v``. For generative targets, those values can be drawn per
  star from a ``distribution``.
* a **model** turning that column into the attenuation curve of the emergent
  light: ``screen`` (default; per-star differential extinction with a
  ``dust_extinction`` law) or ``wg00`` (Witt & Gordon 2000 radiative-transfer
  curves from ``dust_attenuation``; geometry-aware, scattering included).

Because it varies per star, it cannot be folded into the deduplicated spectra
of a point-source table. It is therefore _not_ resolved here into the SED, but
*emitted*: a per-row ``attenuation`` column in the model's native parameter
(``Av`` [mag] for ``screen``, ``tau_V`` for ``wg00``) plus flat ``table.meta``
keys (see :func:`apply_attenuation_to_table`). The simulation engine applies the
curve per row; :meth:`Attenuation.from_table_meta` + :meth:`Attenuation.curve`
are the reference implementation it can import, and the same curve is what a
target uses to set an ``observed``-anchor weight against the attenuated SED.

Normative order: rest SED -> attenuation (source frame) -> [distance modulus]
-> redshift -> extinction (observer frame) -> anchor.
"""

from typing import Literal
from dataclasses import dataclass
from enum import Enum
from collections.abc import Mapping

import numpy as np
from numpy.typing import NDArray
import astropy.units as u
from synphot import SourceSpectrum, SpectralElement
from synphot.models import Empirical1D

from .extinction import (
    ExtinctionError,
    EmptyExtError,
    UnphysExtError,
    ExtinctionLaw,
    PLAIN_MAG,
    as_plain_mag,
)

__all__ = [
    "AttenuationError",
    "AttenuationModel",
    "DustColumn",
    "ColumnDistribution",
    "Attenuation",
    "parse_attenuation",
    "apply_attenuation_to_table",
    "ATTENUATION_COLUMN",
]

# Name of the per-row source-table column carrying the native amount.
ATTENUATION_COLUMN = "attenuation"


class AttenuationError(ExtinctionError):
    """Structured attenuation error; extends the E13-E18 table with E19-E22."""


class EmptyAttenError(AttenuationError, EmptyExtError):
    pass


class UnphysAttenError(AttenuationError, UnphysExtError):
    pass


class AttenuationModel(Enum):
    """How a dust column becomes an attenuation curve (schema ``attenuationModel``)."""

    SCREEN = "screen"
    WG00 = "wg00"

    @classmethod
    def coerce(cls, value: "AttenuationModel | str | None") -> "AttenuationModel":
        if value is None:
            return cls.SCREEN
        if isinstance(value, cls):
            return value
        try:
            return cls(str(value).lower())
        except ValueError as exc:
            allowed = ", ".join(repr(m.value) for m in cls)
            raise AttenuationError(
                f"attenuation model must be one of {allowed}, not {value!r}"
            ) from exc

    @property
    def native_parameter(self) -> str:
        """Name of the amount the model is parameterized by (``Av`` / ``tau_V``)."""
        return "Av" if self is AttenuationModel.SCREEN else "tau_V"

    @property
    def native_unit(self):
        return u.mag if self is AttenuationModel.SCREEN else u.dimensionless_unscaled


@dataclass(frozen=True)
class DustColumn:
    """One dust column in extinction currency.

    ``kind`` is ``"mag"`` (``value`` is A in ``band``), ``"ebv"`` (E(B-V)) or
    ``"tau_v"`` (V-band optical depth). Conversion to a model's native
    parameter happens in :meth:`Attenuation.native_amount`.
    """

    kind: str
    value: float
    band: str = "V"

    def __post_init__(self):
        if self.kind not in ("mag", "ebv", "tau_v"):
            raise EmptyAttenError(f"unknown column kind {self.kind!r}")
        if self.value < 0:
            raise AttenuationError(f"column must be >= 0, got {self.value}")


@dataclass(frozen=True)
class ColumnDistribution:
    """A per-star column PDF for generative targets (e.g. clusters).

    ``kind`` is one of ``"column_density_pdf"`` / ``"lognormal"`` / ``"uniform"``.
    ``params`` are read per kind; column-valued ones are plain numbers or
    ``"N mag"`` strings in ``currency`` (``"mag"``, measured in ``band``; or
    ``"tau_v"``). Sampling is done by :meth:`sample`, seeded by the caller so the
    whole realization is reproducible under the target's single ``rng_seed``.

    Parameters by kind
    ------------------
    * ``uniform``: ``low``, ``high``.
    * ``lognormal``: ``median``, ``sigma`` (width in ln column).
    * ``column_density_pdf``: a lognormal core (``median``, ``sigma``) handing
      over to a power-law tail ``p(N) ∝ N**-slope`` above ``break``, continuous
      at the break (Alves, Lombardi & Lada 2017, A&A 606, L2, who quote it in
      A_K -> ``band: Ks``). The cloud's gas-to-dust conversion is folded into
      the choice of these parameters. Optional ``min``/``max`` bound the support
      (defaults ``median/100`` and ``20*break``); sampling is grid inverse-CDF.

    Zeroth-order model: columns are drawn independently of stellar position
    (deeper embedding toward the centre is a planned extension sharing the RNG).
    """

    kind: Literal["column_density_pdf", "lognormal", "uniform"]
    params: Mapping[str, object]
    currency: Literal["mag", "tau_v"] = "mag"
    band: str = "V"

    # Consider uning StrEnum instead
    _KINDS = ("column_density_pdf", "lognormal", "uniform")
    _CURRENCIES = ("mag", "tau_v")

    def __post_init__(self):
        if self.kind not in self._KINDS:
            raise AttenuationError(
                f"distribution must be one of {self._KINDS}, not {self.kind!r}"
            )
        if self.currency not in self._CURRENCIES:
            raise EmptyAttenError(
                f"currency must be one of {self._CURRENCIES}, "
                f"not {self.currency!r}"
            )

    def _num(self, key: str, default=None) -> float:
        if key not in self.params:
            if default is None:
                raise EmptyAttenError(f"distribution needs param {key!r}")
            return float(default)
        return as_plain_mag(self.params[key], key)

    def sample(self, n: int, rng: np.random.Generator) -> NDArray:
        """Draw ``n`` non-negative columns (in ``currency``) from ``rng``."""
        if n < 0:
            raise ValueError("n must be >= 0")
        if self.kind == "uniform":
            low, high = self._num("low"), self._num("high")
            if not 0 <= low <= high:
                raise UnphysAttenError(f"need 0 <= low <= high, got {low}, {high}")
            return rng.uniform(low, high, size=n)

        if self.kind == "lognormal":
            median, sigma = self._num("median"), self._num("sigma")
            if median <= 0 or sigma <= 0:
                raise UnphysAttenError("median and sigma must be > 0")
            return rng.lognormal(mean=np.log(median), sigma=sigma, size=n)

        return self._sample_column_density_pdf(n, rng)

    def _sample_column_density_pdf(self, n: int, rng):
        median, sigma = self._num("median"), self._num("sigma")
        brk, slope = self._num("break"), self._num("slope")
        if median <= 0 or sigma <= 0 or brk <= 0:
            raise UnphysAttenError("median, sigma, break must be > 0")
        if slope <= 1:
            raise UnphysAttenError(
                f"tail slope must be > 1 (integrable), got {slope}"
            )
        lo = self._num("min", median / 100.0)
        hi = self._num("max", 20.0 * brk)
        if not 0 < lo < hi:
            raise UnphysAttenError("need 0 < min < max")

        grid = np.linspace(lo, hi, 8192)

        def lognormal_pdf(a):
            return (1.0 / (a * sigma * np.sqrt(2 * np.pi))) * np.exp(
                -((np.log(a) - np.log(median)) ** 2) / (2 * sigma**2)
            )

        core = lognormal_pdf(grid)
        # power-law tail beyond the break, scaled for continuity at the break
        tail = lognormal_pdf(brk) * (grid / brk) ** (-slope)
        pdf = np.where(grid <= brk, core, tail)

        cdf = np.cumsum(pdf)
        cdf /= cdf[-1]
        uniq, idx = np.unique(cdf, return_index=True)  # keep interp monotone
        return np.interp(rng.random(n), uniq, grid[idx])


@dataclass(frozen=True)
class Attenuation:
    """Normalized ``attenuation`` attribute: a column (or column PDF) + a model.

    ``law``/``rv`` belong to ``screen``; ``geometry``/``dust_type``/
    ``dust_distribution`` to ``wg00`` (:func:`parse_attenuation` enforces the
    pairing, E19). Frozen and hashable in the fixed-column case.
    """

    model: AttenuationModel
    column: DustColumn | ColumnDistribution
    law: ExtinctionLaw = ExtinctionLaw.CCM89
    rv: float = 3.1
    geometry: Literal["shell", "cloudy", "dusty"] | None = None
    dust_type: Literal["mw", "smc"] = "mw"
    dust_distribution: Literal["homogeneous", "clumpy"] | None = None

    # Consider uning StrEnum instead
    _GEOMETRIES = ("shell", "cloudy", "dusty")
    _DUST_TYPES = ("mw", "smc")
    _DUST_DISTRIBUTIONS = ("homogeneous", "clumpy")

    # WG00 validity (dust_attenuation.radiative_transfer.WG00): tau_V grid 0.25-50.
    _WG00_TAU_MIN =  0.25
    _WG00_TAU_MAX = 50.0

    # Fallback effective wavelengths (Bessell 2005 Johnson-Cousins) used only when
    # no band-wavelength lookup is supplied; V/B are needed to convert A_V / E(B-V)
    # through the wg00 dust curve. Any other band requires the lookup.
    _FALLBACK_BAND_WAVELENGTHS = {"V": 0.551 * u.um, "B": 0.445 * u.um}

    def __post_init__(self):
        if self.model is AttenuationModel.SCREEN:
            if self.rv <= 0:
                raise UnphysAttenError(f"R_V must be > 0, got {self.rv}")
        else:
            if self.geometry not in self._GEOMETRIES:
                raise AttenuationError(
                    f"wg00 needs geometry in {self._GEOMETRIES}, "
                    f"got {self.geometry!r}"
                )
            if self.dust_distribution not in self._DUST_DISTRIBUTIONS:
                raise AttenuationError(
                    f"wg00 needs dust_distribution in {self._DUST_DISTRIBUTIONS}, "
                    f"got {self.dust_distribution!r}",
                )
            if self.dust_type not in self._DUST_TYPES:
                raise AttenuationError(
                    f"wg00 dust_type must be in {self._DUST_TYPES}, "
                    f"got {self.dust_type!r}"
                )

    @property
    def is_sampled(self) -> bool:
        return isinstance(self.column, ColumnDistribution)

    @property
    def native_parameter(self) -> str:
        return self.model.native_parameter

    def _wg00(self, tau_v: float = 1.0):
        # Lazy import: dust_attenuation is an optional, GitHub-only dependency.
        try:
            from dust_attenuation.radiative_transfer import WG00
        except ModuleNotFoundError as exc:
            raise ModuleNotFoundError(
                "attenuation model 'wg00' needs the optional dust_attenuation "
                "package, which is not on PyPI yet: pip install "
                "git+https://github.com/karllark/dust_attenuation.git"
            ) from exc

        return WG00(
            tau_V=tau_v,
            geometry=self.geometry,
            dust_type=self.dust_type,
            dust_distribution=self.dust_distribution,
        )

    def _band_wavelength(self, band: str, lookup) -> u.Quantity:
        if lookup is not None:
            return u.Quantity(lookup(band))
        if band in self._FALLBACK_BAND_WAVELENGTHS:
            return self._FALLBACK_BAND_WAVELENGTHS[band]
        raise EmptyAttenError(
            f"A_{band} needs a band_wavelength_lookup to convert the column"
        )

    def _extinction_per_unit_column(self, band: str, lookup) -> float:
        """A_band per unit native parameter: A/A_V (screen) or A per tau_V (wg00)."""
        if self.model is AttenuationModel.SCREEN:
            if band == "V":
                return 1.0  # laws are V-normalized: no bandpass needed
            lam = self._band_wavelength(band, lookup)
            return float(np.squeeze(self.law.model(self.rv)(lam)))
        # wg00: the model's own dust extinction curve, evaluated for tau_V = 1
        # (exactly linear in tau_V). NB its tabulation gives ~1.109 mag at
        # 0.551 um rather than the nominal 1.086 = 2.5 log10(e); we use the
        # model's curve consistently for every band rather than special-case V.
        lam = self._band_wavelength(band, lookup)
        return float(np.squeeze(self._wg00().get_extinction(lam, 1.0)))

    def _to_native(self, kind: str, values, band: str, lookup):
        vals = np.asarray(values, dtype=float)
        if self.model is AttenuationModel.SCREEN:
            if kind == "mag":
                native = vals / self._extinction_per_unit_column(band, lookup)
            elif kind == "ebv":
                native = self.rv * vals
            else:  # tau_v -> A_V for a pure screen
                native = (2.5 / np.log(10.0)) * vals
        else:
            if kind == "tau_v":
                native = vals
            elif kind == "mag":
                native = vals / self._extinction_per_unit_column(band, lookup)
            else:  # E(B-V) = A_B - A_V, both linear in tau_V
                per_tau = (
                    self._extinction_per_unit_column("B", lookup)
                    - self._extinction_per_unit_column("V", lookup)
                )
                native = vals / per_tau
        self._check_range(native)
        return native

    def _check_range(self, native) -> None:
        arr = np.atleast_1d(native)
        if np.any(arr < 0):
            raise AttenuationError("negative column")
        if self.model is AttenuationModel.WG00 and np.any(arr > self._WG00_TAU_MAX):
            raise AttenuationError(
                f"wg00 tau_V max is {self._WG00_TAU_MAX}, "
                f"got up to {arr.max():.3g}"
            )

    def native_amount(
        self,
        *,
        n: int | None = None,
        rng=None,
        band_wavelength_lookup=None,
    ):
        """The model's native parameter (``Av`` [mag] or ``tau_V``) for the column.

        A fixed column returns a float. A sampled column needs ``n`` and ``rng``
        (E21 otherwise) and returns an array of ``n`` draws converted to the
        native parameter. ``band_wavelength_lookup`` (``band -> Quantity``) is
        needed for non-V bands (screen) or any band without a built-in fallback.
        """
        col = self.column
        if isinstance(col, DustColumn):
            return float(self._to_native(col.kind, col.value, col.band, band_wavelength_lookup))
        if n is None or rng is None:
            raise AttenuationError("sampling needs n and a seeded rng")
        kind = "mag" if col.currency == "mag" else "tau_v"
        return self._to_native(kind, col.sample(n, rng), col.band, band_wavelength_lookup)

    def curve(self, wave: u.Quantity, amount: float) -> "np.ndarray":
        """Attenuation ``A_att(lambda)`` [mag] for one native ``amount``.

        Wavelengths outside the model's validity range are clamped to the
        nearest in-range value (curve held constant beyond validity), the same
        rule :mod:`.extinction` uses. For ``wg00`` below its tabulated
        ``tau_V = 0.25`` the curve is scaled linearly from the 0.25 curve
        (optically-thin limit: attenuation ∝ column).
        """
        wave = u.Quantity(wave).to(u.um)
        if amount == 0:
            return np.zeros(wave.shape)
        if self.model is AttenuationModel.SCREEN:
            model = self.law.model(self.rv)
            lo, hi = model.x_range  # in 1/um
            x = np.clip(1.0 / wave.value, lo, hi)
            return amount * np.asarray(model((1.0 / x) * u.um), dtype=float)

        scale = 1.0
        tau = amount
        if tau < self._WG00_TAU_MIN:
            scale, tau = amount / self._WG00_TAU_MIN, self._WG00_TAU_MIN
        model = self._wg00(tau)
        lo, hi = model.x_range  # in um
        lam = np.clip(wave.value, lo, hi) * u.um
        return scale * np.asarray(model(lam), dtype=float)

    def transmission_element(self, wave: u.Quantity, amount: float) -> SpectralElement:
        """``10**(-0.4 A_att)`` as a ``synphot`` throughput on ``wave``."""
        trans = 10 ** (-0.4 * self.curve(wave, amount))
        return SpectralElement(Empirical1D, points=wave, lookup_table=trans)

    def attenuate(self, spectrum: SourceSpectrum, amount: float) -> SourceSpectrum:
        """Return ``spectrum`` attenuated by the curve for ``amount``."""
        if amount == 0:
            return spectrum
        if spectrum.waveset is None:
            raise ValueError(
                "spectrum has no waveset (analytic model); attenuate an "
                "empirical/resampled spectrum, or use transmission_element on a grid"
            )
        return spectrum * self.transmission_element(spectrum.waveset, amount)

    def table_meta(self) -> dict:
        """Flat ``table.meta`` keys describing the model (FITS-header friendly)."""
        meta = {
            "attenuation_model": self.model.value,
            "attenuation_parameter": self.native_parameter,
        }
        if self.model is AttenuationModel.SCREEN:
            meta["attenuation_law"] = self.law.value
            meta["attenuation_rv"] = self.rv
        else:
            meta["attenuation_geometry"] = self.geometry
            meta["attenuation_dust_type"] = self.dust_type
            meta["attenuation_dust_distribution"] = self.dust_distribution
        return meta

    @classmethod
    def from_table_meta(cls, meta: Mapping) -> "Attenuation | None":
        """Rebuild the model from ``table.meta`` (the engine-side entry point).

        The column is not part of the meta (it is per row), so the returned
        object carries a zero placeholder column; use :meth:`curve` /
        :meth:`transmission_element` with each row's amount.
        """
        model_name = meta.get("attenuation_model")
        if model_name is None:
            return None
        model = AttenuationModel.coerce(model_name)
        placeholder = DustColumn("tau_v" if model is AttenuationModel.WG00 else "mag", 0.0)
        if model is AttenuationModel.SCREEN:
            return cls(
                model,
                placeholder,
                law=ExtinctionLaw.coerce(meta.get("attenuation_law")),
                rv=float(meta.get("attenuation_rv", 3.1)),
            )
        return cls(
            model,
            placeholder,
            geometry=meta.get("attenuation_geometry"),
            dust_type=meta.get("attenuation_dust_type", "mw"),
            dust_distribution=meta.get("attenuation_dust_distribution"),
        )


# --------------------------------------------------------------------------- #
# Parsing                                                                     #
# --------------------------------------------------------------------------- #

_SCREEN_KEYS = {"law", "rv"}
_WG00_KEYS = {"geometry", "dust_type", "dust_distribution"}
_CURRENCY_KEYS = ("value", "ebv", "tau_v", "distribution")
_ALL_KEYS = {"model", "band", "params", "currency", *_CURRENCY_KEYS, *_SCREEN_KEYS, *_WG00_KEYS}


def parse_attenuation(spec: object) -> Attenuation:
    """Normalize the ``attenuation`` input (schema ``attenuation``).

    * bare ``"8.5 mag"`` / number -> ``screen``, A_V, ccm89, R_V 3.1
    * mapping -> ``model`` (default ``screen``) + exactly one currency among
      ``value`` (+``band``) / ``ebv`` / ``tau_v`` / ``distribution`` (+``params``,
      ``band``, ``currency``), plus the model's own keys only (E19).
    """
    if isinstance(spec, str):
        if not PLAIN_MAG.match(spec):
            raise EmptyAttenError(
                f"bare attenuation sugar must be a plain 'N mag' A_V, got {spec!r}"
            )
        return Attenuation(AttenuationModel.SCREEN, DustColumn("mag", float(spec.split()[0])))
    if isinstance(spec, (int, float)) and not isinstance(spec, bool):
        return Attenuation(AttenuationModel.SCREEN, DustColumn("mag", float(spec)))
    if not isinstance(spec, Mapping):
        raise EmptyAttenError(f"attenuation not understood: {spec!r}")

    unknown = set(spec) - _ALL_KEYS
    if unknown:
        raise EmptyAttenError(f"unknown attenuation key(s) {sorted(unknown)}")

    model = AttenuationModel.coerce(spec.get("model"))
    currencies = [k for k in _CURRENCY_KEYS if k in spec]
    if not currencies:
        raise EmptyAttenError("attenuation needs value, ebv, tau_v or distribution")
    if len(currencies) > 1:
        raise EmptyAttenError(f"attenuation sets multiple of {currencies}")
    currency = currencies[0]

    if currency != "distribution" and ("params" in spec or "currency" in spec):
        raise AttenuationError("'params'/'currency' only pair with 'distribution'")
    if "band" in spec and currency in ("ebv", "tau_v"):
        raise AttenuationError(f"'band' does not pair with {currency!r}")

    # model-specific keys
    if model is AttenuationModel.SCREEN and (set(spec) & _WG00_KEYS):
        raise AttenuationError(f"{sorted(set(spec) & _WG00_KEYS)} are wg00-only")
    if model is AttenuationModel.WG00 and (set(spec) & _SCREEN_KEYS):
        raise AttenuationError(f"{sorted(set(spec) & _SCREEN_KEYS)} are screen-only")

    band = str(spec.get("band", "V"))
    if currency == "value":
        val = spec["value"]
        if not (isinstance(val, str) and PLAIN_MAG.match(val)):
            raise EmptyAttenError(
                f"attenuation 'value' must be a plain 'N mag' amount, got {val!r}"
            )
        column = DustColumn("mag", float(val.split()[0]), band=band)
    elif currency == "ebv":
        column = DustColumn("ebv", float(spec["ebv"]))
    elif currency == "tau_v":
        column = DustColumn("tau_v", float(spec["tau_v"]))
    else:
        params = spec.get("params")
        if not isinstance(params, Mapping) or not params:
            raise EmptyAttenError("distribution needs a non-empty 'params'")
        column = ColumnDistribution(
            kind=str(spec["distribution"]),
            params=dict(params),
            currency=str(spec.get("currency", "mag")),
            band=band,
        )

    if model is AttenuationModel.SCREEN:
        return Attenuation(
            model, column,
            law=ExtinctionLaw.coerce(spec.get("law")),
            rv=float(spec.get("rv", 3.1)),
        )
    return Attenuation(
        model, column,
        geometry=spec.get("geometry"),
        dust_type=str(spec.get("dust_type", "mw")),
        dust_distribution=spec.get("dust_distribution"),
    )


# --------------------------------------------------------------------------- #
# Source-table emission (-> simulation engine)                                #
# --------------------------------------------------------------------------- #


def apply_attenuation_to_table(
    tbl,
    attenuation: "Attenuation | None",
    *,
    rng=None,
    band_wavelength_lookup=None,
):
    """Add the per-row ``attenuation`` column and model meta to a source table.

    Contract with the simulation engine:

    * column :data:`ATTENUATION_COLUMN` (float, one value per row) in the
      model's native parameter -- ``Av`` [mag] for ``screen`` (column unit
      ``mag``), ``tau_V`` for ``wg00`` (dimensionless); ``meta["attenuation_parameter"]``
      names it;
    * flat meta keys from :meth:`Attenuation.table_meta` (``attenuation_model``
      and the model's own parameters);
    * absent column = no attenuation; the table is returned untouched
      (byte-identical);
    * the engine multiplies each row's spectrum by
      ``Attenuation.from_table_meta(meta).transmission_element(wave, row_amount)``
      (or its own equivalent), at observed wavelengths in v1.

    A sampled column draws ``len(tbl)`` values from ``rng`` (E21 if missing); a
    fixed column broadcasts (the N=1, zero-variance case). Returns ``tbl``.
    """
    if attenuation is None:
        return tbl
    if attenuation.is_sampled:
        if rng is None:
            raise AttenuationError(
                "sampled attenuation ('distribution') is for generative "
                "targets only and needs an rng_seed"
            )
        amount = attenuation.native_amount(
            n=len(tbl), rng=rng, band_wavelength_lookup=band_wavelength_lookup
        )
    else:
        amount = attenuation.native_amount(band_wavelength_lookup=band_wavelength_lookup)
    tbl[ATTENUATION_COLUMN] = amount
    tbl[ATTENUATION_COLUMN].unit = attenuation.model.native_unit
    tbl.meta.update(attenuation.table_meta())
    return tbl
