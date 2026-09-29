# -*- coding: utf-8 -*-
"""Any target that will be unresolved for all instruments, e.g. ``Star``."""

import warnings
from collections.abc import Sequence, Mapping
from numbers import Integral

import numpy as np
from numpy.typing import ArrayLike
from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
from synphot import SourceSpectrum, units as synphot_units
from synphot.models import ConstFlux1D

from astar_utils import SpectralType
from astar_utils.guard_functions import guard_same_len
from spextra import Spextrum
from scopesim import Source
from scopesim.source.source_fields import TableSourceField

from .typing_utils import POSITION_TYPE, SPECTRUM_TYPE, BRIGHTNESS_TYPE
from .target import Brightness, SpectrumTarget
from .seeding import resolve_seed, stream, new_rng
from .stellar.luminosity import sample_magnitudes
from .stellar.morphology import UniformRectangularMorphology
from .brightness import (
    is_brightness_spec,
    is_single_brightness,
    FromSpectralType,
    LocatorError,
)


class PointSourceTarget(SpectrumTarget):
    """Base class for Point Source Targets."""

    def __init__(
        self,
        position: POSITION_TYPE | None = None,
        spectrum: SPECTRUM_TYPE | None = None,
        brightness: BRIGHTNESS_TYPE | None = None,
        anchor: str | None = None,
    ) -> None:
        if position is not None:
            self.position = position
        if spectrum is not None:
            self.spectrum = spectrum
        if brightness is not None:
            self.brightness = brightness
        if anchor is not None:
            self.anchor = anchor

    def to_source(self, optical_train=None) -> Source:
        """Convert to ScopeSim Source object."""
        source = Source(
            field=TableSourceField(
                self.to_table(), spectra=self.source_spectra()
            )
        )
        return source

    def to_table(self, local_frame=None) -> Table:
        """Convert to table for Source conversion."""
        tbl = self._create_source_table()
        tbl.add_row(self._to_table_row(local_frame))
        return tbl

    def source_spectra(self, start: int = 0) -> dict[int, SourceSpectrum]:
        """Create spectra dict for Source conversion."""
        # TODO: Deal with redshift from position!
        return {start: self.resolve_spectrum(self.spectrum)}

    def _create_source_table(self) -> Table:
        tbl = Table(
            names=["x", "y", "ref", "weight"],
            units={"x": u.arcsec, "y": u.arcsec},
        )

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
        x_arcsec, y_arcsec = self._xy_arcsec_position(
            self.position, local_frame
        )

        # If not given from parent, resolve now
        if spectrum is None:
            spectrum = self.resolve_spectrum(self.spectrum)
            spectrum = self.redshift_spectrum(spectrum, self.position)

        weight = self._anchored_spectrum_scale(spectrum, self.brightness)

        row = {
            "x": x_arcsec,
            "y": y_arcsec,
            "weight": weight,
            "ref": ref,
        }
        return row


class Star(PointSourceTarget):
    """A single star.

    Examples
    --------
    >>> tgt = Star(
    ...     position=(2, 3),  # [arcsec] from center of FOV
    ...     spectrum="A0V",
    ...     brightness=("R", 15),  # [mag]
    ... )

    For more examples, see also
    `the YAML syntax <../yaml_syntax.html#single-stars>`_.

    """

    @classmethod
    def from_spectral_type(
        cls,
        spectrum: SPECTRUM_TYPE,
        *,
        position: POSITION_TYPE | None = None,
        band: str = "V",
        table: str = "mamajek",
    ) -> "Star":
        """Build a star whose brightness is resolved from its spectral type.

        Convenience wrapper around the ``brightness: {from_spectral_type: ...}``
        resolver: the absolute magnitude for ``spectrum`` is looked up in
        ``table`` (band ``M_<band>``) and ``anchor`` is set to ``absolute``, so
        the standard machinery applies the distance modulus from
        ``position.distance`` (required -- a missing distance is E10) and any
        extinction screens. The resolved value and table name/version are
        recorded in :attr:`~.target.SpectrumTarget.brightness_provenance`.

        Examples
        --------
        >>> tgt = Star.from_spectral_type(
        ...     "K5V", position={"distance": 25 * u.pc}
        ... )
        """
        return cls(
            position=position,
            spectrum=spectrum,
            brightness={"from_spectral_type": table, "band": band},
        )


class Binary(PointSourceTarget):
    """Binary star.

    .. todo:: Fix offset definition via distance and physical separation, see
        :issue:`107` (also applies to other target subclasses).

    Examples
    --------
    >>> tgt = Binary(
    ...     position={"distance": 100*u.pc},
    ...     spectra=("F0V", "M2V"),
    ...     brightness=("R", 10),
    ...     contrast=100.,  # F_pri / F_sec
    ...     offset={"separation": .1*u.AU},
    ... )

    For more examples, see also
    `the YAML syntax <../yaml_syntax.html#binaries>`_.

    """

    def __init__(
        self,
        position: POSITION_TYPE | None = None,
        offset: Mapping[str, float | u.Quantity] | None = None,
        spectra: Sequence[SPECTRUM_TYPE] | None = None,
        brightness: BRIGHTNESS_TYPE | Sequence[BRIGHTNESS_TYPE] | None = None,
        contrast: float | None = None,
        anchor: str | None = None,
    ) -> None:
        if position is not None:
            self.position = position
        if offset is not None:
            self.offset = offset
        if anchor is not None:
            self.anchor = anchor

        if spectra is not None:
            self.primary_spectrum, self.secondary_spectrum = spectra

        # Distinguish "two brightness specs" from "one (locator, amount) spec"
        # structurally (see _is_brightness_pair), so string amounts
        # ('R', '15 mag'), the canonical mapping form, and wavelength/frequency
        # locators all route correctly -- the old (str(), Quantity()|Number())
        # shape silently misread every one of those as a pair.
        is_pair = self._is_brightness_pair(brightness)
        match brightness, contrast:
            case None, None:
                pass  # TODO: What to do here?
            case None, _:
                self.contrast = contrast
            case _, None if is_pair:
                primary, secondary = brightness
                self.brightness = primary
                self.brightness_secondary = secondary
            case _, None:
                self.brightness = brightness
            case _, _ if is_pair:
                raise TypeError(
                    "Either supply brightness and contrast or two brightness "
                    "specs, but not both."
                )
            case _, _:
                self.brightness = brightness
                self.contrast = contrast

    @staticmethod
    def _is_brightness_pair(brightness) -> bool:
        """Is ``brightness`` two specs (primary + secondary) or one spec?

        A *pair* is a length-2 sequence whose **both** elements are themselves
        brightness specs -- each a mapping or a ``(locator, amount)`` sequence.
        A *single* spec is a mapping, or a ``(locator, amount)`` pair whose
        elements are a bare locator and amount (``str`` / ``Quantity`` /
        number), never nested specs. The test is purely structural, so it does
        not care whether the amount was written as a number, a ``Quantity`` or
        a string, nor whether the locator is a band or a wavelength/frequency.
        """
        return (
            isinstance(brightness, Sequence)
            and not isinstance(brightness, (str, bytes))
            and len(brightness) == 2
            and all(is_brightness_spec(element) for element in brightness)
        )

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

        .. todo:: Add support for dimensionless Quantity and other numerical
            types in setter type check.
        """
        return self._contrast

    @contrast.setter
    def contrast(self, contrast: float):
        if not isinstance(contrast, float):
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
            return self._anchored_spectrum_scale(
                secondary_spectrum, self.brightness_secondary
            )
        raise ValueError("Either contrast or secondary brightness is needed.")

    def to_table(self, local_frame=None, spectra=None, refs=None) -> Table:
        """Convert to table for Source conversion."""
        tbl = self._create_source_table()

        # TODO: Add support for parent frame as in base class
        try:
            primary_position = self.position
        except AttributeError:
            # Default to (0, 0)
            primary_position = SkyCoord(0 * u.deg, 0 * u.deg, 10 * u.pc)

        local_frame = primary_position.skyoffset_frame()
        x_arcsec_pri, y_arcsec_pri = self._xy_arcsec_position(
            primary_position, local_frame
        )
        x_arcsec_sec, y_arcsec_sec = self._xy_arcsec_position(
            self.resolve_position(primary_position),
            local_frame,
        )

        spectra, (ref_pri, ref_sec) = self._resolve_spectra_refs(spectra, refs)
        # Use primary position for both to deal with systemic radial velocity
        # TODO: Implement instantaneous radial velocities somehow...
        spectra = {
            ref: self.redshift_spectrum(spec, primary_position)
            for ref, spec in spectra.items()
        }

        primary = {
            "x": x_arcsec_pri,
            "y": y_arcsec_pri,
            "weight": self._anchored_spectrum_scale(
                spectra[ref_pri], self.brightness
            ),
            "ref": ref_pri,
        }

        secondary = {
            "x": x_arcsec_sec,
            "y": y_arcsec_sec,
            "weight": self._resolve_secondary_weight(
                spectra[ref_sec], primary["weight"]
            ),
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
    """Exoplanet (point source) with default spectrum of Neptune.

    Examples
    --------
    >>> tgt = Exoplanet(
    ...     position=(0, 0),
    ...     brightness=("V", 20),
    ... )

    """

    def __init__(
        self,
        position: POSITION_TYPE | None = None,
        offset: Mapping[str, float | u.Quantity] | None = None,
        spectrum: SPECTRUM_TYPE | None = None,
        brightness: BRIGHTNESS_TYPE | None = None,
        contrast: float | None = None,
        anchor: str | None = None,
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
        if anchor is not None:
            self.anchor = anchor

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
    """Planetary system with primary and components.

    Examples
    --------
    >>> tgt = PlanetarySystem(
    ...     position=(0, 0),
    ...     primary=Star(
    ...         spectrum="A0V",
    ...         brightness=("R", 15),
    ...     ),
    ...     components=[
    ...         Exoplanet(
    ...             contrast=1e5,
    ...             offset={"separation": 0.5*u.arcsec},
    ...         ),
    ...     ],
    ... )

    For more examples, see also
    `the YAML syntax <../yaml_syntax.html#exoplanetary>`_.

    """

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

    def to_source(self, optical_train=None) -> Source:
        """Convert to ScopeSim Source object."""
        local_frame = self.position.skyoffset_frame()

        # HACK: Should be able to pass this down
        self.primary.position = self.position

        table = self.primary.to_table(local_frame)
        spectra = self.primary.source_spectra()

        for ref, component in enumerate(
            self.components, start=max(spectra) + 1
        ):
            spectrum = component.resolve_spectrum(component.spectrum)
            spectrum = self.redshift_spectrum(spectrum, self.position)

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


class EmptyStarFieldWarning(UserWarning):
    """A randomly generated star field came out with no stars.

    Not an error: a thin field can legitimately draw zero stars. Escalate with
    ``warnings.simplefilter("error", EmptyStarFieldWarning)`` to make it one.
    """


# TODO: Common base class for multi-component targets
class StarField(PointSourceTarget):
    """Multiple stars around a common field center.

    Per-star ``positions`` are offsets from the field center ``position``
    (default: (0, 0)), in arcsec, with +x = East and +y = North -- the same
    convention as a component ``offset`` relative to its parent. They are
    stored vectorized, as one ``(N, 2)`` offset array (:attr:`offsets`), so
    fields of many thousands of stars stay cheap. Moving the field center moves
    the whole field rigidly.

    ``spectra`` is either one spectrum per star, or a single spectrum shared by
    all stars.

    A field may be empty (zero stars); it still converts to a valid (empty)
    ScopeSim Source.

    Fields made by the constructors :meth:`from_grid`, :meth:`from_random` and
    :meth:`from_random_density` record how they were made in
    :attr:`generated_by` (method, parameters and, for random fields, the
    resolved master seed). This is provenance only: the realized stars are
    authoritative, nothing is regenerated from it.

    .. todo:: Default ``role`` to ``"background"`` (fore-/background fields are
        the typical use), once ``role`` is implemented on targets and in the
        engine; see the schema ``role`` definition.

    Examples
    --------
    >>> tgt = StarField(
    ...     positions=[(0, 0), (1, 1)],  # [arcsec] offset from field center
    ...     spectra=["A0V", "G2V"],
    ...     brightnesses=[10, 15],  # [mag]
    ...     band="V",  # default for bare brightness amounts
    ... )

    A single spectrum for all stars, around an explicit field center:

    >>> tgt = StarField(
    ...     position=SkyCoord(150.1 * u.deg, 2.2 * u.deg),
    ...     positions=np.array([(-5, 0), (0, 0), (5, 0)]),
    ...     spectra="A0V",
    ...     brightnesses=["15 mag(AB)", "16 mag(AB)", "17 mag(AB)"],
    ...     band="V",
    ... )

    For more examples, see also
    `the YAML syntax <../yaml_syntax.html#star-field>`_.

    """

    def __init__(
        self,
        positions: Sequence | np.ndarray | u.Quantity | SkyCoord | None = None,
        spectra: SPECTRUM_TYPE | Sequence[SPECTRUM_TYPE] | None = None,
        brightnesses: Sequence[BRIGHTNESS_TYPE] | None = None,
        band: str | None = None,  # TODO: Proper typing
        position: POSITION_TYPE | None = None,
    ) -> None:
        # Length consistency is checked by each setter against the attributes
        # already set, so no up-front guard (which could not know about a
        # broadcast spectrum) is needed.
        self.band = band
        self.generated_by: dict | None = None
        # Center first: absolute per-star positions are converted to offsets
        # from it on assignment.
        if position is not None:
            self.position = position
        if positions is not None:
            self.positions = positions
        if spectra is not None:
            self.spectra = spectra
        if brightnesses is not None:
            self.brightnesses = brightnesses

    @classmethod
    def from_grid(
        cls,
        shape: int | tuple[int, int],
        spacing: float | u.Quantity,
        *,
        spectra: SPECTRUM_TYPE | ArrayLike,
        brightnesses: BRIGHTNESS_TYPE | ArrayLike,
        band: str | None = None,
        position: POSITION_TYPE | None = None,
    ) -> "StarField":
        """Regularly spaced grid of stars, centered on the field center.

        Artificial, but useful for testing (PSF, photometric linearity,
        distortion, ...).

        Parameters
        ----------
        shape : int or (int, int)
            Number of stars as ``(width, height)``, i.e. ``(nx, ny)``. A single
            int gives a square grid.
        spacing : float or Quantity
            Distance between neighboring stars, the same along x and y. A plain
            number is in arcsec.
        spectra : spectrum or array-like of spectra
            One spectrum for all stars, or a flat sequence of ``nx * ny``.
        brightnesses : brightness or array-like of brightnesses
            One brightness for all stars, or a flat sequence of ``nx * ny``.
            Entries are full specs or bare amounts (located at `band`), as for
            :class:`StarField`. A ``(locator, amount)`` pair counts as one
            brightness (see :func:`~.brightness.is_single_brightness`).
        band : str, optional
            Default band for bare brightness amounts.
        position : optional
            Field center, see :class:`StarField`.

        Notes
        -----
        The grid is symmetric about the field center, so there is a star at
        the center only if both counts are odd.

        Flat per-star values run along x first: the first element is the star
        with the smallest x and y, i.e. bottom-left when plotted with
        ``origin="lower"``, followed by the rest of that row with increasing x,
        then the next row up. Since +x is East, "left" here is West -- the
        mirror image of a standard East-left sky display.

        Examples
        --------
        A 10x10 grid, 5 arcsec apart, with a brightness ramp from 10 to 20 mag:

        >>> tgt = StarField.from_grid(
        ...     10, 5 * u.arcsec,
        ...     spectra="A0V",
        ...     brightnesses=np.linspace(10, 20, 100),
        ...     band="V",
        ... )
        """
        n_x, n_y = _parse_grid_shape(shape)
        n_stars = n_x * n_y

        step = u.Quantity(spacing, u.arcsec)
        if not step.isscalar or not step > 0:
            raise ValueError(f"spacing must be a positive scalar, got {spacing!r}")
        step = step.to_value(u.arcsec)

        # indexing="xy" (default): x varies along axis 1, so a C-order ravel
        # runs along x first, starting with the lowest row.
        x_grid, y_grid = np.meshgrid(
            (np.arange(n_x) - (n_x - 1) / 2) * step,
            (np.arange(n_y) - (n_y - 1) / 2) * step,
        )
        offsets = np.column_stack([x_grid.ravel(), y_grid.ravel()]) << u.arcsec

        if not isinstance(spectra, (str, SpectralType, SourceSpectrum)):
            spectra = _flat_per_star(
                spectra, n_stars, "spectra", _is_single_spectrum
            )
        if is_single_brightness(brightnesses):
            brightnesses = [brightnesses] * n_stars
        else:
            brightnesses = _flat_per_star(
                brightnesses, n_stars, "brightnesses", is_single_brightness
            )

        field = cls(
            positions=offsets,
            spectra=spectra,
            brightnesses=brightnesses,
            band=band,
            position=position,
        )
        field.generated_by = {
            "method": "from_grid",
            "params": {"shape": (n_x, n_y), "spacing": step * u.arcsec},
        }
        return field

    @classmethod
    def from_random(
        cls,
        n_stars: int,
        size: float | u.Quantity | tuple,
        *,
        band: str,
        mag_range: tuple[float, float],
        slope: float = 0.3,
        system: str = "Vega",
        spectra: SPECTRUM_TYPE | ArrayLike,
        position: POSITION_TYPE | None = None,
        exclude_radius: float | u.Quantity | None = None,
        rng_seed: int | None = None,
    ) -> "StarField":
        """Exactly `n_stars` random stars, e.g. as fore- or background.

        Positions are uniform over a rectangle centered on the field center,
        magnitudes follow exponential star counts
        (:func:`~.stellar.luminosity.sample_magnitudes`), and each star's
        spectrum is drawn uniformly from `spectra`.

        Parameters
        ----------
        n_stars : int
            Number of stars (exact; 0 gives an empty field).
        size : float, Quantity or (width, height)
            Extent of the rectangle; a scalar gives a square. Plain numbers
            are arcsec; lengths need a field center `position` with a
            distance.
        band : str
            Band of the magnitudes.
        mag_range : (float, float)
            Bright and faint magnitude limits, in `system`.
        slope : float, optional
            Star-count slope ``d log10(N) / dm`` (default 0.3; 0 is uniform in
            magnitude). The realistic value depends on band and galactic
            latitude, so set it explicitly where it matters.
        system : {"Vega", "AB", "ST"}, optional
            Photometric system of the magnitudes (default Vega).
        spectra : spectrum or array-like of spectra
            One spectrum for all stars, or a flat sequence of choices each
            star draws from uniformly. Repeat entries for integer weights.
        position : optional
            Field center, see :class:`StarField`.
        exclude_radius : float or Quantity, optional
            Keep a central circle free of stars, e.g. around a science target.
            The count stays exact.
        rng_seed : int, optional
            Master seed. ``None`` (default) draws fresh entropy, which is
            recorded in :attr:`generated_by`, so every field is reproducible.

        See Also
        --------
        from_random_density : Poisson-distributed count for a given density.

        Examples
        --------
        >>> tgt = StarField.from_random(
        ...     200, (60, 40),
        ...     band="V", mag_range=(15, 22),
        ...     spectra=["G2V", "K5V", "M2V"],
        ...     exclude_radius=5 * u.arcsec,
        ...     rng_seed=42,
        ... )
        """
        master_seed = resolve_seed(rng_seed)
        field = cls(band=band, position=position)
        center = field.resolve_position()

        morphology = UniformRectangularMorphology(
            n_stars, size, exclude_radius,
            rng_seed=stream(master_seed, "morphology"),
        )
        x_arcsec, y_arcsec = morphology.sample(center)
        magnitudes = sample_magnitudes(
            n_stars, mag_range, slope, new_rng(stream(master_seed, "magnitudes"))
        )

        field.positions = np.column_stack([x_arcsec, y_arcsec]) << u.arcsec
        if isinstance(spectra, (str, SpectralType, SourceSpectrum)):
            field.spectra = spectra
        else:
            choices = list(spectra)
            if not choices or not all(_is_single_spectrum(c) for c in choices):
                raise ValueError(
                    "spectra must be a single spectrum or a non-empty flat "
                    "sequence of choices"
                )
            choices = [field._parse_spectrum(choice) for choice in choices]
            picks = new_rng(stream(master_seed, "spectra")).integers(
                len(choices), size=n_stars
            )
            field.spectra = [choices[pick] for pick in picks]
        field.brightnesses = [
            {"band": band, "value": float(magnitude), "system": system}
            for magnitude in magnitudes
        ]

        field.generated_by = {
            "method": "from_random",
            "rng_seed": master_seed,
            "params": {
                "n_stars": n_stars,
                "size": size,
                "band": band,
                "mag_range": tuple(mag_range),
                "slope": slope,
                "system": system,
                "spectra": spectra,
                "exclude_radius": exclude_radius,
            },
        }
        return field

    @classmethod
    def from_random_density(
        cls,
        density: u.Quantity | str,
        size: float | u.Quantity | tuple,
        *,
        band: str,
        mag_range: tuple[float, float],
        slope: float = 0.3,
        system: str = "Vega",
        spectra: SPECTRUM_TYPE | ArrayLike,
        position: POSITION_TYPE | None = None,
        exclude_radius: float | u.Quantity | None = None,
        rng_seed: int | None = None,
    ) -> "StarField":
        """Random stars at a given surface density, e.g. as a background.

        The number of stars is Poisson-distributed with mean `density` times
        the area stars can occupy (the rectangle minus any exclusion circle),
        so `density` describes the rendered field. Everything else is as in
        :meth:`from_random`, which this wraps with the same master seed.

        Parameters
        ----------
        density : Quantity or str
            Surface density per solid angle, e.g. ``3 / u.arcmin**2`` or
            ``"3 / arcmin2"``.
        size, band, mag_range, slope, system, spectra, position, \
        exclude_radius, rng_seed
            See :meth:`from_random`.

        Warns
        -----
        EmptyStarFieldWarning
            If the draw gives zero stars. The (empty) field is still returned.
        """
        master_seed = resolve_seed(rng_seed)
        density = u.Quantity(density)
        try:
            per_arcsec2 = density.to_value(u.arcsec**-2)
        except u.UnitConversionError:
            raise ValueError(
                f"density must be per solid angle (e.g. 3 / arcmin2), got {density}"
            ) from None
        if not (np.isscalar(per_arcsec2) and per_arcsec2 >= 0):
            raise ValueError(f"density must be a non-negative scalar, got {density}")

        center = cls(position=position).resolve_position()
        # Area only depends on the extent, not on the (still unknown) count.
        area = UniformRectangularMorphology(
            0, size, exclude_radius, rng_seed=master_seed
        ).area(center)
        expected = per_arcsec2 * area.to_value(u.arcsec**2)
        n_stars = int(new_rng(stream(master_seed, "count")).poisson(expected))

        if n_stars == 0:
            warnings.warn(
                EmptyStarFieldWarning(
                    f"Random star field drew 0 stars (expected {expected:.3g} "
                    f"in {area.to_value(u.arcmin**2):.3g} arcmin2 at "
                    f"{density}); returning an empty field."
                ),
                stacklevel=2,
            )

        field = cls.from_random(
            n_stars, size,
            band=band, mag_range=mag_range, slope=slope, system=system,
            spectra=spectra, position=position,
            exclude_radius=exclude_radius, rng_seed=master_seed,
        )
        field.generated_by = {
            "method": "from_random_density",
            "rng_seed": master_seed,
            "params": {
                **field.generated_by["params"],
                "density": density,
                "expected_n_stars": expected,
            },
        }
        return field

    @property
    def position(self) -> SkyCoord:
        """Field center; per-star ``positions`` are offsets from it."""
        return self._position

    @position.setter
    def position(self, position: POSITION_TYPE):
        center = self._parse_position(position)
        if not center.isscalar:
            raise ValueError(
                "StarField `position` is the single field center; per-star "
                "positions go in `positions`"
            )
        self._position = center

    @property
    def offsets(self) -> u.Quantity | None:
        """Per-star offsets from the field center, shape ``(N, 2)`` [arcsec].

        Column 0 is x (+East), column 1 is y (+North).
        """
        return getattr(self, "_offsets", None)

    @property
    def positions(self) -> SkyCoord | None:
        """Absolute per-star positions (array-valued SkyCoord).

        Resolved from :attr:`offsets` against the field center (see
        :meth:`~.target.Target.resolve_position`, default (0, 0)).
        """
        offsets = self.offsets
        if offsets is None:
            return None
        center = self.resolve_position()
        local = SkyCoord(
            lon=offsets[:, 0], lat=offsets[:, 1], frame=center.skyoffset_frame()
        )
        return local.transform_to(center.frame.replicate_without_data())

    @positions.setter
    def positions(self, positions):
        offsets = self._parse_offsets(positions)
        try:
            guard_same_len(offsets, self._per_star_spectra, self.brightnesses)
        except ValueError as err:
            raise ValueError(
                "Positions length doesn't match other attributes"
            ) from err
        self._offsets = offsets

    def _parse_offsets(self, positions) -> u.Quantity:
        """Normalize per-star positions into an ``(N, 2)`` arcsec offset array.

        Accepted, fastest first:

        * an ``(N, 2)`` array or nested sequence of plain numbers [arcsec];
        * an ``(N, 2)`` angular Quantity, or a sequence of angular Quantity
          ``(x, y)`` rows;
        * a sequence of ``(x, y)`` pairs or ``{"x": ..., "y": ...}`` mappings,
          each coordinate a number [arcsec] or an angular Quantity;
        * absolute positions as a (array-valued) SkyCoord or a sequence of
          scalar SkyCoords -- converted to offsets from the *current* field
          center, so set ``position`` first.
        """
        if isinstance(positions, (str, bytes)):
            raise TypeError(f"Unknown positions format {positions!r}")
        if isinstance(positions, SkyCoord) or (
            isinstance(positions, Sequence)
            and len(positions) > 0
            and all(isinstance(pos, SkyCoord) for pos in positions)
        ):
            coords = (
                positions if isinstance(positions, SkyCoord)
                else SkyCoord(list(positions))
            )
            if coords.isscalar:
                coords = coords.reshape((1,))
            local = coords.transform_to(self.resolve_position().skyoffset_frame())
            # Plain values: Longitude and Latitude refuse to stack together.
            # (Offset-frame lon wraps at 180 deg, so West offsets are negative.)
            offsets = np.column_stack(
                [local.lon.to_value(u.arcsec), local.lat.to_value(u.arcsec)]
            ) << u.arcsec
        else:
            try:
                # Fast path: numbers, arrays, angular Quantities (or rows of
                # them) -- converted, never unit-stripped. Pairs mixing numbers
                # and Quantities, and mappings, raise and fall through.
                offsets = u.Quantity(positions, u.arcsec)
            except (TypeError, ValueError):
                # (unit errors re-raise from the per-pair parse)
                offsets = u.Quantity(
                    [self._parse_offset_pair(pos) for pos in positions]
                )

        if offsets.size == 0:
            offsets = offsets.reshape((0, 2))  # empty field
        if offsets.ndim != 2 or offsets.shape[1] != 2:
            raise ValueError(
                f"positions must be N (x, y) pairs, got shape {offsets.shape}"
            )
        return offsets

    @staticmethod
    def _parse_offset_pair(position) -> u.Quantity:
        match position:
            case Mapping() if set(position) == {"x", "y"}:
                x_offset, y_offset = position["x"], position["y"]
            case Mapping():
                # Per-star distances (or anything else) have no meaning here;
                # a field-wide distance belongs on the field center.
                raise ValueError(
                    "per-star position mappings take only 'x' and 'y', got "
                    f"{sorted(position)}; put a distance on the field center "
                    "`position`"
                )
            case u.Quantity() if position.shape == (2,):
                x_offset, y_offset = position
            case (x_offset, y_offset):
                pass
            case _:
                raise TypeError(f"Unknown per-star position {position!r}")
        return u.Quantity(
            [u.Quantity(x_offset, u.arcsec), u.Quantity(y_offset, u.arcsec)]
        )

    @property
    def spectra(self) -> list | None:
        """Per-star spectra; a shared spectrum is repeated for every star."""
        common = getattr(self, "_common_spectrum", None)
        if common is not None:
            n_stars = self._n_stars()
            return [common] * (1 if n_stars is None else n_stars)
        return getattr(self, "_spectra", None)

    @spectra.setter
    def spectra(self, spectra: SPECTRUM_TYPE | Sequence[SPECTRUM_TYPE]):
        if isinstance(spectra, (str, SpectralType, SourceSpectrum)):
            # One spectrum for all stars: no length to match.
            self._common_spectrum = self._parse_spectrum(spectra)
            self._spectra = None
            return

        try:
            guard_same_len(self.offsets, spectra, self.brightnesses)
        except ValueError as err:
            raise ValueError(
                "Spectra length doesn't match other attributes"
            ) from err
        self._common_spectrum = None
        self._spectra = [
            self._parse_spectrum(spectrum) for spectrum in spectra
        ]

    @property
    def _per_star_spectra(self) -> list | None:
        """Explicit per-star spectra for length guards (None if shared)."""
        return getattr(self, "_spectra", None)

    def _n_stars(self) -> int | None:
        """Number of stars, from whichever per-star attribute is set."""
        for per_star in (self.offsets, self._per_star_spectra, self.brightnesses):
            if per_star is not None:
                return len(per_star)
        return None

    @property
    def brightnesses(self):
        try:
            return self._brightnesses
        except AttributeError:
            pass  # return None

    @brightnesses.setter
    def brightnesses(self, brightnesses: Sequence[BRIGHTNESS_TYPE]):
        try:
            guard_same_len(self.offsets, self._per_star_spectra, brightnesses)
        except ValueError as err:
            raise ValueError(
                "Brightnesses length doesn't match other attributes"
            ) from err
        self._brightnesses = [
            self._parse_field_brightness(brightness)
            for brightness in brightnesses
        ]

    def _parse_field_brightness(self, brightness) -> Brightness:
        """Parse one per-star entry: a full spec, or a bare amount + `band`.

        Full specs (mapping or ``(locator, amount)`` sequence) are parsed as
        given; anything else -- a number, a scalar Quantity, an amount string
        like ``"15 mag"`` -- is a bare amount located at the field's `band`.
        The split is structural (:func:`~.brightness.is_brightness_spec`), so
        string amounts and mappings are no longer misrouted.
        """
        if not is_brightness_spec(brightness):
            if self.band is None:
                raise LocatorError(
                    f"bare brightness amount {brightness!r} needs a locator: "
                    "set StarField's `band`, or give a full "
                    "(locator, amount) / mapping spec"
                )
            brightness = (self.band, brightness)

        parsed = self._parse_brightness(brightness)
        if isinstance(parsed, FromSpectralType):
            # The resolver lives on the single-target `brightness` property;
            # per-star resolution (and per-star distances) is not wired up.
            raise TypeError(
                "from_spectral_type is not supported for StarField entries"
            )
        return parsed

    def to_source(self, optical_train=None) -> Source:
        """Convert to ScopeSim Source object."""
        # Offsets are already coordinates in the field center's offset frame,
        # which is exactly ScopeSim's local frame: no transform needed.
        # .round(6) is microarcsec, as in `_xy_arcsec_position`.
        x_positions = self.offsets[:, 0].to_value(u.arcsec).round(6)
        y_positions = self.offsets[:, 1].to_value(u.arcsec).round(6)

        # First-seen order (not set order, which follows per-process string
        # hash randomization), so refs are reproducible across runs.
        spectra_ids = {
            spectrum: spectrum_id
            for spectrum_id, spectrum in enumerate(dict.fromkeys(self.spectra))
        }
        resolved_spectra = {
            # TODO: Implement redshift from position.
            spectrum_id: self.resolve_spectrum(spectrum)
            for spectrum, spectrum_id in spectra_ids.items()
        }

        spec_refs = [spectra_ids[spectrum] for spectrum in self.spectra]

        # One batched scale per unique spectrum: the photometry then runs once
        # per (spectrum, band, amount kind) rather than once per star.
        weights = [None] * len(spec_refs)
        for spectrum_id, spectrum in resolved_spectra.items():
            indices = [i for i, ref in enumerate(spec_refs) if ref == spectrum_id]
            scales = self._anchored_spectrum_scales(
                spectrum, [self.brightnesses[i] for i in indices]
            )
            for i, scale in zip(indices, scales):
                weights[i] = scale

        if not resolved_spectra:
            # Empty field: ScopeSim requires a non-empty spectra mapping even
            # when no row references it; a zero-flux placeholder is harmless.
            resolved_spectra = {
                0: SourceSpectrum(ConstFlux1D, amplitude=0 * synphot_units.PHOTLAM)
            }

        # TODO: Refactor...
        table = Table(
            names=["x", "y", "ref", "weight"],
            units={"x": u.arcsec, "y": u.arcsec},
            data={
                "x": x_positions,
                "y": y_positions,
                "ref": np.asarray(spec_refs, dtype=int),
                "weight": np.asarray(weights, dtype=float),
            },
        )

        # TODO: Figure out if those are really needed
        table.meta["x_unit"] = "arcsec"
        table.meta["y_unit"] = "arcsec"

        return Source(field=TableSourceField(table, spectra=resolved_spectra))


def _parse_grid_shape(shape) -> tuple[int, int]:
    """Normalize a grid shape into ``(nx, ny)`` positive ints."""
    def _is_count(value) -> bool:
        return isinstance(value, Integral) and not isinstance(value, bool)

    if _is_count(shape):
        shape = (shape, shape)
    elif not (isinstance(shape, Sequence) and len(shape) == 2):
        raise TypeError(
            f"shape must be an int or (width, height), got {shape!r}"
        )
    if not all(_is_count(n) for n in shape):
        raise TypeError(f"shape entries must be ints, got {shape!r}")
    if not all(n >= 1 for n in shape):
        raise ValueError(f"shape entries must be at least 1, got {shape!r}")
    return int(shape[0]), int(shape[1])


def _is_single_spectrum(value) -> bool:
    return isinstance(value, (str, SpectralType, SourceSpectrum)) or not (
        isinstance(value, (Sequence, np.ndarray))
    )


def _flat_per_star(values, n_stars: int, name: str, is_single) -> list:
    """Validate a flat, length-`n_stars` sequence of single per-star values."""
    values = list(values)
    if not all(is_single(value) for value in values):
        raise ValueError(
            f"{name} must be a single value or a flat sequence (no nested or "
            "2D arrays)"
        )
    if len(values) != n_stars:
        raise ValueError(
            f"{name}: expected {n_stars} per-star values (nx * ny), "
            f"got {len(values)}"
        )
    return values
