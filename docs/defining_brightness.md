# Defining target brightness
All targets with separate spectral information (basically anything except a full datacube) usually also need some information on how to scale the spectrum.
The only exception to this are user-supplied spectra, which are already flux-calibrated to the desired level.
The `brightness` attribute is the single **flux authority** of a target: no other attribute or model parameter may independently imply a flux (see [Amplitude](#amplitude) below).

## Nomenclature
The name _brightness_ was chosen as the least ambiguous of several rejected alternatives and because it matched with the ESO standard.
It is not usually associated with any particular unit or physical type, but rather more of a general "how bright is my source".

### Rejected alternatives
Some of these were used in `ScopeSim_Templates`, with no particular standard.

#### Magnitude
Physical flux units can also be used to specify brightness.
Calling it _magnitude_ would suggest it only accepts magnitudes.

#### Flux
The opposite as with _magnitude_, this would suggest using a magnitude value is _not_ supported.

#### Amplitude
This name seems more fitting for the scaling of normalized profiles, can be confusing for point sources.
Note that the `amplitude` parameter of the underlying `astropy` profile models is **not accepted** in the `params` of extended targets: a surface-brightness-valued `brightness` (see below) *is* the amplitude specification, so accepting both would over-determine the flux.
Supplying `amplitude` in `params` is an error, with a message pointing here.

#### Luminosity
The second-best alternative, but still more associated with specific units, like $L_\odot$.

## Possible ways to specify brightness
Brightness is a two-slot value — **where** and **how much** — mirroring how astronomers speak: "10 mag in R", "3.5 mJy in K", "5 mJy at 230 GHz".

* **Where** (the locator): a band name (`"R"`, `"Ks"`, `"F160W"`), a wavelength (`"656.3 nm"`), or a frequency (`"230 GHz"`).
  Band names start with a letter and contain no whitespace, so the three forms are structurally unambiguous.
* **How much** (the amount): a magnitude or a flux quantity.
  The interpretation is dispatched on the **physical type of the unit** — the unit *is* the declaration, nothing is implicit:

| amount unit (examples) | physical type | interpretation |
|---|---|---|
| `mag`, `mag(AB)`, bare number | magnitude | magnitude in the band |
| `Jy`, `mJy`, `W / (m2 Hz)` | spectral flux density | flux density at the locator |
| `erg / (s cm2 Angstrom)`, `W / m3` | spectral flux density (per wavelength) | flux density at the locator |
| `W / m2` | energy flux | flux **integrated** over the band (or line) |
| any of the above `/ arcsec2` or `/ sr` | surface brightness / radiance | **surface brightness**, see below |

This covers all four brightness pairs of the ESO target specification (5.3.6) with one uniform syntax:

```yaml
brightness: ["R", "15 mag"]                           # band, magnitude
brightness: ["K", "3.5 mJy"]                          # band, spectral flux density
brightness: ["K", "4.2e-18 W / m2"]                   # band, irradiance
brightness: ["656.3 nm", "1.2e-16 erg / (s cm2 nm)"]  # wavelength, spectral flux density (wav)
brightness: ["230 GHz", "5 mJy"]                      # frequency, spectral flux density
```

The fully explicit (canonical) form is a mapping with exactly one locator key:

```yaml
brightness: {band: V, value: 21.4, system: AB}
brightness: {frequency: 230 GHz, value: 5 mJy}
```

### Photometric systems
The photometric system (zero-point reference) is metadata of a *magnitude*, **not** part of the unit:

* A bare number or a plain `mag` quantity means **Vega**, by convention of this framework (note that `astropy`'s plain `u.mag` carries no zero point of its own).
* AB and ST magnitudes use astropy's function-unit string syntax: `("R", "10 mag(AB)")`, `("R", "10 mag(ST)")`; `"mag(Vega)"` is accepted as an explicit alias.
* `system` in the mapping form is only meaningful for magnitude values, and is mutually exclusive with a system already embedded in the value string.

```{warning}
This framework defaults to **Vega**; the ESO target specification defaults to
**AB**. Converters to/from the ESO FITS keywords must therefore always write
the system explicitly and never rely on either default.
```

AB and ST magnitudes require a band exactly like Vega magnitudes do — the system only fixes the reference spectrum, not the bandpass.

## Integrated flux vs. surface brightness
Whether the amount refers to the **total integrated** flux of the target or to a **surface brightness** is read directly from the unit: a per-solid-angle divisor (`/ arcsec2`, `/ sr`) selects surface brightness; its absence selects integrated flux.
The ESO specification's "/arcsec2 implicit for extended object" is thereby made *explicit* — the same number never has two possible meanings.

For surface-brightness magnitudes, the verified astropy spellings are:

```yaml
brightness: ["V", "21.5 mag / arcsec2"]        # Vega surface brightness
brightness: ["V", "21.5 mag(AB / arcsec2)"]    # AB surface brightness
```

(`"mag(AB) / arcsec2"` does **not** parse — astropy function units cannot be divided — and is rejected.)

### What surface brightness means per profile
For non-uniform profiles, "the" surface brightness needs a reference point.
Each morphology class documents its own; a surface-brightness `brightness` sets the profile amplitude at exactly that point:

| target class | SB reference point |
|---|---|
| `Disk` (uniform) | the uniform surface brightness |
| `Sersic` | at $r_\mathrm{eff}$ (astropy `Sersic2D` amplitude convention) |
| `Powerlaw` | at $r_0$ |
| `Flat` (infinite constant) | the constant surface brightness |

### Validity matrix
Not every combination is meaningful; invalid ones are **errors**, never silent reinterpretations:

| | integrated amount | surface-brightness amount |
|---|---|---|
| point source | OK | error (a point has no solid angle) |
| profile with finite analytic total (`Sersic`, `Disk`, truncated `Powerlaw`) | OK (normalized to the **analytic** total, never a grid sum) | OK (sets amplitude at the reference point) |
| profile without finite total (`Flat`, untruncated `Powerlaw`) | error (the integral diverges) | OK (the only valid form; the flux in the field of view is then legitimately FOV-dependent) |

### How a non-integrable profile carries its flux
A profile without a finite analytic total (`Flat`, untruncated `Powerlaw`) can
carry only a *surface brightness* — there is no intrinsic total to normalize to.
The total is instead realized from the **rendered field of view**: the flux
authority is fixed only once the pixel scale and field dimensions are known (i.e.
at render time), giving a total flux of surface brightness times the field solid
angle, $F = \mathrm{SB}\,A_\mathrm{FOV}$ with $A_\mathrm{FOV} = \Omega_\mathrm{pixel}\,N_\mathrm{pixels}$.

Crucially, the *stored representation is the same* as for finite profiles — a
dimensionless weight map plus a flux-calibrated spectrum — **not** an image in
surface-brightness units. (The spectrum must be flux-calibrated: a shape-only,
dimensionless spectrum is not a valid `synphot.SourceSpectrum`, so the flux cannot
be moved into the image while leaving the spectrum dimensionless.) The spectrum
carries the total and the uniform weight map distributes it, so each pixel
receives $\mathrm{SB}\,\Omega_\mathrm{pixel}$.

The only quantity that differs between profiles is the **effective area** a surface
brightness is multiplied by: the analytic `total_flux_factor` for finite profiles
(field-of-view independent), and $A_\mathrm{FOV}$ for non-integrable ones
(field-of-view dependent). Everything else — the reduction to an equivalent
integrated amount, the weight map, the flux-calibrated spectrum — is identical.

Two consequences, both physically correct rather than artefacts:
* the same non-integrable target yields a **larger total flux in a larger field**
  (more solid angle at the same surface brightness);
* per-pixel values scale with $\Omega_\mathrm{pixel}$, while the flux in a fixed
  aperture is **pixel-scale invariant**.

### Interpretation notes
* **Integrated brightness** of an analytic profile means the total of the *analytic* model — not a surface brightness, and not the flux within the rendered field of view.
  Normalization uses the closed-form integral, so profile flux falling outside the simulated detector window is *not* redistributed into it.
* **Observed vs. intrinsic**: whether `brightness` refers to the observed (extincted) or intrinsic value is controlled by the `anchor` attribute; default `observed`. The same rule applies to surface brightnesses. See [defining extinction](defining_extinction.md).

## Absolute magnitudes
Absolute magnitudes need no new syntax: absolute-ness is not a property of the
*number* but of the **frame it is anchored in**, so it is the third value of the
`anchor` attribute (see [defining extinction](defining_extinction.md) for the
other two):

```yaml
brightness: ["V", "4.83 mag"]
anchor: absolute
```

* `anchor: absolute` means: the value refers to the **unextincted** source at
  **10 pc**. The distance modulus from `position.distance` is applied first,
  then extinction screens — so `absolute` implies `intrinsic` semantics
  (extinction dims *and* reddens).
  Absolute magnitudes that are *not* corrected for extinction are deliberately
  not representable; convert before input.
* A missing `distance` is an **error** — no default distance is assumed.
* An absolute anchor with a **surface-brightness** amount is an **error**:
  surface brightness is distance-invariant, so "absolute surface brightness" is
  a category error.
* The rule applies uniformly to all amount kinds (inverse-square scaling), not
  only magnitudes.

## Brightness from spectral type
As an explicit opt-in — never as a fallback — the brightness can be resolved
from the target's spectral type via a lookup table:

```yaml
spectrum: K5V
brightness: {from_spectral_type: mamajek}
```

This follows the same *resolver* pattern as `{from_map: ...}` for extinction:
resolution happens at load time, and the resolved value **plus the table
name/version** are recorded, keeping exports reproducible.
The resolver yields the absolute magnitude for the spectral type and sets
`anchor: absolute`; the ordinary machinery then applies the distance modulus
(a missing `distance` is an error) and any extinction screens.
Combined with a per-star extinction distribution, this makes a single star
consistent with what the generative cluster classes do internally
(isochrone magnitude → distance → screens).

```{note}
The lookup carries population assumptions (age, metallicity, luminosity
class). It is therefore never applied implicitly when `brightness` is merely
missing — deriving flux from assumptions must be a visible decision in the
target definition. In Python, `Star.from_spectral_type(...)` wraps the same
resolver.
```
_To be implemented._

## To be discussed
Multiple brightness values per target (sparse SED sampling), with dynamic
selection of the anchoring band by the consumer? (Cf. the multi-magnitude
groups in the ESO specification.)

Monochromatic AB magnitudes at a stated wavelength (band-free, cf. quasar
continuum conventions)?
