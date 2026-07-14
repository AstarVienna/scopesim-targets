# Defining extinction
Extinction describes the wavelength-dependent dimming and reddening of a target's light by dust along the line of sight.
In ScopeSim-Targets, this always refers to **interstellar** extinction (Galactic foreground, molecular clouds) or **circumstellar** extinction (envelopes around embedded YSOs).
Atmospheric extinction is handled by `ScopeSim` itself as part of the instrument/site model, and extragalactic (host-internal) attenuation is out of scope, consistent with the ESO target specification.

Alongside `position`, `spectrum` and `brightness`, extinction is the fourth core target attribute.
Unlike the other three, it is optional and *inherited*: most targets never specify it directly.

## Nomenclature
The name _extinction_ was chosen for the attribute as the most established and least ambiguous term.

### Rejected alternatives
#### Reddening
Describes only the color effect, not the dimming.
Whether the dimming actually applies to a given target depends on the anchor frame (see below), but the underlying physical quantity always includes it.

#### Attenuation
In the literature, _attenuation_ usually refers to the extragalactic case, where dust and stars are mixed and scattering geometry matters.
That is exactly the case we exclude, so using the term would be misleading.

#### A_V / Av
Ties the concept to one specific band.
Extinction toward embedded regions is more commonly measured (and specified) in the near-infrared, e.g. $A_{K_s}$.

## The screen model
A single extinction _screen_ is described by:
* an **amount**, given either as a magnitude in a band (`value` + `band`, default band `V`), or as a color excess `ebv` = $E(B-V)$,
* a reddening **law** (default `ccm89`, i.e. Cardelli, Clayton & Mathis 1989, matching the ESO specification),
* a total-to-selective ratio **rv** = $R_V$ (default `3.1`; dense-cloud sightlines typically want 4–5).

The law vocabulary follows the [`dust_extinction`](https://dust-extinction.readthedocs.io/) package.
Screens act as pure foreground screens: no differential extinction within a source, no scattering geometry.

### Possible ways to specify a screen
* A bare `astropy.Quantity` in `u.mag` (or the equivalent quantity string in YAML): interpreted as $A_V$ with the default law and $R_V$.
* A dictionary/keyword form with any of `value`, `band`, `ebv`, `law`, `rv`.
  Exactly one of `value` or `ebv` must be given (`band` only pairs with `value`).
* `{from_map: <map_id>}`: resolve the amount from a directional dust map using the target's `position` (including `distance`).
  _Parsed but not yet resolved (see “Implementation status”)._
  Resolution happens at load time; the resolved value **and the map name/version** are stored, so that exports remain reproducible even when maps are revised.

## Inheritance and composition
Extinction cascades through the scene hierarchy, using the same rule as `role`:
the nearest explicitly set value wins, walking up the containment chain, and a `TargetScene`-level screen acts as the shared default for all components.
Two consequences of this rule:

* **Override replaces, it never adds.**
  A component that wants "the scene's foreground *plus* its own envelope" lists both screens explicitly:

  ```yaml
  extinction:
    - 0.3 mag                              # Galactic foreground, restated
    - {value: 8.5 mag, band: Ks, rv: 4.7}  # local envelope
  ```

  Screens in a list are applied multiplicatively (in transmission), in the listed order.
* **An empty list is an explicit opt-out**, e.g. for a field star *in front of* the cloud that the scene-level screen describes:

  ```yaml
  extinction: []
  ```

## Interaction with brightness: anchor frames
Extinction affects both the spectral shape and the flux scale, so its meaning depends on whether the target's `brightness` is an *observed* (i.e. already extincted) or an *intrinsic* value.
This is controlled by the `anchor` attribute:

* `anchor: observed` (default): the reddened spectrum is scaled so that its synthetic photometry in the brightness band matches the given value.
  This is the catalog and ETC convention — extinction then only changes *colors* relative to that band.
* `anchor: intrinsic`: the unextincted spectrum is scaled to the given brightness first, then the screens are applied.
  Extinction then both dims *and* reddens.
  This is the natural setting for generative targets, where brightness comes from an isochrone.
* A third value, `anchor: absolute`, extends `intrinsic` for values referring to 10 pc (the distance modulus is applied before the screens); see [defining brightness](defining_brightness.md).

The normative order of operations is:
rest-frame spectrum → redshift/velocity shift → extinction screens (observer frame) → anchoring (for `observed`).
See [defining spectra](defining_spectra.md) for the rest-frame convention and [defining brightness](defining_brightness.md) for scaling.

## Sampled extinction for generative targets
Generative targets (e.g. the young-cluster classes) may specify extinction as a *distribution* rather than a value, mirroring how their population and morphology parameters work.
A per-star amount is then drawn at expansion time, seeded by the target's `rng_seed`, and the realized values become ordinary per-component screens in the expanded (and exported) form.

```yaml
extinction:
  distribution: column_density_pdf   # lognormal core + power-law tail
  params:
    av_median: 6 mag     # median A_V of the lognormal core
    sigma: 0.45          # width of the core, in ln(A_V)
    av_break: 12 mag     # A_V where the power-law tail takes over
    slope: 2.0           # tail index: p(A_V) proportional to A_V**-slope (>1)
  law: ccm89
  rv: 4.5
```

The `column_density_pdf` form follows the observed shape of molecular-cloud column-density PDFs (Alves, Lombardi & Lada 2017, A&A 606, L2): a lognormal core that hands over to a power-law tail at `av_break`, continuous at the break. The cloud's gas-to-dust conversion is folded into the choice of these $A_V$-space parameters. The simpler `lognormal` (`av_median`, `sigma`) and `uniform` (`low`, `high`) distributions are also available. All amounts accept plain numbers or `"N mag"` strings.

Note that this is a zeroth-order model: extinction is sampled independently of stellar position.
Correlated sampling (e.g. deeper embedding toward the cluster center) is a planned extension.

## Implementation status
The single-screen and screen-list model is **implemented**: bare `A_V` sugar, the canonical mapping (`value`/`band`/`ebv`/`law`/`rv`), multiplicative screen lists, the empty-list opt-out, and the `observed`/`intrinsic`/`absolute` anchor semantics above. Screens resolve to a single total $A_V$ (band amounts are converted via the law at the band's effective wavelength; `ebv` via $A_V = R_V\,E(B-V)$; list amounts sum).

How the reddening is *delivered* depends on the field type, for scaling reasons rather than physical ones:

* **Point sources** keep their deduplicated, intrinsic spectra and carry a per-row `Av` column in the source table, with `extinction_law`, `extinction_rv` and `anchor` in the table `.meta`. ScopeSim applies the true curve per row at render time. This keeps the spectrum set collapsed for large generated fields/clusters (where per-spectrum flux handling would otherwise scale badly), and is backward compatible: a ScopeSim build that ignores the column simply renders without extinction. The `anchor` is handled entirely in the emitted `weight` (the `observed` weight is set against the reddened SED, `intrinsic`/`absolute` against the intrinsic one), so a supporting ScopeSim applies the same true curve regardless of anchor.
* **Extended (profile) sources** have a single SED with no dedup concern, so the screen is applied to that SED in place — self-contained and independent of ScopeSim support.

**v1 restrictions.** One `law` and one `R_V` per source table: a screen *list* must share a single law and $R_V$ (mixed laws/$R_V$ are rejected), which suffices because a single physical scene rarely mixes laws. Per-star (per-row) $A_V$ is supported for generative targets via the sampled-extinction distributions below.

**Sampled extinction** (the `column_density_pdf` / `lognormal` / `uniform` distributions above) is **implemented** for `ZeroAgeCluster`: the cluster's mandatory `rng_seed` seeds a `SeedSequence` that spawns three independent children — for the population (IMF), the morphology (positions) and the extinction draw — so the whole realization is reproducible under one seed. Each star draws its own $A_V$, which becomes a per-row `Av` column with the same `law`/`rv`/`anchor` meta as a fixed screen; a fixed screen on a cluster broadcasts as the zero-variance special case. A cluster's frame is fixed to `intrinsic` (its fluxes are unreddened apparent magnitudes with the distance modulus already applied, so extinction dims them; `observed` has no meaning without a single anchor band), so `anchor` is not a cluster parameter.

**Not yet implemented.** The `role`-style *cascade/inheritance* across a scene hierarchy is not wired — `extinction` is currently a plain local attribute, shared by construction for `Binary`/`StarField` and overridable per `Star`. `{from_map: ...}` is parsed but raises on resolution (needs the directional-map backend).

## To be discussed
Position-correlated extinction sampling for embedded clusters?

Default `from_map` resolver as the suggested value in interactive use (cf. ESO spec 5.3.7)?
