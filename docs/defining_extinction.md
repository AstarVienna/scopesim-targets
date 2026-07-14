# Defining extinction and attenuation
Dust along the way from a source to the observer dims and reddens its light.
ScopeSim-Targets splits this into **two** attributes, following the distinction the [`dust_extinction`](https://dust-extinction.readthedocs.io/) package draws between *extinction* and *attenuation*:

| attribute | what it is | whose property | varies within a component? | where it is applied |
|---|---|---|---|---|
| `extinction` | the **sightline**: Galactic foreground dust between observer and scene | the scene (cascades; per-component override) | no — one screen | in ScopeSim-Targets, in place on the spectrum |
| `attenuation` | the **source**: natal cloud, circumstellar envelope | the component (never inherited) | yes — one column per star | by the simulation engine, from a per-row table column |

*Extinction* in the strict sense is what the pair method measures: a single star behind a dust screen, no geometry.
*Attenuation* is the effect of dust on groups of stars mixed with dust, or on a star inside its own dust shell, where the column differs from star to star and light scattered back *into* the line of sight matters.
An embedded young cluster needs both: the Galactic foreground (extinction) and the natal cloud it sits in (attenuation).

Atmospheric extinction is handled by `ScopeSim` itself as part of the instrument/site model, and extragalactic (host-internal) attenuation of *extended* sources is out of scope, consistent with the ESO target specification.

## Nomenclature
The name _extinction_ was chosen for the foreground attribute as the most established and least ambiguous term; _attenuation_ is the term `dust_extinction` and the radiative-transfer literature use for the source-side case, which is exactly our second attribute.

### Rejected alternatives
#### Reddening
Describes only the color effect, not the dimming.
Whether the dimming actually applies to a given target depends on the anchor frame (see below), but the underlying physical quantity always includes it.

#### Attenuation (for the foreground)
Attenuation refers to dust *mixed with* or *around* the sources, where geometry matters.
That is not the foreground case — but it is precisely the local case, so the term names the second attribute rather than being rejected outright.

#### A_V / Av
Ties the concept to one specific band.
Extinction toward embedded regions is more commonly measured (and specified) in the near-infrared, e.g. $A_{K_s}$ — and the Alves et al. column-density PDFs used for clusters are given in $A_K$.

#### Column density
Both attributes are, physically, dust columns; the literature we draw the cluster distributions from (Alves, Lombardi & Lada 2017) speaks only of column density.
A column alone does not say how it acts on the light, which is why `attenuation` pairs a column with a *model* (below).

## Extinction: the foreground screen
A single extinction _screen_ is described by:
* an **amount**, given either as a magnitude in a band (`value` + `band`, default band `V`), or as a color excess `ebv` = $E(B-V)$, or resolved from a dust map (`from_map`),
* a reddening **law** (default `ccm89`, i.e. Cardelli, Clayton & Mathis 1989, matching the ESO specification),
* a total-to-selective ratio **rv** = $R_V$ (default `3.1`).

The law vocabulary (`ccm89`, `f99`, `g23`) follows `dust_extinction`.
The screen is a pure foreground screen: no differential extinction within a component, no scattering geometry.

### Possible ways to specify the screen
* A bare `astropy.Quantity` in `u.mag` (or the equivalent quantity string in YAML): interpreted as $A_V$ with the default law and $R_V$.
* A dictionary/keyword form with any of `value`, `band`, `ebv`, `law`, `rv`.
  Exactly one of `value` or `ebv` must be given (`band` only pairs with `value`).
* `{from_map: <map_id>}`: resolve the amount from a directional dust map using the target's `position` (including `distance`).
  _Parsed but not yet resolved (see “Implementation status”)._
  Resolution happens at load time; the resolved value **and the map name/version** are stored, so that exports remain reproducible even when maps are revised.

There is deliberately **no list of screens**.
For a fixed law and $R_V$, $A(\lambda)$ is linear in $A_V$, so two screens compose *exactly* into one with the summed amount; and dust local to a source is not a second foreground screen but `attenuation`.

### Inheritance
Extinction cascades through the scene hierarchy, using the same rule as `role`:
the nearest explicitly set value wins, walking up the containment chain, and a `TargetScene`-level screen acts as the shared default for all components.
An override **replaces** the inherited screen.
The explicit opt-out — e.g. for a field star *in front of* the cloud the scene-level screen describes — is a zero screen:

```yaml
extinction: 0 mag
```

### Delivery
Because the screen is uniform over a component, it is resolved in ScopeSim-Targets and applied to the spectrum **in place**, for every field type.
For point-source tables (`StarField`, clusters) this means the *shared, deduplicated* templates are reddened once; no table column is emitted for the foreground.

## Attenuation: dust at the source
`attenuation` has a two-part grammar mirroring how the column-density literature speaks:

1. a dust **column**, in extinction currency — exactly one of
   * `value` + `band` (default `V`): the column as extinction in a band, e.g. the $A_{K_s}$ of a NICER map,
   * `ebv`: the column as $E(B-V)$,
   * `tau_v`: the column as V-band optical depth,
   * `distribution` + `params` (generative targets only): a PDF from which one column per star is drawn, see below;
2. a **model** saying how that column attenuates the emergent light:

| `model` | curve | own keys | validity |
|---|---|---|---|
| `screen` (default) | per-star *differential extinction*: $A(\lambda)$ from a `dust_extinction` law | `law` (default `ccm89`), `rv` (default `3.1`; dense clouds 4–5) | the law's range (`g23`: 0.09–32 µm) |
| `wg00` | Witt & Gordon (2000) radiative-transfer curves from [`dust_attenuation`](https://dust-attenuation.readthedocs.io/): scattering into the line of sight, geometry-aware | `geometry` (`shell` \| `cloudy` \| `dusty`, **required**), `dust_distribution` (`homogeneous` \| `clumpy`, **required**), `dust_type` (`mw` default \| `smc`) | 0.1–3.0 µm, $\tau_V \in [0.25, 50]$ |

`geometry` and `dust_distribution` have no defaults: the geometry *is* the physics, and hiding it behind a default would be exactly the kind of implicit assumption the rest of the spec avoids.
`shell` is the circumstellar envelope / embedded YSO; `cloudy` and `dusty` are stars inside a dust cloud.

The integrated-light attenuation averages (`c00`, `n09`, `l02`) are deliberately **not** in v1: they describe *unresolved* populations, stop at 2.2 µm (Calzetti 2000 is exactly zero there), and have no meaning for a resolved star.
Adding them is an additive, versioned schema change.

```yaml
attenuation: 8.5 mag                                   # sugar: model screen, A_V, ccm89, R_V 3.1
attenuation: {value: 0.97 mag, band: Ks, rv: 4.7}      # screen, column given as A_Ks
attenuation: {model: wg00, geometry: shell, dust_distribution: homogeneous,
              value: 0.97 mag, band: Ks}               # envelope, same column
```

### Column conversion
Each model is parameterised by one **native** amount: $A_V$ for `screen`, $\tau_V$ for `wg00`.
A column given in another currency is converted with the model's *own* dust curve:
`screen` divides a band amount by the law's $A(\lambda_\mathrm{band})/A_V$ (and uses $A_V = R_V\,E(B-V)$, $A_V = 1.086\,\tau_V$);
`wg00` divides by its tabulated extinction per unit $\tau_V$ at the band wavelength (exactly linear in $\tau_V$; $E(B-V)$ via $A_B - A_V$).
Note that the WG00 tabulation gives 1.109 mag of extinction per unit $\tau_V$ at 0.551 µm rather than the nominal 1.086; the conversion uses the model's curve consistently for every band rather than special-casing $V$.

### Scope
Attenuation is a **point-source** attribute (`Star`, `Binary`, `StarField`, the cluster classes).
It is never scene-level and never inherited: it describes the component's own environment.
Profiles (`Sersic`, `Disk`, …) do not accept it — host-internal attenuation of extended sources is out of scope.

### Sampled attenuation for generative targets
Generative targets (e.g. the young-cluster classes) may give the column as a *distribution* rather than a value, mirroring how their population and morphology parameters work.
One column per star is drawn at expansion time, seeded by the target's `rng_seed`, and the realized values become the ordinary per-row column in the expanded (and exported) form.

```yaml
attenuation:
  model: screen
  law: ccm89
  rv: 4.5
  band: Ks                           # the PDF is written in A_K, as Alves et al. do
  distribution: column_density_pdf   # lognormal core + power-law tail
  params:
    median: 0.6 mag    # median A_Ks of the lognormal core
    sigma: 0.45        # width of the core, in ln(A_K)
    break: 1.2 mag     # A_Ks where the power-law tail takes over
    slope: 2.0         # tail index: p(A_K) proportional to A_K**-slope (>1)
```

The `column_density_pdf` form follows the observed shape of molecular-cloud column-density PDFs (Alves, Lombardi & Lada 2017, A&A 606, L2): a lognormal core that hands over to a power-law tail at `break`, continuous at the break.
The cloud's gas-to-dust conversion is folded into the choice of these parameters.
The simpler `lognormal` (`median`, `sigma`) and `uniform` (`low`, `high`) distributions are also available.
Column-valued parameters accept plain numbers or `"N mag"` strings in `band`; with `currency: tau_v` they are dimensionless optical depths instead (the natural currency for `wg00`).

Note that this is a zeroth-order model: the column is sampled independently of stellar position.
Correlated sampling (e.g. deeper embedding toward the cluster center) is a planned extension sharing the same RNG infrastructure.

### Delivery: the source-table contract
Because the column varies per star, attenuation cannot be folded into the deduplicated spectra of a point-source field.
It is therefore **not** resolved into the SED by ScopeSim-Targets but *emitted* for the simulation engine:

* a per-row column **`attenuation`** (float) in the model's native amount — $A_V$ in `mag` for `screen`, dimensionless $\tau_V$ for `wg00`; `meta["attenuation_parameter"]` names it (`Av` / `tau_V`);
* flat `table.meta` keys: `attenuation_model`, plus the model's own parameters — `attenuation_law`, `attenuation_rv` (screen) or `attenuation_geometry`, `attenuation_dust_type`, `attenuation_dust_distribution` (wg00) — and `anchor`;
* an absent column means *no attenuation*; a table without attenuation is byte-identical to before.

The engine multiplies each row's spectrum by $10^{-0.4\,A_\mathrm{att}(\lambda)}$ for that row's amount.
`scopesim_targets.attenuation.Attenuation.from_table_meta(meta).curve(wave, amount)` is the reference implementation it can import; it clamps wavelengths outside the model's validity range to the boundary value (curve held constant, the same rule as for extinction) and, for `wg00` below its tabulated $\tau_V = 0.25$, scales the 0.25 curve linearly (optically-thin limit).
A build that ignores the column simply renders without attenuation.

**v1 limitations.** One model (and one law/$R_V$ or geometry) per source table.
The curve is applied at *observed* wavelengths although the dust is at the source; this is exact for the $z \approx 0$ science cases and a documented approximation otherwise (the engine would shift the curve by $1+z$).

## Interaction with brightness: anchor frames
Dust affects both the spectral shape and the flux scale, so its meaning depends on whether the target's `brightness` is an *observed* (i.e. already dimmed) or an *intrinsic* value.
This is controlled by the `anchor` attribute:

* `anchor: observed` (default): the attenuated and reddened spectrum is scaled so that its synthetic photometry in the brightness band matches the given value.
  This is the catalog and ETC convention — dust then only changes *colors* relative to that band.
* `anchor: intrinsic`: the unextincted spectrum is scaled to the given brightness first, then attenuation and the screen are applied.
  Dust then both dims *and* reddens.
  This is the natural setting for generative targets, where brightness comes from an isochrone.
* A third value, `anchor: absolute`, extends `intrinsic` for values referring to 10 pc (the distance modulus is applied before the screen); see [defining brightness](defining_brightness.md).

For point sources the anchor lives entirely in the emitted `weight`: an `observed` weight is set against the attenuated *and* extincted SED (ScopeSim-Targets evaluates the attenuation curve for that purpose, using the same reference implementation), so the engine applies the same per-row curve regardless of anchor.

The normative order of operations is:
rest-frame spectrum → attenuation (source frame) → distance modulus (`absolute`) → redshift/velocity shift → extinction screen (observer frame) → anchoring (for `observed`).
See [defining spectra](defining_spectra.md) for the rest-frame convention and [defining brightness](defining_brightness.md) for scaling.

## Implementation status
**Extinction** (single foreground screen) is **implemented**: bare `A_V` sugar, `0 mag` opt-out, the canonical mapping (`value`/`band`/`ebv`/`law`/`rv`), resolution to one $A_V$ (band amounts via the law at the band's effective wavelength; `ebv` via $A_V = R_V\,E(B-V)$), and in-place reddening of the spectrum — for profiles and for the shared templates of point-source tables alike.
The former list form and `[]` opt-out are rejected (E17).
`{from_map: ...}` is parsed but raises on resolution (needs the directional-map backend).

**Attenuation** is **implemented** at the module level (`attenuation.py`): both models, all column currencies with conversion through the model's own curve, the three distributions, the table contract above with `from_table_meta` as the engine-side entry point, and the `ZeroAgeCluster` wiring (fixed column broadcasts; sampled column drawn from the attenuation child of the cluster's seed tree, so the whole realization is reproducible under one `rng_seed`).
`Star`, `Binary`, `StarField` and `PlanetarySystem` emit the column (a fixed column broadcast to every row; planets inherit the primary's) and set an `observed` weight against the attenuated, reddened SED; `Exoplanet`-level or profile-level attenuation is rejected (E22).
The `wg00` model needs `dust_attenuation`, which is not yet on PyPI (`pip install git+https://github.com/karllark/dust_attenuation.git`) and is imported lazily.

**Not yet implemented.** The `role`-style *cascade/inheritance* of `extinction` across a scene hierarchy — it is currently a plain local attribute. Position-correlated column sampling. `from_map`.

## To be discussed
Position-correlated attenuation sampling for embedded clusters?

Default `from_map` resolver as the suggested value in interactive use (cf. ESO spec 5.3.7)?

Should the engine shift the attenuation curve by $1+z$ (needs `z`/`radial_velocity` in the table meta), or is the v1 observed-wavelength application sufficient for the foreseeable science cases?
