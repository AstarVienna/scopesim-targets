---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.19.5
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
  language: python
---

# Nebular backgrounds for embedded clusters
Young clusters are often still embedded in the gas and dust of their natal star-forming region.
This page shows how the two background spectra documented in {doc}`background_sed_orion` and {doc}`background_sed_ophiuchus`
become extended targets that can be laid under a simulated cluster.

Both backgrounds are {class}`~scopesim_targets.extended_source.CoredPowerLaw` targets with the *same* geometry;
they differ only in spectrum, surface-brightness level and anchor frame.
A pair of simulations with one or the other therefore isolates the effect of the background itself.
The geometry is matched to the King profile of the cluster, and the distance is taken from the same {class}`~astropy.coordinates.SkyCoord` as the cluster's position.
The flux bookkeeping behind all of this is explained in {doc}`flux_scaling`.

## Setup

```{code-cell} ipython3
:tags: [hide-output]

import os
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord
import pandas as pd
import seaborn as sns
from myst_nb import glue

from scopesim_targets.extended_source import CoredPowerLaw2D

# Folder holding the SED files, build scripts and nebular_backgrounds.py.
BKG_DIR = Path(os.environ.get("SCOPESIM_BACKGROUNDS_DIR", "../backgrounds")).resolve()
sys.path.insert(0, str(BKG_DIR))

from clusterbkg import nebular_backgrounds as nb

# King profile of the embedded reference cluster (embedded_cluster.yaml).
KING = {"r_core": 0.13 * u.pc, "r_tide": 2.0 * u.pc}
DISTANCE = 414 * u.pc
```

The following cell defines helper functions used only on this page.
Hidden by default to reduce clutter, expand if you want to see them:

```{code-cell} ipython3
:tags: [hide-cell]

sns.set_theme(style="darkgrid", palette="Set2")
tbl_sty = [{
    "selector": "th, td",
    "props": [("text-align", "right"), ("padding", "0.2em 0.8em")],
}]

def position(distance):
    """Sky position of the reference scene, at the given distance."""
    return SkyCoord(84.7 * u.deg, -6.9 * u.deg, distance=distance)

def angle(length, distance):
    return (length / distance * u.rad).to(u.arcsec)

def level_mag(brightness):
    """Magnitude of a ("band", "<m> mag / arcsec2") brightness definition."""
    return float(brightness[1].split()[0])

def extent_arcsec(grid):
    half = grid["pixel_scale"].to_value(u.arcsec / u.pixel) * grid["width"] / 2
    return (-half, half, -half, half)
```

## The profile
Around its illuminating star, the surface brightness of the Orion Nebula falls roughly as {math}`D^{-2}`, closer to {math}`D^{-3}` further out {cite:p}`ODell2010`.
The backgrounds use a cored power law, truncated at an outer radius:

```{math}
:label: bkg-profile

S(r) \;=\; S_\mathrm{ref} \left[\frac{1 + (r/r_\mathrm{c})^2}{1 + (r_\mathrm{ref}/r_\mathrm{c})^2}\right]^{-p/2}
\quad (r \le r_\mathrm{out}), \qquad S(r) = 0 \quad (r > r_\mathrm{out}),
```

where {math}`S_\mathrm{ref}` is the surface brightness given in `brightness`, applied at {math}`r_\mathrm{ref}`.
The truncation makes the total finite for any index {math}`p`:

```{math}
:label: bkg-total

P \;=\; \pi\, r_\mathrm{c}^2 \left(1 + \frac{r_\mathrm{ref}^2}{r_\mathrm{c}^2}\right)^{p/2}
\frac{(1 + X)^{1 - p/2} - 1}{1 - p/2}, \qquad X = \left(\frac{r_\mathrm{out}}{r_\mathrm{c}}\right)^2 ,
```

which becomes {math}`\pi r_\mathrm{c}^2 (1 + r_\mathrm{ref}^2/r_\mathrm{c}^2) \ln(1 + X)` at {math}`p = 2`.
So the background does not depend on the field of view, exactly like the King-profile cluster it is matched to.
Its weightmap sums to the fraction of the background's flux that falls inside the field.

The radii are set as follows:

* {math}`r_\mathrm{c}` is the King core radius of the cluster.
* {math}`r_\mathrm{out}` is the King tidal radius.
* {math}`p` = {glue:text}`bkg_index:.0f`, the inner-Huygens value of {cite:t}`ODell2010`.
* {math}`r_\mathrm{ref}` is *not* tied to the cluster. For Orion, it is where the reference level of the spectrum was measured:
  the projected distance of the PDRs4All H II-region aperture from θ¹ Ori C {cite:p}`Peeters2024`, which is {glue:text}`bkg_rref:.3f` pc.

```{code-cell} ipython3
:tags: [remove-cell]

glue("bkg_index", float(nb.DEFAULT_INDEX), display=False)
glue("bkg_rref", float(nb.ORION_R_REF.to_value(u.pc)), display=False)
```

{numref}`fig:bkg-radial` shows both backgrounds against the projected King density of the cluster {cite:p}`King1962`, as a function of physical radius.
Plotted this way, the figure is independent of distance.

```{code-cell} ipython3
---
tags: [hide-input]
mystnb:
  figure:
    name: fig:bkg-radial
    caption: Radial surface-brightness profiles of the two backgrounds (left axis)
      and the projected King density of the reference cluster (right axis), against
      projected radius.
---
r = np.logspace(-2.5, np.log10(2.5), 600)  # pc
shape = CoredPowerLaw2D(r_core=KING["r_core"].to_value(u.pc), index=nb.DEFAULT_INDEX,
                        r_ref=nb.ORION_R_REF.to_value(u.pc), r_out=KING["r_tide"].to_value(u.pc))
prof = shape(r, 0)

x, xt = r / KING["r_core"].value, KING["r_tide"].value / KING["r_core"].value
king = np.clip(1 / np.sqrt(1 + x**2) - 1 / np.sqrt(1 + xt**2), 0, None) ** 2
king /= king[0]

fig, ax = plt.subplots(figsize=(7, 4), layout="compressed")
with np.errstate(divide="ignore"):
    for (name, level), color in zip((("Orion", nb.ORION_BRIGHTNESS), ("Ophiuchus", nb.OPHIUCHUS_BRIGHTNESS)), ("k", "C1")):
        m_ref = level_mag(level)
        ax.plot(r, m_ref - 2.5 * np.log10(prof), color=color, label=f"{name} background ({level[0]})")
ax.invert_yaxis(); ax.set_xscale("log")
ax.set_xlabel("projected radius [pc]"); ax.set_ylabel("surface brightness [mag/arcsec²]")
for rad, label in ((KING["r_core"], "r_core"), (nb.ORION_R_REF, "r_ref"), (KING["r_tide"], "r_tide")):
    ax.axvline(rad.to_value(u.pc), color="0.5", ls=":", lw=1)
    ax.text(rad.to_value(u.pc), 0.98, f" {label}", transform=ax.get_xaxis_transform(), va="top", fontsize=8)
ax2 = ax.twinx()
ax2.plot(r, king, color="C2", ls="--", label="King density (cluster)")
ax2.set_yscale("log"); ax2.set_ylim(1e-4, 2); ax2.set_ylabel("relative stellar surface density"); ax2.grid(False)
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="lower left")
plt.show()
```

## Building the backgrounds
Both functions take the cluster's position and King radii.
All radii must be physical lengths, and the factories reject angles.
The distance comes from the {class}`~astropy.coordinates.SkyCoord`:

```{code-cell} ipython3
pos = position(DISTANCE)
orion = nb.orion_background(pos, KING["r_core"], KING["r_tide"])
oph = nb.ophiuchus_background(pos, KING["r_core"], KING["r_tide"])
```

```{code-cell} ipython3
:tags: [hide-input]

pd.DataFrame({
    "brightness at r_ref": [" ".join(nb.ORION_BRIGHTNESS), " ".join(nb.OPHIUCHUS_BRIGHTNESS)],
    "anchor": [nb.ORION_ANCHOR, nb.OPHIUCHUS_ANCHOR],
    "spectrum file": [nb.ORION_SED, nb.OPHIUCHUS_SED],
    "total P [arcsec²]": [orion.total_flux_factor().to_value(u.arcsec**2),
                          oph.total_flux_factor().to_value(u.arcsec**2)],
}, index=["Orion", "Ophiuchus"]).style.format(
{"total P [arcsec²]": "{:.4g}"}
).set_table_attributes("class='dataframe'"
).set_table_styles(tbl_sty)
```

The two levels differ by {glue:text}`bkg_dlevel:.1f` mag/arcsec², and their anchor frames differ too:

* The Orion spectrum is intrinsic, so a foreground screen dims it once `extinction` is applied.
  The relevant screen is A_V = 0.9--1.9 with R_V = 5.5 {cite:p}`Peeters2024`; pass it via the `extinction` argument.
* The Ophiuchus level was measured as observed, and the foreground there is negligible {cite:p}`Casassus2008`.

```{code-cell} ipython3
:tags: [remove-cell]

glue("bkg_dlevel", level_mag(nb.OPHIUCHUS_BRIGHTNESS) - level_mag(nb.ORION_BRIGHTNESS), display=False)
```

(sec:bkg-distance)=
## Pushing the scene away
Because the radii are lengths, moving the scene only requires a {class}`~astropy.coordinates.SkyCoord` at the new distance.
The radii are converted to angles at construction, so both the cluster and the background need rebuilding, rather than changing the distance in place.
Surface brightness is distance-invariant: the level in mag/arcsec² stays the same, only the angular extent shrinks.
So as the stars of the cluster fade with distance, the background in front of which they must be detected does not.

```{code-cell} ipython3
grid = {"pixel_scale": 0.5 * u.arcsec / u.pixel, "width": 241, "height": 241}
distances = [414 * u.pc, 1656 * u.pc, 4140 * u.pc]
backgrounds = {d: nb.orion_background(position(d), KING["r_core"], KING["r_tide"]) for d in distances}
sources = {d: bkg.to_source(grid) for d, bkg in backgrounds.items()}
```

```{code-cell} ipython3
---
tags: [hide-input]
mystnb:
  figure:
    name: fig:bkg-distance
    caption: The Orion background in a 2′ field at three distances, as surface brightness.
      Circles mark the core (solid), reference (dashed) and tidal (dotted) radii where
      they fall inside the field.
---
omega = (grid["pixel_scale"] * u.pix).to_value(u.arcsec) ** 2
m_ref = level_mag(nb.ORION_BRIGHTNESS)
fig, axes = plt.subplots(1, 3, figsize=(11, 3.9), layout="compressed", sharey=True)
for ax, d in zip(axes, distances):
    w = sources[d].fields[0].data
    rel = w * backgrounds[d].total_flux_factor().to_value(u.arcsec**2) / omega  # = S / S_ref
    with np.errstate(divide="ignore"):
        sb = np.where(rel > 0, m_ref - 2.5 * np.log10(rel), np.nan)
    im = ax.imshow(sb, origin="lower", extent=extent_arcsec(grid), cmap="magma_r", vmin=m_ref - 2, vmax=m_ref + 5)
    for rad, ls in ((KING["r_core"], "-"), (nb.ORION_R_REF, "--"), (KING["r_tide"], ":")):
        ax.add_patch(plt.Circle((0, 0), angle(rad, d).value, fill=False, color="w", ls=ls, lw=1))
    ax.set_title(f"d = {d.to(u.kpc):.3g}", fontsize=9); ax.set_xlabel("x [arcsec]"); ax.grid(False)
axes[0].set_ylabel("y [arcsec]")
fig.colorbar(im, ax=axes, label="Ks surface brightness [mag/arcsec²]")
plt.show()
```

```{code-cell} ipython3
:tags: [hide-input]

rows = {}
for d, bkg in backgrounds.items():
    rows[f"{d.to(u.kpc):.3g}"] = {
        "r_core [″]": angle(KING["r_core"], d).value,
        "r_ref [″]": angle(nb.ORION_R_REF, d).value,
        "r_tide [″]": angle(KING["r_tide"], d).value,
        "total P [arcsec²]": bkg.total_flux_factor().to_value(u.arcsec**2),
        "fraction in field": sources[d].fields[0].data.sum(),
    }
scaling = pd.DataFrame(rows).T
p_ratio = scaling["total P [arcsec²]"].iloc[0] / scaling["total P [arcsec²]"].iloc[-1]
np.testing.assert_allclose(p_ratio, (distances[-1] / distances[0]).decompose().value ** 2, rtol=1e-6)
glue("bkg_frac_near", float(100 * scaling["fraction in field"].iloc[0]), display=False)
glue("bkg_frac_far", float(100 * scaling["fraction in field"].iloc[-1]), display=False)
scaling.style.format({"total P [arcsec²]": "{:.4g}", "fraction in field": "{:.1%}"}, precision=1
).set_table_attributes("class='dataframe'"
).set_table_styles(tbl_sty)
```

At 414 pc the tidal radius lies far outside a 2′ field, and only {glue:text}`bkg_frac_near:.1f` % of the background's total flux falls inside it.
At 4.1 kpc, {glue:text}`bkg_frac_far:.0f` % falls inside.
The total {math}`P` falls exactly as {math}`d^{-2}`, which the cell above asserts.

## Orion vs. Ophiuchus
With the same position and radii, the two backgrounds render to identical weightmaps.
Only the spectrum and its scaling differ, which is what makes the comparison controlled:

```{code-cell} ipython3
w_orion = orion.to_source(grid).fields[0].data
w_oph = oph.to_source(grid).fields[0].data
np.testing.assert_array_equal(w_orion, w_oph)
```

```{admonition} The Ophiuchus geometry is borrowed
:class: note

No radial law has been measured for the diffuse emission of L1688.
The Ophiuchus background therefore uses the Orion geometry, and only its spectrum and level are Ophiuchus-specific.
That level is an upper envelope, measured on the ρ Oph W ridge {cite:p}`Casassus2008`.
Pass `index=0` for a uniform disk out to the tidal radius instead.
```

## Using the backgrounds in a scene
The background goes under the cluster with the same position and the same King radii.
Building the cluster itself is shown elsewhere; schematically:

```python
pos = SkyCoord(84.7 * u.deg, -6.9 * u.deg, distance=414 * u.pc)
king = {"n_stars": 1000, "r_core": 0.13 * u.pc, "r_tide": 2.0 * u.pc}
cluster = ZeroAgeCluster(pos, IMFPopulation, {"n_stars": 1000},
                         KingProfileMorphology, king, rng_seed=1138)
background = nb.orion_background(pos, king["r_core"], king["r_tide"])
```

Two constraints and one current limitation apply:

* **The tidal radius must exceed** {math}`r_\mathrm{ref}`. Otherwise the surface brightness would be anchored where the profile is zero, and the target raises an error.
* **Moving the scene means rebuilding both targets** with a new {class}`~astropy.coordinates.SkyCoord`, as shown in {ref}`sec:bkg-distance`.
* **The profile is always centred on the field.** Extended targets do not yet honour `position` or `offset` in rendering.

## References
```{bibliography}
:filter: docname in docnames
```
