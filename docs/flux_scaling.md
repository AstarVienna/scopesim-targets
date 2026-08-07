---
file_format: mystnb
kernelspec:
  name: python3
numbering:
  heading_1: true
  heading_2: true
---

# Spectra and weights: How flux is stored in the Source object
Every non-cube target stores the same two things: a set of dimensionless _weights_ and a _spectrum_ scaled to a specific flux.
The total flux stored in the source is always the same product of those two:

```{math}
:label: fdef
F \;=\; W F_\mathrm{band}, \qquad F_\mathrm{band} \;=\; \int S(\lambda)\, T(\lambda)\; \mathrm{d}\lambda ,
```

where {math}`S(\lambda)` is the stored spectrum, {math}`T(\lambda)` the bandpass ("filter") throughput,
and {math}`W` is the weight (a single number for a point source, the sum of a 2D map for an extended source).
This page explains the concept behind this and discusses various special cases and consequences.
It assumes you've read the target brightness syntax {doc}`defining_brightness`.
The following pages {doc}`flux_construction` and {doc}`flux_engine` build on the basic concepts discussed here.

## Setup
```{code-cell} ipython3
:tags: [hide-output]

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import pandas as pd
import seaborn as sns
from synphot import Observation
from spextra import Spextrum, Passband

from scopesim_targets.point_source import Star
from scopesim_targets.extended_source import Box, Flat, Gaussian, Sersic
from scopesim_targets.target import FILTER_SYSTEM

flat_spec = Spextrum.flat_spectrum(
    5*u.ABmag, waves=np.linspace(3000, 35000, 6000)*u.AA
)

# 3.55 Jy is a realistic medium-bright star and lands at a round ~1000 photons
# per second per cm2 in V, which keeps every number below legible.
TOTAL_FLUX = 3.55*u.Jy
```

The following cell defines some helper functions used in this page,
but not otherwise required to use ScopeSim and ScopeSim-Targets.
Hidden by default to reduce clutter, expand if you want to see them:

```{code-cell} ipython3
:tags: [hide-cell]
sns.set_theme(style="darkgrid", palette="Set2")
bbox_flx = {"boxstyle": "round", "fc": "white", "ec": "0.5", "alpha": 0.85}

def band_flux(spectrum, band):
    return Observation(spectrum, band).integrate()

def total_flux(source, band):
    field = source.fields[0]  # assuming single source field in these examples

    if hasattr(field, "data"):
        # 2D image (Box, Gaussian, ...)
        source_weight = float(np.nansum(field.data))
    else:
        # point-source table
        source_weight = float(np.sum(field.field["weight"]))

    return source_weight * band_flux(field.spectrum, band)

def extent_arcsec(grid):
    scale = grid["pixel_scale"].to_value(u.arcsec / u.pixel)
    hw = scale * grid["width"] / 2
    hh = scale * grid["height"] / 2
    return (-hw, hw, -hh, hh)

def add_pixel_grid(ax, grid, color="white", lw=.8, ls=":", alpha=.6,
                   force_major_subset=True):
    scale = grid["pixel_scale"].to_value(u.arcsec/u.pixel)
    width = int(round(grid["width"]))
    height = int(round(grid["height"]))
    x_edges = (np.arange(width + 1) - width / 2) * scale
    y_edges = (np.arange(height + 1) - height / 2) * scale
    ax.set_xticks(x_edges, minor=True)
    ax.set_yticks(y_edges, minor=True)
    if force_major_subset:
        ax.set_xticks(x_edges[::max(1, round(len(x_edges) / 6))])
        ax.set_yticks(y_edges[::max(1, round(len(y_edges) / 6))])
    ax.grid(which="both" if force_major_subset else "minor",
            color=color, linewidth=lw, linestyle=ls, alpha=alpha)
    ax.tick_params(which="minor", length=0)
```

## Flux of a spectrum in a filter
The simplest target is a single `Star`.
Building one and pulling it apart shows exactly what any target stores:
a template spectrum, and a single `weight` {math}`w` that scales it to the requested brightness.
Their product is the star's flux-calibrated SED:

```{code-cell} ipython3
star = Star(
    position=(0, 0),
    spectrum=flat_spec,
    brightness=("V", TOTAL_FLUX),
)
src_pnt = star.to_source()

template = src_pnt.fields[0].spectrum          # the shared, un-scaled template
weight = float(src_pnt.fields[0].field["weight"][0])  # the single scale factor
scaled_spectrum = weight * template           # flux-calibrated to 3.55 Jy in V
```

```{code-cell} ipython3
band = Passband(f"{FILTER_SYSTEM.name}/V")
waves = np.linspace(3800, 8200, 2000)*u.AA

S = scaled_spectrum(waves)  # scaled spectrum evaluated on wavelengths
T = band(waves)             # throughput evaluated on wavelengths
FV = Observation(scaled_spectrum, band).integrate()  # Integral under curve
```

```{tip}
:class: margin
We hide all the boring plotting and formatting code by default.
Feel free to expand the cell to see it!
```

```{code-cell} ipython3
:tags: [hide-input]
:label: fig:spectrumflux

fig, ax1 = plt.subplots(figsize=(7, 4))
ax1.plot(waves, S, color="C0", lw=1.4, label=r"spectrum $S(\lambda)$")
ax1.fill_between(
    waves.value, 0, (S * T).value, color="C2", alpha=0.7,
    label=r"transmitted $S(\lambda)\,T(\lambda)$",
)
ax1.set_xlabel(r"wavelength [$\mathrm{\AA}$]")
ax1.set_ylabel(f"flux density [{S.unit.to_string('latex_inline')}]")
ax1.set_xlim((4000, 8000))
ax1.set_ylim((0, 1.4))

ax1.text(
    6500, 0.85, "flat spectrum in ABmag", color="C0", fontsize=10,
    rotation=-9, rotation_mode="anchor", va="bottom",
)

ax2 = ax1.twinx()
ax2.plot(waves, T, color="C1", ls="--", lw=1.5, label=r"passband $T(\lambda)$")
ax2.set_ylabel("throughput")
ax2.set_ylim(0, 1)
ax2.grid(False)

(line_s, line_f), (label_s, label_f) = ax1.get_legend_handles_labels()
(line_t,), (label_t,) = ax2.get_legend_handles_labels()
ax1.legend(
    [line_s, line_t, line_f],
    [label_s, label_t, label_f],
    loc="upper right", fontsize=9,
)

ax1.text(
    5500, 0.25,
    f"F(V) = {FV.to_string(format='latex', formatter='.0f')}",
    ha="center", fontsize=11,
    bbox=bbox_flx,
)
plt.show()
```

The dashed orange line show the transmission curve of the filter {math}`T(\lambda)`, with its y-axis on the right side, between 0 and 1.
Initially, the spectrum is normalized internally via the `effstim` in PHOTLAM units, given a flux in a filter ("bandpass").
Once it reaches the ScopeSim engine (specifically the `FieldOfView` class),
it is integrated over the selected filter (not necessarily the same that it was normalized to),
assuming the simulation is running in imaging mode -- spectroscopy is somewhat more complex.
The result of that integration {math}`F_\mathrm{band} \,=\, F(V)` corresponds directly to the blue shaded area in the plot.
In terms of units, integrating over the wavelength eliminated the {math}`\lambda^{-1}` dimension, so we're left with photons per second per cm².
That remaining area dimension is finally eliminated by multiplying with the telescope's collecting area,
but we will keep it around in these examples to stay instrument-agnostic.
The integrated flux value is of course not the same as the {math}`3.55\;\mathrm{Jy}` we specified in the input,
because by folding in the wavelength, we implicitly converted it using the energy of each photon at any given wavelength (kann man das so sagen??).

```{admonition} Note on the "flat" spectrum shown here
:class: note
"Flat in AB magnitude" means flat in {math}`f_\nu`; converted to the per-wavelength photon density that `integrate()` sums, a constant {math}`f_\nu` becomes a falling {math}`\propto 1/\lambda` power law.
That conversion -- between magnitude systems, {math}`f_\nu`, {math}`f_\lambda` and PHOTLAM -- is worth reading once in
[synphot's units documentation](https://synphot.readthedocs.io/en/latest/synphot/units.html#counts-and-magnitudes);
using a spectrum that is deliberately _not_ flat in PHOTLAM means a wrong reference-wavelength scaling anywhere downstream would actually show up rather than hide in a constant.
```

## Point sources are scaled via `weight`
For a point source, the `weight` is a simple scalar, so {math}`W = w` and eq. {eq}`fdef` becomes:

```{math}
:label: pnttotal
F \;=\; w \cdot F_\mathrm{band}(\text{template})
```

Using the example from above, we can introspect those numbers:
```{code-cell} ipython3
round(weight, 3)
```
```{code-cell} ipython3
band_flux(template, band).round()
```
```{code-cell} ipython3
(weight * band_flux(template, band)).round()
```
```{code-cell} ipython3
FV.round()
```

and finally we can see how we get the originally specified {math}`3.55\;\mathrm{Jy}` using the `.effstim` method.

```{code-cell} ipython3
Observation(scaled_spectrum, band).effstim(u.Jy).round(4)
```

Why keep the scale in a separate `weight` instead of just storing the scaled spectrum?
Because point sources can come in crowds, i.e. star fields, clusters, etc..
A `StarField` of {math}`10^4` stars drawn from 30 spectral templates stores **30** spectra and {math}`10^4` weights --
each row points at a shared template and carries only its own scalar.
Folding the scale into the spectrum would force {math}`10^4` distinct spectra, which explodes the required memory to store that source object.
So for point sources the spectrum stays the shared template and the `weight` carries the flux.

## Extended sources need a weightmap
Unlike a point source, an extended source has a spatial shape, meaning we need to store more information.
To avoid the memory and performance cost of a full 3D data cube, ScopeSim uses a "trick" to combine a 2D image with a 1D spectrum.
The consequence of this is that the spectrum at every point of such an extended source is identical,
meaning physical effects like e.g. rotation leading to line-of-sight velocity shifts or
different chemical signatures or temparture in different regions of the same object cannot be represented with such a source description.
This limitation is usually acceptable for photometric applications, whereas spectroscopic simulations might need to
pay the price of a full 3D source description to accurately simulate their science case.

A practical consequence of having a single spectrum per 2D image is that the "deduplication" argument for point source clusters outlined in the previous section no longer applies.
As a result, it is more useful in this case to _scale the spectrum itself_ to the specified flux directly.
The role of the scalar **weight** {math}`w` for point sources is now taken by the 2D image, which now acts as a **weightmap** {math}`w_{ij}`.
Instead of scaling a shared spectrum to the flux spcified for a point source, it now decides _how much_ of the total flux in the scaled spectrum ends up _where_.
The sum of this weightmap scales the total flux that enters the simulated optical train.
Because we now scale the spectrum itself, we would expect {math}`\sum_{ij} w_{ij} = 1` to always hold.
We'll see in the following section why that's **not** always the case, but if we assume this for a second, eq. {eq}`fdef` becomes:

```{math}
:label: wmaptotal
F \;=\; W \cdot F_\mathrm{band}(\text{scaled spectrum}), \qquad W \;=\; \sum_{ij} w_{ij} \approx 1
```

```{TODO}
Mention concept of "rasterized point sources" and show comparison of point source grid to 2+1D source.
Might be better on the "engine" page, showing the actual ImagePlane, where a point source has to rasterize.
```

It's worth mentioning at this point that multiple "2+1D" sources can be combined into an "overlay" to create a simulation scene, each with a separate spectrum.
This can be used to map different components of the same scene, where each spectrum's flux is distributed at different location in the 2D image.
An example of this that is used in a few ScopeSim notebooks is a spiral galaxy with two stellar populations,
each with an average spectrum and a weightmap of where in the galaxy it is distributed.
The resulting simulated observation is essentially a weighted overlay of both populations' light.

```{TODO}
Link the "spiral two component" example from Templates here.
```

### Weightmap and flux
The example below creates a simple rectangular `Box` target to visualize the concept of extended sources.
We'll use the same total flux and the same simple spectrum as with the point source example above, so all the numbers should be identical.
Since the target itself is a parametrized description, we need to define a `grid` to evaluate it on -- more about that in [](#sec:raster) below.
The `Box` is immediately converted to a {py:class}`~scopesim.Source` object.
The conversion method (also referred to as "concretization" elsewhere in the documentation) take the `grid` as an argument.

```{code-cell} ipython3
grid = {
    "pixel_scale": 1*u.arcsec/u.pixel,
    "width": 11,   # pixels in x-direction in the field-of-view
    "height": 11,  # pixels in y-direction in the field-of-view
}

src_box = Box(
    spectrum=flat_spec,
    brightness=("V", TOTAL_FLUX),
    params={"x_width": 5, "y_width": 3},
).to_source(grid)
```

We can now introspect the data stored in the `Source` object, both the 2D image and the 1D spectrum (or rather the flux stored there).
In the notation used above, these are {math}`w_{ij}` and {math}`F_\mathrm{band} \,=\, F(V)`, respectively.
Combining the two, we get {math}`F` in units of {math}`\mathrm{ph\,s^{-1}\,cm^{-2}}`,
which is what ScopeSim converts to further downstream -- more about that later in {doc}`flux_engine`.

```{code-cell} ipython3
wmap = src_box.fields[0].data  # dimensionless 2D weightmap, sums to ~1
FV_box = band_flux(src_box.fields[0].spectrum, band)
image = wmap * FV_box          # image in ph / s / cm2  (per pixel)
```

The plots below are arguably not very spectacular, but serve the purpose of showing the rasterization,
as well as the unit conversion according to eq. {eq}`wmaptotal` (stricly speaking the result is {math}`F_{ij}` rather than {math}`F`),
with the {math}`F(V)` taken from the spectrum, as shown in fig. [1](#fig:spectrumflux).

```{code-cell} ipython3
:tags: [hide-input]
fig, axes = plt.subplots(
    1, 3, figsize=(9, 3.8), gridspec_kw={"width_ratios": [1, 0.6, 1]},
    layout="compressed",
)

im0 = axes[0].imshow(
    wmap, extent=extent_arcsec(grid),
    origin="lower", cmap="magma",
)
add_pixel_grid(axes[0], grid, force_major_subset=False)
axes[0].set_aspect("equal")
axes[0].grid(False, which="major")
axes[0].set_title("weight map\n(dimensionless)")

axes[1].axis("off")
axes[1].text(0.5, 0.7, "\u00d7", ha="center", va="center", fontsize=20)
axes[1].text(0.5, 0.3, "=", ha="center", va="center", fontsize=20)
axes[1].text(
    0.5, 0.5,
    FV_box.to_string(format="latex", formatter=".0f"),
    ha="center", va="center", fontsize=14,
    bbox=bbox_flx | {"fc": "C2"},
)

im2 = axes[2].imshow(
    image.value, extent=extent_arcsec(grid),
    origin="lower", cmap="plasma",
)
add_pixel_grid(axes[2], grid, force_major_subset=False)
axes[2].set_aspect("equal")
axes[2].grid(False, which="major")
axes[2].set_title(f"flux image\n[{image.unit.to_string('latex_inline')}]")

fig.subplots_adjust(wspace=0.1, right=0.9)
plt.show()
```

Comparing this to the point source of identical total specified brightness above,
we can see that the total flux in the spectra is the same and matches the sum of the flux that ends up in the image.
This shows us that all the scalings are consistent, both the "total flux via float weight" of the point source,
and the "total flux in spectrum, distributed over image" of the extended one.
In other words, the product `weight(s) x band_flux(spectrum)` is constant, even if the scaling shifts.

```{code-cell} ipython3
total_flux(src_box, band).round()
```
```{code-cell} ipython3
total_flux(src_pnt, band).round()
```
```{code-cell} ipython3
image.sum().round()
```

## Brightness profiles
### Exact vs. cutoff profiles
The spatial structure of the weightmap is given by its brightness profile {math}`p(x, y)`, implemented in the abstract `BrightnessProfile` base class.
The exact formula behind this depends on the concrete profile.
In the exapmle above, we've seen the simplest case: a rectangular `Box`.
While obviously not very realistic for an astrophysical object (although "boxy" galaxies exist), it is useful as a simply test case to visualize the general concept.
It is also easy to see why the naïve approach for the total of the weightmap doesn't work in practice:
If we strictly normalized the weightmap to {math}`\sum_{ij} w_{ij} = 1` in all cases,
the spatial extent of the image would need to be large enough to encompass the whole profile,
otherwise the region that actually ends up _within_ the image might end up too bright,
because the cut-off flux is artificially reintroduced into the field-of-view.
For the simple `Box` case, the trival solution would be to always force the image to be large enough,
although this could lead to performance issues for a very large box on a relatively fine grid,
where only a small portion of the box would actually be inside the field-of-view.

The conceptually much more complicated issue is that most realistic profiles will not have a finite spatial extent, the way the `Box` does.
Enter a simple `Gaussian` brightness profile: Mathematically only barely more complex than the `Box`, it technically carries on infinitely.
Of course, solutions such as cutting it off at a very reasonable {math}`5\sigma` might solve this,
but not only is that rather arbitrary, it also leaves the same potential performance issue as described above.
And much worse, many even more realistic profiles such as the well-known `Sersic` have no such solutions.

Thus it is much more useful to normalize the weightmap to the **intrinsic total** of the brightness profile instead.
Consider for example the analytic integrals of the two simple cases we just discussed, `Box` and `Gaussian`,
in both cases we ignore the offsets {math}`x_0` and {math}`y_0` and any angle {math}`\theta` (or pretend we did a coordinate transformation...):

#### Rectangular `Box` profile
The `Box` is defined as constant ({math}`A`) within the rectangle, and zero outside of it.

```{math}
:label: pbox
p(x, y)\:=\:\left\{\begin{array}{rl}
A : & \frac{w_x}{2} \leq x \leq \frac{w_x}{2} \text{ and } \frac{w_y}{2} \leq y \leq \frac{w_y}{2} \\
0 : & \text{else}
\end{array}\right.
```

where {math}`w_x` and {math}`w_y` are the width and height, respectively (using Astropy's naming conventions).
The integral for that is trivial ("length × width"):

```{math}
P\:=\:\iint p(x, y)\,\mathrm{d}x\,\mathrm{d}y \;=\; A\,w_x\,w_y
```

Now because we established previously that the flux of the source is entirely stored in its spectrum,
the constant {math}`A` devides out, along with its unit, or in other words we can set {math}`A=1`.
This simplifies the integral even more to just {math}`P=w_x\,w_y`.

#### 2D `Gaussian` profile
With the same simplifications, the integral of the `Gaussian` brightness profile becomes:

```{math}
P\:=\:\iint A\,e^{-\frac{x^2}{2\sigma_x^2}-\frac{y^2}{2\sigma_y^2}}\,\mathrm{d}x\,\mathrm{d}y
 \:=\: 2\pi A\,\sigma_x\sigma_y \;\;\xrightarrow{A=1}\;\; 2\pi\,\sigma_x\sigma_y
```

Note the physical dimensions here: Because both {math}`w_x` and {math}`w_y` -- or {math}`\sigma_x` and {math}`\sigma_y` respectively -- are lengths
(or rather angles, but usually small, so we can pretend the sky is linear to make things simpler),
the integral {math}`P` becomes an **area** (again, strictly speaking a solid angle),
equivalent to the area a uniform profile of the same total brightness would cover.
This will become important a bit later.

(sec:raster)=
### Rasterization
While the intrinsic brightness profile is {math}`p(x, y)`, the simulations work on a finite grid instead of continous coordinates.
We will call the angular pixel scale ("length") {math}`\varpi` and the solid angle ("area") covered by one pixel {math}`\Omega=\varpi^2`.
The profile can thus be written in terms of these rasterized coordinates {math}`p(x_i, y_j)`.
Remember this is now a dimensionless scale that simply defines the spatial structure of the profile.
In the simple normalized `Box` case, we can assume {math}`p(x_i, y_j) \equiv 1` everywhere _inside_ the box.
The real implementation is slightly more complicated to account for cases where {math}`w_x` and {math}`w_y` are not integer multiples of {math}`\varpi`,
but we can again simplify this.
Note that this limitation is of course no concern for the `Gaussian`, or any other smooth profile.

The aforementioned weights {math}`w_{ij}` are then just the pixel values of the rasterized image.
To scale them correctly to the instrisic total {math}`P` that the spectrum's flux corresponds to,
we have to scale each pixel by its contribution to the total.
This can be expressed in terms of "area":

```{math}
w_{ij} \;=\; \frac{\Omega}{P} \: p(x_i, y_j)
```

```{TODO}
{math}`A_\mathrm{eff}` here??
```

For a finite profile (e.g. `Box`) fully within the field-of-view,
this will result in {math}`\sum_{ij} w_{ij} = 1.0`, ignoring any floating point imprecisions.
For a technically infinite profile (e.g. `Gaussian`), which will always only be partially inside a finite grid,
this must always result in {math}`\sum_{ij} w_{ij} < 1.0`, although the practical difference can be minimal.

In cases where a finite profile overflows the field-of-view, the expected analytic total can still be calculated.
For the `Box`, this can be written as {math}`\text{area}(\text{box} \cap \text{FOV}) / (w_x\,w_y)`.
For the `Gaussian` the integral over a centered window of width {math}`W` and height {math}`H` the map sums to
{math}`\operatorname{erf}\left(\frac{\sqrt{2}\,W}{4\,\sigma_x}\right)\operatorname{erf}\left(\frac{\sqrt{2}\,H}{4\,\sigma_y}\right) < 1`.
We will use these formulae to cross-check what the code does in the examples below.
Since the pixel weights are already scaled to their contriution to the _total_ integral {math}`P` anyway,
the scaling is unaffected by the extent of {math}`p` relative to the field-of-view,
the only practical difference being {math}`\sum_{ij} w_{ij}`.

```{tip}
The spectrum stores the total flux, which always corresponds to the _intrinsic total_.
An image which cuts off parts of the brightness profile still carries the _same_ spectrum.
The sum of the weight map decides how much of the spectrum's defined flux falls into the field-of-view.
```

The example below shows a `Box` and a `Gaussian` with identical total brightness (thus storing the same spectrum),
scaled to (almost, rounded for the `Gaussian`) the same intrinsic {math}`P` "area".
Since the `Box` fits comforably withing the field-of-view, its weightmap sums to exactly one.
The `Gaussian`, being technically infinite, loses about 6 % of its flux outside the window.

```{TODO}
Check analytic Gauss solution and compare to real grid sum!
```


```{code-cell} ipython3
box_grid = {"pixel_scale": 0.5*u.arcsec/u.pixel, "width": 15, "height": 15}

wmap_box = Box(
    spectrum=flat_spec,
    brightness=("V", TOTAL_FLUX),
    params={"x_width": 5, "y_width": 3},
).to_source(box_grid).fields[0].data

wmap_gauss = Gaussian(
    spectrum=flat_spec,
    brightness=("V", TOTAL_FLUX),
    params={"x_stddev": 2, "y_stddev": 1.2},
).to_source(box_grid).fields[0].data
```

```{code-cell} ipython3
:tags: [hide-input]
vmax = max(wmap_box.max(), wmap_gauss.max())
extent = extent_arcsec(box_grid)

fig, axes = plt.subplots(1, 2, figsize=(9, 4), sharey=True, layout="compressed")
for ax, wmap, title in zip(
    axes, (wmap_box, wmap_gauss), ("Box (exact edges)", "Gaussian (truncated wings)")
):
    im = ax.imshow(wmap, origin="lower", extent=extent, cmap="mako", vmax=vmax)
    add_pixel_grid(ax, box_grid, force_major_subset=False)
    ax.set_aspect("equal")
    ax.grid(False, which="major")
    ax.set_xticks(np.arange(*np.ceil(extent)[:2]))
    ax.set_title(f"{title}\n$\\sum w_{{ij}}$ = {wmap.sum():.4f}")
fig.colorbar(im, ax=axes[1])
plt.show()
```

### Non-integrable profiles
Not all brightness profiles have a closed-form analytic total.
A notebale case relevant for e.g. calibration sources in ScopeSim is the `Flat`, which has a constant uniform brightness.
Again assuming a normalized amplitude (the constant value) of {math}`C=1`, the integral for its total diverges:

```{math}
P\:=\:\iint_{-\infty}^{\infty} C\,\mathrm{d}x\,\mathrm{d}y \;\rightarrow\;\infty
```

We can solve this by using the field-of-view that the profile is rasterized onto as the limits for our integration,
in other words the rendered field-of-view closes the otherwise open integral.
This is a very practical solution, because any flux outside of it is lost to the simulation anyway.
Using the same nomenclature as above and calling the result an "effective area" {math}`A_\mathrm{eff}` instead of {math}`P`,
a window of width {math}`W` and height {math}`H` results in integration limits of:

```{math}
a_x = -\frac{W}{2} \quad b_x = \frac{W}{2}\\
a_y = -\frac{H}{2} \quad b_y = \frac{H}{2}
```

so the integral above becomes

```{math}
A_\mathrm{eff}\:=\:\int_{-H/2}^{H/2}\int_{-W/2}^{W/2} C\,\mathrm{d}x\,\mathrm{d}y
```

which we can simplify, because we're just integrating a constant

```{math}
A_\mathrm{eff}\:=\:\int_{0}^{H}\int_{0}^{W} C\,\mathrm{d}x\,\mathrm{d}y \:=\: W\,H
```

```{tip}
:class: margin
A `Flat` field is identical to a `Box` target with dimensions matching the FOV!
```

If you've been following along, you will notice that this is basically identical to a `Box` whose dimensions correspond exactly to the field-of-view!
We can also write the same total using the detector pixel scale and number of pixels {math}`N` on the detector instead of its on-sky dimensions:

```{math}
A_\mathrm{eff} \:=\: \Omega\,N \:=\: A_\mathrm{FOV}.
```

This makes it a bit easier to conceptualize how the "total flux" of a `Flat` is calculated,
i.e. what the brightness stored in the spectrum corresponds to, which is simply what falls into the field-of-view.
It also means the weightmap image of a `Flat` will always sum to 1.
The consequence of this is that the total flux that the spectrum is scaled to depends directly on the field-of-view,
so the same `Flat` target realized for different parameters will produce a differently scaled spectrum.
But that is actually physically correct, because a larger field-of-view will catch more light of an identical flatfield.

### Surface brightness
A detail we ignored so far while deriving the `Flat` concept is how the brightness is initially specified.
Stating the total intrinsic flux -- as we did so far for integrable profiles -- would be inherently tied to the field-of-view dimensions and pixel scale,
which is probelmatic if they're not exactly known at the time of defining the target.
Also, the same target definition could not be reused for different instruments, or even different modes on the same instrument in many cases.
What makes a lot more sense for these kind of targets is to define a **surface brightness** instead.
This is also how well-resolved extended astronomical sources are commonly described in general, and we'll get to how this works even for the integrable cases from before.

Because the total flux {math}`F` is normalized to whatever {math}`A_\mathrm{eff}` is,
a given surface brightness {math}`S` can be converted to the total flux {math}`F` via the formula below.
Dimensionally, this works nicely because the "per solid angle" part of {math}`S` cancels with the unit of {math}`A_\mathrm{eff}`, leaving just a flux unit.

```{math}
F\:=\: S \cdot A_\mathrm{eff}
```

```{note}
For better or worse, astromomers like to give surface brightnesses in units of magnitudes per on-sky "area", e.g. {math}`\frac{\mathrm{mag}}{\mathrm{arcsec}^2}`.
We should pause for a moment to see how that's dealt with:
Because Astropy cannot directly parse such units if they come in a system other than plain "mag", e.g. ABmag,
such quantities can only be given as a string, e.g. `"20 mag(AB) / arcsec2"`.
These are converted to linear units before they even hit the `Target` constructor, using this formula:

\begin{equation}
m = m_\mathrm{SB} - 2.5 \log_{10} \frac{A_\mathrm{eff}}{\Omega_m}
\end{equation}

where {math}`\Omega_m` the solid-angle unit the surface brightness was given in,
e.g. {math}`\mathrm{arcsec}^2` for {math}`\frac{\mathrm{mag}}{\mathrm{arcsec}^2}`.
```

Of course, the closed-form targets we discussed previously also accept a surface brightness instead of a total flux.
We can confirm that this gives the expected results with a few practical, easy-to-calculte examples.

#### Flat scaling on different grids
The example below shows a `Flat` target compared to two `Box` targets, once with the same surface brightness specified,
and once with a total brightness caclulated analytically for the dimensions of the square box.
These targets are then evaluated onto two grids, one exactly matching the dimensions of the boxes,
and one with a factor of 2 finer pixel scale, resulting in a factor of 4 (or rather 0.25) in terms of area (both per pixel and for the whole field-of-view).
All three targets share the same simple spectrum as we used throughout this page.
The table below shows the calculated weightmap sums, flux-in-band for the scaled spectra, and total flux inside the field-of-view.
Again we've hidden the code for the calculations by default, expand to see the details.

```{code-cell} ipython3
targets = {
    "Flat": Flat(
        spectrum=flat_spec,
        brightness=("V", 44.4*u.mJy/u.arcsec**2),
    ),
    "Box SB": Box(
        spectrum=flat_spec,
        brightness=("V", 44.4*u.mJy/u.arcsec**2),
        params={"x_width": 4, "y_width": 4},
    ),
    "Box total": Box(
        spectrum=flat_spec,
        brightness=("V", 710*u.mJy),
        params={"x_width": 4, "y_width": 4},
    ),
}

grids = [
    {"pixel_scale":  0.5*u.arcsec/u.pixel, "width": 8, "height": 8},  # coarse
    {"pixel_scale": 0.25*u.arcsec/u.pixel, "width": 8, "height": 8},  # fine
]
```

```{code-cell} ipython3
:tags: [hide-input]

rows = {}

for grid in grids:
    wmap_sums = {
        (grid["pixel_scale"], name): tg.to_source(grid).fields[0].data.sum()
        for name, tg in targets.items()
    }
    spec_fluxes = {
        (grid["pixel_scale"], name): band_flux(tg.to_source(grid).fields[0].spectrum, band)
        for name, tg in targets.items()
    }
    totals = {
        key: wmap_sums[key] * spec_fluxes[key]
        for key in wmap_sums
    }

    rows.update({
        key: {
            "Weight sum": wmap_sums[key],
            "Spectrum flux": spec_fluxes[key],
            "Total (w × f)": totals[key],
        }
        for key in wmap_sums
    })

df = pd.DataFrame.from_dict(rows, orient="index")
df.index = pd.MultiIndex.from_tuples(df.index, names=["pixel scale", "target"])
df.sort_index(ascending=False).style.format({
    "Weight sum": "{:.3f}",
    "Spectrum flux": "{:.0f}",
    "Total (w × f)": "{:.0f}",
}).set_table_attributes("class='dataframe align-left'").set_properties(**{
    "text-align": "right",
    "padding": "0.2em 0.8em",
})
```

We can see that the two `Box` targets behave identically on both pixel grids.
Both scale their spectrum to the same intrinsic total {math}`A_\mathrm{eff}`,
and in both cases only the central quarter (4 × 4 out of 16 × 16) of the weightmap ends up inside the field-of-view,
correctly resulting in {math}`\sum_{ij} w_{ij} = 0.25` on the finer grid.
The `Flat` on the other hand always sums to {math}`\sum_{ij} w_{ij} = 1.0` and rescales the spectrum instead,
and we can see that's where the factor of 0.25 ends up instead for the finer grid.
Importantly, the totals as per eq. {eq}`wmaptotal` all agree in the end, producing the same factor of 0.25 between the two grids.

#### Surface brightness in mag
To ensure that surface brightnesses given in magnitudes are scaled correctly, we can come up with a simple example.
Here we create a `Box` with the dimensions of {math}`1\,\times\,1^{\prime\prime}` carrying a surface brightness of {math}`10\;\mathrm{\frac{mag\left(AB\right)}{arcsec^{2}}}`.
Because {math}`A_\mathrm{eff} = 1\;\mathrm{arcsec}^2`, it should produce the exact same flux as a point source with {math}`10\;\mathrm{mag\left(AB\right)}`.

```{code-cell} ipython3
grid = {"pixel_scale":  0.1*u.arcsec/u.pixel, "width": 45, "height": 45}

src_box = Box(
    spectrum=flat_spec,
    brightness=("V", "10 mag(AB) / arcsec2"),  # mag SB has to be a string
    params={"x_width": 1, "y_width": 1},  # 1 arcsec2
).to_source(grid)
src_pnt = Star(
    position=(0, 0),
    spectrum=flat_spec,
    brightness=("V", "10 mag(AB)"),  # total mag can be a string
).to_source()

f_box = total_flux(src_box, band)
f_pnt = total_flux(src_pnt, band)

print(f"1 arcsec2 box @ 10 mag(AB)/arcsec2:  {f_box:.3f}")
print(f"point source  @ 10 mag(AB)        :  {f_pnt:.3f}")
print(f"\nratio:  {f_box/f_pnt:.8f}")
```

Next, we can double the side length of the `Box`, increasing its area by a factor of 4.
If we keep the surface brightness the same, we would expect the same factor 4 in the total flux.
This corresponds to {math}`\sim 1.5\;\mathrm{mag}` difference, which we will use for the point soruce.

```{code-cell} ipython3
src_box = Box(
    spectrum=flat_spec,
    brightness=("V", "10 mag(AB) / arcsec2"),
    params={"x_width": 2, "y_width": 2},  # 4 arcsec2
).to_source(grid)
src_pnt = Star(
    position=(0, 0),
    spectrum=flat_spec,
    brightness=("V", 8.5*u.ABmag),  # total mag can also be a quantity
).to_source()

f_box = total_flux(src_box, band)
f_pnt = total_flux(src_pnt, band)

print(f"1 arcsec2 box @ 10 mag(AB)/arcsec2:  {f_box:.3f}")
print(f"point source  @ 10 mag(AB)        :  {f_pnt:.3f}")
print(f"\nratio:  {f_box/f_pnt:.8f}")
```

The {math}`\sim 0.5\,\%` deviation is easily explained because the {math}`8.5\;\mathrm{mag}` were rounded.
So the result is consistent.

### Totals and reference points of Profiles
The table below provides an overview of the currently implemented `BrightnessProfile` subclasses and their analytic integral (if it exists),
as well as the reference point convention for a stated surface brightness.
It can be read as "how total / surface brightness is translated" for a given profile.

| Profile    | Analytic integral {math}`P`                                | Reference point |
|------------|------------------------------------------------------------|-----------------|
| `Box`      | {math}`w_x\,w_y`                                           | _uniform_       |
| `Disk`     | {math}`\pi\,R_0^2`                                         | _uniform_       |
| `Ring`     | {math}`\pi\,\left(r_\mathrm{out}^2-r_\mathrm{in}^2\right)` | _uniform_       |
| `Gaussian` | {math}`2\pi\,\sigma_x\sigma_y`                             | peak value
| `Sersic`   | {math}`` | at {math}`R_\mathrm{eff}` |
| `Flat`     | no                                                         | _uniform_       |

The example below illustrates these different reference points for a `Sersic`, a `Gaussian`, and a `Box` with identical surface brightnesses and spectra.
The {math}`R_\mathrm{eff}` of the `Sersic` was chosen to match the {math}`\sigma_x` of the `Gaussian`
and the height of the `Box` matched its {math}`\sigma_y`.

```{code-cell} ipython3
surface_brightnes = ("V", 44.4*u.mJy/u.arcsec**2)  # same for all 3 targets
fine_grid = {"pixel_scale":  0.2*u.arcsec/u.pixel, "width": 45, "height": 45}

src_sers = Sersic(
    spectrum=flat_spec,
    brightness=surface_brightnes,
    params={
        "r_eff": 2*u.arcsec,
        "n": 1,
        "ellip": 0.4,
    },
).to_source(fine_grid)

src_gauss = Gaussian(
    spectrum=flat_spec,
    brightness=surface_brightnes,
    params={"x_stddev": 2, "y_stddev": 1.2},
).to_source(fine_grid)

src_box = Box(
    spectrum=flat_spec,
    brightness=surface_brightnes,
    params={"x_width": 5, "y_width": 2.4},
).to_source(fine_grid)

flux_sers = band_flux(src_sers.fields[0].spectrum, band)
flux_gauss = band_flux(src_gauss.fields[0].spectrum, band)
flux_box = band_flux(src_box.fields[0].spectrum, band)
```

```{code-cell} ipython3
:tags: [hide-input]
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.5, 4), sharey=True, layout="compressed")
xax = np.linspace(-4.4, 4.4, 45)

ax1.axvline(-2, c="k", ls=":")
ax1.axvline(2, c="k", ls=":")
ax1.axhline(0.5, c="C3", ls="-.")
ax1.plot(xax, src_sers.fields[0].data[22, :] * flux_sers, label="Sersic")
ax1.plot(xax, src_gauss.fields[0].data[22, :] * flux_gauss, label="Gaussian")
ax1.plot(xax, src_box.fields[0].data[22, :] * flux_box, label="Box")
ax1.set_ylim(top=1.1)
ax1.legend(loc="upper right")
ax1.set_title("Cut in x-direction")
ax1.set_xlabel("arcsec")
ax1.set_ylabel("Flux [$\\mathrm{ph\\,s^{-1}\\,cm^{-2}}$]")

ax2.axvline(-1.2, c="k", ls=":")
ax2.axvline(1.2, c="k", ls=":")
ax2.axhline(0.5, c="C3", ls="-.")
ax2.plot(xax, src_sers.fields[0].data[:, 22] * flux_sers, label="Sersic")
ax2.plot(xax, src_gauss.fields[0].data[:, 22] * flux_gauss, label="Gaussian")
ax2.plot(xax, src_box.fields[0].data[:, 22] * flux_box, label="Box")
ax2.set_ylim(top=1.1)
ax2.legend(loc="upper right")
ax2.set_title("Cut in y-direction")
ax2.set_xlabel("arcsec");
```

Dotted lines are {math}`\sigma_x` and {math}`\sigma_y` of the `Gaussian`,
which in x-direction also corresponds to {math}`R_\mathrm{eff}` of the `Sersic`.

The dash-dotted line shows the surface brightness translated to photon flux, matching the constant level of the `Box`.
In the left panel, we can see that the `Sersic` intersects this line exactly at the specified {math}`R_\mathrm{eff}` (dotted vertical lines),
which is the reference point for surface brightness in that profile.
The `Gaussian` on the other hand has its reference point at the peak, and we can see that this is indeed where the surface brightness matches.
Note that this results in much grater total integrated brightness of the `Sersic` profile,
which is correct considing the different conventions for surface brightness reference points.

## Next: building the maps, and surviving an instrument
```{TODO}
Rewrite this section.
```

This page treated the weight map as given and concentrated on the scalar it multiplies.
{doc}`flux_construction` turns that around: it builds the maps for a `Box` and a `Gaussian` in full,
shows that stating a brightness as a total or as a surface brightness produces byte-identical maps,
and works through what a finite window does (and does not) do to a cutoff profile's flux.
{doc}`flux_engine` then pushes the finished sources through a real ScopeSim optical train and checks that the detector reproduces the same totals.
