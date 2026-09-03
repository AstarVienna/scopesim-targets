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
file_format: mystnb
numbering:
  heading_1: true
  heading_2: true
---

# Rasterization and FOV details
While we briefly looked at the rasterization of parametrized targets previously in {doc}`flux_scaling`,
there are some details we skipped there, which will be discussed here:
most importantly, the effect different grid resolutions have on the realized weightmap,
how the underlying machinery deals with odd- and even-numbered pixel grids,
and how profiles without a sharp edge are affected by this.

In normal operation, the rasterization is done by ScopeSim with the grid parameters of the simulation setup,
which we will explore subsequently in {doc}`flux_engine`, so this intermediate step is rarely seen in practice.
It is however still useful to see how the rasterization affects the weight and flux stored in the {class}`~scopesim.Source` object.
In that sense, this page should be seen as "advanced information" for those who are interested.
If you just trust ScopeSim to handle everything consistently, feel free to skip this page.

## Setup
```{code-cell} ipython3
:tags: [hide-output]

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import astropy.units as u
import pandas as pd
import seaborn as sns
from scipy.special import erf
from synphot import Observation
from spextra import Spextrum, Passband

from scopesim_targets.point_source import Star
from scopesim_targets.extended_source import Box, Flat, Gaussian, Sersic
from scopesim_targets.target import FILTER_SYSTEM

flat_spec = Spextrum.flat_spectrum(
    5*u.ABmag, waves=np.linspace(3000, 35000, 6000)*u.AA
)
band = Passband(f"{FILTER_SYSTEM.name}/V")

# 3.55 Jy is a realistic medium-bright star and lands at a round ~1000 photons
# per second per cm2 in V, which keeps every number below legible.
TOTAL_FLUX = 3.55*u.Jy
```

Like on the previous page, some helper functions for this page.
Hidden by default to reduce clutter, expand if you want to see them:

```{code-cell} ipython3
:tags: [hide-cell]

sns.set_theme(style="darkgrid", palette="Set2")
bbox_flx = {"boxstyle": "round", "fc": "white", "ec": "0.5", "alpha": 0.85}
tbl_sty = [{
    "selector": "th, td",
    "props": [("text-align", "right"), ("padding", "0.2em 0.8em")],
}]

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

def get_tick_limits(grid):
    scale = grid["pixel_scale"].value
    offset = np.array([
        *2*[int(not (grid["width"] % 2))],
        *2*[int(not (grid["height"] % 2))],
    ]) * (grid["pixel_scale"] * 0.5*u.pix).to_value(u.arcsec)
    return np.ceil(np.ceil(extent_arcsec(grid)/scale)*scale) + offset

def add_pixel_grid(ax, grid, color="gray", lw=1.2, ls=":", alpha=1.,
                   major_step=False):
    scale = grid["pixel_scale"].to_value(u.arcsec/u.pixel)
    width = int(round(grid["width"]))
    height = int(round(grid["height"]))
    x_edges = (np.arange(width + 1) - width / 2) * scale
    y_edges = (np.arange(height + 1) - height / 2) * scale
    ax.set_xticks(x_edges, minor=True)
    ax.set_yticks(y_edges, minor=True)
    if major_step:
        limits = get_tick_limits(grid)
        ax.set_xticks(np.arange(*limits[:2], step=major_step))
        ax.set_yticks(np.arange(*limits[2:], step=major_step))
    ax.grid(False, which="major")  # reset seaborn's whitegrid
    ax.grid(which="minor", color=color, linewidth=lw, linestyle=ls, alpha=alpha)
    ax.tick_params(which="minor", length=0)

def x_cross_section(image, pixel_scale=1.0):
    """
    Extract the profile from the center outward along +x, handling both
    even- and odd-sized axes. The plot always starts exactly at x=0:
    for an odd axis, the central pixel is shown as a half-width bin
    (only its outward half); for an even axis, the two central pixels
    are folded/averaged as before.

    Returns (edges, values) ready for ax.stairs.
    """
    image = np.asarray(image, dtype=float)
    ny, nx = image.shape

    # --- collapse to a single representative row through y = 0 ---
    if ny % 2 == 1:
        row = image[ny // 2, :]
    else:
        row = (image[ny // 2 - 1, :] + image[ny // 2, :]) / 2.0

    # --- build the outward profile along x ---
    if nx % 2 == 1:
        cx = nx // 2
        values = row[cx:]
        centers = np.arange(values.size) * pixel_scale          # 0, 1, 2, ...
    else:
        cx_r = nx // 2
        cx_l = cx_r - 1
        k = np.arange(cx_r)
        values = (row[cx_r + k] + row[cx_l - k]) / 2.0          # fold about x=0
        centers = (k + 0.5) * pixel_scale                       # 0.5, 1.5, 2.5, ...

    edges = np.concatenate([centers - 0.5 * pixel_scale,
                             [centers[-1] + 0.5 * pixel_scale]])
    edges[0] = max(edges[0], 0.0)   # clip central bin to its outward half only
    return edges, values
```

## Different pixel scales
While briefly explored in {ref}`sec:flatgrids`, this section will expand somewhat on the effects that a different --
and sometimes suboptimal -- pixel grid has on a target's rasterization.
We will again use the {class}`~scopesim_targets.extended_source.Box` and {class}`~scopesim_targets.extended_source.Gaussian` target classes introduced previously as examples,
as both come with their own intricacies: The `Box` has a sharp border,
which means if the pixel grid does not happen to line up with it,
the rasterization must ensure the sum of the weightmap is still consistent.
The `Gaussian` on the other hand does _not_ have a defined border, meaning the rasterization itself is easier,
but the second consequence is that no finite pixel grid can (theoretically) capture its full total.

(sec:boxgrids)=
### `Box` in different grids
The example below shows the same `Box` target realized onto four different pixel grids.
The grids are given in terms of their side length (to ensure an odd number of pixels) at a constant field-of-view, resulting in pixel scales,
which are sometimes but not always integer multiples of the box's side length.
As a consequence, the rasterization has to distribute the remaining weights into the pixels along the perimeter.
The weight which ends up in each pixel naturally scales with the square of the pixel scale, i.e. the "pixel area".
The number in each box is the per-pixel weight of all the pixels which fully land inside the rasterized box.
As shown in the table below the plot, if this weight is divided by the pixel's "area", the resulting number stays constant,
meaning the **flux per solid angle**, i.e. the surface brightness, is identical, independent of the grid parameters.
This is important, because otherwise the meaning of a target's stated brightness would change with the grid.
To confirm this, the table also lists the total flux of the rasterized target,
which is calculated as the spectrum's flux multiplied by the sum of the weightmap.
As discussed on the previous page, the latter only counts the weight _within_ the field-of-view,
but because in this example, the numbers are chosen such that the box comfortably fits,
and because it has a sharp cutoff, no weight falls outside of it here.

```{code-cell} ipython3
fov = 9*u.arcsec
npixs = [9, 15, 25, 45] * u.pixel
scales = (fov / npixs).to(u.arcsec/u.pixel)

grids = [
    {"pixel_scale": scale, "width": int(npix.value), "height": int(npix.value)}
    for scale, npix in zip(scales, npixs)
]

box = Box(
    spectrum=flat_spec,
    brightness=("V", TOTAL_FLUX),
    params={"x_width": 5, "y_width": 3},
)
```

```{code-cell} ipython3
:tags: [hide-input]

fig, axes = plt.subplots(2, 2, figsize=(7, 5), layout="compressed")
lws = [1.2, 1.0, 0.8, 0.6]
extent = extent_arcsec(grids[0])  # identical for all

fluxes = []
maxs = []
for ax, grid, lw in zip(axes.flatten(), grids, lws):
    src = box.to_source(grid)
    fluxes.append(total_flux(src, band))
    wmap = src.fields[0].data
    maxs.append(wmap.max())
    im = ax.imshow(wmap, origin="lower", extent=extent, cmap="magma", vmin=0)
    cbar = fig.colorbar(im, ax=ax, location="right", shrink=1.0, pad=0)
    cbar.formatter.set_powerlimits((-3, -3))
    cbar.update_ticks()
    add_pixel_grid(ax, grid, major_step=2, color="white", lw=lw, alpha=.8)
    ax.set_aspect("equal")
    ax.text(
        0.95, 0.95,
        f"{grid['width']}\u00d7{grid['height']} pix\n"
        f"{grid['pixel_scale'].to_string(format="latex_inline")}",
        transform=ax.transAxes, va="top", ha="right",
        color="black", bbox=bbox_flx, fontsize=10,
    )
    ax.text(0, 0, wmap.max().round(4), ha="center", va="center", fontsize=16)

axes[1, 0].set_xlabel("arcsec")
axes[1, 1].set_xlabel("arcsec")
axes[0, 0].set_ylabel("arcsec")
axes[1, 0].set_ylabel("arcsec")
axes[0, 0].tick_params(labelbottom=False)
axes[0, 1].tick_params(labelbottom=False)
axes[0, 1].tick_params(labelleft=False)
axes[1, 1].tick_params(labelleft=False)
fig.text(
    0.99, 0.5,
    "dimensionless weight",
    fontsize=12,
    rotation=90,
    va="center",
    ha="right",
)
plt.subplots_adjust(wspace=0.5, hspace=0.15)
plt.show()
```

```{code-cell} ipython3
:tags: [hide-input]

pd.DataFrame({
    "npix": npixs,
    "scale": [s for s in scales],  # to get unit in cells
    "weight/area": [m/s.value**2 for s, m in zip(scales, maxs)],
    "total flux": fluxes,
}).style.format({
    "npix": "{:.0f}",
    "scale": "{:.2f}",
    "weight/area": "{:.4f}",
    "total flux": "{:.0f}",
}).set_table_attributes("class='dataframe'"
).set_table_styles(tbl_sty).hide(axis="index")
```

### `Box` in a single pixel
In the next example we use the same `Box` target, but evaluated onto three different grids.
Each has a pixel scale large enough so that the whole box fits well within a single pixel.
However, whether a grid has a _single central pixel_ depends on its parity.
We examine all three possible cases:

* an _even_ number of pixels along both axes, meaning the center falls exactly between 4 pixels
* _even_ in one direction, but _odd_ in the other, so 2 pixels share the center
* an _odd_ number of pixels in both directions, finally resulting in a single central pixel

In the plots below, the center point is always marked with a teal cross.

```{code-cell} ipython3
large_grids = [
    {"pixel_scale": 10*u.arcsec/u.pixel} | npix for npix in [
    {"title": "even grid", "width": 4, "height": 4},
    {"title": "odd/even grid", "width": 5, "height": 4},
    {"title": "odd grid", "width": 5, "height": 5},
]]
```

```{code-cell} ipython3
:tags: [hide-input]

fig, axes = plt.subplots(1, 3, figsize=(8, 4), layout="compressed")
cmap = plt.get_cmap("magma", 5)
sums = []
maxs = []
for ax, large_grid in zip(axes, large_grids):
    src_large = box.to_source(large_grid)
    wmap = src_large.fields[0].data
    sums.append(wmap.sum())
    maxs.append(wmap.max())

    im = ax.imshow(
        wmap, origin="lower", extent=extent_arcsec(large_grid),
        cmap=cmap, vmin=0, vmax=1,
    )
    ax.scatter(0, 0, s=50, c="C0", marker="x")
    add_pixel_grid(ax, large_grid, major_step=10, color="white", alpha=.8)
    ax.set_xlabel("arcsec")
    ax.set_title(large_grid["title"])
    ax.set_aspect("equal")

axes[0].set_ylabel("arcsec")
axes[1].tick_params(labelleft=False)
axes[2].tick_params(labelleft=False)

cbar = fig.colorbar(
    im, ax=axes[-1], location="right", shrink=1.0,
    label="dimensionless weight",
    ticks=[0, 0.2, 0.4, 0.6, 0.8, 1],
)
plt.subplots_adjust(wspace=0.2)
plt.show()
```

```{code-cell} ipython3
:tags: [hide-input]

pd.DataFrame({
    "grid": [f"{g['width']} \u00d7 {g['height']}" for g in large_grids],
    "weight/pix": maxs,
    "weightmap sum": sums,
}).style.format({
    "weight/pix": "{:.3f}",
    "weightmap sum": "{:.1f}",
}).set_table_attributes("class='dataframe'"
).set_table_styles(tbl_sty).hide(axis="index")
```

The results for the three cases are what we would expect:

* In the case of the even grid, the flux is distributed evenly between the 4 central pixels, which thus each store a value of 0.25.
* Similarly with the odd and even grid, the flux is distributed over 2 pixels which get 0.5 each.
* Finally for the all-odd grid, the whole weight ends up in the central pixel, which contains a value of 1.

### `Gaussian` truncation
In cases where a finite profile overflows the field-of-view, the expected analytic total can still be calculated.
For the `Box`, this can be written as {math}`\text{area}(\text{box} \cap \text{FOV}) / (w_x\,w_y)`.
For the `Gaussian`, integrated over a centered window of width {math}`W` and height {math}`H` (both in angular units), the map sums to
{math}`\operatorname{erf}\left(\frac{\sqrt{2}\,W}{4\,\sigma_x}\right)\operatorname{erf}\left(\frac{\sqrt{2}\,H}{4\,\sigma_y}\right) < 1`.
Since the pixel weights are already scaled to their contribution to the _total_ integral {math}`P` anyway,
the scaling is unaffected by the extent of {math}`p` relative to the field-of-view,
the only practical difference being {math}`\sum_{ij} w_{ij}`.
In the example below, we use an asymmetrical `Gaussian` to showcase this effect.

```{code-cell} ipython3
sx, sy = 6.0, 4.0  # arcsec
gauss_grid = {"pixel_scale": 0.5*u.arcsec/u.pixel, "width": 25, "height": 25}

gaussian = Gaussian(
    spectrum=flat_spec,
    brightness=("V", TOTAL_FLUX),
    params={"x_stddev": sx, "y_stddev": sy},
)
src_gauss = gaussian.to_source(gauss_grid)
wmap = src_gauss.fields[0].data

fraction = wmap.sum()
spec_total = Observation(src_gauss.fields[0].spectrum, band).integrate()
in_window = total_flux(src_gauss, band)

Lx = gauss_grid["width"]*u.pix * gauss_grid["pixel_scale"]
Ly = gauss_grid["height"]*u.pix * gauss_grid["pixel_scale"]
analytic = erf((np.sqrt(2)*Lx)/(4*sx*u.arcsec)) * erf((np.sqrt(2)*Ly)/(4*sy*u.arcsec))
```

```{TODO}
rm the *u.pix and *u.arcsec above once units are fully in grid and model
```

```{code-cell} ipython3
:tags: [hide-input]

fig, ax = plt.subplots(figsize=(6.5, 6), layout="compressed")
im = ax.imshow(
    wmap, origin="lower", extent=extent_arcsec(gauss_grid),
    cmap="magma", vmin=.55e-3, vmax=1.65e-3,
)
add_pixel_grid(ax, gauss_grid, major_step=2, color="white", lw=1)
ax.add_patch(Ellipse((0, 0), 2*sx, 2*sy, fill=False, ec="C0", ls="--", lw=1.5))
ax.text(.75*sx, .7*sy, r"1$\sigma$", color="C0", fontsize=14)
ax.set_aspect("equal")
ax.set_xlabel("arcsec")
ax.set_ylabel("arcsec")
cbar = fig.colorbar(im, ax=ax, label="dimensionless weight", shrink=1.0)
cbar.formatter.set_powerlimits((-3, -3))
cbar.update_ticks()
plt.show()
```

```{code-cell} ipython3
:tags: [hide-input]

pd.DataFrame({
    "quantity": ["weightmap sum", "analytic window fraction"],
    "value": [fraction, analytic],
}).style.format({"value": "{:.3f}"}
).set_table_attributes("class='dataframe'"
).set_table_styles(tbl_sty).hide(axis="index")
```

```{code-cell} ipython3
:tags: [hide-input]

pd.DataFrame({
    "quantity": [
        "total flux in spectrum",
        "in-window flux (fraction \u00d7 total)"
    ],
    "value": [spec_total, in_window],
}).style.format({"value": "{:.1f}"}
).set_table_attributes("class='dataframe'"
).set_table_styles(tbl_sty).hide(axis="index")
```

The agreement of the weightmap sum and the analytic solution down to {math}`<0.2\,\%` proves that the scaling works correctly.
As mentioned, the total flux _stored in the spectrum_ is unaffected by truncation,
whereas the flux _inside the field-of-view window_ is the product of the spectrum's flux and the weightmap.

### `Gaussian` in different grids
In the previous example, we used a rather fine grid that sampled the `Gaussian` profile very well.
If the grid is much coarser relative to the profile (as characterized by its {math}`\sigma_x` and {math}`\sigma_y`),
the deviation from the analytic solution increases slightly due to the finite approximation of the sampling.

Below we use a symmetrical `Gaussian` sampled onto different pixel grids, while keeping the field-of-view constant,
similar to the `Box` example shown in {ref}`sec:boxgrids`.

```{code-cell} ipython3
fov = 16.0*u.arcsec
scales = [3.2, 1.6, 0.8, 0.4, 0.2, 0.1] * u.arcsec/u.pixel
npixs = (fov / scales).to_value(u.pixel)

grids = [
    {"pixel_scale": scale, "width": int(npix), "height": int(npix)}
    for scale, npix in zip(scales, npixs)
]

sigma = 4.0
gauss = Gaussian(
    spectrum=flat_spec,
    brightness=("V", TOTAL_FLUX),
    params={"x_stddev": sigma, "y_stddev": sigma},
)

carried = []
profiles = []
for grid in grids:
    wmap = gauss.to_source(grid).fields[0].data
    carried.append(wmap.sum())
    profiles.append(x_cross_section(wmap, grid["pixel_scale"].value))

analytic = erf((np.sqrt(2)*fov.to_value(u.arcsec))/(4*sigma))**2
```

```{code-cell} ipython3
:tags: [hide-input]

fig, ax = plt.subplots(figsize=(5.5, 4), layout="compressed")
ax.plot(scales, carried, marker="o", c="C0")
ax.set_xscale("log")
ax.set_xinverted(True)
ax.set_xticks(scales.value, labels=[str(s.value) for s in scales])
ax.set_ylim(0.9, 0.93)
ax.axhline(analytic, color="gray", ls=":", label="analytic erf fraction")
ax.set_xlabel(
    r"$\varpi_\mathrm{src}$ "
    f"[{scales.unit.to_string(format='latex')}]"
)
ax.set_ylabel("carried fraction  (weightmap sum)")
ax.set_title(f"Gaussian $\\sigma$={sigma}\u2033, fixed {fov} window")
ax.legend()
plt.show()
```

We can see that the weightmap sum approaches the analytic total the finer the grid is.
The y-axis on this plot is zoomed in, so the difference is somewhat exaggerated:
The actual deviation even for the worst grid is only around {math}`1\,\%`, see table below.

As previously discussed in {ref}`sec:boxgrids`, the weight stored in each pixel scales with the pixel's "area", i.e. {math}`w_{ij}\propto\Omega`.
The plot below visualizes this by showing the cross-section of each weightmap, from the center outwards along the x-direction.
Unlike the `Box` example, where the constant weight within the box is directly scaled by the pixel's "area",
the `Gaussian` (or any other profile with spatial structure) scales the entire profile.
Since the pixel scales of the grids used here are always a factor of 2 smaller or larger than the previous one,
the weight in each pixel scales by a factor of 4 between them.
The y-axis in the plot below is scaled logarithmically, which makes this less obvious between individual grids,
but it shows very nicely between the first and last one.
The difference in scale between those is a factor of 32, thus the weight must scale by a factor of 1024,
and indeed we can read from the y-axis that the first and last lines are 3 orders of magnitude apart.

```{code-cell} ipython3
:tags: [hide-input]

fig, ax = plt.subplots(figsize=(7, 4), layout="compressed")
for x in profiles[0][0]:
    ax.axvline(x, c="gray", ls="--", alpha=.5)
ax.axvline(sigma, c="black", ls="-", alpha=.5)
ax.text(sigma+.05, 4e-4, r"1$\sigma$", rotation=90)
ax.set_xticks(profiles[2][0])
for (edges, values), scale in zip(profiles, scales):
    ax.stairs(values, edges, baseline=None, label=scale.to_string(format="latex_inline"))
ax.set_xlabel("x [arcsec]")
ax.set_ylabel("weightmap value")
ax.set_ylim(bottom=1e-5)
ax.set_yscale("log")
ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1))
plt.show()
```

The grid lines along the x-axis match the pixels of the 0.8 arcsec/pix grid.
Gray dashed lines mark the pixels of the coarsest 3.2 arcsec/pix grid, solid gray line is set at {math}`1\,\sigma`.
Note that in this case, because the number of pixels is odd (5 × 5), it starts with a half-pixel in the center.
All subsequent finer grids are even-numbered and thus can start with a full pixel.
We can also draw the same plot, but this time divide each pixel by its "area".
Now the lines fall nicely on top of each other, and we can also observe the reason
for the larger deviation from the analytic weightmap sum on the coarsest grid:
If we "eyeball" the quadrature of each pixel, we can see that the two outer pixels "cancel out" the over- and undercutting of the curve,
whereas the central half-pixel's value is kept at the same peak value, but takes longer to fall off,
thus bringing slightly more total weight into the field-of-view as the others,
where the central peak value is distributed between the central 4 pixels (even grids).

```{code-cell} ipython3
:tags: [hide-input]

fig, ax = plt.subplots(figsize=(7, 4), layout="compressed")
for x in profiles[0][0]:
    ax.axvline(x, c="gray", ls="--", alpha=.5)
ax.axvline(sigma, c="black", ls="-", alpha=.5)
ax.text(sigma+.05, 2e-3, r"1$\sigma$", rotation=90)
ax.set_xticks(profiles[2][0])
for (edges, values), scale in zip(profiles, scales):
    ax.stairs(
        values/scale.value**2, edges, baseline=None,
        label=scale.to_string(format="latex_inline"),
    )
ax.set_xlabel("x [arcsec]")
ax.set_ylabel("weightmap value")
ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1))
plt.show()
```


```{code-cell} ipython3
:tags: [hide-input]

pd.DataFrame({
    "scale [arcsec/px]": scales,
    "px / sigma": sigma / scales,
    "carried fraction": carried,
    "deviation from analytic": carried/analytic - 1,
}).style.format({
    "scale [arcsec/px]": "{:.1f}",
    "px / sigma": "{:.2f}",
    "carried fraction": "{:.4f}",
    "deviation from analytic": "{:.2%}",
}).set_table_attributes("class='dataframe'"
).set_table_styles(tbl_sty).hide(axis="index")
```

```{TODO}
Keep both plots?

Move table further up? Or at least mention it?
```

```{TODO}
Consider linking to the next page(s).
```
