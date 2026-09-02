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

# Extended source flux through an instrument
Now that we have established in {doc}`flux_scaling` how flux is defined and stored in both point sources and extended profiles,
and have seen the details of rasterization in {doc}`flux_construction`,
the next logical step is to see what happens if we push these sources through a full ScopeSim instrument simulation.
To keep things simple enough to understand the effect of the instrument,
we exclude the atmosphere, the point spread function (PSF) and all detector-level effects.
That means we will "intercept" the simulated image at the image plane level,
theoretically right before the simulated photons hit the detector.

The unit we will see in these images is {math}`ph\,s^{-1}`, i.e. photons per second.
This is because the "per area" dimension of the "flux images" we saw previously has been eliminated by multiplying with the telescope's collecting area,
which happens inside the {class}`~scopesim.optics.FieldOfView` class within ScopeSim.
If we included the quantum efficiency effect, this would technically change to electrons per second,
but that is a property of the detector, so it is disabled here.

Throughout this page, we will use two different "pixel scales":

* {math}`\varpi_\mathrm{src}` -- the pixel scale of the **rasterized source** weightmap (the one we saw on the previous pages).
* {math}`\varpi_\mathrm{inst}` -- the pixel scale of the **detector** (set on the optical train).

Under normal circumstances, those two scales would be identical, i.e. {math}`\varpi_\mathrm{src}\equiv\varpi_\mathrm{inst}`.
The whole point of using parametrized models for extended sources is to be able to always evaluate them onto the same pixel grid that the detector uses,
so any rasterization and interpolation side-effects are minimized.
But they are technically independent quantities, and on this page we will explore the relation between them,
and validate that ScopeSim _can_ deal with cases where the two _do not_ match.

## Setup
```{code-cell} ipython3
:tags: [hide-output]

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch
import astropy.units as u
import pandas as pd
import seaborn as sns
from scipy.special import erf
from synphot import Observation
from spextra import Spextrum, Passband

from scopesim import example_simulation

from scopesim_targets.point_source import Star, StarField
from scopesim_targets.extended_source import Box, Flat, Gaussian
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

Some more helper functions to set up and read out the ScopeSim simulations:

```{code-cell} ipython3
def setup_sim(grid, props=None):
    props = props or {}
    scale_inst = grid["pixel_scale"]
    pixel_scale = scale_inst.to_value(u.arcsec/u.pixel)
    plate_scale = (scale_inst / (10*u.um/u.pix)).to_value(u.arcsec/u.mm)
    sim = example_simulation(
        properties={
            "!INST.pixel_scale": pixel_scale,
            "!INST.plate_scale": plate_scale,
            "!DET.width": grid["width"],
            "!DET.height": grid["height"],
            "!OBS.filter_name": "BrGamma",
        } | props,
        ignore_effects=["psf", "atmospheric_radiometry", "qe_curve"],
    )
    return sim

def detected(sim, source):
    sim(source, update=True)
    data = sim.optical_train.image_planes[0].data
    bunit = u.Unit(sim.optical_train.image_planes[0].header["BUNIT"])
    return data.clip(min=1e-3) * bunit
```

And finally, here are the reusable target definitions:

```{code-cell} ipython3
star = Star(
    position=(0, 0),
    spectrum=flat_spec,
    brightness=("V", TOTAL_FLUX),
)

box_params = {"x_width": 6, "y_width": 4}
box_total = Box(
    spectrum=flat_spec,
    brightness=("V", TOTAL_FLUX),
    params=box_params,
)
```

## Point, 2D-integrated, and 2D-surface-brightness agree
To see how the flux "survives" through the optical train, we create three different sources, all with the same total:
A simple point source (`star`), a rectangular box with the flux given as the total (`box_total`),
and another box of identical size, but with the flux specified as a surface brightness.
Below we use the definition of the box's "area" to calculate a surface brightness that produces the same total.

```{TODO}
Replace the manual grid once `.to_source` is in Simulation!

Remove `*u.arcsec` in A_eff once box params are in arcsec.
```

```{code-cell} ipython3
A_eff = box_params["x_width"]*u.arcsec * box_params["y_width"]*u.arcsec
sb = (TOTAL_FLUX / A_eff).to(u.Jy/u.sr)

box_sb = Box(
    spectrum=flat_spec,
    brightness=("V", sb),
    params=box_params,
)

sb
```

Next we set up a very small fake detector with an odd number of pixels so everything lines up nicely.
The filter is deliberately set to the very narrow "BrGamma", to limit the dynamic range on the resulting image.
Each source is run through the exact same simulation,
and in each case we look at the noiseless image just before the detector, as described above.

```{code-cell} ipython3
:tags: [hide-output]

grid = {"pixel_scale": 1*u.arcsec/u.pixel, "width": 9, "height": 9}

sim_9 = setup_sim(grid)
img_point = detected(sim_9, star.to_source())
img_box = detected(sim_9, box_total.to_source(grid))
img_sb = detected(sim_9, box_sb.to_source(grid))

images = {
    "point": img_point,
    "box (total flux)": img_box,
    "box (surface brightness)": img_sb,
}
```

```{code-cell} ipython3
:tags: [hide-input]

fig, axes = plt.subplots(1, 3, figsize=(9, 4), layout="compressed")
for ax, (title, img) in zip(axes, images.items()):
    im = ax.imshow(
        img.value, origin="lower", extent=extent_arcsec(grid),
        cmap="inferno", norm="log", vmin=1, vmax=1e5,
    )
    add_pixel_grid(ax, grid)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])

axes[0].scatter(0, 0, c="red", label="Star Target")
axes[1].add_patch(Rectangle(
    (0 - box_params["x_width"]/2, 0 - box_params["y_width"]/2),
    box_params["x_width"], box_params["y_width"],
    fill=False, ec="red", ls="--", lw=2, label="Box Target",
))
axes[2].add_patch(Rectangle(
    (0 - box_params["x_width"]/2, 0 - box_params["y_width"]/2),
    box_params["x_width"], box_params["y_width"],
    fill=False, ec="red", ls="--", lw=2, label="Box Target",
))
axes[0].legend(loc="upper right")
axes[1].legend(loc="upper right")
axes[2].legend(loc="upper right")

cbar = fig.colorbar(
    im, ax=axes[-1], location="right", shrink=1.0, pad=0.05,
    label="detected flux",
)
plt.subplots_adjust(wspace=0.2)
plt.show()
```

```{code-cell} ipython3
:tags: [hide-input]

pd.DataFrame({
    "representation": list(images.keys()),
    "detected flux": [img.sum() for img in images.values()],
}).style.format({"detected flux": "{:.0f}"}
).set_table_attributes("class='dataframe'"
).set_table_styles(tbl_sty).hide(axis="index")
```

The edges of the box target run exactly through the surrounding pixel centers,
so the flux in the edge and corner pixels is half and a quarter of the central pixels, respectively.
In the point source, the same flux is obviously much more concentrated (ignoring the PSF),
so the value in the central pixel is ~24 times higher than inside the box, in this case.

The final values in {math}`ph\,s^{-1}` we see in the table above are notably larger than the {math}`F(V)` we saw in {doc}`flux_scaling`,
even though we are using the exact same {math}`3.55\;\mathrm{Jy}` as an input.
This is easily explained by the fact that {math}`F(V)` is in {math}`ph\,s^{-1}\,cm^{-2}`, i.e. _per area_,
whereas now we are simulating a telescope with known collecting area,
so this dimension is multiplied out in the simulation, thus increasing the total value.
Additionally, as mentioned above we are using a much narrower filter than the brightness definition,
which was given in the V-band; otherwise the number would be even higher.
We can calculate all of this manually to verify the number produced by ScopeSim:

```{code-cell} ipython3
thru = sim_9.optical_train["filter_wheel"].current_filter.throughput
area = sim_9.settings["!TEL.area"]

area
```

```{code-cell} ipython3
total_flux(star.to_source(), thru) * area << u.ph/u.s
```

which matches the numbers seen in the table, so everything is indeed consistent.

## Points and pixels
```{note}
:class: margin

ScopeSim _can_ place point source centroids with sub-pixel accuracy, but that's a different story.
```
An important concept to understand in the context of this page it that once we are in the detector's pixel grid,
the difference between point sources and extended ones becomes a matter of scale only.
The two subsections below discuss both ends of this.

### An unresolved extended source is a point source
One pixel is now the smallest scale that we can resolve -- meaning a point source is exactly that (ignoring again the PSF).
Likewise, any extended source that is _smaller_ than that scale will look identical to a point source.
The example below shows exactly that: again a point source and box with identical flux,
but the box is only a quarter the length of one detector pixel, or 1/16 the area.
To avoid any rasterization issues, we sample the box on a smaller grid, as shown below.

```{code-cell} ipython3
:tags: [hide-output]

inst_grid = {"pixel_scale": 2.0*u.arcsec/u.pixel, "width": 9, "height": 9}
tiny_grid = {"pixel_scale": 0.5*u.arcsec/u.pixel, "width": 4, "height": 4}

tiny_box_src = Box(
    spectrum=flat_spec,
    brightness=("V", TOTAL_FLUX),
    params={"x_width": 0.5, "y_width": 0.5},
).to_source(tiny_grid)

sim_9 = setup_sim(inst_grid)
img_tiny = detected(sim_9, tiny_box_src)
img_pt = detected(sim_9, star.to_source())

images_unres = {"point": img_pt, "0.5\u2033 box (unresolved)": img_tiny}
```

```{code-cell} ipython3
:tags: [hide-input]

fig = plt.figure(figsize=(8, 6), layout="compressed")
gs = fig.add_gridspec(2, 2, height_ratios=[1, .6])
axes = [
    fig.add_subplot(gs[0, 0]),
    fig.add_subplot(gs[0, 1]),
]

ax_zoom = fig.add_subplot(gs[1, 1])

for ax, (title, img) in zip(axes, images_unres.items()):
    im = ax.imshow(
        img.value, origin="lower", extent=extent_arcsec(inst_grid),
        cmap="inferno", norm="log", vmin=1, vmax=1e5,
    )
    add_pixel_grid(ax, inst_grid)
    ax.set_aspect("equal")
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])

axes[0].scatter(0, 0, c="red", label="Star Target")
axes[0].legend(loc="upper right")
fig.colorbar(
    im, ax=axes[1], location="right", shrink=1.0, pad=0.05,
    label="detected flux",
)

wmap = tiny_box_src.fields[0].data

im_zoom = ax_zoom.imshow(
    wmap, origin="lower", extent=extent_arcsec(tiny_grid),
    cmap=plt.get_cmap("magma", 5), vmin=0, vmax=1,
)

add_pixel_grid(ax_zoom, tiny_grid, major_step=.5)
ax_zoom.set_xlabel("arcsec")
ax_zoom.set_ylabel("arcsec")
ax_zoom.set_aspect("equal")

ax_zoom.set_title("Zoomed central pixel")

ax_zoom.add_patch(Rectangle(
    (-.25, -.25), .5, .5,
    fill=False, ec="red", ls="--", lw=2, label="Box Target",
))

ax_zoom.legend(loc="upper right")

fig.colorbar(
    im_zoom, ax=ax_zoom,
    location="right", shrink=1.0, pad=0.05,
    label="dimensionless weight",
    ticks=[0, .2, .4, .6, .8, 1],
)

fig.add_artist(ConnectionPatch(
    xyA=(-1, -1), coordsA=axes[1].transData,
    xyB=(-1, 1), coordsB=ax_zoom.transData,
    color="C0", ls="--", lw=1.5, alpha=.8,
))
fig.add_artist(ConnectionPatch(
    xyA=(1, -1), coordsA=axes[1].transData,
    xyB=(1, 1), coordsB=ax_zoom.transData,
    color="C0", ls="--", lw=1.5, alpha=.8,
))

plt.subplots_adjust(wspace=0.2, hspace=0.2)
plt.show()
```

```{code-cell} ipython3
:tags: [hide-input]

pd.DataFrame({
    "target": list(images_unres.keys()),
    "detected flux": [img.sum() for img in images_unres.values()],
}).style.format({"detected flux": "{:.0f}"}
).set_table_attributes("class='dataframe'"
).set_table_styles(tbl_sty).hide(axis="index")
```

Both cases place the same flux in the central pixel, matching the identically stated brightness.

### A resolved extended source is a grid of point sources
```{TODO}
Look into bug when no subpix.
```
Perhaps conceptually a bit more surprising but essential to understanding the internals of ScopeSim is the second consequence of pixelization:
If we construct a "cluster" of point sources placed exactly at the centers of an extended source's pixels
and scale the flux correctly, we end up with an image that is indistinguishable from the true extended source.
Indeed this visualizes how ScopeSim handles a rasterized 2D source internally,
by placing the flux in the same bins as the individual point sources.

```{code-cell} ipython3
:tags: [hide-output]

grid = {"pixel_scale": 1*u.arcsec/u.pixel, "width": 8, "height": 8}
sim_8 = setup_sim(grid, props={"!SIM.sub_pixel.flag": True})
src_box = box_total.to_source(grid)

# TODO: make these lists dynamic
pxposx = [-2.5, -1.5, -.5, .5, 1.5, 2.5]
pxposy = [-1.5, -.5, .5, 1.5]
xs, ys = np.meshgrid(pxposx, pxposy)

src_24stars = StarField(
    band="V",
    positions=[(x, y) for x, y in zip(xs.flatten(), ys.flatten())],
    spectra=24*[flat_spec],
    brightnesses=24*[TOTAL_FLUX/24]
).to_source()

img_box = detected(sim_8, src_box)
img_24stars = detected(sim_8, src_24stars)
images24 = {"Box": img_box, "Star grid": img_24stars}
```

```{code-cell} ipython3
:tags: [hide-input]

fig, axes = plt.subplots(1, 2, figsize=(8, 4), layout="compressed")
for ax, (title, img) in zip(axes, images24.items()):
    im = ax.imshow(
        img.value, origin="lower", extent=extent_arcsec(grid),
        cmap="inferno", #norm="log", vmin=1, vmax=1e5,
    )
    # add_pixel_grid(ax, grid, color="gray")
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])

add_pixel_grid(axes[0], grid)
add_pixel_grid(axes[1], grid)

axes[0].add_patch(Rectangle(
    (0 - box_params["x_width"]/2, 0 - box_params["y_width"]/2),
    box_params["x_width"], box_params["y_width"],
    fill=False, ec="red", ls="--", lw=2, label="Box Target",
))
axes[1].scatter(
    src_24stars.fields[0].field["x"],
    src_24stars.fields[0].field["y"],
    c="red", label="Star grid target",
)
axes[0].legend(loc="upper right")
axes[1].legend(loc="upper right")

cbar = fig.colorbar(
    im, ax=axes[-1], location="right", shrink=1.0, pad=0.05,
    label="detected flux",
)
cbar.formatter.set_powerlimits((4, 4))
cbar.update_ticks()
plt.subplots_adjust(wspace=0.2)
plt.show()
```

```{code-cell} ipython3
:tags: [hide-input]

pd.DataFrame({
    "target": list(images24.keys()),
    "detected flux": [img.sum() for img in images24.values()],
}).style.format({"detected flux": "{:.0f}"}
).set_table_attributes("class='dataframe'"
).set_table_styles(tbl_sty).hide(axis="index")
```

As we can see, the total flux in the image again matches in both cases.

## Non-matching pixel scales
In most previous examples (with the exception of the tiny box),
the extended source's brightness profile was evaluated onto the detector's grid,
and indeed this is the "normal" case for parametrized targets.
But there might be circumstances where this is not the case, so for the sake of completeness,
it is worth investigating what happens when {math}`\varpi_\mathrm{src}` and {math}`\varpi_\mathrm{inst}` _do not_ match.

The example below explores both scenarios:

* Fixed {math}`\varpi_\mathrm{inst}` (i.e. the same simulated detector) and varying {math}`\varpi_\mathrm{src}`.
* Fixed {math}`\varpi_\mathrm{src}` and varying {math}`\varpi_\mathrm{inst}` (i.e. the same source on different detectors)

Both over- and undersampled cases are covered.
The detector grid is large enough to cover the box target,
except in the most extreme cases where some flux is lost to interpolation.

```{code-cell} ipython3
:tags: [hide-output]

inst_grid = {"pixel_scale": 0.2*u.arcsec/u.pixel, "width": 55, "height": 55}
sim_55 = setup_sim(inst_grid)

src_scales = [0.05, 0.1, 0.2, 0.4, 0.8, 1.6]*u.arcsec/u.pixel
fluxes_src = []
for scale_src in src_scales:
    src_grid = {"pixel_scale": scale_src, "width": 128, "height": 128}
    img_box = detected(sim_55, box_total.to_source(src_grid))
    # Subtract clipped background for different img sizes
    fluxes_src.append(img_box.sum() - img_box.size*.001*u.ph/u.s)

inst_scales = [0.05, 0.1, 0.2, 0.4, 0.8, 1.6]*u.arcsec/u.pixel
src_grid = {"pixel_scale": 0.2*u.arcsec/u.pixel, "width": 32, "height": 32}
src_box = box_total.to_source(src_grid)

fluxes_inst = []
for scale_inst in inst_scales:
    inst_grid = {"pixel_scale": scale_inst, "width": 127, "height": 127}
    sim_75 = setup_sim(inst_grid)
    img_box = detected(sim_75, src_box)
    # Subtract clipped background for different img sizes
    fluxes_inst.append(img_box.sum() - img_box.size*.001*u.ph/u.s)
```

```{code-cell} ipython3
:tags: [hide-input]

fig, ax = plt.subplots(figsize=(5, 3), layout="compressed")
ax.plot(
    src_scales, [flux.value for flux in fluxes_src],
    marker="o", c="C0", label="change src scale",
)
ax.plot(
    inst_scales, [flux.value for flux in fluxes_inst],
    marker="o", c="C1", label="change inst scale",
)
ax.set_xscale("log")
ax.set_xticks(src_scales.value, labels=[str(s.value) for s in src_scales])
ax.set_xlabel(
    r"$\varpi$ "
    f"[{src_scales.unit.to_string(format='latex')}]"
)
ax.legend()
ax.set_ylabel("detected flux")
plt.show()
```

```{code-cell} ipython3
:tags: [hide-input]

pd.DataFrame({
    "scale": [scale for scale in src_scales],
    "detected flux for changing src scale": fluxes_src,
    "detected flux for changing inst scale": fluxes_inst,
}).style.format({
    "detected flux for changing src scale": "{:.0f}",
    "detected flux for changing inst scale": "{:.0f}",
}
).set_table_attributes("class='dataframe'"
).set_table_styles(tbl_sty).hide(axis="index")
```

The point at {math}`0.2\;\mathrm{\frac{{}^{\prime\prime}}{pix}}` where both lines meet corresponds to the fixed scale in both cases,
meaning under normal conditions, everything is still consistent.
The worst deviation occurs for the smallest detector scale,
where the fixed number of pixels results in a reduced field-of-view,
which -- even though the whole box target still fits inside -- seems to lose some flux to interpolation.
Note however the extremely zoomed-in scale on the y-axis here:
Even the largest discrepancy is actually happening below the {math}`0.01\,\%` scale!
