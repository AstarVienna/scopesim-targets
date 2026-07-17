---
file_format: mystnb
kernelspec:
  name: python3
---

# Extended-Source Flux: What Gets Stored
A target can be handed to ScopeSim as a point source, as a 2D image plus a spectrum, or as a datacube.
However you describe it, the same physical object must carry the same total flux. {doc}`flux_scaling` opened up one half of that statement -- how a spectrum,
seen through a filter, becomes the single scalar {math}`F_\mathrm{band}` that the weights multiply.
This page opens up the *other* half: how a profile's geometry becomes the dimensionless weight map itself,
using {py:class}`~scopesim_targets.extended_source.Box` and {py:class}`~scopesim_targets.extended_source.Gaussian`,
compared against a plain {py:class}`~scopesim_targets.point_source.Star`.
Everything here happens *before* a source ever reaches an optical train -- {doc}`flux_engine` continues the story through a real instrument.

Two things are worth keeping separate in what follows:

* **How the brightness is stated** -- as a *total* (integrated) flux, or as a *surface brightness* (an amount per solid angle, sampled at some reference point of the profile).
  Both describe the same kind of object; ScopeSim must read either one and arrive at the same total.
* **How much of that total lands inside your field of view** -- a uniform box fully contained in the window carries all of it;
  a smooth profile's wings can extend past the edge, so only a *fraction* of the analytic total is actually rendered.

See {doc}`defining_brightness` for the full brightness grammar.

## Setup
```{code-cell} ipython3
:tags: [hide-output]

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from matplotlib.patches import Ellipse
import astropy.units as u
import pandas as pd
import seaborn as sns
from scipy.special import erf
from synphot import Observation
from spextra import Spextrum, Passband

from scopesim_targets.extended_source import Box, Gaussian
from scopesim_targets.point_source import Star
from scopesim_targets.target import FILTER_SYSTEM

sns.set_theme(style="dark", palette="Set2")
```

Every example below shares one spectral shape and one bandpass, so only the *brightness* and the *pixel scale* ever change between representations.

```{code-cell} ipython3
# Flat in AB magnitude (i.e. flat in f_nu) -- deliberately *not* flat in
# PHOTLAM, so a wrong reference-wavelength scaling would actually show up.
flat_spec = Spextrum.flat_spectrum(
    5*u.ABmag, waves=np.linspace(3000, 35000, 6000)*u.AA
)

def passband(name="V"):
    return Passband(f"{FILTER_SYSTEM.name}/{name}")

TOTAL_FLUX = 3.55*u.Jy
BOX_PARAMS = {"x_width": 5, "y_width": 3}  # arcsec
```

Like on the previous page, some helper functions for this page.
Hidden by default to reduce clutter, expand if you want to see them:

```{code-cell} ipython3
:tags: [hide-cell]
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

## Same box, two ways to say how bright it is
A `Box` can be given a total flux directly, or an equivalent surface brightness.
Both describe the same rectangle of light:

```{code-cell} ipython3
:tags: [hide-output]

optical_train = {
    "pixel_scale": 1*u.arcsec/u.pixel, "width": 11, "height": 11,
}

box_total = Box(spectrum=flat_spec, brightness=("V", TOTAL_FLUX), params=BOX_PARAMS)

A_eff = box_total.total_flux_factor()  # analytic area of the box, in arcsec**2
sb_equivalent = (TOTAL_FLUX / A_eff).to(u.Jy/u.sr)
box_sb = Box(spectrum=flat_spec, brightness=("V", sb_equivalent), params=BOX_PARAMS)

src_total = box_total.to_source(optical_train)
src_sb = box_sb.to_source(optical_train)
```

```{code-cell} ipython3
data_total = src_total.fields[0].data
data_sb = src_sb.fields[0].data
vmax = max(data_total.max(), data_sb.max())

fig, axes = plt.subplots(1, 2, figsize=(8, 4), sharey=True)
extent = extent_arcsec(optical_train)
for ax, data, title in zip(
    axes, (data_total, data_sb), ("stated as total flux", "stated as surface brightness")
):
    im = ax.imshow(data, origin="lower", extent=extent, cmap="viridis", vmin=0, vmax=vmax)
    add_pixel_grid(ax, optical_train)
    ax.set_aspect("equal")
    ax.set_title(title)
    ax.set_xlabel("x [arcsec]")
axes[0].set_ylabel("y [arcsec]")

fig.subplots_adjust(right=0.86, wspace=0.3)
pos = axes[1].get_position()
cbar_ax = fig.add_axes([0.89, pos.y0, 0.02, pos.y1 - pos.y0])
fig.colorbar(im, cax=cbar_ax, label=r"dimensionless weight")
plt.show()

print("weight maps identical:",
      np.allclose(src_total.fields[0].data, src_sb.fields[0].data))
print(f"flux (total-flux form):         {source_total_flux(src_total, passband()):.4f}")
print(f"flux (surface-brightness form): {source_total_flux(src_sb, passband()):.4f}")
```

The two images are pixel-for-pixel identical, and both reconstruct the same {math}`3.55\;\mathrm{Jy}`.
That is not a coincidence of this particular example: the weight map only ever encodes the profile's *shape*,
normalized by its analytic area -- it is built the same way regardless of how the brightness was stated.
The *amount* lives entirely in the attached spectrum, see {doc}`flux_scaling` for details on that.
Stating a surface brightness only changes how that one number is derived (surface brightness {math}`\times` area, rather than being given directly);
it never touches the image.

```{note}
Flux has exactly one authority per source. If you supply a `brightness` *and*
an `amplitude` model parameter, both would independently imply a flux --
`Box`/`Gaussian` refuse this rather than silently pick a winner:
```

```{code-cell} ipython3
from scopesim_targets.brightness import BrightnessError

try:
    Box(spectrum=flat_spec, brightness=("V", TOTAL_FLUX),
        params=BOX_PARAMS | {"amplitude": 2.0})
except BrightnessError as err:
    print(f"raised as expected: {err}")
```

## Resolution doesn't matter -- for a sharp-edged box
A box is a *constant* profile with a hard edge, and `Box` computes each pixel's overlap with the rectangle in closed form (an exact pixel integral, not a point-sample).
So its carried flux does not depend at all on the pixel scale of the input image ({math}`\theta_\mathrm{src}`) -- there is no edge quantization to worry about:

```{code-cell} ipython3
fov = 9*u.arcsec
npixs = [9, 15, 25, 45]*u.pixel

fig, axes = plt.subplots(2, 2, figsize=(8, 6))
fluxes = []
scales = []
for ax, npix in zip(axes.flatten(), npixs):
    scale = round((fov / npix).to(u.arcsec/u.pixel), 2)
    scales.append(scale)
    grid = {"pixel_scale": scale, "width": npix.value, "height": npix.value}
    box = Box(spectrum=flat_spec, brightness=("V", TOTAL_FLUX), params=BOX_PARAMS)
    src = box.to_source(grid)
    fluxes.append(source_total_flux(src, passband()).to(u.Jy))
    im = ax.imshow(src.fields[0].data, origin="lower",
                   extent=extent_arcsec(grid), cmap="viridis")
    add_pixel_grid(ax, grid, force_major_subset=True)
    ax.set_aspect("equal")
    ax.text(1.08, 0.5, f"{scale.to_string(format="latex_inline")}\n({npix.value}\u00d7{npix.value})",
            transform=ax.transAxes, va="center", ha="left", fontsize=10, clip_on=False)

axes[1, 0].set_xlabel("arcsec")
axes[1, 1].set_xlabel("arcsec")
axes[0, 0].set_ylabel("arcsec")
axes[1, 0].set_ylabel("arcsec")
for ax in axes[0, :]:
    ax.tick_params(labelbottom=False)
for ax in axes[:, 1]:
    ax.tick_params(labelleft=False)
plt.subplots_adjust(wspace=0.5, hspace=0.15, right=0.82)
plt.show()
```

```{code-cell} ipython3
---
mystnb:
  figure:
    align: left
---
# :tags: [hide-input]
pd.DataFrame({
    "theta_src [arcsec/px]": [s.value for s in scales],
    "npix": [n.value for n in npixs],
    "flux [Jy]": [f.to_value(u.Jy) for f in fluxes],
}).set_index("theta_src [arcsec/px]").style.set_table_attributes(
"class='dataframe align-left'")
```

The rectangle simply gets sampled more finely from left to right;
the recovered flux stays flat at {math}`3.55\;\mathrm{Jy}` to many decimal places at every resolution.

## A Gaussian profile: the image only shows what's in the window
A smooth profile is different: its wings extend to infinity in principle, so whatever finite grid you render it on necessarily misses some of the light.
A `Gaussian` has a trivial closed-form total,

```{math}
\iint A\,e^{-\frac{(x-x_0)^2}{2\sigma_x^2}-\frac{(y-y_0)^2}{2\sigma_y^2}}\,
\mathrm{d}x\,\mathrm{d}y = 2\pi A\,\sigma_x\sigma_y ,
```

and `Gaussian` normalizes the image with that **analytic** total (unit amplitude here, so {math}`A_\mathrm{eff} = 2\pi\sigma_x\sigma_y`), not the sum of the rendered grid.
The stored spectrum therefore always carries the *full* stated brightness, and the image openly reports how much of it is actually inside the window.
Because the Gaussian is smooth (no sharp edge, no central cusp), that window fraction even has its own closed form:
for {math}`theta = 0` and a centered window of half-widths {math}`L_x, L_y` it is the product of error functions
{math}`\operatorname{erf}\!\bigl(L_x/\sqrt{2}\,\sigma_x\bigr)\,\operatorname{erf}\!\bigl(L_y/\sqrt{2}\,\sigma_y\bigr)`,
giving an analytic solution for the map sum.

```{code-cell} ipython3
sx, sy = 6.0, 4.0  # arcsec
gauss_grid = {"pixel_scale": 0.5*u.arcsec/u.pixel, "width": 28, "height": 28}

gaussian = Gaussian(
    spectrum=flat_spec, brightness=("V", TOTAL_FLUX),
    params={"x_stddev": sx, "y_stddev": sy},
)
src_gauss = gaussian.to_source(gauss_grid)

fraction = source_weight(src_gauss)
total_in_spectrum = Observation(src_gauss.fields[0].spectrum, passband()).effstim(u.Jy)
in_window = source_total_flux(src_gauss, passband())

# Independent analytic solution: centered window, half-width L = 28*0.5/2 = 7 arcsec.
L = gauss_grid["width"] * 0.5 / 2
frac_analytic = float(erf(L/(np.sqrt(2)*sx)) * erf(L/(np.sqrt(2)*sy)))

fig, ax = plt.subplots(figsize=(7, 6.5))
im = ax.imshow(
    src_gauss.fields[0].data * 1e3, origin="lower", extent=extent_arcsec(gauss_grid),
    cmap="magma", norm=PowerNorm(gamma=0.5),
)
add_pixel_grid(ax, gauss_grid, force_major_subset=True)
ax.add_patch(Ellipse((0, 0), 2*sx, 2*sy, fill=False, ec="cyan", ls="--", lw=1.5))
ax.text(sx, 0.5, r"  1$\sigma$", color="cyan", fontsize=9)
ax.set_xlabel("arcsec"); ax.set_ylabel("arcsec")
fig.colorbar(im, ax=ax, shrink=0.8, label=r"weight (dimensionless)  [$\times 10^{-3}$]")
plt.show()
```

```{code-cell} ipython3
:tags: [hide-input]
pd.DataFrame({
    "quantity": ["carried fraction (weight-map sum)",
                 "analytic window fraction (erf product)",
                 "flux carried by spectrum (analytic total)",
                 "in-window flux (fraction \u00d7 total)"],
    "value": [fraction, frac_analytic,
              total_in_spectrum.to_value(u.Jy), in_window.to_value(u.Jy)],
    "unit": ["(dimensionless)", "(dimensionless)", "Jy", "Jy"],
}).set_index("quantity")
```

The spectrum carries the full {math}`3.55\;\mathrm{Jy}` exactly, no matter how the window clips the profile, and the measured weight-map sum lands on the analytic erf fraction to a fraction of a percent.
Only ~70% of the flux is actually inside this particular field of view; the rest is in the wings beyond the edge.
This is deliberate: if the image were instead re-normalized so it always summed to 1 (matching the window rather than the profile),
the clipped light would be silently redistributed back among the visible pixels,
and the in-window flux would come out at the full {math}`3.55\;\mathrm{Jy}` regardless of how much of the profile was actually excluded.

## Resolution matters for a smooth profile -- until it's resolved
Unlike the box, a Gaussian's carried flux is a *quadrature approximation* of an integral, so it can in principle depend on the input pixel scale.
Sampled well relative to {math}`\sigma`, though, it barely moves, because the window fraction being measured is a geometric property,
not a discretization artifact -- and we can check that directly against the erf solution:

```{code-cell} ipython3
sig = 6.0
window_fov = 16.0*u.arcsec  # physically fixed window; divides evenly at every scale
theta_src_values = [3.2, 1.6, 0.8, 0.4, 0.2, 0.1]*u.arcsec/u.pixel

carried = []
for theta in theta_src_values:
    npix = int(round((window_fov / theta).to_value(u.pixel)))
    grid = {"pixel_scale": theta, "width": npix, "height": npix}
    g = Gaussian(
        spectrum=flat_spec,
        brightness=("V", TOTAL_FLUX),
        params={"x_stddev": sig, "y_stddev": sig},
    )
    carried.append(source_weight(g.to_source(grid)))

L = window_fov.to_value(u.arcsec) / 2
frac_oracle = float(erf(L/(np.sqrt(2)*sig))**2)

df_sweep = pd.DataFrame({
    "theta_src [arcsec/px]": theta_src_values.value,
    "px_per_sigma": sig / theta_src_values.value,
    "carried fraction": carried,
})

fig, ax = plt.subplots(figsize=(6, 4))
sns.lineplot(data=df_sweep, x="theta_src [arcsec/px]", y="carried fraction",
             marker="o", ax=ax)
ax.set_xscale("log")
ax.axhline(frac_oracle, color="gray", ls=":", label="analytic erf fraction")
ax.set_xlabel(r"$\theta_{src}$ [arcsec/pixel]")
ax.set_ylabel("carried fraction  (weight-map sum)")
ax.set_title(f"Gaussian $\\sigma$={sig}\u2033, fixed {window_fov} window")
ax.legend()
plt.show()
```

```{code-cell} ipython3
:tags: [hide-input]
df_sweep
```

The curve is essentially flat: even at ~2 pixels per {math}`\sigma` the estimate is within about a percent of the analytic fraction, and it converges rapidly as the profile is resolved.
If anything it approaches from slightly *above* -- where a Sersic's central cusp makes coarse pixel-center sampling under-count,
the Gaussian has no cusp, so the midpoint rule is accurate almost immediately.
This is exactly why {doc}`flux_engine` treats the engine's conservation guarantee (*whatever total the input image has, the detector must reproduce it*) as a separate claim from the *template's* quadrature accuracy shown here:
refining the input grid is how you fix the latter; it is not the engine's job.

## Point, box, and Gaussian: the same physical flux
Finally, put all three representations side by side: a `Star` (all its flux in a single number, no image at all),
the fully-contained `Box` from the first section, and a Gaussian compact enough, relative to its window, to be (almost) fully contained too:

```{code-cell} ipython3
star = Star(position=(0, 0), spectrum=flat_spec, brightness=("V", TOTAL_FLUX))

compact_grid = {"pixel_scale": 0.2*u.arcsec/u.pixel, "width": 256, "height": 256}
compact_gauss = Gaussian(
    spectrum=flat_spec,
    brightness=("V", TOTAL_FLUX),
    params={"x_stddev": 2.0, "y_stddev": 2.0},
)

df_compare = pd.DataFrame({
    "representation": ["Star (point)", "Box (contained)", "Gaussian (compact)"],
    "parameters": [
        "point source, no spatial extent",
        f"x_width={BOX_PARAMS['x_width']}\u2033, y_width={BOX_PARAMS['y_width']}\u2033",
        "x_stddev=y_stddev=2\u2033",
    ],
    "flux [Jy]": [
        source_total_flux(star.to_source(), passband()).to_value(u.Jy),
        source_total_flux(box_total.to_source(optical_train), passband()).to_value(u.Jy),
        source_total_flux(compact_gauss.to_source(compact_grid), passband()).to_value(u.Jy),
    ],
})

fig, ax = plt.subplots(figsize=(6, 4))
sns.barplot(data=df_compare, x="representation", y="flux [Jy]",
            hue="representation", palette="Set2", legend=False, width=0.6, ax=ax)
ax.axhline(TOTAL_FLUX.to_value(u.Jy), color="k", ls="--", lw=1, label="stated brightness")
ax.set_ylabel(f"reconstructed total flux [{TOTAL_FLUX.unit}]")
ax.set_ylim(3.50, 3.60)
ax.set_xlabel("")
ax.tick_params(axis="x", labelsize=10)
ax.legend()
for i, f in enumerate(df_compare["flux [Jy]"]):
    ax.annotate(f"{f:.4f}", (i, f), ha="center", va="bottom", fontsize=9)
plt.show()
```

```{code-cell} ipython3
:tags: [hide-input]
df_compare.set_index("representation")
```

All three sit right at {math}`3.55\;\mathrm{Jy}` -- the Gaussian case is very slightly under, because even a compact profile technically has infinite wings;
at this size relative to its window that shortfall is far below what the plot can show.
A point source has no such shortfall: with no spatial extent, there is nothing for a finite window to clip.

## Next: does the detector agree?
Everything on this page happened before a single photon reached an instrument -- it is a statement about the *stored* artifact a target produces.
{doc}`flux_engine` pushes the exact same sources through a real ScopeSim optical train and checks that what actually lands on the image plane still agrees,
across representations and across *two* independent pixel scales -- the source's own grid and the detector's.
