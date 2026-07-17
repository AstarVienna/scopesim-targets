---
file_format: mystnb
kernelspec:
  name: python3
---

# Extended-Source Flux Through a Real Instrument
{doc}`flux_construction` showed what a `Box`, a `Gaussian`, and a `Star` each store: a dimensionless weight map (or table weight) plus a spectrum that carries the absolute flux.
This page pushes those exact same sources through a real ScopeSim `OpticalTrain` and checks that what actually lands on the
image plane still agrees -- across representations, and across **two independent pixel scales** that are easy to conflate:

* {math}`\theta_\mathrm{src}` -- the pixel scale of the *source's own* image (set when you build the `Source`),
* {math}`\theta_\mathrm{inst}` -- the pixel scale of the *detector* (set on the optical train).

The PSF is switched off throughout. Flux conservation across PSF convolution
is a property of ScopeSim's convolution machinery, already covered elsewhere;
switching it off here isolates the bookkeeping this page is actually about -- resampling and detector binning -- from unrelated blurring.

## Setup
```{code-cell} ipython3
:tags: [hide-output]

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import pandas as pd
import seaborn as sns
from spextra import Spextrum

from scopesim import load_example_optical_train
from scopesim_targets.point_source import Star
from scopesim_targets.extended_source import Box

sns.set_theme(style="dark", palette="Set2")

flat_spec = Spextrum.flat_spectrum(
    5*u.ABmag, waves=np.linspace(3000, 35000, 6000)*u.AA
)
TOTAL_FLUX = 3.55*u.Jy
BOX_PARAMS = {"x_width": 5, "y_width": 3}  # arcsec
```

Like on the previous page, some helper functions for this page.
Hidden by default to reduce clutter, expand if you want to see them:

```{code-cell} ipython3
:tags: [hide-cell]
def train(theta_inst, npix):
    pixel_scale = theta_inst.to_value(u.arcsec/u.pixel)
    plate_scale = (theta_inst / (10*u.um/u.pix)).to_value(u.arcsec/u.mm)
    opt = load_example_optical_train(properties={
        "!INST.pixel_scale": pixel_scale,
        "!INST.plate_scale": plate_scale,
        "!DET.width": npix,
        "!DET.height": npix,
    })
    opt["psf"].include = False
    opt["atmospheric_radiometry"].include = False
    return opt

def detected(opt, source):
    opt.observe(source, update=True)
    data = opt.image_planes[0].data
    return data.clip(min=1e-3)

def add_pixel_grid(ax, shape, color="white", lw=.8, ls=":", alpha=.6,
                   force_major_subset=False):
    ny, nx = shape
    x_edges = np.arange(nx + 1) - 0.5
    y_edges = np.arange(ny + 1) - 0.5
    ax.set_xticks(x_edges, minor=True)
    ax.set_yticks(y_edges, minor=True)
    if force_major_subset:
        ax.set_xticks(x_edges[::max(1, round(len(x_edges) / 6))])
        ax.set_yticks(y_edges[::max(1, round(len(y_edges) / 6))])
    ax.grid(which="both" if force_major_subset else "minor",
            color=color, linewidth=lw, linestyle=ls, alpha=alpha)
    ax.tick_params(which="minor", length=0)

def src_grid(theta_src, npix):
    return {"pixel_scale": theta_src, "width": npix, "height": npix}

def point_source():
    return Star(
        position=(0, 0), spectrum=flat_spec, brightness=("V", TOTAL_FLUX)
    ).to_source()

def box_integrated(theta_src, npix, params=BOX_PARAMS):
    return Box(
        spectrum=flat_spec, brightness=("V", TOTAL_FLUX), params=params
    ).to_source(src_grid(theta_src, npix))

def box_surface_brightness(theta_src, npix, params=BOX_PARAMS):
    A_eff = params["x_width"]*u.arcsec * params["y_width"]*u.arcsec
    sb = (TOTAL_FLUX / A_eff).to(u.Jy/u.sr)
    return Box(
        spectrum=flat_spec, brightness=("V", sb), params=params
    ).to_source(src_grid(theta_src, npix))
```

## Point, 2D-integrated, and 2D-surface-brightness agree
With the source grid matched to the detector grid ({math}`\theta_\mathrm{src}\equiv\theta_\mathrm{inst}`, so there is no resampling to confound the picture),
all three representations of the same {math}`3.55\;\mathrm{Jy}` source should land on the detector as the same total:

```{code-cell} ipython3
theta = 0.5*u.arcsec/u.pixel
opt = train(theta, npix=15)

img_point = detected(opt, point_source())
img_box = detected(opt, box_integrated(theta, npix=15))
img_sb = detected(opt, box_surface_brightness(theta, npix=15))

images = {
    "point": img_point,
    "box (total flux)": img_box,
    "box (surface brightness)": img_sb,
}

fig, axes = plt.subplots(1, 3, figsize=(9, 4))
for ax, (title, img) in zip(axes, images.items()):
    im = ax.imshow(img, origin="lower", cmap="inferno", norm="log",
                   vmin=1, vmax=1e6)
    add_pixel_grid(ax, img.shape, force_major_subset=True)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])

fig.subplots_adjust(right=0.86, wspace=0.15)
pos = axes[-1].get_position()
cbar_ax = fig.add_axes([0.89, pos.y0, 0.02, pos.y1 - pos.y0])
fig.colorbar(im, cax=cbar_ax, label="detected flux")
plt.show()
```

```{code-cell} ipython3
:tags: [hide-input]
df_agree = pd.DataFrame({
    "representation": list(images.keys()),
    "detected flux": [float(np.nansum(img)) for img in images.values()],
})
df_agree["flux / point"] = df_agree["detected flux"] / df_agree["detected flux"].iloc[0]
df_agree.set_index("representation")
```

The point source appears as a single bright pixel, the box as a small rectangle spreading the same total counts over many pixels -- and all three totals agree to numerical precision.

## Two independent pixel scales
The remaining question is what happens when {math}`\theta_\mathrm{src}` and {math}`\theta_\mathrm{inst}` *don't* match.
There are two distinct cases, and they exercise different code paths:

* {math}`\theta_\mathrm{src}` finer than {math}`\theta_\mathrm{inst}` -- the engine **downsamples** the input (summation-like);
* {math}`\theta_\mathrm{src}` coarser than {math}`\theta_\mathrm{inst}` -- the engine **interpolates / upsamples** the input, the riskier direction.

A sharp-edged box used to make this subtle to check: its edges rarely land exactly on a pixel boundary,
and a rendering bug meant the detected flux could pick up a small, {math}`\theta_\mathrm{src}`-dependent bias from that misalignment.
With that fixed upstream, the detected flux is simply flat across the whole sweep, whether the engine is downsampling a finer input grid or interpolating a coarser one:

```{code-cell} ipython3
theta_inst = 0.2*u.arcsec/u.pixel
opt = train(theta_inst, 55)

theta_src_values = [0.05, 0.1, 0.2, 0.4, 0.8, 1.6]*u.arcsec/u.pixel
detected_flux = [
    float(np.nansum(detected(opt, box_integrated(theta_src, 55))))
    for theta_src in theta_src_values
]

df_scales = pd.DataFrame({
    "theta_src [arcsec/px]": theta_src_values.value,
    "detected flux": detected_flux,
})

fig, ax = plt.subplots(figsize=(5, 3))
sns.lineplot(data=df_scales, x="theta_src [arcsec/px]", y="detected flux",
             marker="o", ax=ax)
ax.set_xscale("log")
ax.set_xlabel(r"$\theta_{src}$ [arcsec/px]")
ax.set_ylabel("detected flux")

plt.tight_layout()
plt.show()
```

```{code-cell} ipython3
:tags: [hide-input]
df_scales
```

The detected flux stays flat across more than an order of magnitude in {math}`\theta_\mathrm{src}`.
Downsampling a finer input grid and interpolating a coarser one leave the total equally untouched.

## Detector scale sweep, source fixed
Fixing the source (a point, so there is no {math}`\theta_\mathrm{src}` at all) and sweeping only {math}`\theta_\mathrm{inst}` should leave the detected flux untouched:
a point source's photons land in whichever pixel they land in, but the total is unaffected by how finely the detector happens to be divided up.

```{code-cell} ipython3
theta_inst_values = [100, 200, 400, 800]*u.mas/u.pixel
fluxes = []
for theta_inst in theta_inst_values:
    opt = train(theta_inst.to(u.arcsec/u.pixel), npix=25)
    fluxes.append(float(np.nansum(detected(opt, point_source()))))

df_theta_inst = pd.DataFrame({
    "theta_inst [mas/px]": theta_inst_values.value,
    "detected flux": fluxes,
})

fig, ax = plt.subplots(figsize=(5, 3))
sns.lineplot(data=df_theta_inst, x="theta_inst [mas/px]", y="detected flux",
             marker="o", ax=ax)
ax.set_xlabel(r"$\theta_{inst}$ [mas/px]")
ax.set_ylabel("detected flux")
ax.set_title("point source, detector scale swept")
plt.show()
```

```{code-cell} ipython3
:tags: [hide-input]
df_theta_inst
```

## An unresolved extended source is a point source
As a last check: an extended source much smaller than a detector pixel should be indistinguishable from a point source of the same brightness --
the *unresolved* limit of the point/box agreement shown above, this time stressing sub-pixel placement and detector binning rather than the flux conservation itself.
A 0.5 arcsec box, well sampled on its own fine grid, on a large 2 arcsec detector pixel:

```{code-cell} ipython3
theta_inst = 2.0*u.arcsec/u.pixel
opt = train(theta_inst, npix=15)

tiny_box = box_integrated(
    0.5*u.arcsec/u.pixel, npix=15, params={"x_width": 0.5, "y_width": 0.5},
)
img_tiny = detected(opt, tiny_box)
img_pt = detected(opt, point_source())

images_unres = {"point": img_pt, "0.4\u2033 box (unresolved)": img_tiny}

fig, axes = plt.subplots(1, 2, figsize=(7, 4))
for ax, (title, img) in zip(axes, images_unres.items()):
    im = ax.imshow(img, origin="lower", cmap="inferno", norm="log",
                   vmin=1, vmax=1e6)
    add_pixel_grid(ax, img.shape, force_major_subset=True)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])

fig.subplots_adjust(right=0.84, wspace=0.15)
pos = axes[-1].get_position()
cbar_ax = fig.add_axes([0.87, pos.y0, 0.02, pos.y1 - pos.y0])
fig.colorbar(im, cax=cbar_ax, label="detected flux")
plt.show()

print(f"point:          {float(np.nansum(img_pt)):>6.2f}")
print(f"unresolved box: {float(np.nansum(img_tiny)):>6.2f}")
```

Both land in the same single detector pixel, with the same total -- exactly what "unresolved" should mean.

## What isn't covered here
This page (like {doc}`flux_construction`) is 2D-only.
A datacube adds a third axis and one more subtlety -- the pixel solid angle has to be divided out when forming the cube and multiplied back in on collapse --
which is a large enough topic to deserve its own page later.
The imaging invariant demonstrated above (the same physical source gives the same detected flux, regardless of representation or pixel scale) is the same one the cube case extends,
just with {math}`\Omega` and the spectral bin width {math}`\Delta\lambda` both entering exactly once.
