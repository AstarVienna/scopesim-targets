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

# The Ophiuchus background spectrum
This page documents `ophiuchus_background_v1.fits`, the surface-brightness spectrum of the diffuse near-infrared emission of the ρ Ophiuchi cloud L1688.
It is the faint counterpart to the Orion spectrum in {doc}`background_sed_orion`, with the same wavelength grid, units and file layout,
so the two can be swapped in a scene ({doc}`nebular_backgrounds`).

The physics is completely different from Orion.
ρ Oph has no conspicuous H II region: UV light from its hottest stars heats and dissociates the exposed cloud layers, but does not ionize hydrogen {cite:p}`Casassus2008`.
Its diffuse near-infrared light is attributed instead to UV-fluorescent H₂ ro-vibrational lines.
As for Orion, this page reads the shipped file, recomputes its diagnostics from `build_ophiuchus_sed.py`, and asserts that they agree.

## Setup

```{code-cell} ipython3
:tags: [hide-output]

import os
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io import fits
from astropy.table import QTable
import pandas as pd
import seaborn as sns
from IPython.display import Markdown
from myst_nb import glue

# Folder holding the SED files, build scripts and nebular_backgrounds.py.
BKG_DIR = Path(os.environ.get("SCOPESIM_BACKGROUNDS_DIR", "../clusterbkg")).resolve()
sys.path.insert(0, str(BKG_DIR))

import build_orion_nebula_sed as ori
import build_ophiuchus_sed as oph
import nebular_backgrounds as nb

SED_FILE = BKG_DIR / nb.OPHIUCHUS_SED
sed = QTable.read(SED_FILE, hdu="SPECTRUM")
header = fits.getheader(SED_FILE, "SPECTRUM")
orion_sed = QTable.read(BKG_DIR / nb.ORION_SED, hdu="SPECTRUM")
h2_lines = QTable.read(BKG_DIR / "ophiuchus_h2_linelist_v1.ecsv")
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
SVO = Path(__import__("importlib.util").util.find_spec("spextra").origin).parent / "data" / "svo"
BANDS = {"J": "2MASS/2MASS.J", "H": "2MASS/2MASS.H", "Ks": "2MASS/2MASS.Ks"}
WAVE = sed["WAVELENGTH"].value
SB_UNIT = sed["FLUX"].unit

def diag(key):
    """Diagnostic value stored in the FITS header by the build."""
    return header[f"DIAG {key}"]

def cite(key):
    return "design choice" if key is None else f"{{cite:p}}`{key}`"

def md_table(df):
    """DataFrame -> MyST Markdown table (so citation roles inside resolve)."""
    head = "| " + " | ".join(df.columns) + " |"
    rule = "|" + "---|" * len(df.columns)
    body = ["| " + " | ".join(str(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return Markdown("\n".join([head, rule, *body]))

def rebin(lam, flux, r_grid=3000, r_out=300):
    """Flux-conserving rebin to R~r_grid, then Gaussian smoothing to R~r_out."""
    from scipy.ndimage import gaussian_filter1d
    edges = np.exp(np.arange(np.log(lam[0]), np.log(lam[-1]), 1 / r_grid))
    cum = np.concatenate([[0], np.cumsum(0.5 * (flux[1:] + flux[:-1]) * np.diff(lam))])
    binned = np.diff(np.interp(edges, lam, cum)) / np.diff(edges)
    return np.sqrt(edges[1:] * edges[:-1]), gaussian_filter1d(binned, r_grid / r_out / 2.355)

def band_mag(wave, flux, band):
    lam_t, trans = ori.svo_curve(SVO / BANDS[band])
    f0 = ori.p(f"F0_{band}", u.erg / u.s / u.cm**2 / u.AA)
    return -2.5 * np.log10(ori.band_mean_flam(wave, flux, lam_t, trans) / f0)

def band_mjysr(wave, flux, band):
    lam_t, trans = ori.svo_curve(SVO / BANDS[band])
    return ori.band_mean_inu_mjysr(wave, flux, lam_t, trans)[0]
```

## The anchor: diffuse 2MASS emission on ρ Oph W
The only published absolute near-infrared measurement of this emission is the analysis by {cite:t}`Casassus2008` of the 2MASS images of the ρ Oph W photodissociation region.
That region is illuminated by the B2IV+B3IV binary HD 147889, with a UV field of χ ≈ 400 {cite:p}`Habart2003`.
The numbers they report are the inputs for this spectrum:

```{code-cell} ipython3
:tags: [remove-input]

obs = {
    "J": oph.p("I_Ks_obs", u.MJy / u.sr) * oph.p("ratio_J_K").value,
    "H": oph.p("I_Ks_obs", u.MJy / u.sr) * oph.p("ratio_H_K").value,
    "Ks": oph.p("I_Ks_obs", u.MJy / u.sr),
}
noise = oph.p("I_noise_2MASS", u.MJy / u.sr)
glue("oph_noise", float(noise), display=False)
glue("oph_10s1_diluted", float(oph.p("I_10S1_Ks_diluted", u.MJy / u.sr)), display=False)
pd.DataFrame({"observed [MJy/sr]": obs}).T.style.format("{:.2f}").set_table_styles(tbl_sty)
```

These are peak intensities above the local background, with a mosaic noise of {glue:text}`oph_noise:.2f` MJy/sr, so the measurement has S/N ≈ 3 in Ks and less in J and H.
{cite:t}`Casassus2008` also estimate that the H₂ 1-0 S(1) line alone contributes {glue:text}`oph_10s1_diluted:.1f` MJy/sr to Ks, based on the 1-0 S(1) map of {cite:t}`Habart2003`.
They conclude that ro-vibrational H₂ can supply the whole Ks flux.

## The H₂ fluorescence model
Line data come from the {cite:t}`Roueff2019` calculations as shipped with pdrtpy {cite:p}`Pound2023`.
Level populations follow a two-temperature description,

```{math}
:label: oph-h2pop

N(v, J) \;\propto\; g(v, J)\, \exp\!\left[-\frac{E(v,0)}{k T_\mathrm{vib}}\right]
\exp\!\left[-\frac{E(v,J) - E(v,0)}{k T_\mathrm{rot}}\right] ,
```

with line intensities {math}`I_{ul} \propto N_u A_{ul}\, h\nu_{ul}`. The two temperatures are set as follows:

* {math}`T_\mathrm{vib}` is fixed by the radiative-pumping ratio 2-1 S(1)/1-0 S(1) of {cite:t}`BlackvanDishoeck1987`, as quoted by {cite:t}`LePetit2006`.
* {math}`T_\mathrm{rot}` is the 1-0 S excitation temperature of the {cite:t}`Peeters2024` atomic-PDR template. This is an Orion Bar value, as no ρ Oph W excitation diagram exists.

The spectrum is normalized so that 1-0 S(1) alone gives the Casassus et al. Ks level.
The cell below recomputes the model and checks it against the file:

```{code-cell} ipython3
t, e0 = oph.h2_table()
t_rot = oph.p("T_rot", u.K)
t_vib, i10 = oph.solve_t_vib(t, e0, t_rot)
rel = oph.h2_relative_intensities(t, e0, t_vib, t_rot)
rel /= rel[i10]
lam_all = 1e8 / np.asarray(t["sigma"], float)
ks_curve = ori.svo_curve(SVO / BANDS["Ks"])
i_10s1 = oph.p("I_10S1_Ks_diluted", u.MJy / u.sr) / oph.band_mean_inu_of_lines(
    lam_all[[i10]], np.array([1.0]), *ks_curve)
share = (oph.band_mean_inu_of_lines(lam_all[[i10]], np.array([i_10s1]), *ks_curve)
         / oph.band_mean_inu_of_lines(lam_all, rel * i_10s1, *ks_curve))

np.testing.assert_allclose(t_vib, diag("T_vib_K"), rtol=1e-6)
np.testing.assert_allclose(i_10s1, diag("I_10S1_erg_s_cm2_sr"), rtol=1e-6)
np.testing.assert_allclose(share, diag("share_10S1_in_Ks_model"), rtol=1e-6)

glue("oph_tvib", float(t_vib), display=False)
glue("oph_trot", float(t_rot), display=False)
glue("oph_i10s1", float(i_10s1 * 1e5), display=False)
glue("oph_share", float(100 * share), display=False)
glue("oph_share_bvd", float(100 * oph.p("share_10S1_BvD87").value), display=False)
glue("oph_nlines", int(diag("n_lines_H2")), display=False)
glue("oph_fwhm", float(diag("FWHM_H2_kms")), display=False)
```

This gives {math}`T_\mathrm{vib}` = {glue:text}`oph_tvib:.0f` K at {math}`T_\mathrm{rot}` = {glue:text}`oph_trot:.0f` K,
and I(1-0 S(1)) = {glue:text}`oph_i10s1:.2f` × 10⁻⁵ erg s⁻¹ cm⁻² sr⁻¹.
The file keeps {glue:text}`oph_nlines` lines, each with the pure thermal width of H₂ at the H₂-layer gas temperature of {cite:t}`Habart2003`, {glue:text}`oph_fwhm:.1f` km/s.

Two checks of the model against independent values:

* 1-0 S(1) carries {glue:text}`oph_share:.0f` % of the H₂ flux in Ks,
  against {glue:text}`oph_share_bvd:.1f` % in the {cite:t}`BlackvanDishoeck1987` model quoted by {cite:t}`Casassus2008`.
* The K-band line ratios are compared below with the extinction-corrected JWST ratios of the {cite:t}`Peeters2024` H II-region template (their Table 3).

```{code-cell} ipython3
:tags: [hide-input]

# Peeters et al. (2024), Table 3, H II-region template, extinction-corrected
# (foreground formalism), in 1e-5 erg s-1 cm-2 sr-1:
p24 = {"1-0 S(0)": 4.92, "1-0 S(1)": 13.59, "1-0 S(2)": 5.32, "1-0 S(3)": 10.32,
       "1-0 Q(3)": 10.05, "1-0 Q(4)": 3.02, "1-0 Q(5)": 5.83, "1-0 Q(7)": 3.28}
model = {str(n): float(v) for n, v in zip(h2_lines["transition"], h2_lines["ratio_to_10S1"])}
pd.DataFrame({
    "model / 1-0 S(1)": {k: model[k] for k in p24},
    "Peeters et al. 2024 / 1-0 S(1)": {k: v / p24["1-0 S(1)"] for k, v in p24.items()},
}).style.format("{:.2f}").set_table_styles(tbl_sty)
```

The ratios match to about 10 % for most lines.
The para line 1-0 S(0) is about 40 % low, consistent with an ortho/para ratio below the value of 3 built into the statistical weights;
{cite:t}`Habart2003` find ≈ 1 in the v = 0 levels of ρ Oph W.

## The spectrum

```{code-cell} ipython3
---
tags: [hide-input]
mystnb:
  figure:
    name: fig:oph-sed
    caption: Ophiuchus background spectrum. Top, the fluorescent H₂ lines and the
      scattered-light upper limit at native sampling; bottom, compared with the Orion
      spectrum at R ~ 300.
---
fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(9, 7), sharex=True, layout="compressed")
ax0.plot(WAVE / 1e4, sed["H2_LINES"].value, color="C1", lw=0.4, label="H2_LINES = FLUX")
ax0.plot(WAVE / 1e4, sed["SCATTERED_UL"].value, color="C3", lw=1, label="SCATTERED_UL (upper limit)")
ax0.set_yscale("log"); ax0.set_ylim(1e-21, 3e-15); ax0.legend(fontsize=8, loc="upper right")
ax0.set_ylabel(f"[{SB_UNIT.to_string('latex_inline')}]")
ax0.set_title("native sampling", fontsize=9)

for table, col, color, label in ((orion_sed, "FLUX", "k", "Orion FLUX (intrinsic)"),
                                 (sed, "FLUX", "C1", "Ophiuchus FLUX (observed)"),
                                 (sed, "SCATTERED_UL", "C3", "Ophiuchus SCATTERED_UL")):
    cen, smooth = rebin(table["WAVELENGTH"].value, table[col].value)
    ax1.plot(cen / 1e4, smooth, color=color, lw=1, ls="--" if "UL" in col else "-", label=label)
for band, (lo, hi) in {"J": (1.11, 1.36), "H": (1.50, 1.80), "Ks": (1.99, 2.31), "L'": (3.43, 4.13)}.items():
    ax1.axvspan(lo, hi, color="0.92", zorder=0)
    ax1.text((lo + hi) / 2, 2e-14, band, ha="center")
ax1.set_yscale("log"); ax1.set_ylim(1e-20, 1e-13); ax1.legend(fontsize=8, loc="upper right", ncol=3)
ax1.set_xlabel("vacuum wavelength [µm]"); ax1.set_ylabel(ax0.get_ylabel())
ax1.set_title("Orion vs. Ophiuchus, flux-conserving rebin to R ~ 300", fontsize=9)
plt.show()
```

{numref}`fig:oph-sed` shows the result: a pure line forest, with no continuum,
about two orders of magnitude below Orion across J to Ks.

## Broadband surface brightness
Synthetic 2MASS photometry of `FLUX`, as for Orion {cite:p}`Rodrigo2012,Skrutskie2006,Cohen2003`, compared with the measurement.
The deviation uses the Ks noise as representative for all three bands.

```{code-cell} ipython3
:tags: [hide-input]

rows = {}
for band in BANDS:
    model_i = band_mjysr(WAVE, sed["FLUX"].value, band)
    model_m = band_mag(WAVE, sed["FLUX"].value, band)
    rows[band] = {
        "observed [MJy/sr]": obs[band],
        "model [MJy/sr]": model_i,
        "deviation": f"{(model_i - obs[band]) / noise:+.1f} σ",
        "model [mag/arcsec²]": model_m,
        "observed [mag/arcsec²]": model_m - 2.5 * np.log10(obs[band] / model_i),
        "Orion [mag/arcsec²]": band_mag(orion_sed["WAVELENGTH"].value, orion_sed["FLUX"].value, band),
        "SCATTERED_UL [MJy/sr]": band_mjysr(WAVE, sed["SCATTERED_UL"].value, band),
    }
phot = pd.DataFrame(rows).T
glue("oph_ks", float(phot.loc["Ks", "model [mag/arcsec²]"]), display=False)
for band in BANDS:
    glue(f"oph_sig_{band}", float((phot.loc[band, "model [MJy/sr]"] - obs[band]) / noise), display=False)
glue("oph_ks_obs", float(phot.loc["Ks", "observed [mag/arcsec²]"]), display=False)
glue("oph_dks", float(phot.loc["Ks", "model [mag/arcsec²]"] - phot.loc["Ks", "Orion [mag/arcsec²]"]), display=False)

# The target definition must use the level of the shipped SED:
ks_target = float(nb.OPHIUCHUS_BRIGHTNESS[1].split()[0])
assert abs(phot.loc["Ks", "model [mag/arcsec²]"] - ks_target) < 0.005, \
    f"nebular_backgrounds.OPHIUCHUS_BRIGHTNESS ({ks_target}) is out of date"

phot.style.format(precision=2).set_table_styles(tbl_sty)
```

The model gives Ks = {glue:text}`oph_ks:.2f` mag/arcsec², against {glue:text}`oph_ks_obs:.2f` from 2MASS,
which puts it {glue:text}`oph_dks:.1f` mag/arcsec² below Orion.

```{admonition} The J–H colour tension
:class: warning

Pure fluorescent H₂ deviates by {glue:text}`oph_sig_J:+.1f` σ in J and {glue:text}`oph_sig_Ks:+.1f` σ in Ks, but by {glue:text}`oph_sig_H:+.1f` σ in H.
UV-pumped cascades emit more in J than in H {cite:p}`FosterGoodman2006`, so fluorescent H₂ is intrinsically *blue* in J–H, while the quoted ρ Oph W colours are *red*.
{cite:t}`Casassus2008` note that these colours are also inconsistent with scattered light, and discuss a continuum from very small grains.
No such red continuum has a referenced level, so none is added; treat H as uncertain by about a factor of 3.
```

### Scattered starlight
`SCATTERED_UL` is a grey-albedo scattering of the HD 147889 binary, modelled as the two blackbodies of the {cite:t}`Casassus2008` TLUSTY fit.
It is normalized so that scattered light *alone* reproduces the observed J intensity, which it cannot exceed.
Together with the colour argument above, that makes it a strict upper limit, and it is not part of `FLUX`.

## Strongest lines

```{code-cell} ipython3
:tags: [hide-input]

top = np.argsort(np.asarray(h2_lines["intensity"].value))[::-1][:10]
pd.DataFrame({
    "λ_vac [µm]": h2_lines["lambda_vac"][top].to_value(u.um),
    "I [erg s⁻¹ cm⁻² arcsec⁻²]": h2_lines["intensity"][top].value,
    "I / I(1-0 S(1))": h2_lines["ratio_to_10S1"][top],
}, index=list(h2_lines["transition"][top])).style.format(
    {"λ_vac [µm]": "{:.5f}", "I [erg s⁻¹ cm⁻² arcsec⁻²]": "{:.3e}", "I / I(1-0 S(1))": "{:.3f}"}
).set_table_styles(tbl_sty)
```

## Parameters
Every number entering the build, read from `build_ophiuchus_sed.PARAMS`:

```{code-cell} ipython3
---
tags: [remove-input]
mystnb:
  markdown_format: myst
---
md_table(pd.DataFrame(
    [(f"`{k}`", f"{getattr(v, 'value', v):.6g} {un.to_string() or ''}", cite(ref), note)
     for k, (v, un, ref, note) in oph.PARAMS.items()],
    columns=["parameter", "value", "source", "where / what"],
))
```

## Limitations
* **The anchor is weak.** At S/N ≈ 3 in Ks and less in J and H, treat the absolute level as uncertain by about a factor of 2.
* **It holds for one position.** {cite:t}`Casassus2008` detect extended diffuse 2MASS emission only on the ρ Oph W ridge, besides compact nebulosity around HD 147889, S1 and SR 3.
  The 31 GHz emission east of SR 3 is undetectable in 2MASS, so cloud-interior fields are fainter than the noise level.
  Use this level as an upper envelope.
* **L/M is a strong lower limit.** PAH emission and the dust continuum are omitted; the ambient 3.29 µm emission at the cloud periphery is PAH-dominated {cite:p}`Dartois2025`.
* **The excitation parameters come from Orion Bar data** and a generic fluorescence model, not from ρ Oph W itself.
* **The faint H II region around HD 147889** is not included.

## Rebuilding the file
This page never writes data. To regenerate the FITS file and line list, run the build script; it needs pdrtpy and PyNeb, as does this page.

```bash
python build_ophiuchus_sed.py
```

## References
```{bibliography}
:filter: docname in docnames
```
