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

# The Orion Nebula background spectrum
This page documents `orion_nebula_background_v1.fits`, the intrinsic (unreddened) surface-brightness spectrum of the Orion Nebula (M42) H II region.
It is the bright, line- and continuum-rich counterpart to the faint Ophiuchus spectrum in {doc}`background_sed_ophiuchus`;
how both are turned into extended background targets for embedded-cluster simulations is shown in {doc}`nebular_backgrounds`.

The file covers 0.97--5.40 µm in vacuum wavelengths, in units of {math}`\mathrm{erg\,s^{-1}\,cm^{-2}\,\AA^{-1}\,arcsec^{-2}}`.
It is produced by `build_orion_nebula_sed.py`, which keeps every input number in a `PARAMS` dictionary together with its source.
This page does not rebuild the file. It reads the shipped file, recomputes its key diagnostics from the build code,
and **asserts that the two agree**, so the documentation fails to build if the data and the code drift apart.

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
import pyneb as pn
from dust_extinction.parameter_averages import G23
from IPython.display import Markdown
from myst_nb import glue

# Folder holding the SED files, build scripts and nebular_backgrounds.py.
BKG_DIR = Path(os.environ.get("SCOPESIM_BACKGROUNDS_DIR", "../clusterbkg")).resolve()
sys.path.insert(0, str(BKG_DIR))

import build_orion_nebula_sed as ori
import nebular_backgrounds as nb

SED_FILE = BKG_DIR / nb.ORION_SED
sed = QTable.read(SED_FILE, hdu="SPECTRUM")
header = fits.getheader(SED_FILE, "SPECTRUM")
lines = QTable.read(BKG_DIR / "orion_nebula_linelist_v1.ecsv")
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
COMPONENTS = ("H_LINES", "HE_LINES", "NEB_CONT", "PAH")
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

def band_mag(flux, band):
    lam_t, trans = ori.svo_curve(SVO / BANDS[band])
    f0 = ori.p(f"F0_{band}", u.erg / u.s / u.cm**2 / u.AA)
    return -2.5 * np.log10(ori.band_mean_flam(WAVE, flux, lam_t, trans) / f0)

def band_fraction(component, band):
    lam_t, trans = ori.svo_curve(SVO / BANDS[band])
    total = ori.band_mean_flam(WAVE, sed["FLUX"].value, lam_t, trans)
    return ori.band_mean_flam(WAVE, sed[component].value, lam_t, trans) / total
```

## What the spectrum contains
The Orion Nebula is a blister H II region ionized mainly by θ¹ Ori C.
In the near-infrared its light is dominated by the ionized gas, so the spectrum is built from the following components:

* **H I recombination lines:** Case B emissivities of {cite:t}`Storey1995`, computed with PyNeb {cite:p}`Luridiana2015`
  for the Paschen through n→8 series, with upper levels up to n = 50.
* **He I recombination lines:** emissivities of {cite:t}`Porter2012`, corrected in {cite:t}`Porter2013`,
  with the He⁺ abundance of {cite:t}`Esteban2004`.
  The 2.058 µm line is the exception: it is tied to the 1.70 µm line with the ratio observed by
  {cite:t}`Peeters2024`, which confirms {cite:t}`Marconi1998`, because Case B overpredicts it.
* **Nebular continuum:** free--bound {cite:p}`Ercolano2006`, free--free {cite:p}`Storey1991` and two-photon emission.
* **PAH emission at 3.3--3.6 µm:** the band profiles and intensities that {cite:t}`Peeters2024` fit with PAHFIT {cite:p}`Smith2007`.
  Following {cite:t}`Chown2024`, this emission comes from the face-on PDR *behind* the ionized gas.
* **Scattered starlight** (column `SCATTERED_UL`, not included in `FLUX`): an upper limit, discussed in {ref}`sec:orion-scattered`.

All gas components scale with the emission measure.
For a line or continuum emissivity {math}`\varepsilon_\lambda = 4\pi j_\lambda / (n_\mathrm{e} n_\mathrm{p})`,

```{math}
:label: orion-em

I_\lambda \;=\; \frac{\varepsilon_\lambda(T_\mathrm{e}, n_\mathrm{e})}{4\pi}\,\mathrm{EM},
\qquad \mathrm{EM} \;=\; \int n_\mathrm{e}\, n_\mathrm{p}\, \mathrm{d}l ,
```

evaluated at the electron temperature, density and He⁺ abundance of {cite:t}`Esteban2004`.
Line profiles are Gaussians with the pure thermal width at {math}`T_\mathrm{e}`.

```{code-cell} ipython3
fwhm = pd.DataFrame({
    "FWHM [km/s]": [diag("FWHM_HI_kms"), diag("FWHM_HeI_kms")],
    "lines": [int(diag("n_lines_HI")), int(diag("n_lines_HeI"))],
}, index=["H I", "He I"])
fwhm.style.format({"FWHM [km/s]": "{:.1f}"}).set_table_styles(tbl_sty)
```

## Absolute level
{cite:t}`Peeters2024` derive an emission measure from Brγ for their H II-region template, assuming Case B at their own adopted temperature and density.
Inverting eq. {eq}`orion-em` with *their* conditions recovers the dereddened Brγ intensity they measured.
The same Brγ at our adopted conditions then fixes the emission measure for all other lines and the continuum.
The cell below recomputes this from the build code and checks it against the file:

```{code-cell} ipython3
h1 = pn.RecAtom("H", 1)
em_eff, i_brg = ori.emission_measure_at_adopted_conditions(h1)
em_eff_pc = em_eff / (1 * u.pc * u.cm**-6).to_value(u.cm**-5)

np.testing.assert_allclose(em_eff_pc, diag("EM_eff_pc_cm6"), rtol=1e-6)
np.testing.assert_allclose(i_brg, diag("I_Brg_ref_erg_s_cm2_sr"), rtol=1e-6)

glue("orion_em_ref", float(ori.p("EM_ref").value / 1e6), display=False)
glue("orion_em_eff", float(em_eff_pc / 1e6), display=False)
glue("orion_i_brg", float(i_brg * 1e3), display=False)
print(f"EM = {em_eff_pc:.3e} pc cm-6,  I(Brg) = {i_brg:.3e} erg s-1 cm-2 sr-1")
```

The published emission measure, {glue:text}`orion_em_ref:.1f` × 10⁶ pc cm⁻⁶, corresponds to
I(Brγ) = {glue:text}`orion_i_brg:.2f` × 10⁻³ erg s⁻¹ cm⁻² sr⁻¹.
At our adopted conditions, the same Brγ needs EM = {glue:text}`orion_em_eff:.2f` × 10⁶ pc cm⁻⁶.

This anchor can be checked independently against the `orion` template bundled with speXtra (J. R. Walsh, 1996).
That template's dereddened total Hβ comes from 84.2 GHz free--free emission {cite:p}`Gordon1986`.
Spread over the 5′ Huygens Region {cite:p}`ODell2010`, it gives a mean Hβ surface brightness that can be compared with ours:

```{code-cell} ipython3
i_hb = h1.getEmissivity(ori.p("T_e", u.K), ori.p("n_e", u.cm**-3), lev_i=4, lev_j=2) * em_eff / (4 * np.pi)
i_hb_arcsec2 = i_hb / ori.SR_TO_ARCSEC2
area = np.pi * (ori.p("D_Huygens", u.arcsec) / 2) ** 2
i_hb_walsh = ori.p("F_Hb_total", u.erg / u.s / u.cm**2) / area

np.testing.assert_allclose(i_hb_arcsec2, diag("I_Hb_erg_s_cm2_arcsec2"), rtol=1e-6)
glue("orion_hb_model", float(i_hb_arcsec2 * 1e12), display=False)
glue("orion_hb_walsh", float(i_hb_walsh * 1e12), display=False)
glue("orion_hb_dev", float(100 * (i_hb_arcsec2 / i_hb_walsh - 1)), display=False)
glue("orion_ew_obs", float(ori.p("EW_obs_Hb", u.AA)), display=False)
```

The model gives {glue:text}`orion_hb_model:.2f` × 10⁻¹² erg s⁻¹ cm⁻² arcsec⁻², the radio route {glue:text}`orion_hb_walsh:.2f` × 10⁻¹²:
two independent anchors agree to {glue:text}`orion_hb_dev:+.0f` %.

## Validating the nebular continuum
{cite:t}`ODell2010` quote the dust-free equivalent width of Hβ for two densities at 9000 K, i.e. the ratio of the Hβ line to the pure nebular continuum.
Recomputing it with the continuum used here tests the continuum emissivities directly:

```{code-cell} ipython3
cont = pn.Continuum()
odell = {200: 1600.0, 6000: 1800.0}  # EW(Hbeta) [A] at 9000 K, O'Dell & Harris (2010), Sect. 5.3.2
ew = {n: 1 / cont.get_continuum(tem=9000.0, den=float(n), He1_H=ori.p("y_He+").value,
                                 wl=np.array([4861.0]), HI_label="4_2")[0] for n in odell}
pd.DataFrame({
    "O'Dell & Harris 2010 [Å]": odell,
    "this continuum [Å]": ew,
    "difference": {n: f"{100 * (ew[n] / odell[n] - 1):+.0f} %" for n in odell},
}).rename_axis("n_e [cm⁻³]").style.format(precision=0).set_table_styles(tbl_sty)
```

## The spectrum

```{code-cell} ipython3
---
tags: [hide-input]
mystnb:
  figure:
    name: fig:orion-sed
    caption: Orion Nebula background spectrum. Top, the components at native sampling;
      bottom, the total rebinned to R ~ 300, with the 2MASS and L′ bands shaded.
---
fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(9, 7), sharex=True, layout="compressed")
um = WAVE / 1e4
for col, color in (("NEB_CONT", "C0"), ("PAH", "C2"), ("SCATTERED_UL", "C3")):
    ax0.plot(um, sed[col].value, color=color, lw=1, label=col)
ax0.plot(um, sed["FLUX"].value, color="k", lw=0.4, label="FLUX (total)")
ax0.set_yscale("log"); ax0.set_ylim(1e-19, 3e-15); ax0.legend(fontsize=8, ncol=2)
ax0.set_ylabel(f"[{SB_UNIT.to_string('latex_inline')}]")
ax0.set_title("native sampling", fontsize=9)

for col, color, label in (("FLUX", "k", "FLUX"), ("NEB_CONT", "C0", "nebular continuum")):
    cen, smooth = rebin(WAVE, sed[col].value)
    ax1.plot(cen / 1e4, smooth, color=color, lw=1, label=label)
for band, (lo, hi) in {"J": (1.11, 1.36), "H": (1.50, 1.80), "Ks": (1.99, 2.31), "L'": (3.43, 4.13)}.items():
    ax1.axvspan(lo, hi, color="0.9", zorder=0)
    ax1.text((lo + hi) / 2, 2e-15, band, ha="center")
ax1.set_yscale("log"); ax1.set_ylim(5e-18, 5e-15); ax1.legend(fontsize=8, loc="lower left")
ax1.set_xlabel("vacuum wavelength [µm]"); ax1.set_ylabel(ax0.get_ylabel())
ax1.set_title("flux-conserving rebin to R ~ 300", fontsize=9)
plt.show()
```

{numref}`fig:orion-sed` shows the components at native sampling (top) and the total at a resolution typical of broadband work (bottom).
The steps in the nebular continuum are the hydrogen recombination edges at 1.46, 2.28, 3.28 and 4.47 µm,
where the free--bound contribution of one more level switches off.
Just redward of each edge, the converging high-n lines of the same series blend into a pseudo-continuum.
At R ~ 300 they leave a ripple wherever the line spacing is close to the smoothing width.

## Broadband surface brightness
The table gives the synthetic 2MASS photometry of `FLUX`, using the response curves bundled with speXtra
{cite:p}`Rodrigo2012,Skrutskie2006` and the zero points of {cite:t}`Cohen2003`.
"Reddened" applies the foreground screen A_V = 0.9--1.9 of {cite:t}`Peeters2024`, using the G23 curve with R_V = 5.5.
The component columns give each component's share of the band flux; the scattered-light upper limit would come *on top*.

```{code-cell} ipython3
:tags: [hide-input]

ext = G23(Rv=5.5)
in_range = (WAVE > 1e4 / ext.x_range[1]) & (WAVE < 1e4 / ext.x_range[0])

def reddened_mag(av, band):
    att = np.ones_like(WAVE)
    att[in_range] = ext.extinguish(WAVE[in_range] * u.AA, Av=av)
    return band_mag(sed["FLUX"].value * att, band)

rows = {}
for band in BANDS:
    lam_t, trans = ori.svo_curve(SVO / BANDS[band])
    rows[band] = {
        "mag/arcsec² (intrinsic)": band_mag(sed["FLUX"].value, band),
        "reddened": f"{reddened_mag(0.9, band):.2f} – {reddened_mag(1.9, band):.2f}",
        "MJy/sr": ori.band_mean_inu_mjysr(WAVE, sed["FLUX"].value, lam_t, trans)[0],
        **{c: band_fraction(c, band) for c in COMPONENTS},
        "+ SCATTERED_UL": band_fraction("SCATTERED_UL", band),
    }
phot = pd.DataFrame(rows).T
glue("orion_ks", float(phot.loc["Ks", "mag/arcsec² (intrinsic)"]), display=False)
glue("orion_ks_cont", float(100 * phot.loc["Ks", "NEB_CONT"]), display=False)

# The target definition must use the level of the shipped SED:
ks_target = float(nb.ORION_BRIGHTNESS[1].split()[0])
assert abs(phot.loc["Ks", "mag/arcsec² (intrinsic)"] - ks_target) < 0.005, \
    f"nebular_backgrounds.ORION_BRIGHTNESS ({ks_target}) is out of date"

pct = {c: "{:.0%}" for c in (*COMPONENTS, "+ SCATTERED_UL")}
phot.style.format({"mag/arcsec² (intrinsic)": "{:.2f}", "MJy/sr": "{:.0f}", **pct}).set_table_styles(tbl_sty)
```

In Ks, the nebular continuum carries {glue:text}`orion_ks_cont:.0f` % of the total of {glue:text}`orion_ks:.2f` mag/arcsec²; the lines are a minority even there.
The L/M band is a lower limit (see {ref}`sec:orion-limits`).

(sec:orion-scattered)=
### Scattered starlight
The observed EW(Hβ) of {glue:text}`orion_ew_obs:.0f` Å near θ¹ Ori C {cite:p}`ODell2010` is much smaller than the dust-free value,
because Trapezium light is scattered by dust in the PDR behind the nebula.
`SCATTERED_UL` models it as a θ¹ Ori C blackbody {cite:p}`SimonDiaz2006` scattered with a grey albedo,
normalized to reproduce that observed EW.
{cite:t}`ODell2010` find the scattered light is bluer than the illuminating star, so a grey albedo *over*-predicts the near-infrared part.
The true contribution lies between zero and the last column of the table above, which is why it is not part of `FLUX`.

## Key lines

```{code-cell} ipython3
:tags: [hide-input]

key = ["Pag", "HeI10830.0", "Pab", "Br10", "HeI17003.0", "HeI20581.0", "Brg", "Brb", "Pfg", "Bra", "Pfb"]
sel = lines[[n in key for n in lines["name"]]]
pd.DataFrame({
    "transition": [t or "—" for t in sel["transition"]],
    "λ_vac [µm]": sel["lambda_vac"].to_value(u.um),
    "I [erg s⁻¹ cm⁻² arcsec⁻²]": sel["intensity"].value,
    "I / I(Brγ)": sel["ratio_to_Brg"],
}, index=list(sel["name"])).style.format(
    {"λ_vac [µm]": "{:.5f}", "I [erg s⁻¹ cm⁻² arcsec⁻²]": "{:.3e}", "I / I(Brγ)": "{:.3f}"}
).set_table_styles(tbl_sty)
```

## Parameters
Every number entering the build, read from `build_orion_nebula_sed.PARAMS`:

```{code-cell} ipython3
---
tags: [remove-input]
mystnb:
  markdown_format: myst
---
md_table(pd.DataFrame(
    [(f"`{k}`", f"{getattr(v, 'value', v):.6g} {un.to_string() or ''}", cite(ref), note)
     for k, (v, un, ref, note) in ori.PARAMS.items()],
    columns=["parameter", "value", "source", "where / what"],
))
```

The PAH bands {cite:p}`Peeters2024`, as positions and widths (Table H.1) and H II-template intensities (Table H.2):

```{code-cell} ipython3
:tags: [remove-input]

pd.DataFrame(ori.PAH_BANDS, columns=["band", "λ₀ [µm]", "FWHM [µm]", "I [10⁻³ erg s⁻¹ cm⁻² sr⁻¹]"]
             ).set_index("band").style.set_table_styles(tbl_sty)
```

(sec:orion-limits)=
## Limitations
* **Metal lines are omitted:** [Fe II] 1.257/1.644 µm, [S II] 1.03 µm, O I 1.129/1.317 µm and [C I] are all present in Orion spectra
  {cite:p}`Peeters2024,Walmsley2000`. Their broadband share is expected to be minor but is not quantified here.
* **H₂ from the background PDR is omitted.** It matters for narrow-band H₂ filters, and reaches about half of F212N at the dissociation fronts {cite:p}`PDRs4AllXIII`.
* **The L/M dust continuum is omitted**, so L/M values are lower limits.
* **He I lines beyond 2.16 µm are missing**, because they are not in the Porter tables shipped with PyNeb.
* **He I 1.083 µm** uses the Case B emissivity, although {cite:t}`Blagrave2007` show radiative-transfer deviations for He I in Orion.
* **H I upper levels stop at n = 50**, which leaves a small dip just redward of each recombination edge.
* **The anchor inversion uses Storey & Hummer emissivities**, while {cite:t}`Peeters2024` used those of {cite:t}`ProzeskySmits2018`. The difference is expected to be at the percent level.
* **Beyond 5.27 µm the spectrum is physics only**, without a JWST data anchor.

## Rebuilding the file
This page never writes data.
To regenerate the FITS file and line list from `PARAMS`, run the build script in the backgrounds folder.
It needs PyNeb at build time; this page needs it too, for the checks above.

```bash
python build_orion_nebula_sed.py
```

## References
```{bibliography}
:filter: docname in docnames
```
