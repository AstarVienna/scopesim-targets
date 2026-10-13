"""
Build a reference-traceable 0.97-5.40 um surface-brightness SED of the diffuse
near-IR emission of the rho Ophiuchi cloud (L1688), anchored on the rho Oph W
photodissociation region, for use as the "faint background" counterpart of the
Orion Nebula SED (build_orion_nebula_sed.py). Same grid, units and file layout.

Physics: rho Oph has no conspicuous H II region (Casassus2008); its diffuse
near-IR light at rho Oph W is attributed to UV-fluorescent H2 ro-vibrational
lines (Casassus2008 Sect. 8.4.2, based on the H2 1-0 S(1) map of Habart2003).
FLUX = fluorescent H2 lines. A grey-albedo scattered-starlight column is an
upper limit, as for Orion. PAH/dust emission (L/M) is NOT included (no
referenced level was available) -> L/M values are lower limits.

Every number lives in PARAMS with a reference key (see references.bib).
Unit: erg s-1 cm-2 Angstrom-1 arcsec-2, vacuum wavelengths; distance-invariant.
Unlike the Orion file, the anchor is an OBSERVED (not dereddened) intensity.

Build-time dependencies: numpy, scipy, astropy, pdrtpy (Roueff2019 H2 data).
"""

from pathlib import Path

import numpy as np
import astropy.units as u
import astropy.constants as const
from astropy.io import fits
from astropy.table import QTable
from scipy.optimize import brentq
from pdrtpy.molecule import H2

import build_orion_nebula_sed as orion  # shared helpers: grid, profiles, photometry

VERSION = "1.0.0"
OUTDIR = Path(__file__).resolve().parent

REFS = {
    "Casassus2008": "Casassus et al. 2008, MNRAS 391, 1075",
    "Habart2003": "Habart et al. 2003, A&A 397, 623",
    "Roueff2019": "Roueff et al. 2019, A&A 630, A58 (H2 transition data)",
    "Pound2023": "Pound & Wolfire 2023, AJ 165, 25 (pdrtpy, ships Roueff2019)",
    "BlackvanDishoeck1987": "Black & van Dishoeck 1987, ApJ 322, 412",
    "LePetit2006": "Le Petit et al. 2006, ApJS 164, 506 (quotes BvD87 2-1/1-0)",
    "Peeters2024": "Peeters et al. 2024, A&A 685, A74 (PDRs4All III)",
    "Cohen2003": "Cohen, Wheaton & Megeath 2003, AJ 126, 1090 (2MASS ZPs)",
    "Skrutskie2006": "Skrutskie et al. 2006, AJ 131, 1163 (2MASS)",
    "Rodrigo2012": "Rodrigo, Solano & Bayo 2012, IVOA WD (SVO FPS)",
    "Meija2016": "Meija et al. 2016, Pure Appl. Chem. 88, 265 (atomic weights)",
    "Schlegel1998": "Schlegel, Finkbeiner & Davis 1998, ApJ 500, 525 (via Casassus2008)",
    "Tiesinga2021": "Tiesinga et al. 2021, RvMP 93, 025010 (CODATA 2018)",
}

PARAMS = {
    # ---- absolute anchor: 2MASS diffuse emission on rho Oph W ----------------
    "I_Ks_obs": (
        1.0,
        u.MJy / u.sr,
        "Casassus2008",
        "Sect. 8.4.2: peak Ks intensity above background ~1 MJy/sr",
    ),
    "I_noise_2MASS": (
        0.35,
        u.MJy / u.sr,
        "Casassus2008",
        "Sect. 8.4.2: noise of the 2MASS Ks mosaic",
    ),
    "ratio_J_K": (
        1.0 / 2.7,
        u.dimensionless_unscaled,
        "Casassus2008",
        "Sect. 8.4.1: diffuse J:H:K ~ 1:2.2:2.7",
    ),
    "ratio_H_K": (
        2.2 / 2.7,
        u.dimensionless_unscaled,
        "Casassus2008",
        "Sect. 8.4.1: diffuse J:H:K ~ 1:2.2:2.7",
    ),
    "I_10S1_Ks_diluted": (
        0.2,
        u.MJy / u.sr,
        "Casassus2008",
        "Sect. 8.4.2: H2 1-0 S(1) diluted in Ks, from the "
        "Habart2003 1-0 S(1) map",
    ),
    "share_10S1_BvD87": (
        0.218,
        u.dimensionless_unscaled,
        "Casassus2008",
        "Sect. 8.4.2: 1-0 S(1) = 21.8% of H2 in Ks "
        "(BlackvanDishoeck1987 Table 3); validation only",
    ),
    # ---- H2 fluorescent level populations (two-temperature description) -----
    "R_21S1_10S1": (
        0.54,
        u.dimensionless_unscaled,
        "LePetit2006",
        "Sect. 9.2: 2-1 S(1)/1-0 S(1) = 0.54 for radiative pumping "
        "(BlackvanDishoeck1987); fixes T_vib",
    ),
    "T_rot": (
        1598.0,
        u.K,
        "Peeters2024",
        "Table 4: 1-0 S excitation temperature, atomic-PDR template "
        "(foreground formalism), 1598 +/- 109 K",
    ),
    "line_rel_min": (
        1e-3,
        u.dimensionless_unscaled,
        None,
        "numerical choice: keep lines >= 1e-3 x I(1-0 S(1))",
    ),
    # ---- line width -----------------------------------------------------------
    "T_gas": (
        330.0,
        u.K,
        "Habart2003",
        "abstract: T_gas = 330 K in the H2-emitting layer (300-345 K)",
    ),
    "m_H": (1.008, u.u, "Meija2016", "standard atomic weight of H"),
    # ---- scattered starlight (upper limit): HD 147889 binary ------------------
    "T_eff_1": (
        23000.0,
        u.K,
        "Casassus2008",
        "Sect. 8.1.1: TLUSTY fit, primary",
    ),
    "T_eff_2": (
        20000.0,
        u.K,
        "Casassus2008",
        "Sect. 8.1.1: TLUSTY fit, secondary",
    ),
    "area_ratio": (
        1.3,
        u.dimensionless_unscaled,
        "Casassus2008",
        "Sect. 8.1.1: surface-area ratio of the two components",
    ),
    # ---- context (not used numerically) --------------------------------------
    "A_V_foreground": (
        0.1,
        u.mag,
        "Schlegel1998",
        "upper limit outside the dust clouds, as quoted in "
        "Casassus2008 Sect. 8.4.1",
    ),
    "distance": (
        135.0,
        u.pc,
        "Habart2003",
        "HD 147889 parallax; context only",
    ),
    # ---- spectral range: identical to the Orion SED ---------------------------
    "lam_min": orion.PARAMS["lam_min"],
    "lam_max": orion.PARAMS["lam_max"],
    "F0_J": orion.PARAMS["F0_J"],
    "F0_H": orion.PARAMS["F0_H"],
    "F0_Ks": orion.PARAMS["F0_Ks"],
}

SR_TO_ARCSEC2 = orion.SR_TO_ARCSEC2


def p(name, unit=None):
    val, un, _, _ = PARAMS[name]
    q = val * un
    return q.to_value(unit) if unit is not None else q


# --------------------------------------------------------------------------
# Fluorescent H2
# --------------------------------------------------------------------------
def h2_table():
    """Roueff2019 transitions (via pdrtpy) with upper-level vibrational
    ground energies E(v, J=0) attached."""
    t = H2().transition_data
    e_lev = {(int(r["vu"]), int(r["Ju"])): float(r["Tu"]) for r in t}
    e_lev[(0, 0)] = 0.0
    e0 = np.array([e_lev[(int(v), 0)] for v in t["vu"]])
    return t, e0


def h2_relative_intensities(t, e0, t_vib, t_rot):
    """Line intensities (arbitrary scale) for populations
    N(v,J) = g(v,J) exp(-E(v,0)/T_vib) exp(-[E(v,J)-E(v,0)]/T_rot);
    I_ul proportional to N_u A_ul sigma_ul. g includes the nuclear-spin
    weight in the Roueff2019 table (i.e. ortho/para = 3)."""
    return (
        np.asarray(t["gu"], float)
        * np.exp(-e0 / t_vib)
        * np.exp(-(np.asarray(t["Tu"], float) - e0) / t_rot)
        * np.asarray(t["A"], float)
        * np.asarray(t["sigma"], float)
    )


def solve_t_vib(t, e0, t_rot):
    names = np.array(t["Line"])
    i10, i21 = [np.where(names == n)[0][0] for n in ("H210S1", "H221S1")]
    target = p("R_21S1_10S1").value
    f = lambda tv: (
        h2_relative_intensities(t, e0, tv, t_rot)[i21]
        / h2_relative_intensities(t, e0, tv, t_rot)[i10]
        - target
    )
    return brentq(f, 1e3, 5e4), i10


def band_mean_inu_of_lines(lam_aa, inten_sr, lam_t, trans):
    """Photon-weighted band-mean I_nu [MJy/sr] of delta lines (erg s-1 cm-2
    sr-1), same convention as the synthetic photometry of the grid."""
    t_at = np.interp(lam_aa, lam_t, trans, left=0, right=0)
    flam = np.sum(inten_sr * lam_aa * t_at) / np.trapezoid(
        lam_t * trans, lam_t
    )
    piv = np.sqrt(
        np.trapezoid(lam_t * trans, lam_t) / np.trapezoid(trans / lam_t, lam_t)
    )
    return flam * piv**2 / const.c.to_value(u.AA / u.s) / 1e-23 / 1e6


def main():
    lam_lo, lam_hi = p("lam_min", u.AA), p("lam_max", u.AA)
    svo = (
        Path(
            __import__("importlib.util").util.find_spec("spextra").origin
        ).parent
        / "data"
        / "svo"
    )
    curves = {
        k: orion.svo_curve(svo / f"2MASS/2MASS.{k}") for k in ("J", "H", "Ks")
    }

    # ---- H2 line list, normalised to the Casassus2008 1-0 S(1) level --------
    t, e0 = h2_table()
    t_rot = p("T_rot", u.K)
    t_vib, i10 = solve_t_vib(t, e0, t_rot)
    rel = h2_relative_intensities(t, e0, t_vib, t_rot)
    rel /= rel[i10]
    lam_all = 1e8 / np.asarray(t["sigma"], float)  # vacuum Angstrom
    one = band_mean_inu_of_lines(
        np.array([lam_all[i10]]), np.array([1.0]), *curves["Ks"]
    )
    i_10s1 = p("I_10S1_Ks_diluted", u.MJy / u.sr) / one  # erg s-1 cm-2 sr-1
    keep = (
        (lam_all >= lam_lo)
        & (lam_all <= lam_hi)
        & (rel >= p("line_rel_min").value)
    )
    share_ks = band_mean_inu_of_lines(
        lam_all[[i10]], np.array([i_10s1]), *curves["Ks"]
    ) / band_mean_inu_of_lines(lam_all, rel * i_10s1, *curves["Ks"])
    rows = [
        (
            "H2",
            str(t["Line"][i]),
            str(t["Transition"][i]),
            lam_all[i],
            rel[i] * i_10s1,
            "Roueff2019",
        )
        for i in np.where(keep)[0]
    ]

    # pure thermal Doppler FWHM of H2 at the H2-layer gas temperature
    fwhm = {
        "H2": np.sqrt(
            8 * np.log(2) * const.k_B * p("T_gas") / (2 * p("m_H"))
        ).to_value(u.km / u.s)
    }
    grid = orion.build_grid([(r[3], fwhm["H2"]) for r in rows], lam_lo, lam_hi)
    comp = {"H2_LINES": orion.gaussian_lines(grid, rows, fwhm) / SR_TO_ARCSEC2}

    # ---- scattered starlight upper limit: grey albedo, HD 147889 binary ------
    shape = orion.planck_lambda(grid, p("T_eff_1", u.K)) * p(
        "area_ratio"
    ).value + orion.planck_lambda(grid, p("T_eff_2", u.K))
    # normalised so that scattered light ALONE reproduces the observed J
    # intensity (Casassus2008): it cannot exceed that, hence an upper limit.
    lt, tr = curves["J"]
    unit_inu, _ = orion.band_mean_inu_mjysr(
        grid, shape, lt, tr
    )  # shape read as per arcsec2
    i_j_obs = p("I_Ks_obs", u.MJy / u.sr) * p("ratio_J_K").value
    comp["SCATTERED_UL"] = (
        shape * i_j_obs / unit_inu
    )  # erg s-1 cm-2 A-1 arcsec-2
    total = comp["H2_LINES"].copy()

    # ---- diagnostics / synthetic photometry ----------------------------------
    diag = {
        "T_vib_K": t_vib,
        "T_rot_K": t_rot,
        "I_10S1_erg_s_cm2_sr": i_10s1,
        "share_10S1_in_Ks_model": share_ks,
        "share_10S1_in_Ks_BvD87": p("share_10S1_BvD87").value,
        "FWHM_H2_kms": fwhm["H2"],
        "n_lines_H2": len(rows),
    }
    obs = {
        "J": p("I_Ks_obs", u.MJy / u.sr) * p("ratio_J_K").value,
        "H": p("I_Ks_obs", u.MJy / u.sr) * p("ratio_H_K").value,
        "Ks": p("I_Ks_obs", u.MJy / u.sr),
    }
    phot = {}
    for band in ("J", "H", "Ks"):
        lt, tr = curves[band]
        f0 = p(f"F0_{band}", u.erg / u.s / u.cm**2 / u.AA)
        e = {
            "mag_arcsec2_Vega": -2.5
            * np.log10(orion.band_mean_flam(grid, total, lt, tr) / f0)
        }
        e["MJy_sr"], e["lam_pivot_A"] = orion.band_mean_inu_mjysr(
            grid, total, lt, tr
        )
        e["obs_MJy_sr"] = obs[band]
        e["obs_mag_arcsec2_Vega"] = e["mag_arcsec2_Vega"] - 2.5 * np.log10(
            obs[band] / e["MJy_sr"]
        )
        e["dev_sigma"] = (e["MJy_sr"] - obs[band]) / p(
            "I_noise_2MASS", u.MJy / u.sr
        )
        e["scattered_UL_MJy_sr"], _ = orion.band_mean_inu_mjysr(
            grid, comp["SCATTERED_UL"], lt, tr
        )
        phot[band] = e
    lt, tr = orion.svo_curve(svo / "Gemini/NIRI.Lprime-G0207w")
    phot["Lp_NIRI"] = {
        "MJy_sr": orion.band_mean_inu_mjysr(grid, total, lt, tr)[0]
    }

    # ---- write products (same layout as the Orion file) ----------------------
    sb_unit = u.erg / u.s / u.cm**2 / u.AA / u.arcsec**2
    cols = [
        fits.Column("WAVELENGTH", "D", "Angstrom", array=grid),
        fits.Column("FLUX", "E", sb_unit.to_string(), array=total),
    ]
    cols += [
        fits.Column(k, "E", sb_unit.to_string(), array=v)
        for k, v in comp.items()
    ]
    hdu = fits.BinTableHDU.from_columns(cols, name="SPECTRUM")
    h = hdu.header
    h["OBJECT"] = "rho Ophiuchi (L1688) diffuse emission, rho Oph W"
    h["VERSION"] = VERSION
    h["WAVEAIR"] = (False, "wavelengths are in vacuum")
    h["SBUNIT"] = ("arcsec2", "FLUX is a surface brightness per arcsec2")
    h["DISTINV"] = (True, "surface brightness: no distance modulus")
    h["REDDENED"] = (True, "anchored on OBSERVED intensities (not dereddened)")
    h["FLUXDEF"] = "FLUX = H2_LINES"
    h["COMPS"] = "H2_LINES"
    h["SCATUL"] = "SCATTERED_UL is an upper limit, not in FLUX"
    h["LMINCOMP"] = (True, "L/M lower limit: PAH and dust continuum omitted")
    for i, (k, (val, un, ref, note)) in enumerate(PARAMS.items()):
        h[f"HIERARCH PAR{i:02d} NAME"] = k
        h[f"HIERARCH PAR{i:02d} VALUE"] = float(getattr(val, "value", val))
        h[f"HIERARCH PAR{i:02d} UNIT"] = un.to_string() or "1"
        h[f"HIERARCH PAR{i:02d} REF"] = ref or "none"
    for i, (k, v) in enumerate(REFS.items()):
        h[f"HIERARCH REF{i:02d}"] = f"{k}: {v}"
    for k, v in diag.items():
        h[f"HIERARCH DIAG {k}"] = float(v)
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(
        OUTDIR / "ophiuchus_background_v1.fits", overwrite=True
    )

    lines = QTable(
        rows=[
            (
                s,
                n,
                tr_,
                lam * u.AA,
                i / SR_TO_ARCSEC2 * u.erg / u.s / u.cm**2 / u.arcsec**2,
                ref,
            )
            for s, n, tr_, lam, i, ref in rows
        ],
        names=(
            "species",
            "name",
            "transition",
            "lambda_vac",
            "intensity",
            "data_ref",
        ),
    )
    lines["ratio_to_10S1"] = lines["intensity"].value / (
        i_10s1 / SR_TO_ARCSEC2
    )
    lines.sort("lambda_vac")
    lines.meta.update(
        {
            "T_vib_K": float(t_vib),
            "T_rot_K": float(t_rot),
            "anchor": "Casassus2008 1-0 S(1) diluted in 2MASS Ks = 0.2 MJy/sr",
        }
    )
    lines.write(OUTDIR / "ophiuchus_h2_linelist_v1.ecsv", overwrite=True)
    return grid, total, comp, diag, phot, lines


if __name__ == "__main__":
    import json

    grid, total, comp, diag, phot, lines = main()
    print(json.dumps({"diag": diag, "phot": phot}, indent=1, default=float))
    top = np.argsort(lines["ratio_to_10S1"])[::-1][:12]
    print(lines[top]["transition", "lambda_vac", "ratio_to_10S1"])
