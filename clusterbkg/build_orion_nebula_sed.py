"""
Build a reference-traceable 0.97-5.40 um surface-brightness SED of the Orion
Nebula (M42) H II region, for use as the diffuse background of a mock embedded
cluster in ScopeSim-Targets.

Every number used below lives in PARAMS with a key into REFS (and into
references.bib). Nothing numeric is hard-coded inside the physics functions.

Output intensity unit: erg s-1 cm-2 Angstrom-1 arcsec-2 (vacuum wavelengths).
This is a surface brightness, i.e. distance-INVARIANT: do not apply a distance
modulus to it. The file carries an absolute scale so the loader can rescale it
to the requested brightness; the reference level is the Orion Bar H II-region
aperture of Peeters et al. (2024).

Build-time dependencies: numpy, astropy, pyneb (Luridiana et al. 2015).
Runtime (loading the product): astropy (+ synphot / speXtra if desired).
"""

from pathlib import Path

import numpy as np
import astropy.units as u
import astropy.constants as const
from astropy.io import fits
from astropy.io.votable import parse_single_table
from astropy.table import QTable
import pyneb as pn

VERSION = "1.0.0"
OUTDIR = Path(__file__).resolve().parent

# --------------------------------------------------------------------------
# References (keys match references.bib)
# --------------------------------------------------------------------------
REFS = {
    "Esteban2004": "Esteban et al. 2004, MNRAS 355, 229",
    "Peeters2024": "Peeters et al. 2024, A&A 685, A74 (PDRs4All III)",
    "ODell2010": "O'Dell & Harris 2010, AJ 140, 985",
    "SimonDiaz2006": "Simon-Diaz et al. 2006, A&A 448, 351",
    "Storey1995": "Storey & Hummer 1995, MNRAS 272, 41 (H I Case B)",
    "Porter2012": "Porter et al. 2012, MNRAS 425, L28 (He I)",
    "Porter2013": "Porter et al. 2013, MNRAS 433, L89 (He I corrigendum)",
    "Ercolano2006": "Ercolano & Storey 2006, MNRAS 372, 1875 (free-bound)",
    "Storey1991": "Storey & Hummer 1991, CoPhC 66, 129 (free-free Gaunt)",
    "Luridiana2015": "Luridiana, Morisset & Shaw 2015, A&A 573, A42 (PyNeb)",
    "Smith2007": "Smith et al. 2007, ApJ 656, 770 (PAHFIT Drude profile)",
    "Cohen2003": "Cohen, Wheaton & Megeath 2003, AJ 126, 1090 (2MASS ZPs)",
    "Skrutskie2006": "Skrutskie et al. 2006, AJ 131, 1163 (2MASS)",
    "Rodrigo2012": "Rodrigo, Solano & Bayo 2012, IVOA WD (SVO FPS)",
    "Morton2000": "Morton 2000, ApJS 130, 403 (air-to-vacuum)",
    "Tiesinga2021": "Tiesinga et al. 2021, RvMP 93, 025010 (CODATA 2018)",
    "Tsivilev": "Tsivilev et al., Astronomy Reports (Orion A RRLs, 13 mm)",
    "Gordon1986": "Gordon 1986, ApJ 308, 288 (via speXtra orion template)",
    "Marconi1998": "Marconi et al. 1998, A&A 330, 696",
    "Meija2016": "Meija et al. 2016, Pure Appl. Chem. 88, 265 (atomic weights)",
}

# --------------------------------------------------------------------------
# Parameters: (value, unit, reference key, where/what)
# --------------------------------------------------------------------------
PARAMS = {
    # physical conditions of the ionised gas (single self-consistent source)
    "T_e": (
        8730.0,
        u.K,
        "Esteban2004",
        "T(O II+III) = 8730 +/- 320 K, as adopted in Sect. 5",
    ),
    "n_e": (
        8900.0,
        u.cm**-3,
        "Esteban2004",
        "n_e = 8900 +/- 200 cm-3, as adopted in Sect. 5",
    ),
    "y_He+": (
        0.0874,
        u.dimensionless_unscaled,
        "Esteban2004",
        "Sect. 5/Table 4: He+/H+ = 0.0874 +/- 0.0006",
    ),
    "y_He++": (
        0.0,
        u.dimensionless_unscaled,
        "Tsivilev",
        "RRLs give He++/H+ < 7e-4 in Orion A; set to 0",
    ),
    # empirical correction of the radiative-transfer-sensitive singlet line
    "R_HeI_2058_1700": (
        8.5,
        u.dimensionless_unscaled,
        "Peeters2024",
        "Sect. 6.2: observed He I 2.058/1.70 um = 7-10 "
        "(confirming Marconi1998); mid-range adopted",
    ),
    # absolute anchor: Brgamma-derived emission measure, H II template
    "EM_ref": (
        3.3e6,
        u.pc * u.cm**-6,
        "Peeters2024",
        "Table 1: EM(Brgamma), H II-region template",
    ),
    "T_e_P24": (
        1.0e4,
        u.K,
        "Peeters2024",
        "Sect. 6.1: Case B T_e adopted for H I analysis",
    ),
    "n_e_P24": (
        1.0e3,
        u.cm**-3,
        "Peeters2024",
        "Sect. 6.1: Case B n_e adopted for H I analysis",
    ),
    # scattered starlight (upper limit, see README)
    "EW_obs_Hb": (
        423.0,
        u.AA,
        "ODell2010",
        "Sect. 5.1.4: M4 group EW(Hb) = 423 +/- 6 A, similar "
        "to samples closest to theta1 Ori C",
    ),
    "T_eff_th1C": (
        39000.0,
        u.K,
        "SimonDiaz2006",
        "theta1 Ori C: T_eff = 39000 +/- 1000 K",
    ),
    # spectral range: NIRSpec template range 0.97-5.27 um, extended to M band
    "lam_min": (0.97, u.um, "Peeters2024", "NIRSpec template blue limit"),
    "lam_max": (
        5.40,
        u.um,
        None,
        "design choice: covers M band; no data anchor > 5.27 um",
    ),
    "m_He": (4.002602, u.u, "Meija2016", "standard atomic weight of He"),
    "n_upper_max": (
        50,
        u.dimensionless_unscaled,
        "Storey1995",
        "highest upper level used from the SH95 tables",
    ),
    # 2MASS zero points (in-band, photon-weighted), for the sanity check
    "F0_J": (
        3.129e-13,
        u.W / u.cm**2 / u.um,
        "Cohen2003",
        "2MASS J zero point",
    ),
    "F0_H": (
        1.133e-13,
        u.W / u.cm**2 / u.um,
        "Cohen2003",
        "2MASS H zero point",
    ),
    "F0_Ks": (
        4.283e-14,
        u.W / u.cm**2 / u.um,
        "Cohen2003",
        "2MASS Ks zero point",
    ),
    # cross-check of the absolute level via the speXtra orion template
    "F_Hb_total": (
        1.40e-7,
        u.erg / u.s / u.cm**2,
        "Gordon1986",
        "dereddened total Hb of the speXtra 'orion' template "
        "(84.2 GHz free-free scaled to Hb)",
    ),
    "D_Huygens": (
        5.0,
        u.arcmin,
        "ODell2010",
        "Sect. 1: Huygens Region ~0.7 pc = 5 arcmin",
    ),
}

# AIB (PAH) 3.3-3.6 um complex, PAHFIT decomposition of the H II-region
# template: position and FWHM from Peeters2024 Table H.1, integrated
# intensity [1e-3 erg s-1 cm-2 sr-1] from Peeters2024 Table H.2.
PAH_BANDS = [  # name, lambda0 [um], FWHM [um], P [1e-3 erg s-1 cm-2 sr-1]
    ("3.29", 3.291, 0.03762, 10.58),
    ("3.39", 3.395, 0.00995, 0.30),
    ("3.40", 3.405, 0.02691, 1.01),
    ("3.42", 3.4253, 0.015, 0.15),
    ("3.46", 3.464, 0.07012, 1.61),
    ("3.51", 3.516, 0.0271, 0.39),
    ("3.56", 3.561, 0.02, 0.14),
]
PAH_REF = "Peeters2024 Tables H.1 (profiles) and H.2 (H II template)"

SERIES = {3: "Pa", 4: "Br", 5: "Pf", 6: "Hu", 7: "n7", 8: "n8"}
GREEK = "abgd"  # alpha..delta; higher members named by upper level (e.g. Br10)

SR_TO_ARCSEC2 = (1 * u.sr).to_value(u.arcsec**2)


def p(name, unit=None):
    """Return a parameter value (optionally converted)."""
    val, un, _, _ = PARAMS[name]
    q = val * un
    return q.to_value(unit) if unit is not None else q


# --------------------------------------------------------------------------
# Atomic physics helpers
# --------------------------------------------------------------------------
def hydrogen_vacuum_wavelength(n_up, n_lo):
    """Vacuum wavelength [Angstrom] from the Rydberg formula with the
    hydrogen reduced-mass correction; constants are CODATA 2018
    (Tiesinga2021) as shipped by astropy."""
    r_h = const.Ryd / (1 + const.m_e / const.m_p)
    wavenumber = r_h * (1 / n_lo**2 - 1 / n_up**2)
    return (1 / wavenumber).to_value(u.AA)


def air_to_vacuum(lam_air):
    """Air -> vacuum [Angstrom], refractive index of Morton (2000),
    solved iteratively."""
    lam_vac = np.array(lam_air, dtype=float)
    for _ in range(4):
        s2 = (1e4 / lam_vac) ** 2
        n = 1 + 8.34254e-5 + 2.406147e-2 / (130 - s2) + 1.5998e-4 / (38.9 - s2)
        lam_vac = lam_air * n
    return lam_vac


def thermal_fwhm_kms(mass):
    """Pure thermal Doppler FWHM at T_e [km/s]: sqrt(8 ln2 k T / m)."""
    t_e = p("T_e")
    return np.sqrt(8 * np.log(2) * const.k_B * t_e / mass).to_value(u.km / u.s)


# --------------------------------------------------------------------------
# Absolute anchor
# --------------------------------------------------------------------------
def emission_measure_at_adopted_conditions(h1):
    """Invert Peeters2024's EM(Brgamma) with their stated Case B conditions
    to recover the dereddened Brgamma intensity, then re-express it as the
    EM that gives the same Brgamma at the adopted (Esteban2004) conditions.
    Convention: I = eps * EM / (4 pi), eps = 4 pi j / (n_e n_p)."""
    em_ref = p("EM_ref", u.cm**-5)
    eps_p24 = h1.getEmissivity(
        p("T_e_P24", u.K), p("n_e_P24", u.cm**-3), lev_i=7, lev_j=4
    )
    i_brg = eps_p24 * em_ref / (4 * np.pi)  # erg s-1 cm-2 sr-1
    eps_adopt = h1.getEmissivity(
        p("T_e", u.K), p("n_e", u.cm**-3), lev_i=7, lev_j=4
    )
    em_eff = 4 * np.pi * i_brg / eps_adopt  # cm-5
    return em_eff, i_brg


# --------------------------------------------------------------------------
# Line lists
# --------------------------------------------------------------------------
def hydrogen_lines(h1, em_eff, lam_lo, lam_hi):
    t_e, n_e = p("T_e", u.K), p("n_e", u.cm**-3)
    rows = []
    for n_lo in SERIES:
        for n_up in range(n_lo + 1, int(p("n_upper_max").value) + 1):
            lam = hydrogen_vacuum_wavelength(n_up, n_lo)
            if not lam_lo <= lam <= lam_hi:
                continue
            eps = h1.getEmissivity(t_e, n_e, lev_i=n_up, lev_j=n_lo)
            dn = n_up - n_lo
            name = (
                f"{SERIES[n_lo]}{GREEK[dn - 1]}"
                if dn <= len(GREEK)
                else f"{SERIES[n_lo]}{n_up}"
            )  # e.g. Brg, Br10
            rows.append(
                (
                    "H I",
                    name,
                    f"{n_up}-{n_lo}",
                    lam,
                    eps * em_eff / (4 * np.pi),
                    "Storey1995",
                )
            )
    return rows


def helium_lines(he1, em_eff, lam_lo, lam_hi):
    t_e, n_e = p("T_e", u.K), p("n_e", u.cm**-3)
    y = p("y_He+").value
    rows = []
    for label in he1.labels:
        lam_air = float(label)
        if lam_air < 2000:  # PyNeb labels: air above 2000 A
            continue
        lam = float(air_to_vacuum(lam_air))
        if not lam_lo <= lam <= lam_hi:
            continue
        eps = he1.getEmissivity(t_e, n_e, label=label)
        rows.append(
            [
                "He I",
                f"HeI{label}",
                "",
                lam,
                y * eps * em_eff / (4 * np.pi),
                "Porter2012,Porter2013",
            ]
        )
    # He I 2.058 um (2p1P-2s1S) is set by 584 A resonance-line transfer and is
    # not well described by Case B; tie it to the Case-B-safe 1.70 um triplet
    # line (same upper level as 4471 A) with the observed ratio.
    i_1700 = [r[4] for r in rows if r[1] == "HeI17003.0"][0]
    for r in rows:
        if r[1] == "HeI20581.0":
            r[4] = p("R_HeI_2058_1700").value * i_1700
            r[5] = "Porter2012 x Peeters2024 (obs. 2.058/1.70)"
    return [tuple(r) for r in rows]


# --------------------------------------------------------------------------
# Spectral components
# --------------------------------------------------------------------------
def build_grid(lines, lam_lo, lam_hi, r_base=2.0e4, n_sig=6, per_sig=10):
    """Log grid at R_base plus dense sampling (sigma/per_sig) around lines.
    Grid spacing is a numerical choice, not a physical parameter."""
    base = np.exp(np.arange(np.log(lam_lo), np.log(lam_hi), 1 / r_base))
    extra = []
    for lam0, fwhm_kms in lines:
        sig = lam0 * fwhm_kms / const.c.to_value(u.km / u.s) / 2.3548
        extra.append(
            lam0 + np.linspace(-n_sig, n_sig, 2 * n_sig * per_sig + 1) * sig
        )
    grid = np.unique(np.concatenate([base, *extra, [lam_lo, lam_hi]]))
    return grid[(grid >= lam_lo) & (grid <= lam_hi)]


def gaussian_lines(grid, rows, fwhm_by_species):
    """Sum of Gaussians; each renormalised on the grid so the integral equals
    the tabulated line intensity exactly."""
    spec = np.zeros_like(grid)
    c_kms = const.c.to_value(u.km / u.s)
    for species, _, _, lam0, inten, _ in rows:
        sig = lam0 * fwhm_by_species[species] / c_kms / 2.3548
        m = np.abs(grid - lam0) < 8 * sig
        prof = np.exp(-0.5 * ((grid[m] - lam0) / sig) ** 2)
        spec[m] += inten * prof / np.trapezoid(prof, grid[m])
    return spec  # erg s-1 cm-2 sr-1 A-1


def nebular_continuum(cont, grid, i_hb):
    """H I + He I free-bound (Ercolano2006), free-free (Storey1991) and
    two-photon continuum from PyNeb, normalised to I(Hbeta)."""
    per_hb = cont.get_continuum(
        tem=p("T_e", u.K),
        den=p("n_e", u.cm**-3),
        He1_H=p("y_He+").value,
        He2_H=p("y_He++").value,
        wl=grid,
        HI_label="4_2",
    )
    return per_hb * i_hb  # erg s-1 cm-2 sr-1 A-1


def pah_complex(grid):
    """Drude profiles (PAHFIT definition, Smith2007) in I_nu, each normalised to its
    integrated intensity P (Peeters2024 Table H.2), returned as I_lambda."""
    lam_um = grid * 1e-4
    nu = const.c.to_value(u.um / u.s) / lam_um
    order = np.argsort(nu)
    i_lam = np.zeros_like(grid)
    for _, lam0, fwhm, power in PAH_BANDS:
        gam = fwhm / lam0
        prof_nu = gam**2 / ((lam_um / lam0 - lam0 / lam_um) ** 2 + gam**2)
        norm = np.trapezoid(prof_nu[order], nu[order])
        i_nu = power * 1e-3 * prof_nu / norm  # erg s-1 cm-2 sr-1 Hz-1
        i_lam += i_nu * nu / grid  # -> per Angstrom
    return i_lam


def planck_lambda(grid, t_eff):
    x = (
        (const.h * const.c / (grid * u.AA * const.k_B * t_eff * u.K))
        .decompose()
        .value
    )
    return 1.0 / (grid**5 * np.expm1(x))  # arbitrary scale


def scattered_starlight(grid, i_hb, i_cont_hb):
    """Upper limit: grey-albedo backscatter of a theta1 Ori C blackbody,
    normalised so that total continuum at Hbeta reproduces EW_obs(Hb).
    Grey albedo over-predicts the red/NIR part because ODell2010 find the
    scattered light is bluer than the illuminating star."""
    i_sca_hb = i_hb / p("EW_obs_Hb", u.AA) - i_cont_hb  # at 4861 A
    shape = planck_lambda(grid, p("T_eff_th1C", u.K))
    shape_hb = planck_lambda(np.array([4861.33]), p("T_eff_th1C", u.K))[0]
    return i_sca_hb * shape / shape_hb, i_sca_hb


# --------------------------------------------------------------------------
# Synthetic photometry (sanity check only)
# --------------------------------------------------------------------------
def svo_curve(path):
    tab = parse_single_table(str(path)).to_table()
    return np.asarray(tab["Wavelength"], float), np.asarray(
        tab["Transmission"], float
    )


def band_mean_flam(grid, flam, lam_t, trans):
    """Photon-weighted in-band mean F_lambda (the Cohen2003 convention)."""
    t = np.interp(grid, lam_t, trans, left=0, right=0)
    return np.trapezoid(flam * grid * t, grid) / np.trapezoid(grid * t, grid)


def band_mean_inu_mjysr(grid, flam_arcsec2, lam_t, trans):
    fl = band_mean_flam(grid, flam_arcsec2, lam_t, trans)
    t = np.interp(grid, lam_t, trans, left=0, right=0)
    lam_piv = np.sqrt(
        np.trapezoid(grid * t, grid) / np.trapezoid(t / grid, grid)
    )
    fnu = (
        fl * lam_piv**2 / const.c.to_value(u.AA / u.s)
    )  # erg s-1 cm-2 Hz-1 arcsec-2
    return fnu * SR_TO_ARCSEC2 / 1e-23 / 1e6, lam_piv  # MJy sr-1


# --------------------------------------------------------------------------
def main():
    lam_lo, lam_hi = p("lam_min", u.AA), p("lam_max", u.AA)
    h1, he1, cont = pn.RecAtom("H", 1), pn.RecAtom("He", 1), pn.Continuum()

    em_eff, i_brg_ref = emission_measure_at_adopted_conditions(h1)
    t_e, n_e = p("T_e", u.K), p("n_e", u.cm**-3)
    i_hb = h1.getEmissivity(t_e, n_e, lev_i=4, lev_j=2) * em_eff / (4 * np.pi)

    rows = hydrogen_lines(h1, em_eff, lam_lo, lam_hi) + helium_lines(
        he1, em_eff, lam_lo, lam_hi
    )
    fwhm = {
        "H I": thermal_fwhm_kms(const.m_p + const.m_e),
        "He I": thermal_fwhm_kms(p("m_He")),
    }
    grid = build_grid([(r[3], fwhm[r[0]]) for r in rows], lam_lo, lam_hi)

    comp = {  # erg s-1 cm-2 sr-1 A-1
        "H_LINES": gaussian_lines(
            grid, [r for r in rows if r[0] == "H I"], fwhm
        ),
        "HE_LINES": gaussian_lines(
            grid, [r for r in rows if r[0] == "He I"], fwhm
        ),
        "NEB_CONT": nebular_continuum(cont, grid, i_hb),
        "PAH": pah_complex(grid),
    }
    i_cont_hb = nebular_continuum(cont, np.array([4861.33]), i_hb)[0]
    comp["SCATTERED_UL"], i_sca_hb = scattered_starlight(grid, i_hb, i_cont_hb)
    comp = {k: v / SR_TO_ARCSEC2 for k, v in comp.items()}  # -> per arcsec2
    total = comp["H_LINES"] + comp["HE_LINES"] + comp["NEB_CONT"] + comp["PAH"]

    # ---------------- diagnostics -------------------------------------------
    diag = {
        "EM_eff_pc_cm6": em_eff / (1 * u.pc * u.cm**-6).to_value(u.cm**-5),
        "I_Brg_ref_erg_s_cm2_sr": i_brg_ref,
        "I_Hb_erg_s_cm2_arcsec2": i_hb / SR_TO_ARCSEC2,
        "EW_Hb_dustfree_A": i_hb / i_cont_hb,
        "scatter_to_nebular_cont_at_Hb": i_sca_hb / i_cont_hb,
        "FWHM_HI_kms": fwhm["H I"],
        "FWHM_HeI_kms": fwhm["He I"],
        "n_lines_HI": sum(r[0] == "H I" for r in rows),
        "n_lines_HeI": sum(r[0] == "He I" for r in rows),
    }
    # cross-check vs the speXtra 'orion' template (Walsh 1996 / Gordon 1986)
    area = np.pi * (p("D_Huygens", u.arcsec) / 2) ** 2
    diag["I_Hb_spextra_huygens_mean"] = (
        p("F_Hb_total", u.erg / u.s / u.cm**2) / area
    )

    svo = (
        Path(
            __import__("importlib.util").util.find_spec("spextra").origin
        ).parent
        / "data"
        / "svo"
    )
    phot = {}
    for band, fname, zp in (
        ("J", "2MASS/2MASS.J", "F0_J"),
        ("H", "2MASS/2MASS.H", "F0_H"),
        ("Ks", "2MASS/2MASS.Ks", "F0_Ks"),
    ):
        lt, tr = svo_curve(svo / fname)
        f0 = p(zp, u.erg / u.s / u.cm**2 / u.AA)
        entry = {
            "mag_arcsec2_Vega": -2.5
            * np.log10(band_mean_flam(grid, total, lt, tr) / f0)
        }
        entry["MJy_sr"], entry["lam_pivot_A"] = band_mean_inu_mjysr(
            grid, total, lt, tr
        )
        for k, v in comp.items():
            entry[f"frac_{k}"] = band_mean_flam(
                grid, v, lt, tr
            ) / band_mean_flam(grid, total, lt, tr)
        phot[band] = entry
    lt, tr = svo_curve(svo / "Gemini/NIRI.Lprime-G0207w")
    entry = {}
    entry["MJy_sr"], entry["lam_pivot_A"] = band_mean_inu_mjysr(
        grid, total, lt, tr
    )
    for k, v in comp.items():
        entry[f"frac_{k}"] = band_mean_flam(grid, v, lt, tr) / band_mean_flam(
            grid, total, lt, tr
        )
    phot["Lp_NIRI"] = entry

    # ---------------- write products -----------------------------------------
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
    h["OBJECT"] = "Orion Nebula (M42) H II region background"
    h["VERSION"] = VERSION
    h["WAVEAIR"] = (False, "wavelengths are in vacuum")
    h["SBUNIT"] = ("arcsec2", "FLUX is a surface brightness per arcsec2")
    h["DISTINV"] = (True, "surface brightness: no distance modulus")
    h["REDDENED"] = (False, "intrinsic; apply foreground extinction")
    h["FLUXDEF"] = "FLUX = H_LINES + HE_LINES + NEB_CONT + PAH"
    h["SCATUL"] = "SCATTERED_UL is an upper limit, not in FLUX"
    for i, (k, (val, un, ref, note)) in enumerate(PARAMS.items()):
        h[f"HIERARCH PAR{i:02d} NAME"] = k
        h[f"HIERARCH PAR{i:02d} VALUE"] = float(getattr(val, "value", val))
        h[f"HIERARCH PAR{i:02d} UNIT"] = un.to_string() or "1"
        h[f"HIERARCH PAR{i:02d} REF"] = ref or "none"
    for i, (k, v) in enumerate(REFS.items()):
        h[f"HIERARCH REF{i:02d}"] = f"{k}: {v}"
    h["HIERARCH PAHREF"] = PAH_REF
    for k, v in diag.items():
        h[f"HIERARCH DIAG {k}"] = float(v)
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(
        OUTDIR / "orion_nebula_background_v1.fits", overwrite=True
    )

    lines = QTable(
        rows=[
            (
                s,
                n,
                t,
                lam * u.AA,
                i / SR_TO_ARCSEC2 * u.erg / u.s / u.cm**2 / u.arcsec**2,
                ref,
            )
            for s, n, t, lam, i, ref in rows
        ],
        names=(
            "species",
            "name",
            "transition",
            "lambda_vac",
            "intensity",
            "emissivity_ref",
        ),
    )
    brg = [r[4] for r in rows if r[1] == "Brg"][0] / SR_TO_ARCSEC2
    lines["ratio_to_Brg"] = lines["intensity"].value / brg
    lines.sort("lambda_vac")
    lines.meta.update(
        {
            "T_e": p("T_e").value,
            "n_e": p("n_e").value,
            "y_He+": p("y_He+").value,
            "conditions_ref": "Esteban2004",
            "anchor": "Peeters2024 EM(Brg) of H II template",
            "EM_eff_pc_cm6": float(diag["EM_eff_pc_cm6"]),
        }
    )
    lines.write(OUTDIR / "orion_nebula_linelist_v1.ecsv", overwrite=True)
    return grid, total, comp, diag, phot, lines


if __name__ == "__main__":
    import json

    grid, total, comp, diag, phot, lines = main()
    print(json.dumps({"diag": diag, "phot": phot}, indent=1, default=float))
    print(
        lines[
            [
                i
                for i, n in enumerate(lines["name"])
                if n
                in (
                    "Pab",
                    "Pag",
                    "Brg",
                    "Bra",
                    "Brb",
                    "Pfb",
                    "Pfg",
                    "HeI10830.0",
                    "HeI20581.0",
                    "HeI17003.0",
                )
            ]
        ]
    )
