"""
Runtime loader for the ScopeSim-Targets background SEDs
(orion_nebula_background_v1.fits, ophiuchus_background_v1.fits).
Supersedes load_orion_background.py (same behaviour, any file).

Each FITS file stores a surface brightness in erg s-1 cm-2 A-1 arcsec-2
(vacuum wavelengths) in the SPECTRUM extension: FLUX plus its components.
synphot has no per-solid-angle flux unit, so the spectrum is returned as FLAM
*per arcsec2*; the solid angle is kept in the metadata. speXtra's Spextrum is
a synphot SourceSpectrum subclass and can wrap it.

Surface brightness is distance-invariant: rescale to the requested surface
brightness, never apply a distance modulus. Check meta["REDDENED"]: the Orion
file is intrinsic (apply foreground extinction), the Ophiuchus file is anchored
on observed intensities.
"""

import re

import astropy.units as u
from astropy.table import QTable
from synphot import SourceSpectrum
from synphot.models import Empirical1D


def flux_components(meta):
    """Components summed into FLUX, from COMPS or else the FLUXDEF header card."""
    if "COMPS" in meta:
        return tuple(c.strip() for c in meta["COMPS"].split(","))
    return tuple(re.findall(r"[A-Z0-9_]+", meta["FLUXDEF"].split("=", 1)[1]))


def load_background(path, components=None, extra=()):
    """Return (SourceSpectrum per arcsec2, meta).

    components : tuple of column names; default = those defining FLUX.
    extra      : additional columns to add, e.g. ("SCATTERED_UL",) for the
                 scattered-light upper limit.
    The spectrum is zero outside the tabulated range (0.97-5.40 um); callers
    must refuse bands not fully covered rather than rely on it."""
    tab = QTable.read(path, hdu="SPECTRUM")
    comps = tuple(components or flux_components(tab.meta)) + tuple(extra)
    unknown = [c for c in comps if c not in tab.colnames or c == "WAVELENGTH"]
    if unknown:
        raise ValueError(
            f"unknown components {unknown}; available: {tab.colnames[1:]}"
        )
    sb = sum(tab[c] for c in comps)
    flam_per_arcsec2 = (sb * u.arcsec**2).to(u.erg / u.s / u.cm**2 / u.AA)
    spec = SourceSpectrum(
        Empirical1D,
        points=tab["WAVELENGTH"],
        lookup_table=flam_per_arcsec2,
        keep_neg=False,
    )
    meta = dict(tab.meta)
    meta.update({
        "solid_angle": 1*u.arcsec**2,
        "components": comps,
        "wave_range": (tab["WAVELENGTH"].min(), tab["WAVELENGTH"].max()),
    })
    return spec, meta
