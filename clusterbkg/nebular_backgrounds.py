# -*- coding: utf-8 -*-
"""Orion- and Ophiuchus-based diffuse backgrounds for the embedded-cluster runs.

Both backgrounds are :class:`~scopesim_targets.extended_source.CoredPowerLaw`
targets with the *same* geometry; they differ only in spectrum, surface
brightness level and anchor frame, so a pair of simulations isolates the
effect of the background itself.

* Radii are physical lengths, meant to be matched to the cluster's King
  profile: ``r_core`` -> core radius of the nebular profile, ``r_tide`` ->
  outer truncation ``r_out``. The truncation makes the total finite for any
  index, so the background does not depend on the field of view.
* The distance is taken from ``position`` -- pass the *same* SkyCoord as for
  the cluster. Pushing the scene away means building both with a SkyCoord at
  the new distance (radii are converted to angles at construction).
* Surface brightness is distance-invariant: only the angular size changes with
  distance, never the brightness per arcsec2.

Example (reference scene: King r_core = 0.13 pc, r_tide = 2.0 pc, 414 pc)::

    pos = SkyCoord(84.7 * u.deg, -6.9 * u.deg, distance=414 * u.pc)
    king = {"n_stars": 1000, "r_core": 0.13 * u.pc, "r_tide": 2.0 * u.pc}
    cluster = ZeroAgeCluster(pos, IMFPopulation, {"n_stars": 1000},
                             KingProfileMorphology, king, rng_seed=1138)
    orion = orion_background(pos, king["r_core"], king["r_tide"])
    oph = ophiuchus_background(pos, king["r_core"], king["r_tide"])

References: see references.bib next to this file (keys in brackets).
"""
from pathlib import Path

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import QTable
from synphot.models import Empirical1D

from scopesim_targets.extended_source import CoredPowerLaw

from load_background import flux_components

SED_DIR = Path(__file__).resolve().parent
ORION_SED = "orion_nebula_background_v1.fits"
OPHIUCHUS_SED = "ophiuchus_background_v1.fits"

# ---- shared geometry ------------------------------------------------------------
# S(Hbeta) falls roughly as D^-2 near theta1 Ori C, closer to D^-3 further out
# [ODell2010, Sect. 5]. D^-2 is the inner-Huygens value.
DEFAULT_INDEX = 2.0

# Reference point of the Orion surface brightness: projected distance of the
# PDRs4All H II-region aperture from theta1 Ori C, 111.4 arcsec [Peeters2024,
# Table A.2], expressed as a length at 414 pc (the reference-scene distance;
# TODO: confirm it is also the distance adopted in [Peeters2024]). As a length
# it stays attached to the nebula when the scene is moved.
ORION_R_REF = (111.4 * u.arcsec).to_value(u.rad) * 414 * u.pc

# ---- surface-brightness anchors -----------------------------------------------
# Orion: 2MASS Ks synthetic photometry of the intrinsic SED v1 [README.md,
# "Broadband results"; 2MASS RSR via SVO/speXtra, Cohen2003 zero points].
# Intrinsic frame: a foreground screen (A_V 0.9-1.9, G23, R_V = 5.5
# [Peeters2024, Sect. 6.1]) dims it when extinction is applied.
ORION_BRIGHTNESS = ("K", "13.72 mag / arcsec2")
ORION_ANCHOR = "intrinsic"

# Ophiuchus: 2MASS Ks synthetic photometry of the H2 model [README_ophiuchus.md,
# "Anchor and validation"]; the 2MASS measurement is 18.64 [Casassus2008,
# Sect. 8.4.2]. Observed frame (foreground A_V < 0.1 [Schlegel1998, via
# Casassus2008]). rho Oph W level: an upper envelope for cloud-interior fields.
OPHIUCHUS_BRIGHTNESS = ("K", "19.06 mag / arcsec2")
OPHIUCHUS_ANCHOR = "observed"


def _require_length(name: str, value) -> None:
    if not (
        isinstance(value, u.Quantity) and value.unit.physical_type == "length"
    ):
        raise ValueError(
            f"{name} must be a physical length (e.g. pc) so it scales with the "
            f"distance in position; got {value!r}"
        )


def load_sed(filename: str, components=None, sed_dir: Path | None = None):
    """Background SED as a speXtra Spextrum, flux per arcsec2.

    ``components`` defaults to those defining ``FLUX`` in the file (e.g. add
    ``"SCATTERED_UL"`` for the scattered-light upper limit).
    """
    from spextra import (
        Spextrum,
    )  # lazy: speXtra may fetch its database on import

    tab = QTable.read(Path(sed_dir or SED_DIR) / filename, hdu="SPECTRUM")
    comps = tuple(components or flux_components(tab.meta))
    sb = sum(tab[c] for c in comps)
    # Surface brightness -> flux of one arcsec2 (synphot has no per-sr units);
    # the brightness then rescales it, so only the shape and a finite absolute
    # scale matter here.
    flam = (sb * u.arcsec**2).to(u.erg / u.s / u.cm**2 / u.AA)
    return Spextrum(
        modelclass=Empirical1D,
        points=tab["WAVELENGTH"],
        lookup_table=flam,
        keep_neg=False,
    )


def nebular_background(
    position: SkyCoord,
    spectrum,
    brightness,
    anchor: str,
    r_core: u.Quantity,
    r_tide: u.Quantity,
    *,
    index: float = DEFAULT_INDEX,
    r_ref: u.Quantity = ORION_R_REF,
    extinction=None,
) -> CoredPowerLaw:
    """Generic King-matched cored power-law background (all radii lengths)."""
    for name, value in (
        ("r_core", r_core),
        ("r_tide", r_tide),
        ("r_ref", r_ref),
    ):
        _require_length(name, value)
    return CoredPowerLaw(
        position=position,
        spectrum=spectrum,
        brightness=brightness,
        anchor=anchor,
        extinction=extinction,
        params={
            "r_core": r_core,
            "index": index,
            "r_ref": r_ref,
            "r_out": r_tide,
        },
    )


def orion_background(
    position: SkyCoord,
    r_core: u.Quantity,
    r_tide: u.Quantity,
    *,
    index: float = DEFAULT_INDEX,
    r_ref: u.Quantity = ORION_R_REF,
    extinction=None,
    components=None,
) -> CoredPowerLaw:
    """Orion-Nebula-like H II-region background (intrinsic frame).

    ``r_core``/``r_tide``: the cluster's King radii (lengths). The surface
    brightness ``ORION_BRIGHTNESS`` applies at ``r_ref``, which must lie
    inside ``r_tide``. Pass ``extinction`` for the foreground screen once
    screens are applied in the pipeline.
    """
    return nebular_background(
        position,
        load_sed(ORION_SED, components),
        ORION_BRIGHTNESS,
        ORION_ANCHOR,
        r_core,
        r_tide,
        index=index,
        r_ref=r_ref,
        extinction=extinction,
    )


def ophiuchus_background(
    position: SkyCoord,
    r_core: u.Quantity,
    r_tide: u.Quantity,
    *,
    index: float = DEFAULT_INDEX,
    r_ref: u.Quantity = ORION_R_REF,
    components=None,
) -> CoredPowerLaw:
    """rho-Oph-W-like fluorescent-H2 background (observed frame).

    Same geometry as :func:`orion_background` by default, for a controlled
    comparison: no radial law has been measured for the diffuse emission of
    L1688, so the shape is the Orion one and only the spectrum and the level
    (``OPHIUCHUS_BRIGHTNESS`` at ``r_ref``) are Ophiuchus-specific. Use
    ``index=0`` for a uniform disk out to ``r_tide`` instead.
    """
    return nebular_background(
        position,
        load_sed(OPHIUCHUS_SED, components),
        OPHIUCHUS_BRIGHTNESS,
        OPHIUCHUS_ANCHOR,
        r_core,
        r_tide,
        index=index,
        r_ref=r_ref,
    )
