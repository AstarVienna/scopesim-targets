# -*- coding: utf-8 -*-
"""Referenced random textures for the nebular backgrounds.

The texture statistics come from the surface-brightness analysis of the inner
Orion Nebula by Arthur, Medina & Henney (2016) [Arthur2016]: a log-normal
brightness PDF and a broken power-law power spectrum. Two values are not used
as published but forward-calibrated, because the published estimators are not
the texture's parameters:

* the 2D slope ``beta`` is chosen so that 1D "slit" power spectra of the
  texture reproduce the measured 1D slope in the measured window;
* ``sigma`` is chosen so that the textured *profile*, binned and analysed like
  the observed maps, reproduces the measured PDF width. The smooth profile
  already contributes part of that width, so ``sigma`` depends on ``r_core``.

The angular scales are converted to lengths at the reference-scene distance,
the same convention as ``ORION_R_REF``; the texture then scales with the scene.

References: see references.bib (keys in brackets).
"""
import warnings
from functools import lru_cache

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from scipy.optimize import brentq

from scopesim_targets.extended_source import CoredPowerLaw2D
from scopesim_targets.texture import LognormalTexture

from .nebular_backgrounds import DEFAULT_INDEX

# ---- measured statistics [Arthur2016] -------------------------------------------
# Break scales of the surface-brightness power spectra (Sect. 3.1, 4.4): flat
# (uncorrelated) above 22", a shallow power law between 22" and 8", very steep
# below 8".
ORION_OUTER_ANGLE = 22 * u.arcsec
ORION_INNER_ANGLE = 8 * u.arcsec
# 1D (slit) slopes of the H-alpha thick-slice, i.e. surface-brightness, power
# spectra: regime II (Table 2) and regime III (Table 3).
ORION_GAMMA_T2, ORION_GAMMA_T2_ERR = -2.79, 0.15
ORION_GAMMA_T3 = -4.26
# Width of the log-normal surface-brightness PDF, 0.45-0.6 for all lines
# (Sect. 4.3, Fig. 15); the mid-range is adopted.
ORION_SIGMA_OBS_RANGE = (0.45, 0.6)
ORION_SIGMA_OBS = float(np.mean(ORION_SIGMA_OBS_RANGE))
# How that width was measured (Sect. 4.3, Fig. 15): maps binned 16 x 16 pixels
# of 0.54" (Fig. 8 caption), normalized by the mean, PDF weighted by the
# surface brightness, log-normal fitted where the density exceeds 0.3 x max.
ORION_BIN = 16 * 0.54 * u.arcsec
PDF_CORE_FRACTION = 0.3
# Observed field and position of theta1 Ori C, J2000 (footnote 1).
ORION_FIELD = {"ra": ("5h35m10.4s", "5h35m23.2s"), "dec": ("-5d21m36s", "-5d26m10s")}
THETA1_ORI_C = ("5h35m16.5s", "-5d23m23s")

# ---- texture parameters ------------------------------------------------------------
# Angles expressed as lengths at the reference-scene distance (cf. ORION_R_REF).
SCENE_DISTANCE = 414 * u.pc
# Small-scale tail: isotropic 1D -> 2D conversion, beta_2D = 1 - gamma_1D.
ORION_BETA_INNER = 1 - ORION_GAMMA_T3
# Mid-range slope, forward-calibrated by calibrate_beta(); the plain conversion
# 1 - ORION_GAMMA_T2 = 3.79 would give a 1D slope of -3.0 in the narrow 8-22"
# window, because of the adjacent steep tail.
ORION_BETA = 3.53

# Numerical choices: seeds used for calibrations (statistics, not realizations)
# and histogram bins of the PDF estimator.
CALIBRATION_SEEDS = tuple(range(8))
PDF_BINS = 40


def _angle_to_length(angle: u.Quantity) -> u.Quantity:
    return (angle.to(u.rad).value * SCENE_DISTANCE).to(u.pc)


def slit_slope(texture: LognormalTexture, k_lo, k_hi) -> float:
    """Slope of the mean 1D power spectrum of ``ln f`` along rows ("slits"),
    fitted between ``k_lo`` and ``k_hi`` (in cycles per texture unit)."""
    lf = texture.log_field
    dx = texture.spacing.to_value(texture.unit)
    ps = np.mean(np.abs(np.fft.rfft(lf, axis=1)) ** 2, axis=0)
    k = np.fft.rfftfreq(lf.shape[1], dx)
    m = (k >= k_lo) & (k <= k_hi)
    return float(np.polyfit(np.log(k[m]), np.log(ps[m]), 1)[0])


def calibrate_beta(seeds=CALIBRATION_SEEDS, n_grid: int = 1024) -> float:
    """2D slope whose 1D slit spectra reproduce ORION_GAMMA_T2 over the
    measured window (1/22" to 1/8"), given the measured small-scale tail."""
    k_lo = 1 / ORION_OUTER_ANGLE.to_value(u.arcsec)
    k_hi = 1 / ORION_INNER_ANGLE.to_value(u.arcsec)

    def mean_slope(beta):
        return np.mean([slit_slope(LognormalTexture(
            sigma=1.0, beta=beta, outer_scale=ORION_OUTER_ANGLE,
            inner_scale=ORION_INNER_ANGLE, beta_inner=ORION_BETA_INNER,
            rng_seed=s, n_grid=n_grid), k_lo, k_hi) for s in seeds])

    return brentq(lambda b: mean_slope(b) - ORION_GAMMA_T2, 2.5, 4.5, xtol=1e-3)


def field_grid():
    """Bin centres of the observed field, offsets from theta1 Ori C [arcsec]."""
    centre = SkyCoord(*THETA1_ORI_C)
    corners = SkyCoord(list(ORION_FIELD["ra"]), list(ORION_FIELD["dec"]))
    dra, ddec = centre.spherical_offsets_to(corners)
    (x0, x1), (y0, y1) = sorted(dra.to_value(u.arcsec)), sorted(ddec.to_value(u.arcsec))
    b = ORION_BIN.to_value(u.arcsec)
    return np.arange(x0 + b / 2, x1, b), np.arange(y0 + b / 2, y1, b)


def pdf_core_width(surface_brightness) -> float:
    """Log-normal width of a binned map, measured as in [Arthur2016]."""
    s = np.asarray(surface_brightness, float).ravel()
    s = s / s.mean()
    hist, edges = np.histogram(np.log(s), bins=PDF_BINS, weights=s, density=True)
    centres = 0.5 * (edges[1:] + edges[:-1])
    core = hist > PDF_CORE_FRACTION * hist.max()
    curvature = np.polyfit(centres[core], np.log(hist[core]), 2)[0] if core.sum() >= 3 else 0.0
    if curvature < 0:
        return float(np.sqrt(-1 / (2 * curvature)))
    # No peaked core to fit (far outside the measured regime, e.g. at the
    # edge of a solver bracket): fall back to the brightness-weighted std.
    ln_s = np.log(s)
    mean = np.sum(s * ln_s) / s.sum()
    return float(np.sqrt(np.sum(s * (ln_s - mean) ** 2) / s.sum()))


def _profile_on_field(r_core: u.Quantity, index: float):
    xs, ys = field_grid()
    rc = (r_core / SCENE_DISTANCE * u.rad).to_value(u.arcsec)
    x, y = np.meshgrid(xs, ys)
    return CoredPowerLaw2D(r_core=rc, index=index)(x, y), xs, ys


def profile_width(r_core: u.Quantity, index: float = DEFAULT_INDEX) -> float:
    """PDF width of the smooth profile alone over the observed field."""
    return pdf_core_width(_profile_on_field(r_core, index)[0])


def minimum_r_core(index: float = DEFAULT_INDEX, sigma_obs: float = ORION_SIGMA_OBS) -> u.Quantity:
    """Smallest nebular core for which the smooth profile alone does not
    exceed the measured PDF width (a larger core leaves room for texture)."""
    return brentq(lambda r: profile_width(r * u.pc, index) - sigma_obs, 0.05, 1.0) * u.pc


@lru_cache(maxsize=64)
def _calibrate_sigma_cached(r_core_pc: float, index: float, sigma_obs: float, seeds: tuple) -> float:
    profile, xs, ys = _profile_on_field(r_core_pc * u.pc, index)
    if pdf_core_width(profile) >= sigma_obs:
        return 0.0
    x, y = _angle_to_length(xs * u.arcsec), _angle_to_length(ys * u.arcsec)
    pix = _angle_to_length(ORION_BIN)

    def forward(sigma):
        return np.mean([pdf_core_width(profile * orion_texture(
            r_core_pc * u.pc, s, index=index, sigma=sigma, n_grid=256).sample(x, y, pixel_size=pix))
            for s in seeds])

    return brentq(lambda sg: forward(sg) - sigma_obs, 0.0, 1.2, xtol=2e-3)


def calibrate_sigma(r_core: u.Quantity, index: float = DEFAULT_INDEX,
                    sigma_obs: float = ORION_SIGMA_OBS, seeds=CALIBRATION_SEEDS) -> float:
    """Texture ``sigma`` for which profile x texture, analysed like the observed
    maps, reproduces ``sigma_obs``. Zero (with a warning) if the smooth profile
    alone is already wider, i.e. for cores below :func:`minimum_r_core`."""
    sigma = _calibrate_sigma_cached(float(r_core.to_value(u.pc)), float(index),
                                    float(sigma_obs), tuple(seeds))
    if sigma == 0.0:
        warnings.warn(
            f"the smooth profile with r_core = {r_core:.3g} is already more structured "
            f"than the Orion Nebula [Arthur2016] (PDF width {profile_width(r_core, index):.2f} "
            f">= {sigma_obs:.2f}); no texture added. Cores >= "
            f"{minimum_r_core(index, sigma_obs):.3f} leave room for texture.", stacklevel=2)
    return sigma


def orion_texture(r_core: u.Quantity, rng_seed: int, *, index: float = DEFAULT_INDEX,
                  sigma: float | None = None, n_grid: int = 2048) -> LognormalTexture:
    """Orion-Nebula-like texture for a background with nebular core ``r_core``.

    ``sigma`` defaults to :func:`calibrate_sigma` for this core and index (it
    does not depend on the seed or on the scene distance). The same texture
    can be passed to :func:`~.nebular_backgrounds.ophiuchus_background` for a
    controlled comparison; no Ophiuchus-specific statistics exist.
    """
    if sigma is None:
        sigma = calibrate_sigma(r_core, index)
    return LognormalTexture(
        sigma=sigma, beta=ORION_BETA, beta_inner=ORION_BETA_INNER,
        outer_scale=_angle_to_length(ORION_OUTER_ANGLE),
        inner_scale=_angle_to_length(ORION_INNER_ANGLE),
        rng_seed=rng_seed, n_grid=n_grid)
