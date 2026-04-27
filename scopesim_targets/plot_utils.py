# -*- coding: utf-8 -*-
"""Plotting convenience functinality.

Example usage
-------------
>>> from matplotlib import pyplot as plt
>>> from scopesim_targets.spectral_classes import StellarParameters
>>> from scopesim_targets.plot_utils import spec_classes_axis, figure_factory_hrd
>>>
>>> x = np.linspace(2000, 15000, 400)
>>> y = np.sin(x / 1000)  # just to have something to plot...
>>>
>>> plt.style.use("seaborn-v0_8")  # -talk, -poster, -paper
>>> fig, ax = figure_factory_hrd()
>>> ax.plot(x, y)  # doctest: +ELLIPSIS
[<matplotlib.lines.Line2D object at ...>]
>>>
>>> stp = StellarParameters()
>>> spclss = list(stp.group_spectral_classes())
>>> spec_classes_axis(ax, spclss)

"""

from collections.abc import Generator

import numpy as np
from matplotlib import axes
from matplotlib import pyplot as plt
from matplotlib.patches import Circle

from .spectral_classes import SpectralClass, TeffRange, teff_range_overlap


BBOX_SPEC_TICK_LABELS = {
    "boxstyle": "round",
    "ec": "black",
    "fc": "gray",
    "alpha": 0.8,  # 1
}

SPEC_AXSPAN_ALPHA = 1  # 0.5

SPEC_TICK_PARAMS = {
    "major": {
        "length": 0,  # REQUIRED: no tick under label
        "pad": 8,
    },
    "minor": {
        "length": 6,
        "width": 1.0,
    },
}


def figure_factory(nrows=1, ncols=1, **kwargs):
    """Generate default fig and ax, to easily modify later."""
    fig, ax = plt.subplots(nrows, ncols, **kwargs)
    return fig, ax


def figure_factory_hrd(**kwargs):
    """Wrap ``figure_factory`` for single axes with inverted x-axis."""
    fig, ax = figure_factory(**kwargs)
    ax.invert_xaxis()  # HRD-style
    return fig, ax


def _shade_spec_classes(ax: axes.Axes, spec_classes: list[SpectralClass]):
    """Shade areas of plot covered by spectral class Teff range."""
    teff_range = TeffRange(*sorted(ax.get_xlim()))

    for spcls in spec_classes:
        if not spcls.is_in_range(teff_range):
            continue

        ax.axvspan(
            *teff_range_overlap(spcls.teff_range, teff_range),
            color=spcls.color,
            alpha=SPEC_AXSPAN_ALPHA,
            zorder=-10,
        )


def _get_midpoints(
    spec_classes: list[SpectralClass],
    teff_range: TeffRange,
) -> Generator[float]:
    """Collect spec class range midpoints, considering globals plot limits."""
    for spcls in spec_classes:
        try:
            yield sum(teff_range_overlap(spcls.teff_range, teff_range)) / 2
        except ValueError:  # no overlap, i.e. out of range
            yield spcls.midpoint


def _get_xticks(
    spec_classes: list[SpectralClass],
    teff_range: TeffRange,
) -> np.ndarray:
    """Collect min and max points, remove duplicates, clip to globals."""
    xticks = np.unique(np.clip(
        [spcls.teff_range for spcls in spec_classes],
        a_min=teff_range.min,
        a_max=teff_range.max,
    ))
    return xticks


def spec_classes_axis(
    ax: axes.Axes,
    spec_classes: list[SpectralClass],
    shade_area: bool = True,
) -> None:
    """Mark spectral classes in plot and optionally shade regions."""
    teff_range = TeffRange(*sorted(ax.get_xlim()))

    if shade_area:
        _shade_spec_classes(ax, spec_classes)

    ax_top = ax.secondary_xaxis("top")

    ax_top.set_xticks(
        list(_get_midpoints(spec_classes, teff_range)),
        [spcls.name for spcls in spec_classes],
        weight="semibold",
    )

    # Can't set directly
    for tick, spcls in zip(ax_top.get_xticklabels(), spec_classes):
        if spcls.color is not None:
            tick.set_color(spcls.color)
        tick.set_bbox(BBOX_SPEC_TICK_LABELS)

    ax_top.set_xticks(_get_xticks(spec_classes, teff_range), minor=True)

    for which in ("major", "minor"):
        ax_top.tick_params(
            axis="x",
            which=which,
            **SPEC_TICK_PARAMS[which]
        )


def draw_circle(
    ax,
    center: tuple[float, float] = (0, 0),
    radius: float = 10,
    edgecolor: str | None = None,
    facecolor: str | None = None,
    linewidth: float | None = None,
    linestyle: str | None = None,
    fill: bool = False,
    label: str = "",
) -> None:
    """Draw a circle in data coordinates."""
    ax.add_artist(Circle(
        center,
        radius,
        edgecolor=edgecolor,
        facecolor=facecolor,
        linewidth=linewidth,
        linestyle=linestyle,
        fill=fill,
        label=label,
    ))
