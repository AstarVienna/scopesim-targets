# -*- coding: utf-8 -*-
"""Various initial mass function (IMF) laws as scipy distributions."""

import numpy as np
from scipy.stats import rv_continuous


class BrokenPowerlaw(rv_continuous):
    def _pdf(self, x, a0, a1, a2, a3, m0, m1, m2, m3):
        m = [m0, m1, m2, m3]
        alpha = [a0, a1, a2, a3]
        mass_ranges = [
            (m[0] < x) & (x <= m[1]),
            (m[1] < x) & (x <= m[2]),
            (m[2] < x) & (x <= m[3]),
            (m[3] < x)
        ]

        def _prod(n):
            return np.prod(
                [(m[i] / m[i-1])**-alpha[i-1] for i in range(2, n+1)],
                axis=0,
            )

        powerlaws = [
            (x / m[1])**-alpha[0],  # n=0
            (x / m[1])**-alpha[1],  # n=1
            _prod(2) * (x / m[2])**-alpha[2],  # n=2
            _prod(3) * (x / m[3])**-alpha[3],  # n=3
        ]

        return np.select(mass_ranges, powerlaws, default=0.0) * 4.2


class LogNormal(rv_continuous):
    def _pdf(self, x, mo, delta, beta):
        return x**-delta * np.exp(-(mo / x)**beta) * 15


# FIXME: These get auto-generated docstrings from the base class, so they should
#        show up on RTD, but don't. See "attributes" in the module template...
# TODO: Consider renaming these to something more generic if the PDFs are used
#       by other laws as well. Find out more about upper and lower bounds
#       of each IMF before that!! Can those ("support") be set in the "freezing"
#       step afterwards??
kroupa_gen = BrokenPowerlaw(name="Kroupa02", a=0.005, b=150, momtype=0)
chabrier_gen = LogNormal(name="Chabrier01", a=0.005, b=150)


def load_default_imfs() -> dict[str, rv_continuous]:
    """Return "frozen" IMF distributions with default values."""
    imfs = {}

    # Kroupa 2002
    # 2002Sci...295...82K
    kroupa_02_params = {
        "a0": .3, "a1": 1.3, "a2": 2.3, "a3": 2.7,
        "m0": .01, "m1": .08, "m2": .5, "m3": 1.,
    }
    imfs["kroupa02"] = kroupa_gen(**kroupa_02_params)

    # Chabrier 2001
    chabrier_01_params = {
        "mo": 716.4, "delta": 3.3, "beta": .25,
    }
    imfs["chabrier01"] = chabrier_gen(**chabrier_01_params)

    return imfs


DEFAULT_IMFS = load_default_imfs()
