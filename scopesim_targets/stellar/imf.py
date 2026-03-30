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


# Kroupa 2002
# 2002Sci...295...82K
params = {
    "a0": .3, "a1": 1.3, "a2": 2.3, "a3": 2.7,
    "m0": .01, "m1": .08, "m2": .5, "m3": 1.,
}
kroupa_gen = BrokenPowerlaw(name="Kroupa02", a=0.005, b=150, momtype=0)
kroupa = kroupa_gen(**params)

# Chabrier 2001
chabrier_gen = LogNormal(name="Chabrier01", a=0.005, b=150)
chabrier = chabrier_gen(mo=716.4, delta=3.3, beta=.25)
