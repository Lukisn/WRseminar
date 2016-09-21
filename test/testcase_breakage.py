#!/usr/bin/env python3

"""
Testcase for pure breakage.

Taken from ...
"""

import numpy as np
import matplotlib.pyplot as plt
from math import exp
from WR.grid import Grid
from WR.methods import FixedPivot, CellAverage


def step(x, L):
    if x < L:
        return 1
    else:
        return 0


def num_delta(x, L, eps=1e-1):
    assert eps > 0

    min = L - eps
    max = L + eps

    if min <= x <= max:
        return 1
    else:
        return 0


def main():
    # TODO: tweak initial delta function representation!
    # initial NDF:
    L = 1

    def f(x, L=L):
        return num_delta(x, L)

    START, END = 0, 1e1
    SECTIONS = 100
    FACTOR = 1.25

    initial_ndf = Grid.create_geometric_end(
        start=START, end=END, sections=SECTIONS, factor=FACTOR, func=f
    )

    # Breakage functions:
    def Gamma(v):
        return v

    def beta(v1, v2):
        return 2 / v2

    # Simulation:
    T0 = 0
    TEND = 1
    STEPS = 10
    EVERY = None
    ORDER = 1

    # setup methods and do simulation:
    # Fixed Pivot Method:
    fp = FixedPivot(
        initial=initial_ndf,
        bre=True, bre_freq=Gamma, child=beta,
        agg=False, gro=False, nuc=False
    )
    fp.simulate(
        start_time=T0, end_time=TEND, steps=STEPS,
        write_every=EVERY, max_order=ORDER
    )
    # Cell Average Technique:
    ca = CellAverage(
        initial=initial_ndf,
        bre=True, bre_freq=Gamma, child=beta,
        agg= False, gro=False, nuc=False
    )
    ca.simulate(
        start_time=T0, end_time=TEND, steps=STEPS,
        write_every=EVERY, max_order=ORDER
    )

    # analytic solution:
    def n(t, x, L=L):
        brace = num_delta(x, L) + (2 * t + t**2 * (L - x)) * step(x, L)
        return exp(-t * x) * brace

    # plot comparison:
    ana_x = initial_ndf.pivots()
    ana_y = []
    for x in ana_x:
        ana_y.append(n(TEND, x))
    plt.plot(ana_x, ana_y, "-", label="analytic")

    ini_x = initial_ndf.pivots()
    ini_y = initial_ndf.densities()
    plt.plot(ini_x, ini_y, label="initial")

    fp_x = fp.result_ndfs[TEND].pivots()
    fp_y = fp.result_ndfs[TEND].densities()
    plt.plot(fp_x, fp_y, ".-", label="fixed pivot")

    ca_x = ca.result_ndfs[TEND].pivots()
    ca_y = ca.result_ndfs[TEND].densities()
    plt.plot(ca_x, ca_y, ".-", label="cell average")

    plt.xlim(1e-5, 1e1)
    #plt.ylim(1e-10, 1e5)
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(loc="best", fontsize="small")
    plt.grid()
    plt.show()


if __name__ == "__main__":
    main()
