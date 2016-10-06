#!/usr/bin/env python3

"""
Demo case for coupled breakage and aggregation.

Taken from Yuan paper case 12.
"""

import matplotlib.pyplot as plt
from math import exp
from WR.cmdline import step, hstep
from WR.grid import Grid
from WR.methods import FixedPivot, CellAverage

# TODO: implement case 12!
def num_delta(x, eps=1e-1):
    assert eps > 0

    min = 0 - eps
    max = 0 + eps
    if min <= x <= max:
        return 1
    else:
        return 0


def main():

    # initial NDF:
    def f(x):
        return num_delta(x - 1)

    START, END = 0, 2
    SECTIONS = 100
    FACTOR = 1.25

    initial_ndf = Grid.create_geometric_end(
        start=START, end=END, sec=SECTIONS, fact=FACTOR, func=f
    )
    #initial_ndf = Grid.create_uniform(
    #    start=START, end=END, sections=SECTIONS, func=f
    #)

    # Breakage functions:
    def Gamma(v):
        return v ** 2

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
    def n(t, x):
        brace = num_delta(x - 1) + 2 * t * hstep(1 - x)
        return exp(-t * x ** 2) * brace

    # plot comparison:
    ana_x = initial_ndf.pivots()
    ana_y = []
    for x in ana_x:
        ana_y.append(n(TEND, x))
    plt.plot(ana_x, ana_y, "-", label="analytic")

    ini_x = initial_ndf.pivots()
    ini_y = initial_ndf.densities()
    plt.plot(ini_x, ini_y, ".-", label="initial")

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
