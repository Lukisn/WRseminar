#!/usr/bin/env python3

"""
Testcase for pure growth.

Taken from Yuan paper Case 1.
"""

import numpy as np
import matplotlib.pyplot as plt
from math import exp
from WR.grid import Grid
from WR.methods import FixedPivot, CellAverage


# TODO: find test case!
def main():
    # initial NDF:
    def f(x):
        res = 60 * x ** 2 * (1 - x) ** 3
        if res < 0:
            return 0
        else:
            return res

    START, END = 0, 1
    SECTIONS = 100
    FACTOR = 1.25

    initial_ndf = Grid.create_geometric_end(
        start=START, end=END, sections=SECTIONS, factor=FACTOR, func=f
    )

    # Growth function:
    def Q(v):
        return -0.5

    # Simulation:
    T0, TEND = 0, 0.5
    STEPS = 50
    EVERY = None
    ORDER = 1

    # setup methods and do simulation:
    # Fixed Pivot Method:
    fp = FixedPivot(
        initial=initial_ndf,
        gro=True, gro_rate=Q,
        bre=False, agg=False, nuc=False
    )
    fp.simulate(
        start_time=T0, end_time=TEND, steps=STEPS,
        write_every=EVERY, max_order=ORDER
    )
    # Cell Average Technique:
    ca = CellAverage(
        initial=initial_ndf,
        gro=True, gro_rate=Q,
        bre=False, agg=False, nuc=False
    )
    ca.simulate(
        start_time=T0, end_time=TEND, steps=STEPS,
        write_every=EVERY, max_order=ORDER
    )

    # analytic solution:
    def n(t, x):
        return max(
            60 * (x + t/2) ** 2 * (1 - (x + t/2)) ** 3,
            0
        )

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

    plt.xlim(0, 1)
    plt.ylim(0, 5)
    plt.xscale("linear")
    plt.yscale("linear")
    plt.legend(loc="best", fontsize="small")
    plt.grid()
    plt.show()


if __name__ == "__main__":
    main()
