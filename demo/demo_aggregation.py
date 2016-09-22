#!/usr/bin/env python3

"""
Demo case for pure aggregation (coagulation).

Taken from Yuan paper case 9.
"""

import numpy as np
import matplotlib.pyplot as plt
from math import exp
from WR.grid import Grid
from WR.methods import FixedPivot, CellAverage

# TODO: implement test case!
def main():
    # initial NDF:
    N0 = 1
    v0 = 1

    def f(x, N0=N0, v0=v0):
        return (3 * N0 / v0) * x ** 2 * exp(-x ** 3 / v0)

    START, END = 0, 1e5
    SECTIONS = 100
    FACTOR = 1.25

    initial_ndf = Grid.create_geometric_end(
        start=START, end=END, sections=SECTIONS, factor=FACTOR, func=f
    )

    # Breakage function:
    Q0 = 1

    def Q(v1, v2, Q0=Q0):
        return Q0

    # Simulation:
    T0 = 0
    TEND = 1
    STEPS = 10
    EVERY = None
    ORDER = 1

    # setup method and do simulation:
    # Fixed Pivot Method:
    fp = FixedPivot(
        initial=initial_ndf,
        agg=True, agg_freq=Q,
        bre=False, gro=False, nuc=False
    )
    fp.simulate(
        start_time=T0, end_time=TEND, steps=STEPS,
        write_every=EVERY, max_order=ORDER
    )
    # Cell Average Technique:
    ca = CellAverage(
        initial=initial_ndf,
        agg=True, agg_freq=Q,
        bre=False, gro=False, nuc=False
    )
    ca.simulate(
        start_time=T0, end_time=TEND, steps=STEPS,
        write_every=EVERY, max_order=ORDER
    )

    # analytic solution
    def n(t, x, N0=N0, v0=v0, Q0=Q0):
        Ta = N0 * Q0 * t
        first = ((12 * N0 * x**2) / (v0 * (Ta + 2)**2))
        inner = (-2 * x**3) / (v0 * (Ta + 2))
        return first * np.exp(inner)

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

    plt.xlim(1e-5, 1e5)
    plt.ylim(1e-10, 1e5)
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(loc="best", fontsize="small")
    plt.grid()
    plt.show()


if __name__ == "__main__":
    main()
