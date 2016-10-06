#!/usr/bin/env python3

"""
Demo case for pure growth
"""

import numpy as np
import matplotlib.pyplot as plt
from math import exp, sqrt, pi
from WR.grid import Grid
from WR.methods import FixedPivot, CellAverage


# TODO: check moments!
# TODO: find out why ca has negative numeric diffusion!?!?!!!
def main():

    # initial NDF:
    def f(x):
        sig = 0.5
        mu = 2
        return 1 / (sqrt(2 * pi * sig**2)) * exp(-(x - mu)**2 / (2 * sig**2))

    START, END = 0, 10
    SECTIONS = 100
    FACTOR = 1.1

    initial_ndf = Grid.create_geometric_end(
        start=START, end=END, sec=SECTIONS, fact=FACTOR, func=f
    )

    # Growth function:
    def Q(v):
        return 1

    # Simulation:
    T0, TEND = 0, 1
    STEPS = 10
    EVERY = 1
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
        sig = 0.5
        mu = 2 + t
        return 1 / (sqrt(2 * pi * sig ** 2)) * exp(
            -(x - mu) ** 2 / (2 * sig ** 2))

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

    #plt.xlim(0, 1e3)
    #plt.ylim(0, 5)
    plt.xscale("linear")
    plt.yscale("linear")
    plt.legend(loc="best", fontsize="small")
    plt.grid()
    plt.show()

    # plot moments comparison:
    times = sorted(fp.result_moments)
    moments = {"fp": {}, "ca": {}}
    for method in moments.keys():
        for order in range(ORDER + 1):
            moments[method][order] = []
            for time in times:
                moments[method][order].append(fp.result_moments[time][order])
    print(moments)

    plt.xlabel("time")
    plt.ylabel("moment")
    for method in moments.keys():
        if method == "ca":
            symbol = ".-"
        else:  # method == "fp"
            symbol = "x-"
        for order in range(ORDER + 1):
            plt.plot(times, moments[method][order], symbol,
                     label="{}-moment{}".format(method, order))

    plt.legend(loc="best", fontsize="small")
    plt.grid()
    plt.show()


if __name__ == "__main__":
    main()
