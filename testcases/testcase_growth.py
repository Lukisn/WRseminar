#!/usr/bin/env python3

"""
Demo case for pure growth
(positive = condensation / negative = evaporation).

Taken from Yuan paper Case 4, (5) (both condensation).
(Also possible Cases 1, 2, 3 (evaporation).)
"""

import numpy as np
import matplotlib.pyplot as plt
from math import exp, sqrt
from WR.grid import Grid
from WR.methods import FixedPivot, CellAverage


# TODO: implement first testcase!
def main():

    # initial NDF:
    def f(x):
        #p = 2
        #q = 8
        #return (2**p / 7**q) * (x - 1)**q * (15 - x)**q
        return 6 * x**3 * exp(-x)

    START, END = 0, 50
    SECTIONS = 100
    FACTOR = 1.25

    initial_ndf = Grid.create_geometric_end(
        start=START, end=END, sections=SECTIONS, factor=FACTOR, func=f
    )

    # Growth function:
    def G(v):
        #K = 0.78
        #return K / v
        return v / 2

    # Simulation:
    T0, TEND = 0, 1
    STEPS = 10
    EVERY = 1
    ORDER = 1

    # setup methods and do simulation:
    # Fixed Pivot Method:
    fp = FixedPivot(
        initial=initial_ndf,
        gro=True, gro_rate=G,
        bre=False, agg=False, nuc=False
    )
    fp.simulate(
        start_time=T0, end_time=TEND, steps=STEPS,
        write_every=EVERY, max_order=ORDER
    )
    # Cell Average Technique:
    ca = CellAverage(
        initial=initial_ndf,
        gro=True, gro_rate=G,
        bre=False, agg=False, nuc=False
    )
    ca.simulate(
        start_time=T0, end_time=TEND, steps=STEPS,
        write_every=EVERY, max_order=ORDER
    )

    # analytic solution:
    def n(t, x):
        #K = 0.78
        #a = sqrt(1 + 2*K*t)
        #b = sqrt(15**2 + 2*K*t)
        #if a < x < b:
        #    root = sqrt(x**2 - 2*K*t)
        #    return f(root)* x / root
        #else:
        #    return 0
        num = (x * exp(-t / 2))**3 * exp(-x * exp(-t / 2))
        den = 6 * exp(t / 2)
        return num / den

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
