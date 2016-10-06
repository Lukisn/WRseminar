#!/usr/bin/env python3

"""
Demo case for pure growth

Taken from Yuan paper Case 4.
"""
import matplotlib.pyplot as plt
from math import exp, sqrt
from WR.grid import Grid
from WR.methods import FixedPivot, CellAverage


# TODO: find problem in the calculation!
def main():
    """Main function.
    """

    # PROBLEM FUNCTIONS: ------------------------------------------------------

    # initial NDF:
    def f(x):
        return 6 * x**3 * exp(-x)

    # Growth function:
    def G(v):
        return v / 2

    # analytic solution:
    def n(t, x):
        num = (x * exp(-t / 2)) ** 3 * exp(-x * exp(-t / 2))
        den = 6 * exp(t / 2)
        return num / den

    # CONSTANTS: --------------------------------------------------------------

    # Grid:
    START, END = 0, 2e1
    SECTIONS = 100
    FACTOR = 1.2

    # Simulation:
    T0, TEND = 0, 1
    STEPS = 100
    EVERY = 1
    ORDER = 1

    # Plotting:
    XSCALE, YSCALE = "linear", "linear"
    XMIN, XMAX = 0, 2e1
    YMIN, YMAX = 0, 1

    # SIMULATION: -------------------------------------------------------------

    # initial NDF:
    initial_ndf = Grid.create_geometric_end(
        start=START, end=END, sec=SECTIONS, fact=FACTOR, func=f
    )

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

    # PLOTTING: ---------------------------------------------------------------

    # plot comparison:
    plt.subplot(211)
    plt.ylabel("NDF")
    ana_x, ana_y = initial_ndf.pivots(), []
    for x in ana_x:
        ana_y.append(n(TEND, x))
    plt.plot(ana_x, ana_y, "-", label="analytic")
    ini_x, ini_y = initial_ndf.pivots(), initial_ndf.densities()
    plt.plot(ini_x, ini_y, ".-", label="initial")
    fp_x, fp_y = fp.result_ndfs[TEND].pivots(), fp.result_ndfs[
        TEND].densities()
    plt.plot(fp_x, fp_y, "x-", label="fixed pivot")
    ca_x, ca_y = ca.result_ndfs[TEND].pivots(), ca.result_ndfs[
        TEND].densities()
    plt.plot(ca_x, ca_y, ".-", label="cell average")
    plt.xlim(XMIN, XMAX)
    plt.ylim(YMIN, YMAX)
    plt.xscale(XSCALE)
    plt.yscale(YSCALE)
    plt.legend(loc="best", fontsize="small")
    plt.grid()

    # calculate errors:
    fp_err_y, ca_err_y = [], []
    for i, x in enumerate(fp_x):
        err = fp_y[i] - n(TEND, x)
        fp_err_y.append(err)
    for i, x in enumerate(ca_x):
        err = ca_y[i] - n(TEND, x)
        ca_err_y.append(err)

    # plot errors:
    plt.subplot(212)
    plt.xlabel("size")
    plt.ylabel("error")
    plt.plot(fp_x, fp_err_y, "x-", label="fixed pivot")
    plt.plot(ca_x, ca_err_y, ".-", label="cell average")
    plt.xlim(XMIN, XMAX)
    # plt.ylim(YMIN, YMAX)
    plt.xscale(XSCALE)
    # plt.yscale(YSCALE)
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
