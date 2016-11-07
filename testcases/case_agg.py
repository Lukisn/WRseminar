#!/usr/bin/env python3

"""
Demo case for pure aggregation (coagulation).

Taken from Yuan paper case 9.
"""
import matplotlib.pyplot as plt
from math import exp, sqrt
from scipy.special import iv  # modified bessel function
from WR.grid import Grid
from WR.methods import FixedPivot, CellAverage


def main():
    """Main function.
    """

    # PROBLEM FUNCTIONS: ------------------------------------------------------

    # initial NDF:
    def f(x):
        return exp(-x)

    # Aggregation function:
    def Q(v1, v2):
        return v1 + v2

    # analytic solution
    def n(t, x):
        ratio = (exp(-t - 2 * x + x * exp(-t))) / (x * sqrt(1 - exp(-t)))
        bessel = iv(1, 2 * x * sqrt(1 - exp(-t)))
        return ratio * bessel

    # CONSTANTS: --------------------------------------------------------------

    # Grid:
    START, END = 0, 1e5
    SECTIONS = 100
    FACTOR = 1.2

    # Simulation:
    T0, TEND = 0, 1
    STEPS = 10
    EVERY = 1
    ORDER = 1

    # Plotting:
    XSCALE, YSCALE = "log", "log"
    XMIN, XMAX = 1e-5, 1e5
    YMIN, YMAX = 1e-10, 1e5

    # SIMULATION: -------------------------------------------------------------

    # initial NDF:
    initial_ndf = Grid.create_geometric_end(
        start=START, end=END, sec=SECTIONS, fact=FACTOR, func=f
    )

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

    # PLOTTING: ---------------------------------------------------------------

    # gather NDF data:
    ini_x, ini_y = initial_ndf.pivots(), initial_ndf.densities()
    ana_x, ana_y = initial_ndf.pivots(), []
    for x in ana_x:
        ana_y.append(n(TEND, x))
    fp_x, fp_y = fp._result_ndfs[TEND].pivots(), fp._result_ndfs[
        TEND].densities()
    ca_x, ca_y = ca._result_ndfs[TEND].pivots(), ca._result_ndfs[
        TEND].densities()

    # calculate errors:
    fp_err_y, ca_err_y = [], []
    for i, x in enumerate(fp_x):
        err = fp_y[i] - n(TEND, x)
        fp_err_y.append(err)
    for i, x in enumerate(ca_x):
        err = ca_y[i] - n(TEND, x)
        ca_err_y.append(err)

    # gather moment data:
    times = sorted(fp._result_moments)  # == sorted(ca.result_moments)
    fp_moment0, fp_moment1 = [], []
    ca_moment0, ca_moment1 = [], []
    for time in times:
        fp_moment0.append(fp._result_moments[time][0])
        fp_moment1.append(fp._result_moments[time][1])
        ca_moment0.append(ca._result_moments[time][0])
        ca_moment1.append(ca._result_moments[time][1])

    # plot NDF comparison and errors:
    # upper subplot: NDF:
    plt.subplot(211)
    # plt.xlabel("size")
    plt.ylabel("NDF")
    plt.plot(ana_x, ana_y, "y-", lw=3, label="analytic")
    plt.plot(ini_x, ini_y, "g.-", lw=2, label="initial")
    plt.plot(fp_x, fp_y, "bx-", label="fixed pivot")
    plt.plot(ca_x, ca_y, "r.-", lw=2, label="cell average")
    plt.xlim(XMIN, XMAX)
    plt.ylim(YMIN, YMAX)
    plt.xscale(XSCALE)
    plt.yscale(YSCALE)
    plt.legend(loc="best", fontsize="small")
    plt.grid()
    # lower subplot: errors:
    plt.subplot(212)
    plt.xlabel("size")
    plt.ylabel("error")
    plt.plot(fp_x, fp_err_y, "bx-", label="fixed pivot")
    plt.plot(ca_x, ca_err_y, "r.-", lw=2, label="cell average")
    plt.xlim(XMIN, XMAX)
    # plt.ylim(YMIN, YMAX)  # comment out = auto
    plt.xscale(XSCALE)
    # plt.yscale(YSCALE)  # comment out  = auto
    plt.legend(loc="best", fontsize="small")
    plt.grid()

    plt.savefig("agg_ndf.eps")
    plt.show()

    # plot Moments comparison:
    plt.xlabel("time")
    plt.ylabel("moment")
    plt.plot(times, fp_moment0, "bx-", label="fp_m0")
    plt.plot(times, fp_moment1, "b.-", label="fp_m1")
    plt.plot(times, ca_moment0, "rx-", lw=2, label="ca_m0")
    plt.plot(times, ca_moment1, "r.-", lw=2, label="ca_m1")
    plt.legend(loc="best", fontsize="small")
    plt.grid()

    plt.savefig("agg_mom.eps")
    plt.show()

if __name__ == "__main__":
    main()
