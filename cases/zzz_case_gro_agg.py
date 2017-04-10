#!/usr/bin/env python3
"""
DEPRECATED CASE
Demo case for negative growth (evaporation) and aggregation (coalescence).

Taken from Yuan paper Case 13.
"""
# standard library imports:
import os
from math import exp
# third party imports:
import matplotlib.pyplot as plt
# application imports:
from sectional.grid import Grid
from sectional.methods import FixedPivot, CellAverage


# TODO: find problem in the calculation!
def main():
    """Main function."""

    # PROBLEM FUNCTIONS: ------------------------------------------------------

    KE = 1
    KC = 0.05

    # initial NDF:
    def f(x):
        return x**2 * exp(-x)

    # Growth function:
    def G(v):
        return -KE + v**(1/3)

    # Aggregation function:
    def Q(v1, v2):
        v1 = v1**(1/3)
        v2 = v2**(1/3)
        return KC * (v1 + v2) * (1/v1 + 1/v2)

    # CONSTANTS: --------------------------------------------------------------

    # Grid:
    START, END = 1e-3, 1e3
    SECTIONS = 100
    FACTOR = 1.1

    # Simulation:
    T0, TEND = 0, 0.1
    STEPS = 10
    EVERY = 1
    ORDER = 1

    # Plotting:
    XSCALE, YSCALE = "log", "log"
    XMIN, XMAX = 1e-3, 1e3
    YMIN, YMAX = 1e-10, 1e2

    # File output:
    WRITE_DATA_FILES = False
    WRITE_PLOT_FILES = True
    FOLDER = os.path.abspath("./results/")
    if not os.path.exists(FOLDER):
        os.makedirs(FOLDER)

    # SIMULATION: -------------------------------------------------------------

    # initial NDF:
    initial_ndf = Grid.create_geometric_end(
        start=START, end=END, sec=SECTIONS, fact=FACTOR, func=f
    )
    if WRITE_DATA_FILES:
        initial_ndf.to_file(os.path.join(FOLDER, "gro_agg_initial_ndf.dat"))

    # Fixed Pivot Method:
    fp = FixedPivot(
        initial=initial_ndf,
        gro=True, gro_rate=G,
        agg=True, agg_freq=Q,
        bre=False, nuc=False
    )
    fp.simulate(start_time=T0, end_time=TEND, steps=STEPS, write_every=EVERY)
    if WRITE_DATA_FILES:
        fp.moments_to_file(
            os.path.join(FOLDER, "gro_agg_fp_moments.dat"), max_order=ORDER
        )
        fp.ndf_to_files(os.path.join(FOLDER, "gro_agg_fp_ndf.dat"))

    # Cell Average Technique:
    ca = CellAverage(
        initial=initial_ndf,
        gro=True, gro_rate=G,
        agg=True, agg_freq=Q,
        bre=False, nuc=False
    )
    ca.simulate(start_time=T0, end_time=TEND, steps=STEPS, write_every=EVERY)
    ca.moments_to_file("results/gro_agg_ca_moments.dat", max_order=ORDER)
    ca.ndf_to_files("results/gro_agg_ca_ndf.dat")

    # PLOTTING: ---------------------------------------------------------------

    # gather NDF data:
    ini_x, ini_y = initial_ndf.pivots(), initial_ndf.densities()
    """
    ana_x, ana_y = initial_ndf.pivots(), []
    for x in ana_x:
        ana_y.append(n(TEND, x))
    """
    fp_x, fp_y = fp._result_ndfs[TEND].pivots(), fp._result_ndfs[
        TEND].densities()
    ca_x, ca_y = ca._result_ndfs[TEND].pivots(), ca._result_ndfs[
        TEND].densities()
    """
    # calculate errors:
    fp_err_y, ca_err_y = [], []
    for i, x in enumerate(fp_x):
        err = fp_y[i] - n(TEND, x)
        fp_err_y.append(err)
    for i, x in enumerate(ca_x):
        err = ca_y[i] - n(TEND, x)
        ca_err_y.append(err)
    """
    # gather moment data:
    times = sorted(fp._result_ndfs.keys())
    fp_moment0, fp_moment1 = [], []
    ca_moment0, ca_moment1 = [], []
    for time in times:
        fp_moment0.append(fp._result_ndfs[time].moment(0))
        fp_moment1.append(fp._result_ndfs[time].moment(1))
        ca_moment0.append(ca._result_ndfs[time].moment(0))
        ca_moment1.append(ca._result_ndfs[time].moment(1))

    # plot NDF comparison and errors:
    # upper subplot: NDF:
    plt.subplot(211)
    # plt.xlabel("size")
    plt.ylabel("NDF")
    #plt.plot(ana_x, ana_y, "y-", lw=3, label="analytic")
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
    #plt.plot(fp_x, fp_err_y, "bx-", label="fixed pivot")
    #plt.plot(ca_x, ca_err_y, "r.-", lw=2, label="cell average")
    plt.xlim(XMIN, XMAX)
    # plt.ylim(YMIN, YMAX)  # comment out = auto
    plt.xscale(XSCALE)
    # plt.yscale(YSCALE)  # comment out  = auto
    plt.legend(loc="best", fontsize="small")
    plt.grid()

    if WRITE_PLOT_FILES:
        plt.savefig(os.path.join(FOLDER, "gro_agg_ndf.eps"))
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

    if WRITE_PLOT_FILES:
        plt.savefig(os.path.join(FOLDER, "gro_agg_mom.eps"))
    plt.show()


if __name__ == "__main__":
    main()