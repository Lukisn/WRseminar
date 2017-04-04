#!/usr/bin/env python3
"""
Demo case for pure growth
"""
# standard library imports:
import os
# third party imports:
import matplotlib.pyplot as plt
from scipy import inf
from scipy.integrate import quad
# application imports:
from sectional.functions import hstep
from sectional.grid import Grid
from sectional.methods import FixedPivot, CellAverage
from util.plotting import *  # plotting styles


def main():
    """Main function."""

    # PROBLEM FUNCTIONS: ------------------------------------------------------

    def f(x):  # initial NDF
        return hstep(x - 0.1) - hstep(x - 0.6)

    def G(v):  # growth function
        return 1

    def n(t, x):  # analytic solution
        return hstep(x - (0.1 + t)) - hstep(x - (0.6 + t))

    # CONSTANTS: --------------------------------------------------------------

    # Grid:
    START, END = 0, 2
    SECTIONS = 100
    FACTOR = 1.0

    # Simulation:
    T0, TEND = 0, 1
    STEPS = 100
    TIME_STEP = (TEND - T0) / STEPS
    EVERY = 1
    ORDER = 1

    # Plotting:
    XSCALE, YSCALE = "linear", "linear"
    XMIN, XMAX = 0, 2
    YMIN, YMAX = 0, 1.5

    # File output:
    WRITE_DATA_FILES = True
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
        initial_ndf.to_file(os.path.join(FOLDER, "growth_initial_ndf.dat"))

    # Fixed Pivot Method:
    fp = FixedPivot(
        initial=initial_ndf,
        gro=True, gro_rate=G,
        bre=False, agg=False, nuc=False
    )
    fp.simulate(start_time=T0, end_time=TEND, steps=STEPS, write_every=EVERY)
    if WRITE_DATA_FILES:
        fp.moments_to_file(
            os.path.join(FOLDER, "growth_fp_moments.dat"), max_order=ORDER
        )
        fp.ndf_to_files(os.path.join(FOLDER, "growth_fp_ndf.dat"))

    # Cell Average Technique:
    ca = CellAverage(
        initial=initial_ndf,
        gro=True, gro_rate=G,
        bre=False, agg=False, nuc=False
    )
    ca.simulate(start_time=T0, end_time=TEND, steps=STEPS, write_every=EVERY)
    if WRITE_DATA_FILES:
        ca.moments_to_file(
            os.path.join(FOLDER, "growth_ca_moments.dat"), max_order=ORDER
        )
        ca.ndf_to_files(os.path.join(FOLDER, "growth_ca_ndf.dat"))

    # CALCULATIONS FOR PLOTTING: ----------------------------------------------

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
    times = sorted(fp._result_ndfs.keys())
    fp_moment0, fp_moment1 = [], []
    ca_moment0, ca_moment1 = [], []
    for time in times:
        fp_moment0.append(fp._result_ndfs[time].moment(0))
        fp_moment1.append(fp._result_ndfs[time].moment(1))
        ca_moment0.append(ca._result_ndfs[time].moment(0))
        ca_moment1.append(ca._result_ndfs[time].moment(1))

    # PLOTTING: ---------------------------------------------------------------

    # plot NDF comparison and errors:
    fig, (upper, lower) = plt.subplots(2, 1, sharex="all")
    # upper subplot - NDF:
    upper.set_title(
        "analytical and discrete NDFs at TEND = {} s, STEP = {} s".format(
            TEND, TIME_STEP
        )
    )
    upper.set_ylabel("NDF")
    upper.plot(ana_x, ana_y, label="ana", **ana_style)
    upper.plot(ini_x, ini_y, label="ini", **initial_style)
    upper.plot(fp_x, fp_y, label="FP", **fp_style)
    upper.plot(ca_x, ca_y, label="CA", **ca_style)
    upper.set_xlim(XMIN, XMAX)
    upper.set_ylim(YMIN, YMAX)
    upper.set_xscale(XSCALE)
    upper.set_yscale(YSCALE)
    upper.legend(**legend_style)
    upper.grid()
    # lower subplot: errors:
    lower.set_xlabel("particle size")
    lower.set_ylabel("relative error")
    lower.plot(fp_x, fp_err_y, label="FP", **fp_style)
    lower.plot(ca_x, ca_err_y, label="CA", **ca_style)
    lower.set_xlim(XMIN, XMAX)
    lower.set_xscale(XSCALE)
    lower.legend(**legend_style)
    lower.grid()

    # tighten layout and show:
    fig.tight_layout()
    if WRITE_PLOT_FILES:
        plt.savefig(os.path.join(FOLDER, "growth_ndf.eps"))
    plt.show()

    # plot Moments comparison:
    plt.xlabel("time in s")
    plt.ylabel("moment")
    plt.plot(times, fp_moment0, label="FP 0", **fp_style0)
    plt.plot(times, fp_moment1, label="FP 1", **fp_style1)
    plt.plot(times, ca_moment0, label="CA 0", **ca_style0)
    plt.plot(times, ca_moment1, label="CA_1", **ca_style1)
    plt.legend(**legend_style)
    plt.grid()

    if WRITE_PLOT_FILES:
        plt.savefig(os.path.join(FOLDER, "growth_mom.eps"))
    plt.show()


if __name__ == "__main__":
    main()
