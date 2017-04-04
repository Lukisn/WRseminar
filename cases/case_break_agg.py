#!/usr/bin/env python3
"""
Demo case for coupled breakage and aggregation.

Taken from Yuan paper case 12.
"""
# standard library imports:
import os
from math import exp, tanh
# third party imports:
import matplotlib.pyplot as plt
from scipy.integrate import quad
# application imports:
from sectional.grid import Grid
from sectional.methods import FixedPivot, CellAverage
from util.plotting import *  # plotting styles


def main():
    """Main function."""

    # PROBLEM FUNCTIONS: ------------------------------------------------------

    PHI = 2

    def phi(t):
        comp = tanh(PHI * t / 2)
        return PHI * (1 + PHI * comp) / (PHI + comp)

    def f(x):  # initial NDF:
        return exp(-x)

    def Gamma(v):  # Breakage frequency
        return PHI**2 * v / 2

    def beta(v1, v2):  # child distribution
        return 2 / v2

    def Q(v1, v2):  # aggregation frequency
        return 1

    def n(t, x):  # analytic solution
        phit = phi(t)
        return phit**2 * exp(-phit * x)

    # CONSTANTS: --------------------------------------------------------------

    # NDF Grid:
    START, END = 0, 1e2
    SECTIONS = 100
    FACTOR = 1.1

    # Simulation:
    T0, TEND = 0, 1
    STEPS = 10  # = (TEND - T0) / 0.1  # 0.1 = TIME_STEP
    TIME_STEP = (TEND - T0) / STEPS
    EVERY = 1
    ORDER = 1

    # Plotting:
    XSCALE, YSCALE = "log", "linear"  # "log" or "linear"
    XMIN, XMAX = 1e-5, 1e3
    YMIN, YMAX = 1e-5, 5

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
        initial_ndf.to_file(os.path.join(FOLDER, "break_agg_initial_ndf.dat"))

    # Fixed Pivot Method:
    fp = FixedPivot(
        initial=initial_ndf,
        bre=True, bre_freq=Gamma, child=beta,
        agg=True, agg_freq=Q,
        gro=False, nuc=False
    )
    fp.simulate(start_time=T0, end_time=TEND, steps=STEPS, write_every=EVERY)
    if WRITE_DATA_FILES:
        fp.moments_to_file(
            os.path.join(FOLDER, "break_agg_fp_moments.dat"), max_order=ORDER
        )
        fp.ndf_to_files(os.path.join(FOLDER, "break_agg_fp_ndf.dat"))

    # Cell Average Technique:
    ca = CellAverage(
        initial=initial_ndf,
        bre=True, bre_freq=Gamma, child=beta,
        agg=True, agg_freq=Q,
        gro=False, nuc=False
    )
    ca.simulate(start_time=T0, end_time=TEND, steps=STEPS, write_every=EVERY)
    if WRITE_DATA_FILES:
        ca.moments_to_file(
            os.path.join(FOLDER, "break_agg_ca_moments.dat"), max_order=ORDER
        )
        ca.ndf_to_files(os.path.join(FOLDER, "break_agg_ca_ndf.dat"))

    # CALCULATIONS FOR PLOTTING: ----------------------------------------------

    # gather NDF data:
    ini_x, ini_y = initial_ndf.pivots(), initial_ndf.densities()
    ana_x, ana_y = initial_ndf.pivots(), []
    for x in ana_x:
        ana_y.append(n(TEND, x))
    fp_x, fp_y = fp._result_ndfs[TEND].pivots(), fp._result_ndfs[TEND].densities()
    ca_x, ca_y = ca._result_ndfs[TEND].pivots(), ca._result_ndfs[TEND].densities()

    # calculate errors:
    fp_err_y, ca_err_y = [], []
    for i, x in enumerate(fp_x):
        ana = n(TEND, x)
        try:
            err = (fp_y[i] - ana) / ana
        except ZeroDivisionError:
            err = None  # matplotlib handles this by not plotting anything
        fp_err_y.append(err)
    for i, x in enumerate(ca_x):
        ana = n(TEND, x)
        try:
            err = (ca_y[i] - ana) / ana
        except ZeroDivisionError:
            err = None  # matplotlib handles this by not plotting anything
        ca_err_y.append(err)

    # gather moment data:
    def moment(order, time, end=END):
        """Calculate analytical moment of `order` at `time` numerically.

        :param order: order of analytic moment.
        :param time: time of evaluation of the analytic solution.
        :param end: upper integration boundary.
        :return: analytical moment of `order` at `time`.
        """
        def integrand(v):
            return v ** order * n(time, v)
        if time == 0:
            mom = 1
        else:
            mom, *_ = quad(integrand, 0, end)
        return mom

    times = sorted(fp._result_ndfs.keys())
    ana_moment0, ana_moment1 = [], []
    fp_moment0, fp_moment1 = [], []
    ca_moment0, ca_moment1 = [], []
    for time in times:
        ana_moment0.append(moment(0, time))
        ana_moment1.append(moment(1, time))
        fp_moment0.append(fp._result_ndfs[time].moment(0))
        fp_moment1.append(fp._result_ndfs[time].moment(1))
        ca_moment0.append(ca._result_ndfs[time].moment(0))
        ca_moment1.append(ca._result_ndfs[time].moment(1))

    # calculate moment errors:
    fp_mom_err_y0, ca_mom_err_y0 = [], []
    fp_mom_err_y1, ca_mom_err_y1 = [], []
    for i, t in enumerate(times):
        try:
            err0 = (fp_moment0[i] - ana_moment0[i]) / ana_moment0[i]
        except ZeroDivisionError:
            err0 = None  # matplotlib handles this by not plotting anything
        try:
            err1 = (fp_moment1[i] - ana_moment1[i]) / ana_moment1[i]
        except ZeroDivisionError:
            err1 = None  # matplotlib handles this by not plotting anything
        fp_mom_err_y0.append(err0)
        fp_mom_err_y1.append(err1)

    for i, t in enumerate(times):
        try:
            err0 = (ca_moment0[i] - ana_moment0[i]) / ana_moment0[i]
        except ZeroDivisionError:
            err0 = None  # matplotlib handles this by not plotting anything
        try:
            err1 = (ca_moment1[i] - ana_moment1[i]) / ana_moment1[i]
        except ZeroDivisionError:
            err1 = None  # matplotlib handles this by not plotting anything
        ca_mom_err_y0.append(err0)
        ca_mom_err_y1.append(err1)

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
    # lower subplot - errors:
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
        plt.savefig(os.path.join(FOLDER, "break_agg_ndf.eps"))
    plt.show()

    # plot moment comparison and errors:
    fig, (upper, lower) = plt.subplots(2, 1, sharex="all")
    # upper subplot - moments:
    upper.set_title("moments over time. STEP = {} s".format(TIME_STEP))
    upper.set_ylabel("moment")
    upper.plot(times, ana_moment0, label="ana 0", **ana_style0)
    upper.plot(times, ana_moment1, label="ana 1", **ana_style1)
    upper.plot(times, fp_moment0, label="FP 0", **fp_style0)
    upper.plot(times, fp_moment1, label="FP 1", **fp_style1)
    upper.plot(times, ca_moment0, label="CA 0", **ca_style0)
    upper.plot(times, ca_moment1, label="CA 1", **ca_style1)
    upper.legend(**legend_style)
    upper.grid()
    # lower subplot - errors:
    lower.set_xlabel("time in s")
    lower.set_ylabel("relative error")
    lower.plot(times, fp_mom_err_y0, label="FP 0", **fp_style0)
    lower.plot(times, fp_mom_err_y1, label="FP 1", **fp_style1)
    lower.plot(times, ca_mom_err_y0, label="CA 0", **ca_style0)
    lower.plot(times, ca_mom_err_y1, label="CA 1", **ca_style1)
    lower.legend(**legend_style)
    lower.grid()
    # tighten layout and show:
    fig.tight_layout()
    if WRITE_PLOT_FILES:
        fig.savefig(os.path.join(FOLDER, "break_agg_mom.eps"))
    plt.show()


if __name__ == "__main__":
    main()
