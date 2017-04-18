#!/usr/bin/env python3
"""Utility function and style definitions for plotting the results."""
import os
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy import inf

MARKER_SIZE = 5

NO_LINE = 0
SMALL_LINE = 0.75
MEDIUM_LINE = 1
THICK_LINE = 2

ANALYTIC_COLOR = "orange"
INITIAL_COLOR = "green"
FIXED_PIVOT_COLOR = "blue"
CELL_AVERAGE_COLOR = "red"

legend_style = {
    "loc": "best", "fontsize": "small", "numpoints": 1,
    "framealpha": 0.75, "fancybox": True
}

ana_style = {
    "linewidth": THICK_LINE, "linestyle": "solid", "color": ANALYTIC_COLOR,
    "marker": None, "markersize": MARKER_SIZE
}
ana_style0 = {
    "linewidth": THICK_LINE, "linestyle": "solid", "color": ANALYTIC_COLOR,
    "marker": "o", "markersize": MARKER_SIZE, "markerfacecolor": "None"
}
ana_style1 = {
    "linewidth": THICK_LINE, "linestyle": "solid", "color": ANALYTIC_COLOR,
    "marker": "s", "markersize": MARKER_SIZE, "markerfacecolor": "None"
}

initial_style = {
    "linewidth": MEDIUM_LINE, "linestyle": "dotted", "color": INITIAL_COLOR,
    "marker": ".", "markersize": MARKER_SIZE, "markerfacecolor": "None"
}

ca_style = {
    "linewidth": NO_LINE, "linestyle": "dashed", "color": CELL_AVERAGE_COLOR,
    "marker": "o", "markersize": MARKER_SIZE, "markerfacecolor": "None"
}
ca_style0 = {
    "linewidth": SMALL_LINE, "linestyle": "dotted", "color": CELL_AVERAGE_COLOR,
    "marker": "^", "markersize": MARKER_SIZE, "markerfacecolor": "None"
}
ca_style1 = {
    "linewidth": SMALL_LINE, "linestyle": "dotted", "color": CELL_AVERAGE_COLOR,
    "marker": "v", "markersize": MARKER_SIZE, "markerfacecolor": "None"
}

fp_style = {
    "linewidth": NO_LINE, "linestyle": "dashed", "color": FIXED_PIVOT_COLOR,
    "marker": "x", "markersize": MARKER_SIZE, "markerfacecolor": "None"
}
fp_style0 = {
    "linewidth": SMALL_LINE, "linestyle": "dotted",  "color": FIXED_PIVOT_COLOR,
    "marker": "v", "markersize": MARKER_SIZE, "markerfacecolor": "None"
}
fp_style1 = {
    "linewidth": SMALL_LINE, "linestyle": "dotted", "color": FIXED_PIVOT_COLOR,
    "marker": "^", "markersize": MARKER_SIZE, "markerfacecolor": "None"
}


def moment_inf(n, order, time):
    """Calculate analytical moment of `order` at `time` numerically.

    :param order: order of analytic moment.
    :param time: time of evaluation of the analytic solution.
    :return: analytical moment of `order` at `time`.
    """

    def integrand(v):
        return v ** order * n(time, v)

    if time == 0:
        mom = 1
    else:
        mom, *_ = quad(integrand, 0, inf)
    return mom


def moment_end(n, order, time, end):
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
        mom, *_ = quad(integrand, 0, end, maxp1=100)
    return mom


def plot_results(initial_ndf, n, fp, ca,
                 END, TEND, TIME_STEP, XMIN, XMAX, YMIN, YMAX, XSCALE, YSCALE,
                 YMIN_ERR, YMAX_ERR, YMIN_MOM_ERR, YMAX_MOM_ERR,
                 WRITE_PLOT_FILES, FOLDER, mom_type, prefix):
    """Plot resulting NDFs and Moments and compare with analytic solution."""
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

    # TODO: find reason for end/inf diversity
    # gather moment data:
    times = sorted(fp._result_ndfs.keys())
    ana_moment0, ana_moment1 = [], []
    fp_moment0, fp_moment1 = [], []
    ca_moment0, ca_moment1 = [], []
    for time in times:
        if mom_type == "inf":
            ana_moment0.append(moment_inf(n, 0, time))
            ana_moment1.append(moment_inf(n, 1, time))
        elif mom_type == "end":
            ana_moment0.append(moment_end(n, 0, time, END))
            ana_moment1.append(moment_end(n, 1, time, END))
        else:
            msg = "Unknown mom_type '{}'. Use mom_type 'inf' or 'end'!"
            raise ValueError(msg.format(mom_type))
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

    # don't plot moment and errors for t = 0:
    # ana_moment0[0], ana_moment1[0] = None, None
    # fp_moment0[0], fp_moment1[0] = None, None
    # ca_moment0[0], ca_moment1[0] = None, None
    # fp_mom_err_y0[0], fp_mom_err_y1[0] = None, None
    # ca_mom_err_y0[0], ca_mom_err_y1[0] = None, None

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
    lower.set_ylim(YMIN_ERR, YMAX_ERR)
    lower.set_xscale(XSCALE)
    lower.legend(**legend_style)
    lower.grid()

    # tighten layout and show:
    fig.tight_layout()
    if WRITE_PLOT_FILES:
        plt.savefig(os.path.join(FOLDER, "{}_ndf.eps".format(prefix)))
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
    lower.set_ylim(YMIN_MOM_ERR, YMAX_MOM_ERR)
    lower.legend(**legend_style)
    lower.grid()
    # tighten layout and show:
    fig.tight_layout()
    if WRITE_PLOT_FILES:
        fig.savefig(os.path.join(FOLDER, "{}_mom.eps".format(prefix)))
    plt.show()
