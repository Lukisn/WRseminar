#!/usr/bin/env python3
"""Utility functions and style definitions for plotting the results."""
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad, quadrature
from scipy import inf


MARKER_SIZE = 5
ANALYTIC_MARKER_SIZE = MARKER_SIZE + 2

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
    "marker": None, "markersize": ANALYTIC_MARKER_SIZE
}
ana_style0 = {
    "linewidth": THICK_LINE, "linestyle": "solid", "color": ANALYTIC_COLOR,
    "marker": "o", "markersize": ANALYTIC_MARKER_SIZE,
    "markerfacecolor": "None"
}
ana_style1 = {
    "linewidth": THICK_LINE, "linestyle": "solid", "color": ANALYTIC_COLOR,
    "marker": "s", "markersize": ANALYTIC_MARKER_SIZE,
    "markerfacecolor": "None"
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
    "marker": "o", "markersize": MARKER_SIZE, "markerfacecolor": "None"
}
ca_style1 = {
    "linewidth": SMALL_LINE, "linestyle": "dotted", "color": CELL_AVERAGE_COLOR,
    "marker": "s", "markersize": MARKER_SIZE, "markerfacecolor": "None"
}

fp_style = {
    "linewidth": NO_LINE, "linestyle": "dashed", "color": FIXED_PIVOT_COLOR,
    "marker": "x", "markersize": MARKER_SIZE, "markerfacecolor": "None"
}
fp_style0 = {
    "linewidth": SMALL_LINE, "linestyle": "dotted",  "color": FIXED_PIVOT_COLOR,
    "marker": "h", "markersize": MARKER_SIZE, "markerfacecolor": "None"
}
fp_style1 = {
    "linewidth": SMALL_LINE, "linestyle": "dotted", "color": FIXED_PIVOT_COLOR,
    "marker": "D", "markersize": MARKER_SIZE, "markerfacecolor": "None"
}


def moment_inf(n, order, time, zero_default=None):
    """Calculate analytical moment of `order` at `time` numerically.

    :param order: order of analytic moment.
    :param time: time of evaluation of the analytic solution.
    :return: analytical moment of `order` at `time`.
    """
    def integrand(v):
        return v ** order * n(time, v)
    if time == 0 and zero_default is not None:
        return zero_default
    mom, *_ = quad(integrand, 0, inf)  # , maxp1=100)
    return mom


def moment_end(n, order, time, end, zero_default=None):
    """Calculate analytical moment of `order` at `time` numerically.

    :param order: order of analytic moment.
    :param time: time of evaluation of the analytic solution.
    :param end: upper integration boundary.
    :return: analytical moment of `order` at `time`.
    """
    def integrand(v):
        return v ** order * n(time, v)
    if time == 0 and zero_default is not None:
        return zero_default
    mom, *_ = quad(integrand, 0, end)  # , maxp1=100
    return mom


def plot_results(initial, solution, fp, ca,
                 ndf_end, t_end, time_step,
                 xscale, yscale,
                 xlim_ndf, ylim_ndf, ylim_ndf_err,
                 ylim_mom, ylim_mom_err,
                 write_plot_files, output_folder, prefix,
                 mom_type="inf", zero_default=None, show_plots=True):
    """Plot resulting NDFs and Moments and compare with analytic solution.
    
    :param initial: initial NDF
    :param solution: analytic solution function
    :param fp: parameterized fixed pivot method object
    :param ca: parameterized cell average method object
    :param ndf_end: upper end value of the discrete NDF
    :param t_end: end time of the simulation
    :param time_step: time step used in the simulation
    :param xscale: x-axis scale to use ("linear" / "log")
    :param yscale: y-axis scale to use ("linear" / "log")
    :param xlim_ndf: x-axis limits for NDF plots (x_min, x_max)
    :param ylim_ndf: y-axis limits for NDF plot (y_min, y_max)
    :param ylim_ndf_err: y-axis limits for NDF error plot
    :param ylim_mom: y-axis limits for moment plot
    :param ylim_mom_err: y-axis limits for moment error plot
    :param write_plot_files: boolean flag for toggling output of plot files
    :param output_folder: directory to write output files to
    :param mom_type: type of moment calculation method to use ("inf" / "end")
    :param prefix: file name prefix for saving files
    :param show_plots: boolean flag for toggling plt.show() calls on and off
    :return: None
    """
    # CALCULATIONS FOR PLOTTING: ----------------------------------------------

    # gather NDF data:
    ini_x, ini_y = initial.pivots(), initial.densities()
    ana_x, ana_y = initial.pivots(), []
    for x in ana_x:
        ana_y.append(solution(t_end, x))
    fp_x, fp_y = fp._result_ndfs[t_end].pivots(), fp._result_ndfs[
        t_end].densities()
    ca_x, ca_y = ca._result_ndfs[t_end].pivots(), ca._result_ndfs[
        t_end].densities()

    # calculate errors:
    fp_err_y, ca_err_y = [], []
    for i, x in enumerate(fp_x):
        ana = solution(t_end, x)
        try:
            err = (fp_y[i] - ana) / ana
        except ZeroDivisionError:
            err = None  # matplotlib handles this by not plotting anything
        fp_err_y.append(err)
    for i, x in enumerate(ca_x):
        ana = solution(t_end, x)
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
            ana_moment0.append(moment_inf(solution, 0, time, zero_default))
            ana_moment1.append(moment_inf(solution, 1, time, zero_default))
        elif mom_type == "end":
            ana_moment0.append(moment_end(solution, 0, time, ndf_end, zero_default))
            ana_moment1.append(moment_end(solution, 1, time, ndf_end, zero_default))
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

    # configure matplotlib to use latex rendering:
    plt.rc("text", usetex=True)

    # plot NDF comparison and errors:
    fig, (upper, lower) = plt.subplots(2, 1, sharex="all")
    # upper subplot - NDF:
    upper.set_title(
        "NDFs bei $t_{{end}} = {:.2f}$ s, $\Delta t = {:.2f}$ s".format(
            t_end, time_step
        )
    )
    upper.set_ylabel("NDF $n(t_{end}, v)$")
    upper.plot(ana_x, ana_y, label="ana", **ana_style)
    upper.plot(ini_x, ini_y, label="ini", **initial_style)
    upper.plot(fp_x, fp_y, label="FP", **fp_style)
    upper.plot(ca_x, ca_y, label="CA", **ca_style)
    # upper.set_xlim(XMIN, XMAX)  # not needed for shared x-axis!?
    upper.set_ylim(ylim_ndf)
    upper.set_xscale(xscale)
    upper.set_yscale(yscale)
    upper.legend(**legend_style)
    upper.grid()
    # lower subplot - errors:
    lower.set_xlabel("Interne Koordinate $v$")
    lower.set_ylabel("Relativer Fehler")
    lower.plot(fp_x, fp_err_y, label="FP", **fp_style)
    lower.plot(ca_x, ca_err_y, label="CA", **ca_style)
    lower.set_xlim(xlim_ndf)
    lower.set_ylim(ylim_ndf_err)
    lower.set_xscale(xscale)
    lower.legend(**legend_style)
    lower.grid()

    # tighten layout and show:
    fig.tight_layout()
    if write_plot_files:
        plt.savefig(os.path.join(output_folder, "{}_ndf.eps".format(prefix)))
    if show_plots:
        plt.show()

    # plot moment comparison and errors:
    fig, (upper, lower) = plt.subplots(2, 1, sharex="all")
    # upper subplot - moments:
    upper.set_title("Zeitverlauf der Momente, $\Delta t = {}$ s".format(
        time_step
    ))
    # plot 0th moment on primary y axis:
    upper.set_ylabel("0. Moment $\mu_0$ (Anzahl)")
    upper.plot(times, ana_moment0, label="ana 0", **ana_style0)
    upper.plot(times, fp_moment0, label="FP 0", **fp_style0)
    upper.plot(times, ca_moment0, label="CA 0", **ca_style0)
    upper.set_ylim(ylim_mom)
    # plot 1st moment on second y axis:
    upper1 = upper.twinx()
    upper1.set_ylabel("1. Moment $\mu_1$ (Masse)")
    upper1.plot(times, ana_moment1, label="ana 1", **ana_style1)
    upper1.plot(times, fp_moment1, label="FP 1", **fp_style1)
    upper1.plot(times, ca_moment1, label="CA 1", **ca_style1)
    upper1.set_ylim(ylim_mom)
    # plot common legend:
    lines0, labels0 = upper.get_legend_handles_labels()
    lines1, labels1 = upper1.get_legend_handles_labels()
    lines, labels = lines0 + lines1, labels0 + labels1
    upper1.legend(lines, labels, **legend_style)
    upper.grid()
    # lower subplot - errors:
    lower.set_xlabel("Zeit $t$ in s")
    lower.set_ylabel("Relativer Fehler")
    lower.plot(times, fp_mom_err_y0, label="FP 0", **fp_style0)
    lower.plot(times, fp_mom_err_y1, label="FP 1", **fp_style1)
    lower.plot(times, ca_mom_err_y0, label="CA 0", **ca_style0)
    lower.plot(times, ca_mom_err_y1, label="CA 1", **ca_style1)
    lower.set_ylim(ylim_mom_err)
    lower.legend(**legend_style)
    lower.grid()

    # tighten layout and show:
    fig.tight_layout()
    if write_plot_files:
        fig.savefig(os.path.join(output_folder, "{}_mom.eps".format(prefix)))
    if show_plots:
        plt.show()
