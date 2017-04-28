#!/usr/bin/env python3
"""
Demo case for pure aggregation (coagulation).

Taken from Yuan paper case 9.
"""
# standard library imports:
import os
from math import exp, sqrt
# third party imports:
import matplotlib.pyplot as plt
from scipy import inf
from scipy.integrate import quad
from scipy.special import iv  # modified bessel function
# application imports:
from sectional.grid import Grid
from sectional.methods import FixedPivot, CellAverage
from util.plotting import *  # plotting styles


def main(show_plots=True):
    """Main function."""

    # PROBLEM FUNCTIONS: ------------------------------------------------------

    def f(x):  # initial NDF
        return exp(-x)

    def Q(v1, v2):   # aggregation function
        return v1 + v2

    def n(t, x):   # analytic solution
        if t == 0 or x == 0:
            return 0
        else:
            ratio = exp(-t - 2 * x + x * exp(-t)) / (x * sqrt(1 - exp(-t)))
            if ratio == 0:
                return 0
            else:
                bessel = iv(1, 2 * x * sqrt(1 - exp(-t)))
                result = ratio * bessel
                return result

    # CONSTANTS: --------------------------------------------------------------

    # NDF Grid:
    START, END = 0, 1e5
    SECTIONS = 100  # A: 100, B: 50, C: 25, D: 100
    FACTOR = 1.2  # A: 1.2, B: 1.4, C: 1.8, D: 1.2

    # Simulation:
    T0, TEND = 0, 1
    STEPS = 10  # A-C: 10, D: 100
    TIME_STEP = (TEND - T0) / STEPS
    EVERY = 1
    ORDER = 1

    # Plotting:
    XSCALE, YSCALE = "log", "log"  # or "linear"
    XLIM_NDF = 1e-5, 1e5
    YLIM_NDF = 1e-10, 1e5
    YLIM_NDF_ERR = -1.1, 1.1
    YLIM_MOM = 0.3, 1.1
    YLIM_MOM_ERR = -0.125, 0.075

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
        initial_ndf.to_file(os.path.join(FOLDER, "agg_initial_ndf.dat"))

    # Fixed Pivot Method:
    fp = FixedPivot(
        initial=initial_ndf,
        agg=True, agg_freq=Q,
        bre=False, gro=False, nuc=False
    )
    fp.simulate(start_time=T0, end_time=TEND, steps=STEPS, write_every=EVERY)
    if WRITE_DATA_FILES:
        fp.moments_to_file(os.path.join(
            FOLDER, "agg_fp_moments.dat"), max_order=ORDER
        )
        fp.ndf_to_files(os.path.join(FOLDER, "agg_fp_ndf.dat"))

    # Cell Average Technique:
    ca = CellAverage(
        initial=initial_ndf,
        agg=True, agg_freq=Q,
        bre=False, gro=False, nuc=False
    )
    ca.simulate(start_time=T0, end_time=TEND, steps=STEPS, write_every=EVERY)
    if WRITE_DATA_FILES:
        ca.moments_to_file(
            os.path.join(FOLDER, "agg_ca_moments.dat"), max_order=ORDER
        )
        ca.ndf_to_files(os.path.join(FOLDER, "agg_ca_ndf.dat"))

    # PLOTTING: ---------------------------------------------------------------

    plot_results(initial=initial_ndf, solution=n, fp=fp, ca=ca,
                 ndf_end=END, t_end=TEND, time_step=TIME_STEP,
                 xscale=XSCALE, yscale=YSCALE,
                 xlim_ndf=XLIM_NDF,
                 ylim_ndf=YLIM_NDF, ylim_ndf_err=YLIM_NDF_ERR,
                 ylim_mom=YLIM_MOM, ylim_mom_err=YLIM_MOM_ERR,
                 write_plot_files=WRITE_PLOT_FILES, output_folder=FOLDER,
                 prefix="agg", mom_type="inf", zero_default=1,
                 show_plots=show_plots)


if __name__ == "__main__":
    main()
