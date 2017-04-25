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
    SECTIONS = 100  # A: 100, B: 100, C: 50, D: 25
    FACTOR = 1.2  # A: 1.2, B: 1.2, C: 1.4, D: 1.8

    # Simulation:
    T0, TEND = 0, 1
    STEPS = 10  # A: 100, B-D: 10
    TIME_STEP = (TEND - T0) / STEPS
    EVERY = 1
    ORDER = 1

    # Plotting:
    XSCALE, YSCALE = "log", "log"  # or "linear"
    XMIN, XMAX = 1e-5, 1e5
    YMIN, YMAX = 1e-10, 1e5

    YMIN_ERR, YMAX_ERR = -1.1, 1.1
    YMIN_MOM_ERR, YMAX_MOM_ERR = -0.125, 0.075

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
    plot_results(initial_ndf, n, fp, ca,
                 END, TEND, TIME_STEP, XMIN, XMAX, YMIN, YMAX, XSCALE, YSCALE,
                 YMIN_ERR, YMAX_ERR, YMIN_MOM_ERR, YMAX_MOM_ERR,
                 WRITE_PLOT_FILES, FOLDER, mom_type="inf", prefix="agg",
                 show_plots=show_plots)


if __name__ == "__main__":
    main()
