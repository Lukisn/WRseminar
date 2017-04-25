#!/usr/bin/env python3
"""
Demo case for pure growth

Taken from Yuan paper Case 1.
"""
# standard library imports:
import os
# third party imports:
import matplotlib.pyplot as plt
from scipy.integrate import quad
# application imports:
from sectional.grid import Grid
from sectional.methods import FixedPivot, CellAverage
from util.plotting import *  # plotting styles


def main(show_plots=True):
    """Main function."""

    # PROBLEM FUNCTIONS: ------------------------------------------------------

    def f(x):  # initial NDF
        return max(60 * x**2 * (1 - x)**3, 0)

    def G(v):  # growth function
        #return -1 / 2
        return 1 / 2

    def n(t, x):  # analytic solution
        # return max(60 * (x - t / 2) ** 2 * (1 - (x - t / 2)) ** 3, 0)
        if x > t/2:
            return max(60 * (x - t / 2)**2 * (1 - (x - t / 2))**3, 0)
        else:
            return 0

    # CONSTANTS: --------------------------------------------------------------

    # Grid:
    START, END = 0, 2
    SECTIONS = 100
    FACTOR = 1.1

    # Simulation:
    T0, TEND = 0, 1
    STEPS = 100
    TIME_STEP = (TEND - T0) / STEPS
    EVERY = 1
    ORDER = 1

    # Plotting:
    XSCALE, YSCALE = "linear", "linear"
    XLIM_NDF= 0, 2
    YLIM_NDF = 0, 2.5
    YLIM_NDF_ERR = -1.1, 1.1
    YLIM_MOM = -0.1, 2.1
    YLIM_MOM_ERR = -0.005, 0.005

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
        initial_ndf.to_file(os.path.join(FOLDER, "growth1_initial_ndf.dat"))

    # Fixed Pivot Method:
    fp = FixedPivot(
        initial=initial_ndf,
        gro=True, gro_rate=G,
        bre=False, agg=False, nuc=False
    )
    fp.simulate(start_time=T0, end_time=TEND, steps=STEPS, write_every=EVERY)
    if WRITE_DATA_FILES:
        fp.moments_to_file(
            os.path.join(FOLDER, "growth1_fp_moments.dat"), max_order=ORDER
        )
        fp.ndf_to_files(os.path.join(FOLDER, "growth1_fp_ndf.dat"))

    # Cell Average Technique:
    ca = CellAverage(
        initial=initial_ndf,
        gro=True, gro_rate=G,
        bre=False, agg=False, nuc=False
    )
    ca.simulate(start_time=T0, end_time=TEND, steps=STEPS, write_every=EVERY)
    if WRITE_DATA_FILES:
        ca.moments_to_file(
            os.path.join(FOLDER, "growth1_ca_moments.dat"), max_order=ORDER
        )
        ca.ndf_to_files(os.path.join(FOLDER, "growth1_ca_ndf.dat"))

    # PLOTTING: ---------------------------------------------------------------

    plot_results(initial=initial_ndf, solution=n, fp=fp, ca=ca,
                 ndf_end=END, t_end=TEND, time_step=TIME_STEP,
                 xscale=XSCALE, yscale=YSCALE,
                 xlim_ndf=XLIM_NDF,
                 ylim_ndf=YLIM_NDF, ylim_ndf_err=YLIM_NDF_ERR,
                 ylim_mom=YLIM_MOM, ylim_mom_err=YLIM_MOM_ERR,
                 write_plot_files=WRITE_PLOT_FILES, output_folder=FOLDER,
                 mom_type="inf", prefix="growth1", show_plots=show_plots)


if __name__ == "__main__":
    main()
