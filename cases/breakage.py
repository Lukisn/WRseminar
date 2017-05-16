#!/usr/bin/env python3
"""
Demo case for pure breakage.

Taken from Yuan paper case 8.
"""
# standard library imports:
import os
from math import exp
# third party imports:
import matplotlib.pyplot as plt
from scipy.integrate import quad
# application imports:
from sectional.functions import hstep, dirac_norm
from sectional.grid import Grid
from sectional.methods import FixedPivot, CellAverage
from util.plotting import *  # plotting styles


def main(show_plots=True):
    """Main function."""

    # PROBLEM FUNCTIONS: ------------------------------------------------------

    def f(x):   # initial NDF
        return dirac_norm(x - 1)

    def Gamma(v):  # breakage frequency
        return v ** 2

    def beta(v1, v2):  # child distribution
        return 2 / v2 * hstep(v2 - v1)

    def n(t, x):  # analytic solution:
        brace = dirac_norm(x - 1) + 2 * t * hstep(1 - x)
        return exp(-t * x ** 2) * brace

    # CONSTANTS: --------------------------------------------------------------

    # NDF Grid:
    START, END = 0, 10
    SECTIONS = 100  # A: 100, B: 50, C: 25, D: 100
    FACTOR = 1.1  # A-D: 1.1

    # Simulation:
    T0, TEND = 0, 10
    STEPS = 10  # A-C: 10, D: 100
    TIME_STEP = (TEND - T0) / STEPS
    EVERY = 1
    ORDER = 1

    # Plotting:
    XSCALE, YSCALE = "log", "log"  # or "linear"
    XLIM_NDF = 1e-5, 1e1  # (x_min, x_max)
    YLIM_NDF = 1e-5, 1e3
    YLIM_NDF_ERR = -1.1, 0.25
    YLIM_MOM = 0.5, 6.5
    YLIM_MOM_ERR = -0.65, 0.25

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
        initial_ndf.to_file(os.path.join(FOLDER, "break_initial_ndf.dat"))

    # Fixed Pivot Method:
    fp = FixedPivot(
        initial=initial_ndf,
        bre=True, bre_freq=Gamma, child=beta,
        agg=False, gro=False, nuc=False
    )
    fp.simulate(start_time=T0, end_time=TEND, steps=STEPS, write_every=EVERY)
    if WRITE_DATA_FILES:
        fp.moments_to_file(
            os.path.join(FOLDER, "break_fp_moments.dat"), max_order=ORDER
        )
        fp.ndf_to_files(os.path.join(FOLDER, "break_fp_ndf.dat"))

    # Cell Average Technique:
    ca = CellAverage(
        initial=initial_ndf,
        bre=True, bre_freq=Gamma, child=beta,
        agg=False, gro=False, nuc=False
    )
    ca.simulate(start_time=T0, end_time=TEND, steps=STEPS, write_every=EVERY)
    if WRITE_DATA_FILES:
        ca.moments_to_file(
            os.path.join(FOLDER, "break_ca_moments.dat"), max_order=ORDER
        )
        ca.ndf_to_files(os.path.join(FOLDER, "break_ca_ndf.dat"))

    # PLOTTING: ---------------------------------------------------------------
    plot_results(initial=initial_ndf, solution=n, fp=fp, ca=ca,
                 ndf_end=END, t_end=TEND, time_step=TIME_STEP,
                 xscale=XSCALE, yscale=YSCALE,
                 xlim_ndf=XLIM_NDF,
                 ylim_ndf=YLIM_NDF, ylim_ndf_err=YLIM_NDF_ERR,
                 ylim_mom=YLIM_MOM, ylim_mom_err=YLIM_MOM_ERR,
                 write_plot_files=WRITE_PLOT_FILES, output_folder=FOLDER,
                 prefix="break", mom_type="end", zero_default=1,
                 show_plots=show_plots)


if __name__ == "__main__":
    main()
