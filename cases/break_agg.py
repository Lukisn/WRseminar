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


def main(show_plots=True):
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
    START, END = 0, 1e2  # = 100 ;)
    SECTIONS = 100  # A/B: 100, C1/C2: 50, D1/D2: 25
    FACTOR = 1.1  # A/B: 1.1, C1/C2: 1.2, D1/D2: 1.4

    # Simulation:
    T0, TEND = 0, 1
    STEPS = 10  # A/C1/D1: 100, B/C2/D2: 10
    TIME_STEP = (TEND - T0) / STEPS
    EVERY = 1
    ORDER = 1

    # Plotting:
    XSCALE, YSCALE = "log", "linear"  # "log" or "linear"
    XLIM_NDF = 1e-5, 1e3
    YLIM_NDF = 1e-5, 5
    YLIM_NDF_ERR = -1.1, 1.1
    YLIM_MOM = 0.7, 1.9
    YLIM_MOM_ERR = -0.55, 0.15

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

    # PLOTTING: ---------------------------------------------------------------

    plot_results(initial=initial_ndf, solution=n, fp=fp, ca=ca,
                 ndf_end=END, t_end=TEND, time_step=TIME_STEP,
                 xscale=XSCALE, yscale=YSCALE,
                 xlim_ndf=XLIM_NDF,
                 ylim_ndf=YLIM_NDF, ylim_ndf_err=YLIM_NDF_ERR,
                 ylim_mom=YLIM_MOM, ylim_mom_err=YLIM_MOM_ERR,
                 write_plot_files=WRITE_PLOT_FILES, output_folder=FOLDER,
                 prefix="break_agg", mom_type="inf", zero_default=None,
                 show_plots=show_plots)


if __name__ == "__main__":
    main()
