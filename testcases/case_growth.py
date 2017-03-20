#!/usr/bin/env python3

"""
Demo case for pure growth
"""
import matplotlib.pyplot as plt
from sectional.functions import hstep
from sectional.grid import Grid
from sectional.methods import FixedPivot, CellAverage


def main():
    """Main function.
    """

    # PROBLEM FUNCTIONS: ------------------------------------------------------

    # initial NDF:
    def f(x):
        return hstep(x - 0.1) - hstep(x - 0.6)

    # Growth function:
    def G(v):
        return 1

    # analytic solution:
    def n(t, x):
        return hstep(x - (0.1 + t)) - hstep(x - (0.6 + t))

    # CONSTANTS: --------------------------------------------------------------

    # Grid:
    START, END = 0, 2
    SECTIONS = 100
    FACTOR = 1.0

    # Simulation:
    T0, TEND = 0, 1
    STEPS = 100
    EVERY = 1
    ORDER = 1

    # Plotting:
    XSCALE, YSCALE = "linear", "linear"
    XMIN, XMAX = 0, 2
    YMIN, YMAX = 0, 1.5

    # SIMULATION: -------------------------------------------------------------

    # initial NDF:
    initial_ndf = Grid.create_geometric_end(
        start=START, end=END, sec=SECTIONS, fact=FACTOR, func=f
    )
    initial_ndf.to_file("results/growth_initial_ndf.dat")

    # Fixed Pivot Method:
    fp = FixedPivot(
        initial=initial_ndf,
        gro=True, gro_rate=G,
        bre=False, agg=False, nuc=False
    )
    fp.simulate(start_time=T0, end_time=TEND, steps=STEPS, write_every=EVERY)
    fp.moments_to_file("results/growth_fp_moments.dat", max_order=ORDER)
    fp.ndf_to_files("results/growth_fp_ndf.dat")

    # Cell Average Technique:
    ca = CellAverage(
        initial=initial_ndf,
        gro=True, gro_rate=G,
        bre=False, agg=False, nuc=False
    )
    ca.simulate(start_time=T0, end_time=TEND, steps=STEPS, write_every=EVERY)
    ca.moments_to_file("results/growth_ca_moments.dat", max_order=ORDER)
    ca.ndf_to_files("results/growth_ca_ndf.dat")

    # PLOTTING: ---------------------------------------------------------------

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

    # plot NDF comparison and errors:
    # upper subplot: NDF:
    plt.subplot(211)
    # plt.xlabel("size")
    plt.ylabel("NDF")
    plt.plot(ana_x, ana_y, "y-", lw=3, label="analytic")
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
    plt.plot(fp_x, fp_err_y, "bx-", label="fixed pivot")
    plt.plot(ca_x, ca_err_y, "r.-", lw=2, label="cell average")
    plt.xlim(XMIN, XMAX)
    # plt.ylim(YMIN, YMAX)  # comment out = auto
    plt.xscale(XSCALE)
    # plt.yscale(YSCALE)  # comment out  = auto
    plt.legend(loc="best", fontsize="small")
    plt.grid()

    plt.savefig("results/growth_ndf.eps")
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

    plt.savefig("results/growth_mom.eps")
    plt.show()


if __name__ == "__main__":
    main()
