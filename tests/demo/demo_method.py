#!/usr/bin/env python3
# TODO: find errors
# this demo module uses old apis and fails!!!

"""
Module demonstrating the basic usage of the provided classes for the
calculation methods.
"""
from collections import Iterable
from math import exp

import matplotlib.pyplot as plt

from sectional.grid import Grid
from sectional.methods import FixedPivot, CellAverage

# CONSTANTS: ------------------------------------------------------------------

# GRID:
START, END = 0, 1000
SECTIONS = 50
FACTOR = 1.3

# SIMULATION:
END_TIME = 1
STEPS = 100

# PLOTTING:
XSCALE, YSCALE = "log", "log"
XMIN, XMAX = 1e-5, 1e5
YMIN, YMAX = 1e-5, 1e5
EVERY = 1
ORDER = 2


# UTILITY FUNCTIONS: ----------------------------------------------------------

def f(v, N0=1, v0=1):
    """initial number density function.
    """
    return (N0 / v0) * (v / v0) * exp(-v / v0)


def Gamma(v):
    """breakage frequency function.
    """
    return v ** 2


def beta(v1, v2):
    """child number function.
    """
    return 2 / v2


def Q(x1, x2):
    """aggregation frequency function.
    """
    return 1


def G(v):
    """growth rate for particles of size v.
    """
    return 0.1


def S(v):
    """rate of nucleation of particles of size v.
    """
    if v <= 0.0001:
        return 1
    else:
        return 0


def plot_initial_and_current(methods):
    """Plot initial and current NDFs of the given methods.
    """
    if not isinstance(methods, Iterable):
        methods = [methods]

    plt.xlabel("particle volume")
    plt.ylabel("number density")
    for method in methods:
        name = method.__class__.__name__
        ini_pivs, ini_dens = method._initial.pivots(), method._initial.densities()
        curr_pivs, curr_dens = method._current.pivots(), method._current.densities()
        plt.plot(ini_pivs, ini_dens, ".-", label="{} ini".format(name))
        plt.plot(curr_pivs, curr_dens, ".-", label="{} curr".format(name))
    plt.xscale(XSCALE)
    plt.yscale(YSCALE)
    plt.xlim(XMIN, XMAX)
    plt.ylim(YMIN, YMAX)
    plt.legend(loc="best", fontsize="small")
    plt.grid()
    plt.show()


def plot_moments_over_time(methods, max_order):
    """Plot development of the moments over time of the given methods.
    """
    if not isinstance(methods, Iterable):
        methods = [methods]

    plt.xlabel("time")
    plt.ylabel("moment")
    for method in methods:
        name = method.__class__.__name__
        times = sorted(method.result_moments)
        moments = {}
        for order in range(max_order + 1):
            moments[order] = []
            for time in times:
                moments[order].append(method.result_moments[time][order])
        for order in range(max_order + 1):
            plt.plot(times, moments[order], ".-", label="{} mom{}".format(name, order))
    plt.legend(loc="best", fontsize="small")
    plt.grid()
    plt.show()


# DEMO FUNCTIONS: -------------------------------------------------------------

def demo_zero():
    """Demo basic behaviour without specifying the keyword arguments.
    """
    print("Demoing zero behaviour.")
    ini = Grid.create_geometric_end(
        start=START, end=END, sec=SECTIONS, fact=FACTOR, func=f
    )
    fp = FixedPivot(initial=ini)
    ca = CellAverage(initial=ini)
    fp.simulate(
        end_time=END_TIME, steps=STEPS, write_every=EVERY  #, max_order=ORDER
    )
    ca.simulate(
        end_time=END_TIME, steps=STEPS, write_every=EVERY  #, max_order=ORDER
    )
    plot_initial_and_current([fp, ca])
    plot_moments_over_time([fp, ca], max_order=ORDER)


def demo_breakage():
    """Demo pure breakage.
    """
    print("Demoing pure breakage.")
    ini = Grid.create_geometric_end(
        start=START, end=END, sec=SECTIONS, fact=FACTOR, func=f
    )
    fp = FixedPivot(
        initial=ini,
        bre=True, bre_freq=Gamma, child=beta,
        agg=False, gro=False, nuc=False
    )
    ca = CellAverage(
        initial=ini,
        bre=True, bre_freq=Gamma, child=beta,
        agg=False, gro=False, nuc=False
    )
    fp.simulate(
        end_time=END_TIME, steps=STEPS, write_every=EVERY  #, max_order=ORDER
    )
    ca.simulate(
        end_time=END_TIME, steps=STEPS, write_every=EVERY  #, max_order=ORDER
    )
    plot_initial_and_current([fp, ca])
    plot_moments_over_time([fp, ca], max_order=ORDER)


def demo_aggregation():
    """Demo pure aggregation.
    """
    print("Demoing pure aggregation.")
    ini = Grid.create_geometric_end(
        start=START, end=END, sec=SECTIONS, fact=FACTOR, func=f
    )
    fp = FixedPivot(
        initial=ini,
        agg=True, agg_freq=Q,
        bre=False, gro=False, nuc=False
    )
    ca = CellAverage(
        initial=ini,
        agg=True, agg_freq=Q,
        bre=False, gro=False, nuc=False
    )
    fp.simulate(
        end_time=END_TIME, steps=STEPS, write_every=EVERY  #, max_order=ORDER
    )
    ca.simulate(
        end_time=END_TIME, steps=STEPS, write_every=EVERY  #, max_order=ORDER
    )
    plot_initial_and_current([fp, ca])
    plot_moments_over_time([fp, ca], max_order=ORDER)


def demo_breakage_aggregation():
    """Demo combined breakage and aggregation.
    """
    print("Demoing combined breakage and aggregation.")
    ini = Grid.create_geometric_end(
        start=START, end=END, sec=SECTIONS, fact=FACTOR, func=f
    )
    fp = FixedPivot(
        initial=ini,
        bre=True, bre_freq=Gamma, child=beta,
        agg=True, agg_freq=Q,
        gro=False, nuc=False
    )
    ca = CellAverage(
        initial=ini,
        bre=True, bre_freq=Gamma, child=beta,
        agg=True, agg_freq=Q,
        gro=False, nuc=False
    )
    fp.simulate(
        end_time=END_TIME, steps=STEPS, write_every=EVERY  #, max_order=ORDER
    )
    ca.simulate(
        end_time=END_TIME, steps=STEPS, write_every=EVERY  #, max_order=ORDER
    )
    plot_initial_and_current([fp, ca])
    plot_moments_over_time([fp, ca], max_order=ORDER)


def demo_growth():
    """Demo pure growth.
    """
    print("Demoing pure growth.")
    ini = Grid.create_geometric_end(
        start=START, end=END, sec=SECTIONS, fact=FACTOR, func=f
    )
    fp = FixedPivot(
        initial=ini,
        gro=True, gro_rate=G,
        bre=False, agg=False, nuc=False
    )
    ca = CellAverage(
        initial=ini,
        gro=True, gro_rate=G,
        bre=False, agg=False, nuc=False
    )
    fp.simulate(
        end_time=END_TIME, steps=STEPS, write_every=EVERY  #, max_order=ORDER
    )
    ca.simulate(
        end_time=END_TIME, steps=STEPS, write_every=EVERY  #, max_order=ORDER
    )
    plot_initial_and_current([fp, ca])
    plot_moments_over_time([fp, ca], max_order=ORDER)


def demo_nucleation():
    """Demo pure nucleation.
    """
    print("Demoing pure nucleation.")
    ini = Grid.create_geometric_end(
        start=START, end=END, sec=SECTIONS, fact=FACTOR, func=f
    )
    fp = FixedPivot(
        initial=ini,
        nuc=True, nuc_rate=S,
        bre=False, agg=False, gro=False
    )
    ca = CellAverage(
        initial=ini,
        nuc=True, nuc_rate=S,
        bre=False, agg=False, gro=False
    )
    fp.simulate(
        end_time=END_TIME, steps=STEPS, write_every=EVERY  #, max_order=ORDER
    )
    ca.simulate(
        end_time=END_TIME, steps=STEPS, write_every=EVERY  #, max_order=ORDER
    )
    plot_initial_and_current([fp, ca])
    plot_moments_over_time([fp, ca], max_order=ORDER)


def demo_growth_nucleation():
    """Demo combined growth and nucleation.
    """
    print("Demoing combined growth and nucleation.")
    ini = Grid.create_geometric_end(
        start=START, end=END, sec=SECTIONS, fact=FACTOR, func=f
    )
    fp = FixedPivot(
        initial=ini,
        gro=True, gro_rate=G,
        nuc=True, nuc_rate=S,
        bre=False, agg=False
    )
    ca = CellAverage(
        initial=ini,
        gro=True, gro_rate=G,
        nuc=True, nuc_rate=S,
        bre=False, agg=False
    )
    fp.simulate(
        end_time=END_TIME, steps=STEPS, write_every=EVERY  #, max_order=ORDER
    )
    ca.simulate(
        end_time=END_TIME, steps=STEPS, write_every=EVERY  #, max_order=ORDER
    )
    plot_initial_and_current([fp, ca])
    plot_moments_over_time([fp, ca], max_order=ORDER)


def demo_all():
    """Demo combined breakage, aggregation, growth and nucleation.
    """
    print("Demoing combined breakage, aggregation, growth and nucleation.")
    ini = Grid.create_geometric_end(
        start=START, end=END, sec=SECTIONS, fact=FACTOR, func=f
    )
    fp = FixedPivot(
        initial=ini,
        bre=True, bre_freq=Gamma, child=beta,
        agg=True, agg_freq=Q,
        gro=True, gro_rate=G,
        nuc=True, nuc_rate=S
    )
    ca = CellAverage(
        initial=ini,
        bre=True, bre_freq=Gamma, child=beta,
        agg=True, agg_freq=Q,
        gro=True, gro_rate=G,
        nuc=True, nuc_rate=S
    )
    fp.simulate(
        end_time=END_TIME, steps=STEPS, write_every=EVERY  #, max_order=ORDER
    )
    ca.simulate(
        end_time=END_TIME, steps=STEPS, write_every=EVERY  #, max_order=ORDER
    )
    plot_initial_and_current([fp, ca])
    plot_moments_over_time([fp, ca], max_order=ORDER)


# MAIN: -----------------------------------------------------------------------

if __name__ == "__main__":
    demo_zero()
    demo_breakage()
    demo_aggregation()
    demo_breakage_aggregation()
    demo_growth()
    demo_nucleation()
    demo_growth_nucleation()
    demo_all()
