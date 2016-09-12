#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from WR.grid import Grid
from WR.methods import FixedPivot, CellAverage

MIN = 0
MAX = 1000
SECTIONS = 50
FACTOR = 1.3

END_TIME = 1
STEPS = 100

#METHOD = "fixed pivot"
METHOD = "cell average"


def f(v, N0=1, v0=1):
    """initial number density function.
    """
    return (N0 / v0) * (v / v0) * np.exp(-v / v0)
    #return (N0 / v0) * np.exp(-v / v0)


def gamma(v):
    """breakage frequency function.
    """
    return v ** 2


def beta(v1, v2):
    """child number function.
    """
    return 2/v2  # v2/2


def Q(x1, x2):
    """aggregation frequency function.
    """
    return x1 + x2


def G(v):
    """growth rate for particles of size v.
    """
    return 0.1


def S(v):
    """rate of nucleation of particles of size v.
    """
    if v <= 0.01:
        return 1
    else:
        return 0


def plot_initial_and_current(method):
    plt.plot(method._initial.pivots(),
             method._initial.densities(), "g.-", label="initial")
    plt.plot(method._current.pivots(),
             method._current.densities(), "r.-", label="current")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(1.e-5, 1.e5)
    plt.ylim(1.e-5, 1.e5)
    plt.legend(loc="best", fontsize="small")
    plt.grid()
    plt.show()


def demo_zero():
    """Demo basic behaviour without specifying the keyword arguments.
    """
    print("Demoing zero behaviour.")
    ini = Grid.create_geometric_end(
        start=MIN, end=MAX, sections=SECTIONS, factor=FACTOR, func=f
    )
    if METHOD == "fixed pivot":
        method = FixedPivot(initial=ini)
    elif METHOD == "cell average":
        method = CellAverage(initial=ini)
    else:
        raise ValueError("unknown method '{}'!".format(METHOD))
    method.simulate(end_time=END_TIME, steps=STEPS)
    plot_initial_and_current(method)


def demo_breakage():
    """Demo pure breakage.
    """
    print("Demoing pure breakage.")
    ini = Grid.create_geometric_end(
        start=MIN, end=MAX, sections=SECTIONS, factor=FACTOR, func=f
    )
    if METHOD == "fixed pivot":
        method = FixedPivot(
            initial=ini,
            bre=True, bre_freq=gamma, child=beta,
            agg=False, gro=False, nuc=False
        )
    elif METHOD == "cell average":
        method = CellAverage(
            initial=ini,
            bre=True, bre_freq=gamma, child=beta,
            agg=False, gro=False, nuc=False
        )
    else:
        raise ValueError("unknown method '{}'!".format(METHOD))
    method.simulate(end_time=END_TIME, steps=STEPS)
    plot_initial_and_current(method)


def demo_aggregation():
    """Demo pure aggregation.
    """
    print("Demoing pure aggregation.")
    ini = Grid.create_geometric_end(
        start=MIN, end=MAX, sections=SECTIONS, factor=FACTOR, func=f
    )
    if METHOD == "fixed pivot":
        method = FixedPivot(
            initial=ini,
            agg=True, agg_freq=Q,
            bre=False, gro=False, nuc=False
        )
    elif METHOD == "cell average":
        method = CellAverage(
            initial=ini,
            agg=True, agg_freq=Q,
            bre=False, gro=False, nuc=False
        )
    else:
        raise ValueError("unknown method '{}'!".format(METHOD))
    method.simulate(end_time=END_TIME, steps=STEPS)
    plot_initial_and_current(method)


def demo_breakage_aggregation():
    """Demo combined breakage and aggregation.
    """
    print("Demoing combined breakage and aggregation.")
    ini = Grid.create_geometric_end(
        start=MIN, end=MAX, sections=SECTIONS, factor=FACTOR, func=f
    )
    if METHOD == "fixed pivot":
        method = FixedPivot(
            initial=ini,
            bre=True, bre_freq=gamma, child=beta,
            agg=True, agg_freq=Q,
            gro=False, nuc=False
        )
    elif METHOD == "cell average":
        method = CellAverage(
            initial=ini,
            bre=True, bre_freq=gamma, child=beta,
            agg=True, agg_freq=Q,
            gro=False, nuc=False
        )
    else:
        raise ValueError("unknown method '{}'!".format(METHOD))
    method.simulate(end_time=END_TIME, steps=STEPS)
    plot_initial_and_current(method)


def demo_growth():
    """Demo pure growth.
    """
    print("Demoing pure growth.")
    ini = Grid.create_geometric_end(
        start=MIN, end=MAX, sections=SECTIONS, factor=FACTOR, func=f
    )
    if METHOD == "fixed pivot":
        method = FixedPivot(
            initial=ini,
            gro=True, gro_rate=G,
            bre=False, agg=False, nuc=False
        )
    elif METHOD == "cell average":
        method = CellAverage(
            initial=ini,
            gro=True, gro_rate=G,
            bre=False, agg=False, nuc=False
        )
    else:
        raise ValueError("unknown method '{}'!".format(METHOD))
    method.simulate(end_time=END_TIME, steps=STEPS)
    plot_initial_and_current(method)


def demo_nucleation():
    """Demo pure nucleation.
    """
    print("Demoing pure nucleation.")
    ini = Grid.create_geometric_end(
        start=MIN, end=MAX, sections=SECTIONS, factor=FACTOR, func=f
    )
    if METHOD == "fixed pivot":
        method = FixedPivot(
            initial=ini,
            nuc=True, nuc_rate=S,
            bre=False, agg=False, gro=False
        )
    elif METHOD == "cell average":
        method = CellAverage(
            initial=ini,
            nuc=True, nuc_rate=S,
            bre=False, agg=False, gro=False
        )
    else:
        raise ValueError("unknown method '{}'!".format(METHOD))
    method.simulate(end_time=END_TIME, steps=STEPS)
    plot_initial_and_current(method)


def demo_growth_nucleation():
    """Demo combined growth and nucleation.
    """
    print("Demoing combined growth and nucleation.")
    ini = Grid.create_geometric_end(
        start=MIN, end=MAX, sections=SECTIONS, factor=FACTOR, func=f
    )
    if METHOD == "fixed pivot":
        method = FixedPivot(
            initial=ini,
            gro=True, gro_rate=G,
            nuc=True, nuc_rate=S,
            bre=False, agg=False
        )
    elif METHOD == "cell average":
        method = CellAverage(
            initial=ini,
            gro=True, gro_rate=G,
            nuc=True, nuc_rate=S,
            bre=False, agg=False
        )
    else:
        raise ValueError("unknown method '{}'!".format(METHOD))
    method.simulate(end_time=END_TIME, steps=STEPS)
    plot_initial_and_current(method)


def demo_all():
    """Demo combined breakage, aggregation, growth and nucleation.
    """
    print("Demoing combined breakage, aggregation, growth and nucleation.")
    ini = Grid.create_geometric_end(
        start=MIN, end=MAX, sections=SECTIONS, factor=FACTOR, func=f
    )
    if METHOD == "fixed pivot":
        method = FixedPivot(
            initial=ini,
            bre=True, bre_freq=gamma, child=beta,
            agg=True, agg_freq=Q,
            gro=True, gro_rate=G,
            nuc=True, nuc_rate=S
        )
    elif METHOD == "cell average":
        method = CellAverage(
            initial=ini,
            bre=True, bre_freq=gamma, child=beta,
            agg=True, agg_freq=Q,
            gro=True, gro_rate=G,
            nuc=True, nuc_rate=S
        )
    else:
        raise ValueError("unknown method '{}'!".format(METHOD))
    method.simulate(end_time=END_TIME, steps=STEPS)
    plot_initial_and_current(method)


if __name__ == "__main__":
    print("demoing {} method".format(METHOD.upper()))
    demo_zero()
    demo_breakage()
    demo_aggregation()
    demo_breakage_aggregation()
    demo_growth()
    demo_nucleation()
    demo_growth_nucleation()
    demo_all()
