#!/usr/bin/env python3

from WR.grid import Grid
from WR.methods import FixedPivot
import numpy as np
import matplotlib.pyplot as plt

MIN = 0
MAX = 1000
SECTIONS = 100
FACTOR = 1.1

END_TIME = 1
STEPS = 10


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
    return 1


def G(v):
    """growth rate for particles of size v.
    """
    return 1


def S(v):
    """rate of nucleation of particles of size v.
    """
    if v <= 1:
        return 1
    else:
        return 0


def plot_initial_and_current(method):
    plt.plot(method.initial_ndf.pivots(),
             method.initial_ndf.particle_densities(), "g-", label="ini")
    plt.plot(method.current_ndf.pivots(),
             method.current_ndf.particle_densities(), "ro", label="cur")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(1.e-5, 1.e5)
    plt.ylim(1.e-5, 1.e5)
    plt.legend()
    plt.grid()
    plt.show()


def demo_breakage():
    """Demo pure breakage.
    """
    print("Demoing pure breakage.")

    ini = Grid.create_geometric_end(start=MIN, end=MAX, sections=SECTIONS,
                                    factor=FACTOR, func=f)
    #ini._plot()

    method = FixedPivot(initial_ndf=ini,
                        breakage=True,
                        break_freq=gamma, child_number=beta,
                        aggregation=False, growth=False, nucleation=False)
    method.simulate(end_time=END_TIME, steps=STEPS)
    #method.current_ndf._plot()

    plot_initial_and_current(method)

    print("done.\n")


def demo_aggregation():
    """Demo pure aggregation.
    """
    print("Demoing pure aggregation.")

    ini = Grid.create_geometric_end(start=MIN, end=MAX, sections=SECTIONS,
                                    factor=FACTOR, func=f)
    #ini._plot()

    method = FixedPivot(initial_ndf=ini,
                        aggregation=True,
                        agg_freq=Q,
                        breakage=False, growth=False, nucleation=False)
    method.simulate(end_time=END_TIME, steps=STEPS)
    #method.current_ndf._plot()

    plot_initial_and_current(method)

    print("done.\n")


def demo_breakage_aggregation():
    """Demo combined breakage and aggregation.
    """
    print("Demoing combined breakage and aggregation.")

    ini = Grid.create_geometric_end(start=MIN, end=MAX, sections=SECTIONS,
                                    factor=FACTOR, func=f)
    #ini._plot()

    method = FixedPivot(initial_ndf=ini,
                        breakage=True,
                        break_freq=gamma, child_number=beta,
                        aggregation=True,
                        agg_freq=Q,
                        growth=False, nucleation=False)
    method.simulate(end_time=END_TIME, steps=STEPS)
    #method.current_ndf._plot()

    plot_initial_and_current(method)

    print("done.\n")


def demo_growth():
    """Demo pure growth.
    """
    print("Demoing pure growth.")

    ini = Grid.create_geometric_end(start=MIN, end=MAX, sections=SECTIONS,
                                    factor=FACTOR, func=f)
    # ini._plot()

    method = FixedPivot(initial_ndf=ini,
                        growth=True,
                        gro_rate=G,
                        breakage=False, aggregation=False, nucleation=False)
    method.simulate(end_time=END_TIME, steps=STEPS)
    # method.current_ndf._plot()

    plot_initial_and_current(method)

    print("done.\n")


def demo_nucleation():
    """Demo pure nucleation.
    """
    print("Demoing pure nucleation.")

    ini = Grid.create_geometric_end(start=MIN, end=MAX, sections=SECTIONS,
                                    factor=FACTOR, func=f)
    # ini._plot()

    method = FixedPivot(initial_ndf=ini,
                        nucleation=True,
                        nuc_rate=S,
                        breakage=False, aggregation=False, growth=False)
    method.simulate(end_time=END_TIME, steps=STEPS)
    # method.current_ndf._plot()

    plot_initial_and_current(method)

    print("done.\n")

# TODO: write code for combinations
def demo_growth_nucleation():
    """Demo combined growth and nucleation.
    """
    print("Demoing combined growth and nucleation.")
    print("done.\n")


def demo_all():
    """Demo combined breakage, aggregation, growth and nucleation.
    """
    print("Demoing combined breakage, aggregation, growth and nucleation.")
    print("done.\n")


if __name__ == "__main__":
    #demo_breakage()
    #demo_aggregation()
    #demo_breakage_aggregation()
    demo_growth()
    #demo_nucleation()
    #demo_growth_nucleation()
    #demo_all()
