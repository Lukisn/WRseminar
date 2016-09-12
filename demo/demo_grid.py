#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from WR.grid import Grid


def f(v, N0=1, v0=1):
    """initial number density function.
    """
    return (N0/v0) * (v/v0) * np.exp(-v/v0)


def demo_creation():
    """Demo creation of discrete NDFs from a given function.
    """
    pivots = np.linspace(0, 10, 1000)
    densities = f(pivots)

    start= 0
    end = 20
    steps = 20
    factor = 1.35
    corr = True

    uni = Grid.create_uniform(start, end, steps, f, corr)
    geo_step = Grid.create_geometric_step(start, end, steps, factor, f, corr)
    geo_max = Grid.create_geometric_end(start, end, steps, factor, f, corr)

    plt.plot(pivots, densities, "y-", lw=2, label="ana dens")
    plt.plot(uni.pivots(), uni.particles(), "ro", label="uni part")
    plt.plot(uni.pivots(), uni.densities(), "r-", label="uni dens")
    plt.plot(geo_step.pivots(), geo_step.particles(), "go", label="geo step part")
    plt.plot(geo_step.pivots(), geo_step.densities(), "g-", label="geo step dens")
    plt.plot(geo_max.pivots(), geo_max.particles(), "bo", label="geo max part")
    plt.plot(geo_max.pivots(), geo_max.densities(), "b-", label="geo max dens")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(loc="best", fontsize="small")
    plt.grid()
    plt.show()


# TODO: revisit this when implementation of coarsen and refine is done!
def demo_manipulation():
    """Demo manipulationf of the NDFs.
    """
    raise NotImplementedError

    pivots = np.linspace(0, 10, 1000)
    densities = f(pivots)

    start = 0
    end = 20
    steps = 13
    factor = 1.3
    corr = True

    ini = Grid.create_geometric_end(start, end, steps, factor, f, corr)
    man = Grid.create_geometric_end(start, end, steps, factor, f, corr)

    man.refine(2)
    man.coarsen(2)

    plt.plot(pivots, densities, "y-", lw=2, label="ana dens")
    plt.plot(ini.pivots(), ini.particles(), "ro", label="ini part")
    plt.plot(ini.pivots(), ini.densities(), "r-", label="ini dens")
    plt.plot(man.pivots(), man.particles(), "go", label="man part")
    plt.plot(man.pivots(), man.densities(), "g-", label="man dens")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    demo_creation()
    #demo_manipulation()
