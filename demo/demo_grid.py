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

    plt.plot(uni.pivots(), uni.particles(), "rx", label="uni part")
    plt.plot(uni.pivots(), uni.densities(), "ro", label="uni dens")

    plt.plot(geo_step.pivots(), geo_step.particles(), "gx", label="geo step part")
    plt.plot(geo_step.pivots(), geo_step.densities(), "go", label="geo step dens")

    plt.plot(geo_max.pivots(), geo_max.particles(), "bx", label="geo max part")
    plt.plot(geo_max.pivots(), geo_max.densities(), "bo", label="geo max dens")

    plt.xscale("log")
    plt.yscale("log")
    plt.legend(loc="best", fontsize="small")
    plt.grid()
    plt.show()


if __name__ == "__main__":
    demo_creation()
