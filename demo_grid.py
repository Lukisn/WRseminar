#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from WR.grid import Grid


def demo():

    # initial number density function:
    def f(v, N0=1, v0=1):
        return (N0/v0) * (v/v0) * np.exp(-v/v0)

    start= 0
    end = 20
    steps = 42
    factor = 1.35
    corr = True

    uni = Grid.create_uniform(start, end, steps, f, corr)
    geo_step = Grid.create_geometric_step(start, end, steps, factor, f, corr)
    geo_max = Grid.create_geometric_end(start, end, steps, factor, f, corr)

    analytic_pivots = np.linspace(0, 10, 1000)
    analytic_densities = f(analytic_pivots)

    print("boundaries:")
    print("uniform  ", uni.boundaries())
    print("geom_step", geo_step.boundaries())
    print("geom_max ", geo_max.boundaries())

    plt.plot(analytic_pivots, analytic_densities, "y", lw=2, label="ana dens")
    plt.plot(uni.pivots(), uni.particles(), "ro", label="uni part")
    plt.plot(uni.pivots(), uni.particle_densities(), "r-", label="uni dens")
    plt.plot(geo_step.pivots(), geo_step.particles(), "go", label="geo step part")
    plt.plot(geo_step.pivots(), geo_step.particle_densities(), "g-", label="geo step dens")
    plt.plot(geo_max.pivots(), geo_max.particles(), "bo", label="geo max part")
    plt.plot(geo_max.pivots(), geo_max.particle_densities(), "b-", label="geo max dens")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    demo()
