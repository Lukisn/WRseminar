#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from WR.grid import Grid


def demo():

    # initial number density function:
    def f(v, N0=1, v0=1):
        return (N0/v0) * (v/v0) * np.exp(-v/v0)

    uniform = Grid.create_uniform(0, 20, 13, func=f)
    geom_step = Grid.create_geometric_step(0, 20, 13, 1.15, func=f)
    geom_max = Grid.create_geometric_maximum(0, 20, 13, 1.15, func=f)

    analytic_pivots = np.linspace(0, 10, 1000)
    analytic_densities = f(analytic_pivots)

    print("boundaries:")
    print("uniform  ", uniform.boundaries())
    print("geom_step", geom_step.boundaries())
    print("geom_max ", geom_max.boundaries())

    plt.plot(analytic_pivots, analytic_densities, "k-.", label="ana dens")
    plt.plot(uniform.pivots(), uniform.particles(), "ro", label="uni part")
    plt.plot(uniform.pivots(), uniform.particle_densities(), "r-", label="uni dens")
    plt.plot(geom_step.pivots(), geom_step.particles(), "go", label="geo step part")
    plt.plot(geom_step.pivots(), geom_step.particle_densities(), "g-", label="geo step dens")
    plt.plot(geom_max.pivots(), geom_max.particles(), "bo", label="geo max part")
    plt.plot(geom_max.pivots(), geom_max.particle_densities(), "b-", label="geo max dens")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    demo()
