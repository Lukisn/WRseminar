#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from WR.grid import Grid


def demo():

    # initial number density function:
    def f(v, N0=1, v0=1):
        return (N0/v0) * (v/v0) * np.exp(-v/v0)

    uniform = Grid.create_uniform(0, 20, 40, func=f)
    geometric = Grid.create_geometric(0, 20, 40, 1.15, func=f)

    analytic_pivots = np.linspace(0, 10, 1000)
    analytic_densities = f(analytic_pivots)
    plt.plot(analytic_pivots, analytic_densities)

    plt.plot(uniform.pivots(), uniform.particles(), "ro")
    plt.plot(uniform.pivots(), uniform.particle_densities(), "r-")

    plt.plot(geometric.pivots(), geometric.particles(), "go")
    plt.plot(geometric.pivots(), geometric.particle_densities(), "g-")

    plt.show()

if __name__ == "__main__":
    demo()
