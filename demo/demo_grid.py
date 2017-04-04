#!/usr/bin/env python3
"""Module demonstrating the basic usage of the Grid class."""
from math import exp
import matplotlib.pyplot as plt
import numpy as np
from sectional.grid import Grid


def f(v, N0=1, v0=1):
    """initial number density function."""
    return (N0/v0) * (v/v0) * exp(-v/v0)
    # return 1 * v ** 2
vf = np.vectorize(f)


def demo_creation():
    """Demo creation of discrete NDFs from a given function."""
    # constants:
    START, END = 1e-2, 1e1  # discrete size range
    SECTIONS = 10  # number of sections
    FACTOR = 1.5  # factor for neighboring sections sizes (geometric grid)
    CORR = True  # flag for correction of last section
    XSCALE = "log"  # log/linear - for all x scales
    YSCALE = "log"  # log/linear - for NDF y scale
    SCALEFACTOR = 1

    # analytic function:
    pivots = np.linspace(START, END, 1000)
    densities = vf(pivots)

    # generate grids:
    grids = {
        "uni": Grid.create_uniform(start=START, end=END, sec=SECTIONS, func=f, corr=CORR),
        #"step": Grid.create_geometric_step(start=START, end=END, uni_sec=SECTIONS, fact=FACTOR, func=f, corr=CORR),
        "end": Grid.create_geometric_end(start=START, end=END, sec=SECTIONS, fact=FACTOR, func=f, corr=CORR)
    }

    # calculate errors:
    errors = {}  # {"uni": [], "step": [], "end": []}
    rel_errors = {}
    for key in grids.keys():
        errors[key] = []
        rel_errors[key] = []
        for x, y in zip(grids[key].pivots(), grids[key].densities()):
            err = y - f(x)
            rel_err = abs(err / f(x))
            errors[key].append(err)
            rel_errors[key].append(rel_err)
    print(errors)

    # plot data:
    plt.subplot(311)
    plt.plot(pivots, densities, "y.-", label="analytic")
    for key, grid in grids.items():
        plt.plot(grids[key].pivots(), grids[key].densities(), ".", label=key)
        #plt.plot(grids[key].pivots(), grids[key].particles(), "x", label="{} particles".format(key))
    plt.xlim(START / SCALEFACTOR, END * SCALEFACTOR)
    plt.xscale(XSCALE)
    plt.yscale(YSCALE)
    plt.legend(loc="best", fontsize="small")
    plt.grid()

    # plot errors:
    plt.subplot(312)
    for key, errs in errors.items():
        plt.plot(grids[key].pivots(), errs, ".-", label=key)
    plt.xlim(START / SCALEFACTOR, END * SCALEFACTOR)
    plt.xscale(XSCALE)
    #plt.ylim(-1e-3, 1e-3)
    plt.legend(loc="best", fontsize="small")
    plt.grid()

    # plot relative errors:
    plt.subplot(313)
    for key, errs in rel_errors.items():
        plt.plot(grids[key].pivots(), errs, ".-", label=key)
    plt.xlim(START / SCALEFACTOR, END * SCALEFACTOR)
    plt.xscale(XSCALE)
    #plt.ylim(0, 1)
    plt.legend(loc="best", fontsize="small")
    plt.grid()

    plt.show()


if __name__ == "__main__":
    demo_creation()
