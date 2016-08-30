#!/usr/bin/env python3

from WR.grid import Grid
import matplotlib.pyplot as plt


def demo():
    """Show very basic grid usage (factory method creation, coarsening,
    refining) and plot the results for visualization.
    """

    # create uniform grid:
    uni = Grid.create_uniform(minimum=0, size=1, amount=10)
    uni.print_info()
    xbe0 = uni.boundary_list()
    ybe0 = [0] * len(xbe0)
    xpe0 = uni.pivot_list()
    ype0 = [0] * len(xpe0)
    # coarsen
    uni.coarsen_range(start=5, end=9)
    uni.print_info()
    xbe1 = uni.boundary_list()
    ybe1 = [1] * len(xbe1)
    xpe1 = uni.pivot_list()
    ype1 = [1] * len(xpe1)
    # refine
    uni.refine_range(start=3, end=5)
    uni.print_info()
    xbe2 = uni.boundary_list()
    ybe2 = [2] * len(xbe2)
    xpe2 = uni.pivot_list()
    ype2 = [2] * len(xpe2)

    # create non-equidistant grid
    geo = Grid.create_geometric(minimum=0, initial_size=1, factor=1.1, amount=10)
    geo.print_info()
    xbf0 = geo.boundary_list()
    ybf0 = [0] * len(xbf0)
    xpf0 = geo.pivot_list()
    ypf0 = [0] * len(xpf0)
    # coarsen
    geo.coarsen_range(start=5, end=9)
    geo.print_info()
    xbf1 = geo.boundary_list()
    ybf1 = [1] * len(xbf1)
    xpf1 = geo.pivot_list()
    ypf1 = [1] * len(xpf1)
    # refine
    geo.refine_range(start=3, end=5)
    geo.print_info()
    xbf2 = geo.boundary_list()
    ybf2 = [2] * len(xbf2)
    xpf2 = geo.pivot_list()
    ypf2 = [2] * len(xpf2)

    # plot grids
    plt.title("equidistant grid")
    plt.plot(xbe0, ybe0, "rx", label="initial boundary")
    plt.plot(xpe0, ype0, "ro", label="initial pivot")
    plt.plot(xbe1, ybe1, "gx", label="coarsened boundary")
    plt.plot(xpe1, ype1, "go", label="coarsened pivot")
    plt.plot(xbe2, ybe2, "bx", label="refined boundary")
    plt.plot(xpe2, ype2, "bo", label="refined pivot")
    plt.xlim(-1, 11)
    plt.ylim(-1, 5)
    plt.legend()
    plt.grid()
    plt.show()

    plt.title("factor grid")
    plt.plot(xbf0, ybf0, "rx", label="initial boundary")
    plt.plot(xpf0, ypf0, "ro", label="initial pivot")
    plt.plot(xbf1, ybf1, "gx", label="coarsened boundary")
    plt.plot(xpf1, ypf1, "go", label="coarsened pivot")
    plt.plot(xbf2, ybf2, "bx", label="refined boundary")
    plt.plot(xpf2, ypf2, "bo", label="refined pivot")
    plt.xlim(-1, 17)
    plt.ylim(-1, 5)
    plt.legend()
    plt.grid()
    plt.show()


if __name__ == "__main__":
    demo()
