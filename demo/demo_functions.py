#!/usr/bin/env python3

"""
Module demonstrating the basic usage of the predefined functions.
"""
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
from WR.functions import pairwise
from WR.functions import kronecker, zero
from WR.functions import step, hstep
from WR.functions import dirac_norm, dirac_simple, dirac_rect


def demo_utilities():
    """Demo utility functions: pairwise.
    """
    l = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    print("list: {}".format(l))
    print("iterating")
    for index, pair in enumerate(pairwise(l)):
        print(index, pair, l[index], l[index + 1])
    print("iterating in reversed order")
    for index, pair in enumerate(pairwise(reversed(l))):
        rev_index = -1 - index
        print(rev_index, pair, l[rev_index], l[rev_index - 1])


def demo_basic_math():
    """Demo basic math functions: kronecker, zero.
    """
    for i in range(1, 4):  # i = 1...3
        for j in range(1, 4):  # j = 1...3
            print("i,j=({},{})".format(i, j), end="")
            if kronecker(i, j) == 1:
                print(", diagonal")
            else:
                print()
    print(zero(), zero(1), zero(1, 2), zero(x=42), zero(1, 2, x=3))


def demo_step_functions():
    """Demo step functions: step, hstep.
    """
    # constants:
    START, END = -1, 1  # plot interval boundaries
    NUM = 101  # number of samples
    POS = 0.5  # position of the step
    FACTOR = 1  # factor for widening the integration boundaries
    LOWER, UPPER = FACTOR * START, FACTOR * END  # integration boundaries

    # generate discrete function data:
    xs = np.linspace(start=START, stop=END, num=NUM)
    ys = {"step": [], "hstep": []}
    for x in xs:
        ys["step"].append(step(x - POS))  # rising step at POS
        # ys["step"].append(1 - step(x - POS))  # falling step at POS
        ys["hstep"].append(hstep(x - POS))  # rising step at POS
        # ys["hstep"].append(1 - hstep(x - POS))  # falling step at POS

    # calculate integrals:
    integrals = {
        "step": quad(lambda x: step(x - POS), LOWER, UPPER),  # rising step
        # "step": quad(lambda x: 1 - step(x), LOWER, UPPER),  # falling step
        "hstep": quad(lambda x: hstep(x - POS), LOWER, UPPER),  # rising step
        # "hstep": quad(lambda x: 1 - hstep(x), LOWER, UPPER)  # falling step
    }
    for key, integral in integrals.items():
        print("{}: {}".format(key, integral))

    # plot functions:
    for key, value in ys.items():
        plt.plot(xs, ys[key], ".-", label="{}: {}".format(key, integrals[key]))
    plt.legend(loc="best", fontsize="small")
    plt.grid()
    plt.show()


def demo_dirac_functions():
    """Demo numeric Dirac delta functions: norm, simple, rect"""
    # constants:
    START, END = -1, 1  # plot interval boundaries
    NUM = 101  # number of samples
    POS = 0  # position of dirac peak
    PARAM = 1e-1  # dirac peak size parameter
    FACTOR = 10000  # factor for widening the integration boundaries
    LOWER, UPPER = FACTOR * START, FACTOR * END  # integration boundaries

    # generate discrete function data:
    xs = np.linspace(start=START, stop=END, num=NUM)
    ys = {"rect": [], "norm": [], "simple": []}
    for x in xs:
        ys["rect"].append(dirac_rect(x - POS, a=PARAM))
        ys["norm"].append(dirac_norm(x - POS, a=PARAM))
        ys["simple"].append(dirac_simple(x - POS, a=PARAM))

    # calculate integrals:
    integrals = {
        "rect": quad(dirac_rect, LOWER, UPPER, args=(PARAM,)),
        "norm": quad(dirac_norm, LOWER, UPPER, args=(PARAM,)),
        "simple": quad(dirac_simple, LOWER, UPPER, args=(PARAM,))
    }
    for key, integral in integrals.items():
        print("{}: {}".format(key, integral))

    # plot functions:
    for key, value in ys.items():
        plt.plot(xs, ys[key], ".-", label="{}: {}".format(key, integrals[key]))
    plt.legend(loc="best", fontsize="small")
    plt.grid()
    plt.show()


if __name__ == "__main__":
    #demo_utilities()
    demo_basic_math()
    #demo_step_functions()
    #demo_dirac_functions()
