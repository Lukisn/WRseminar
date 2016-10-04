#!/usr/bin/env python3

from math import sqrt, exp, pi
from scipy.integrate import quad, quadrature
import numpy as np
import matplotlib.pyplot as plt


def step(x):
    if x >= 0:
        return 1
    else:  # x < 0
        return 0
vstep = np.vectorize(step)


def hstep(x):
    if x == 0:
        return 0.5
    elif x > 0:
        return 1
    else:  # x < 0
        return 0
vhstep = np.vectorize(hstep)


def dirac(x):
    if x == 0:
        return 1
    else:  # x != 0
        return 0
vdirac = np.vectorize(dirac)


def num_dirac_rect(x, a=1e-1):
    if 0 - a/2 <= x <= 0 + a/2:
        return 1
    else:  # x < 0 - a/2 or x > 0 + a/2
        return 0
vnum_dirac_rect = np.vectorize(num_dirac_rect)


def num_dirac_norm(x, a=1e-1):
    return 1 / (a * sqrt(pi)) * exp(-x**2 / a**2)
vnum_dirac_norm = np.vectorize(num_dirac_norm)


def main():
    # plot functions:
    xs = np.linspace(start=-2, stop=2, num=1000)
    ys = {"step": [], "hstep": [], "dirac": [],
          "num_dirac_rect": [], "num_dirac_norm": []}
    for x in xs:
        ys["step"].append(step(x))
        ys["hstep"].append(hstep(x))
        ys["dirac"].append(dirac(x))
        ys["num_dirac_rect"].append(num_dirac_rect(x))
        ys["num_dirac_norm"].append(num_dirac_norm(x))

    for key, value in ys.items():
        if "dirac" in key:
            plt.plot(xs, ys[key], ".-", label=key)

    plt.legend(loc="best", fontsize="small")
    plt.grid()
    plt.show()

    # calculate integrals:
    lower = -2
    upper = 2
    Is_quad = {
        "step": quad(step, lower, upper),
        "hstep": quad(hstep, lower, upper),
        "dirac": quad(dirac, lower, upper),
        "num_dirac_rect": quad(num_dirac_rect, lower, upper),
        "num_dirac_norm": quad(num_dirac_norm, lower, upper)
    }
    '''
    Is_quadrature = {
        "step": quadrature(vstep, lower, upper),
        "hstep": quadrature(vhstep, lower, upper),
        "dirac": quadrature(vdirac, lower, upper),
        "num_dirac_rect": quadrature(vnum_dirac_rect, lower, upper),
        "num_dirac_norm": quadrature(vnum_dirac_norm, lower, upper)
    }
    '''
    for key, I_quad in Is_quad.items():
        print("{k}: {v}".format(k=key, v=I_quad))


if __name__ == "__main__":
    main()
