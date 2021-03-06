#!/usr/bin/env python3
"""
Module implementing useful functions for use in the entire project.

This module contains general utility functions, mathematical functions and
shortcuts and also functions for numerical representation of mathematical
objects.
"""
from itertools import tee
from math import sqrt, exp, pi, sin


def pairwise(iterable):
    """Pairwise iterator recipe.

    s -> (s0,s1), (s1,s2), (s2, s3), ...
    taken from: https://docs.python.org/3/library/itertools.html

    :param iterable: iterable object.
    :raises TypeError: if object is not iterable.
    :return: chained pairs of iterable objects elements.
    :rtype: zip object. Use list() or tuple() on it to build a list or tuple.
    """
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def kronecker(i, j):
    """Function representing the mathematical Kronecker delta symbol.

    :param i: first argument.
    :param j: second argument.
    :return: 1 if i == j, 0 otherwise.
    :rtype: int.
    """
    if i == j:
        return 1
    else:
        return 0


def zero(*args, **kwargs):
    """Function always returning zero.

    :param args: optional arguments.
    :param kwargs: optional keyword arguments.
    :return: always zero.
    """
    return 0


def step(x):
    """Simple step function.
    
    This is the simple step function where the return value is
    0 for x < 0 and 1 for x >= 0.
    
    This function can also be used to shift the step and reverse the edge:
    - rising step:  ``step(x - POS)``
    - falling step: ``1 - step(x - POS)``

    :param x: function argument.
    :return: result.
    """
    if x < 0:
        return 0
    else:  # x >= 0
        return 1


def hstep(x):
    """Heavyside step function.
    
    This is the heavyside step function whe the return value is
    0 for x < 0, 0.5 for x = 0 and 1 for x > 0.
    
    This function can also be used to shift the step and reverse the edge:
    - rising step:  hstep(x - POS)
    - falling step: 1 - hstep(x - POS)

    :param x: function argument.
    :return: result.
    """
    if x == 0:
        return 0.5
    elif x > 0:
        return 1
    else:  # x < 0
        return 0


def dirac_norm(x, a=1e-3):
    """Numeric Dirac delta function using the gaussian normal distribution.

    :param x: function argument.
    :param a: square root of the variance of the normal distribution.
    :return: numeric approximation of the Dirac delta function.
    """
    return 1 / (a * sqrt(pi)) * exp(-x**2 / a**2)


def dirac_rect(x, a=1e-3):
    """Numeric Dirac delta function using a rectangle function (two steps).

    :param x: function argument.
    :param a: half width of the rectangle step.
    :return: numeric approximation of the Dirac delta function.
    """
    peak = 1 / (2 * a)
    if x == 0 - a or x == 0 + a:
        return peak / 2
    elif 0 - a < x < 0 + a:
        return peak
    else:  # x < 0 - a/2 or x > 0 + a/2
        return 0


def dirac_simple(x, a=1e-3):
    """Numeric Dirac delta function using a simple fraction approach.

    :param x: function argument.
    :param a: function parameter.
    :return: numeric approximation of the Dirac delta function.
    """
    return a / (pi * (x**2 + a**2))


def dirac_sin(x, a=1e-3):
    """Numeric Dirac delta function using a sinus approach.

    This function contains a singularity at x = 0!

    :param x: function argument.
    :param a: function parameter.
    :return: numeric approximation of the Dirac delta function.
    """
    return sin(x / a) / (pi * x)
