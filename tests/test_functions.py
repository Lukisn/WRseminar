#!/usr/bin/env python3
"""Test module for the functions defined in sectional.functions."""
# standard library imports:
import unittest
# third party imports:
from scipy.integrate import quad
# application imports:
from sectional.functions import pairwise, kronecker, zero, step, hstep
from sectional.functions import dirac_norm, dirac_rect, dirac_simple, dirac_sin


class TestPairwise(unittest.TestCase):

    def test_iterable(self):
        l = [1, 2, 3, 4]
        pw = list(pairwise(l))
        self.assertEqual([(1, 2), (2, 3), (3, 4)], pw)

    def test_noniterable(self):
        self.assertRaises(TypeError, pairwise, 1)
        self.assertRaises(TypeError, pairwise, 2.1)


class TestKronecker(unittest.TestCase):

    def test_same_inputs(self):
        self.assertEqual(1, kronecker(0, 0))
        self.assertEqual(1, kronecker(1, 1))
        self.assertEqual(1, kronecker(-1, -1))
        self.assertEqual(1, kronecker(2.1, 2.1))
        self.assertEqual(1, kronecker(2.0, 2))

    def test_different_inputs(self):
        self.assertEqual(0, kronecker(0, 1))
        self.assertEqual(0, kronecker(1, 0))
        self.assertEqual(0, kronecker(1, -1))
        self.assertEqual(0, kronecker(2, 2.000000001))


class TestZero(unittest.TestCase):

    def test_no_input(self):
        self.assertEqual(0, zero())

    def test_input_args(self):
        self.assertEqual(0, zero(0))
        self.assertEqual(0, zero(0, 1.1, "2", []))

    def test_input_kwargs(self):
        self.assertEqual(0, zero(kw=1))
        self.assertEqual(0, zero(kw1=1, kw2=2.1, kw3=["3"]))


class TestStep(unittest.TestCase):

    def test_values(self):
        self.assertEqual(0, step(-1))
        self.assertEqual(1, step(0))
        self.assertEqual(1, step(1))

    def test_integral(self):
        result, *rest = quad(lambda x: step(x), -1, 1)
        self.assertAlmostEqual(1, result, places=15)


class TestHStep(unittest.TestCase):

    def test_values(self):
        self.assertEqual(0, hstep(-1))
        self.assertEqual(0.5, hstep(0))
        self.assertEqual(1, hstep(1))

    def test_integral(self):
        result, *rest = quad(lambda x: hstep(x), -1, 1)
        self.assertAlmostEqual(1, result, places=15)


class TestDirac(unittest.TestCase):

    def test_dirac_norm(self):
        result, *rest = quad(
            lambda x: dirac_norm(x), -1, 1, points=[0, 0-1e-3, 0+1e-3]
        )
        self.assertAlmostEqual(1, result, places=15)

    def test_dirac_rect(self):
        result, *rest = quad(
            lambda x: dirac_rect(x), -1, 1, points=[0, 0-1e-3, 0+1e-3]
        )
        self.assertAlmostEqual(1, result, places=15)

    def test_dirac_simple(self):
        result, *rest = quad(
            lambda x: dirac_simple(x), -1, 1, points=[0]
        )
        self.assertAlmostEqual(1, result, places=2)

    def test_dirac_sin(self):
        result, *rest = quad(
            lambda x: dirac_sin(x), -1, 1, points=[0]
        )
        self.assertAlmostEqual(1, result, places=2)
