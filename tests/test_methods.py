#!/usr/bin/env python3
"""Test module for method classes in sectional.methods."""
# standard library imports:
import unittest
# application imports:
from sectional.methods import Method


class TestMethod(unittest.TestCase):

    def setUp(self):
        self.method = Method(initial=[1, 2, 3])

    def test_creation(self):
        self.assertIsNot(self.method._initial, self.method._current)
        self.assertIsNot(self.method._initial, self.method._previous)
        self.assertIsNot(self.method._current, self.method._previous)

    def test_do_time_step(self):
        with self.assertRaises(NotImplementedError):
            self.method.do_time_step(1)

    def test_simulate(self):
        with self.assertRaises(AssertionError):
            self.method.simulate(end_time=-1, steps=10)
        with self.assertRaises(AssertionError):
            self.method.simulate(end_time=1, steps=-10)

    def test_moments_to_file(self):
        pass

    def test_ndfs_to_files(self):
        pass
