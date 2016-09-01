#!/usr/bin/env python3

import unittest
from WR.grid import Section, Grid

# TODO: checkout Property Based Testing
# https://hypothesis.readthedocs.io/en/latest/index.html


class TestSection(unittest.TestCase):

    def test_creation(self):
        s = Section(0, 1)
        print(s)
        s = Section(-1.0, 1)
        s = Section(-2.3, -1.01)
        self.assertRaises(AssertionError, Section, 1, 0)

    def test_minmax_getters_setters(self):
        s = Section(0, 1)
        self.assertEquals(s.start, 0)
        self.assertEquals(s.end, 1)
        self.assertEquals(s.pivot, (s.end + s.start) / 2)
        self.assertEquals(s.pivot, 0.5)
        s.end = 3
        s.start = 2
        self.assertEquals(s.start, 2)
        self.assertEquals(s.end, 3)
        self.assertEquals(s.pivot, (s.end + s.start) / 2)
        self.assertEquals(s.pivot, 2.5)
        with self.assertRaises(AssertionError):
            s.end = 1
        with self.assertRaises(AssertionError):
            s.start = 4

    def test_range_getters_setters(self):
        s = Section(0, 1)
        self.assertEquals(s.range, (0, 1))
        s.range = (-2, 3)
        self.assertEquals(s.range, (-2, 3))
        with self.assertRaises(IndexError):
            s.range = (-2,)
        with self.assertRaises(TypeError):
            s.range = 2
        with self.assertRaises(AssertionError):
            s.range = (4, -2)
        s.range = (0, 1)
        s.range = (2, 3)
        s.range = (-2, -1)

    def test_particle_getters_setters(self):
        s = Section(0, 1, particles=100)
        self.assertEquals(s.particles, 100)
        self.assertEquals(s.particle_density, 100/s.size)
        s.particles = 200
        self.assertEquals(s.particles, 200)
        self.assertEquals(s.particle_density, 200 / s.size)
        with self.assertRaises(AssertionError):
            s.particles = -1


class TestGrid(unittest.TestCase):

    def test_creation(self):
        g = Grid(0, 1)
        g._info()
        self.assertEquals(g.start, 0)
        self.assertEquals(g.end, 1)
        self.assertRaises(AssertionError, Grid, 1, 0)

    def test_factory_creation(self):
        # uniform:
        u = Grid.create_uniform(start=0, end=10, steps=3)
        self.assertEquals(u.start, 0)
        self.assertEquals(u.end, 10)
        with self.assertRaises(AssertionError):
            Grid.create_uniform(start=0, end=-10, steps=1)
        self.assertRaises(AssertionError, Grid.create_uniform, 0, 1, steps=-1)
        # geometric step:
        g = Grid.create_geometric_step(start=0, end=10, factor=1, steps=3)
        self.assertEquals(u.start, 0)
        self.assertEquals(u.end, 10)
        with self.assertRaises(AssertionError):
            Grid.create_geometric_step(0, end=-10, factor=1, steps=10)
        with self.assertRaises(AssertionError):
            Grid.create_geometric_step(0, end=10, factor=-1, steps=10)
        with self.assertRaises(AssertionError):
            Grid.create_geometric_step(0, end=10, factor=1, steps=-10)
        # geometric maximum:
        g = Grid.create_geometric_end(start=0, end=10, factor=1.1, steps=3)
        self.assertEquals(u.start, 0)
        self.assertEquals(u.end, 10)
        with self.assertRaises(AssertionError):
            Grid.create_geometric_end(0, end=-10, factor=1, steps=10)
        with self.assertRaises(AssertionError):
            Grid.create_geometric_end(0, end=10, factor=-1, steps=10)
        with self.assertRaises(AssertionError):
            Grid.create_geometric_end(0, end=10, factor=1, steps=-10)

    def test_add_remove(self):
        g = Grid(4, 5)
        self.assertEquals(g.start, 4)
        self.assertEquals(g.end, 5)
        g.add_left(1)
        self.assertEquals(g.start, 3)
        self.assertEquals(g.end, 5)
        g.add_right(2)
        self.assertEquals(g.start, 3)
        self.assertEquals(g.end, 7)
        self.assertRaises(AssertionError, g.add_left, -1)
        self.assertRaises(AssertionError, g.add_right, -1)
        self.assertRaises(AssertionError, g.add_left, 1, -1)
        self.assertRaises(AssertionError, g.add_right, 1, -1)
        g.remove_left()
        self.assertEquals(g.start, 4)
        self.assertEquals(g.end, 7)
        g.remove_right()
        self.assertEquals(g.start, 4)
        self.assertEquals(g.end, 5)
        self.assertEquals(len(g), 1)
        self.assertRaises(IndexError, g.remove_left)
        self.assertRaises(IndexError, g.remove_right)

    def test_coarsening_refining(self):
        g = Grid.create_uniform(start=1, end=10, steps=9)
        self.assertEquals(g.start, 1)
        self.assertEquals(g.end, 10)
        self.assertEquals(len(g), 9)
        # no range/no effect:
        g.coarsen_range(1, 1)
        g.coarsen_range(8, 8)
        self.assertEquals(g.start, 1)
        self.assertEquals(g.end, 10)
        self.assertEquals(len(g), 9)
        # effect:
        g.coarsen_range(3, 7)
        self.assertEquals(g.start, 1)
        self.assertEquals(g.end, 10)
        self.assertEquals(len(g), 5)
        self.assertEquals(g._sections[3].start, 4)
        self.assertEquals(g._sections[3].end, 9)
        self.assertEquals(g._sections[3].pivot, 6.5)
        g.refine_range(3, 3)
        self.assertEquals(g.start, 1)
        self.assertEquals(g.end, 10)
        self.assertEquals(len(g), 6)
        g.refine_range(3, 5)
        self.assertEquals(g.start, 1)
        self.assertEquals(g.end, 10)
        self.assertEquals(len(g), 9)
        self.assertRaises(AssertionError, g.coarsen_range, 1, 0)
        self.assertRaises(AssertionError, g.coarsen_range, -1, 1)
        self.assertRaises(AssertionError, g.coarsen_range, 1, 99)
        self.assertRaises(AssertionError, g.refine_range, 1, 0)
        self.assertRaises(AssertionError, g.refine_range, -1, 1)
        self.assertRaises(AssertionError, g.refine_range, 1, 99)


if __name__ == "__main__":
    unittest.main()
