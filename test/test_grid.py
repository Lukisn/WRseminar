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
        self.assertRaises(ValueError, Section, 1, 0)

    def test_minmax_getters_setters(self):
        s = Section(0, 1)
        self.assertEquals(s.minimum, 0)
        self.assertEquals(s.maximum, 1)
        self.assertEquals(s.pivot, (s.maximum + s.minimum) / 2)
        self.assertEquals(s.pivot, 0.5)
        s.maximum = 3
        s.minimum = 2
        self.assertEquals(s.minimum, 2)
        self.assertEquals(s.maximum, 3)
        self.assertEquals(s.pivot, (s.maximum + s.minimum) / 2)
        self.assertEquals(s.pivot, 2.5)
        with self.assertRaises(ValueError):
            s.maximum = 1
        with self.assertRaises(ValueError):
            s.minimum = 4

    def test_range_getters_setters(self):
        s = Section(0, 1)
        self.assertEquals(s.range, (0, 1))
        s.range = (-2, 3)
        self.assertEquals(s.range, (-2, 3))
        with self.assertRaises(IndexError):
            s.range = (-2,)
        with self.assertRaises(TypeError):
            s.range = 2
        with self.assertRaises(ValueError):
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
        with self.assertRaises(ValueError):
            s.particles = -1


class TestGrid(unittest.TestCase):

    def test_creation(self):
        g = Grid(0, 1)
        g.info()
        self.assertEquals(g.minimum, 0)
        self.assertEquals(g.maximum, 1)
        self.assertRaises(ValueError, Grid, 1, 0)

    def test_factory_creation(self):
        # uniform:
        u = Grid.create_uniform(minimum=0, maximum=10, amount=3)
        self.assertEquals(u.minimum, 0)
        self.assertEquals(u.maximum, 10)
        with self.assertRaises(ValueError):
            Grid.create_uniform(minimum=0, maximum=-10, amount=1)
        self.assertRaises(ValueError, Grid.create_uniform, 0, 1, amount=-1)
        # geometric step:
        g = Grid.create_geometric_step(minimum=0, maximum=10, factor=1, amount=3)
        self.assertEquals(u.minimum, 0)
        self.assertEquals(u.maximum, 10)
        with self.assertRaises(ValueError):
            Grid.create_geometric_step(0, maximum=-10, factor=1, amount=10)
        with self.assertRaises(ValueError):
            Grid.create_geometric_step(0, maximum=10, factor=-1, amount=10)
        with self.assertRaises(ValueError):
            Grid.create_geometric_step(0, maximum=10, factor=1, amount=-10)
        # geometric maximum:
        g = Grid.create_geometric_maximum(minimum=0, maximum=10, factor=1.1, amount=3)
        self.assertEquals(u.minimum, 0)
        self.assertEquals(u.maximum, 10)
        with self.assertRaises(ValueError):
            Grid.create_geometric_maximum(0, maximum=-10, factor=1, amount=10)
        with self.assertRaises(ValueError):
            Grid.create_geometric_maximum(0, maximum=10, factor=-1, amount=10)
        with self.assertRaises(ValueError):
            Grid.create_geometric_maximum(0, maximum=10, factor=1, amount=-10)

    def test_add_remove(self):
        g = Grid(4, 5)
        self.assertEquals(g.minimum, 4)
        self.assertEquals(g.maximum, 5)
        g.add_left(1)
        self.assertEquals(g.minimum, 3)
        self.assertEquals(g.maximum, 5)
        g.add_right(2)
        self.assertEquals(g.minimum, 3)
        self.assertEquals(g.maximum, 7)
        self.assertRaises(ValueError, g.add_left, -1)
        self.assertRaises(ValueError, g.add_right, -1)
        self.assertRaises(ValueError, g.add_left, 1, -1)
        self.assertRaises(ValueError, g.add_right, 1, -1)
        g.remove_left()
        self.assertEquals(g.minimum, 4)
        self.assertEquals(g.maximum, 7)
        g.remove_right()
        self.assertEquals(g.minimum, 4)
        self.assertEquals(g.maximum, 5)
        self.assertEquals(len(g), 1)
        self.assertRaises(IndexError, g.remove_left)
        self.assertRaises(IndexError, g.remove_right)

    def test_coarsening_refining(self):
        g = Grid.create_uniform(minimum=1, maximum=10, amount=9)
        self.assertEquals(g.minimum, 1)
        self.assertEquals(g.maximum, 10)
        self.assertEquals(len(g), 9)
        # no range/no effect:
        g.coarsen(1, 1)
        g.coarsen(8, 8)
        self.assertEquals(g.minimum, 1)
        self.assertEquals(g.maximum, 10)
        self.assertEquals(len(g), 9)
        # effect:
        g.coarsen(3, 7)
        self.assertEquals(g.minimum, 1)
        self.assertEquals(g.maximum, 10)
        self.assertEquals(len(g), 5)
        self.assertEquals(g._sections[3].minimum, 4)
        self.assertEquals(g._sections[3].maximum, 9)
        self.assertEquals(g._sections[3].pivot, 6.5)
        g.refine(3, 3)
        self.assertEquals(g.minimum, 1)
        self.assertEquals(g.maximum, 10)
        self.assertEquals(len(g), 6)
        g.refine(3, 5)
        self.assertEquals(g.minimum, 1)
        self.assertEquals(g.maximum, 10)
        self.assertEquals(len(g), 9)
        self.assertRaises(ValueError, g.coarsen, 1, 0)
        self.assertRaises(ValueError, g.coarsen, -1, 1)
        self.assertRaises(ValueError, g.coarsen, 1, 99)
        self.assertRaises(ValueError, g.refine, 1, 0)
        self.assertRaises(ValueError, g.refine, -1, 1)
        self.assertRaises(ValueError, g.refine, 1, 99)


if __name__ == "__main__":
    unittest.main()
