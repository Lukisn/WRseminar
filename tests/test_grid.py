#!/usr/bin/env python3

import unittest

from sectional.grid import find_initial_step, Section, Grid


class TestFunctions(unittest.TestCase):

    def test_find_initial_step(self):
        with self.assertRaises(ValueError):
            find_initial_step(start=1, end=0, steps=1, factor=1)  # start > end
        with self.assertRaises(ValueError):
            find_initial_step(start=0, end=1, steps=-1, factor=1)  # steps <= 0
        with self.assertRaises(ValueError):
            find_initial_step(start=0, end=1, steps=1, factor=0.1)  # factor < 1
        with self.assertRaises(ValueError):
            find_initial_step(start=0, end=1, steps=1, factor=1, max_err=-1)  # max_err <= 0
        with self.assertRaises(ValueError):
            find_initial_step(start=0, end=1, steps=1, factor=1, max_iter=-1)  # max_iter <= 0


class TestSection(unittest.TestCase):

    def test_creation(self):
        s = Section(start=0, end=1, particles=0)  # default section
        print(s)
        Section(start=-1.0, end=1, particles=10)  # section including zero
        Section(start=-2.3, end=-1.01, particles=1.001)  # negative section
        with self.assertRaises(ValueError):
            Section(start=1, end=0)  # start > end
        with self.assertRaises(ValueError):
            Section(start=0, end=1, particles=-1)  # negative particle number

    def test_start_end(self):
        s = Section(start=0, end=1)
        self.assertEquals(s.start, 0)
        self.assertEquals(s.end, 1)
        s.end = 3
        s.start = 2
        self.assertEquals(s.start, 2)
        self.assertEquals(s.end, 3)
        with self.assertRaises(ValueError):
            s.end = 1  # end < start
        with self.assertRaises(ValueError):
            s.start = 4  # start > end

    def test_interval(self):
        s = Section(start=0, end=1)  # default interval
        self.assertEquals(s.interval, (0, 1))
        s.interval = (-2, 3)  # interval including zero
        self.assertEquals(s.interval, (-2, 3))
        s.interval = (-2, -1)  # negative interval
        with self.assertRaises(TypeError):
            s.interval = 2
        with self.assertRaises(IndexError):
            s.interval = (-2,)
        with self.assertRaises(ValueError):
            s.interval = (4, -2)  # start > end

    def test_pivot(self):
        s = Section(0, 1)
        self.assertEquals(s.pivot, (s.end + s.start) / 2)
        self.assertEquals(s.pivot, 0.5)
        s.end = 4
        s.start = 2
        self.assertEquals(s.pivot, (s.end + s.start) / 2)
        self.assertEquals(s.pivot, 3)

    def test_particle(self):
        s = Section(start=0, end=1, particles=100)
        self.assertEquals(s.particles, 100)
        self.assertEquals(s.density, 100 / s.size)
        s.particles = 200
        self.assertEquals(s.particles, 200)
        self.assertEquals(s.density, 200 / s.size)
        with self.assertRaises(ValueError):
            s.particles = -1  # negative particle number


class TestGrid(unittest.TestCase):

    def test_creation(self):
        g = Grid()  # default grid
        if not g._is_seamless():
            raise RuntimeError("found seams!")
        g = Grid(start=0, end=1)  # default grid
        if not g._is_seamless():
            raise RuntimeError("found seams!")
        g._info()
        Grid(start=-1, end=2)  # grid including zero
        Grid(start=-3, end=-1)  # negative grid
        self.assertEquals(g.start, 0)
        self.assertEquals(g.end, 1)
        with self.assertRaises(ValueError):
            Grid(start=1, end=0)  # start > end
        with self.assertRaises(ValueError):
            Grid(start=0, end=1, particles=-1)  # particles < 0

    def test_section(self):
        g = Grid.create_uniform(start=0, end=1, sec=10)
        g._info()
        if not g._is_seamless():
            raise RuntimeError("found seams!")
        self.assertEquals(len(g), 10)
        _ = g[0]  # first section
        _ = g[-1]  # last section
        with self.assertRaises(IndexError):
            _ = g[len(g)]  # right outside index range

    def test_creation_uniform(self):
        u = Grid.create_uniform(start=0, end=10, sec=3)
        if not u._is_seamless():
            raise RuntimeError("found seams!")
        self.assertEquals(u.start, 0)
        self.assertEquals(u.end, 10)
        self.assertEquals(len(u), 3)
        with self.assertRaises(ValueError):
            Grid.create_uniform(start=0, end=-10, sec=1)  # start > end
        with self.assertRaises(ValueError):
            Grid.create_uniform(start=0, end=1, sec=-1)  # sections < 0

    def test_creation_geometric(self):
        g = Grid.create_geometric(start=0, end=10, ini_sec=1, fact=1.1)
        if not g._is_seamless():
            raise RuntimeError("found seams!")
        self.assertEqual(g.start, 0)
        self.assertEqual(g.end, 10)
        with self.assertRaises(ValueError):
            Grid.create_geometric(start=0, end=-10, ini_sec=1, fact=1)  # start > end
        with self.assertRaises(ValueError):
            Grid.create_geometric(start=0, end=1, ini_sec=-1, fact=1)  # initial_step < 0
        with self.assertRaises(ValueError):
            Grid.create_geometric(start=0, end=1, ini_sec=1, fact=0.1)  # factor < 1

    def test_creation_geometric_step(self):
        g = Grid.create_geometric_step(start=0, end=10, fact=1, uni_sec=3)
        if not g._is_seamless():
            raise RuntimeError("found seams!")
        self.assertEquals(g.start, 0)
        self.assertEquals(g.end, 10)
        with self.assertRaises(ValueError):  # start > end
            Grid.create_geometric_step(start=0, end=-10, fact=1, uni_sec=10)
        with self.assertRaises(ValueError):  # negative factor
            Grid.create_geometric_step(start=0, end=10, fact=-1, uni_sec=10)
        with self.assertRaises(ValueError):  # negative sections
            Grid.create_geometric_step(start=0, end=10, fact=1, uni_sec=-10)

    def test_creation_geometric_end(self):
        g = Grid.create_geometric_end(start=0, end=10, fact=1.1, sec=3)
        if not g._is_seamless():
            raise RuntimeError("found seams!")
        self.assertEquals(g.start, 0)
        self.assertEquals(g.end, 10)
        self.assertEqual(len(g), 3)
        with self.assertRaises(ValueError):  # start > end
            Grid.create_geometric_end(start=0, end=-10, fact=1, sec=10)
        with self.assertRaises(ValueError):  # negative factor
            Grid.create_geometric_end(start=0, end=10, fact=-1, sec=10)
        with self.assertRaises(ValueError):  # negative sections
            Grid.create_geometric_end(start=0, end=10, fact=1, sec=-10)

    def test_manipulation(self):
        g = Grid(start=4, end=5)
        if not g._is_seamless():
            raise RuntimeError("found seams!")
        self.assertEquals(g.start, 4)
        self.assertEquals(g.end, 5)

        # adding:
        g.add_left(1)
        self.assertEquals(g.start, 3)
        self.assertEquals(g.end, 5)
        g.add_right(2)
        self.assertEquals(g.start, 3)
        self.assertEquals(g.end, 7)
        if not g._is_seamless():
            raise RuntimeError("found seams!")
        with self.assertRaises(ValueError):
            g.add_left(size=-1)  # negative size
        with self.assertRaises(ValueError):
            g.add_right(size=-1)  # negative size
        with self.assertRaises(ValueError):
            g.add_left(size=1, particles=-1)  # negative number of particles
        with self.assertRaises(ValueError):
            g.add_right(size=1, particles=-1)  # negative number of particles

        # removing:
        g.remove_left()
        self.assertEquals(g.start, 4)
        self.assertEquals(g.end, 7)
        g.remove_right()
        self.assertEquals(g.start, 4)
        self.assertEquals(g.end, 5)
        self.assertEquals(len(g), 1)
        if not g._is_seamless():
            raise RuntimeError("found seams!")
        with self.assertRaises(RuntimeError):
            g.remove_left()  # last section
        with self.assertRaises(RuntimeError):
            g.remove_right()  # last section

    def test_moment(self):
        u = Grid.create_uniform(start=0, end=1, sec=10)  # default zero
        if not u._is_seamless():
            raise RuntimeError("found seams!")
        for order in range(10):
            self.assertEqual(u.moment(order=order), 0)
        g = Grid.create_geometric_end(start=0, end=1, sec=10, fact=1.1)  # default zero
        if not g._is_seamless():
            raise RuntimeError("found seams!")
        for order in range(10):
            self.assertEqual(g.moment(order=order), 0)
        with self.assertRaises(ValueError):
            u.moment(order=-1)  # oder < 0

    def test_lists(self):
        g = Grid.create_uniform(start=0, end=1, sec=10)  # default zero
        if not g._is_seamless():
            raise RuntimeError("found seams!")
        b = g.boundaries()
        self.assertEqual(b[0], g.start)
        self.assertEqual(b[-1], g.end)
        pi = g.pivots()
        self.assertEqual(len(pi) + 1, len(b))
        pa = g.particles()
        d = g.densities()
        self.assertEqual(len(pi), len(pa))
        self.assertEqual(len(pi), len(d))

    def test_seams(self):
        u = Grid.create_uniform(start=0, end=1, sec=10)
        u._sections[1].start += 0.01  # add a little offset to make a seam
        seams = u._find_seams()
        print(seams)
        if not seams:
            raise RuntimeError("no seams!")
