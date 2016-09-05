#!/usr/bin/env python3

"""
module implementing the core classes used for representing a discrete NDF.

These classes are the top level Grid class and the lower level Section class.
A Grid is a connected amount of sections. The sections are defined by their
lower and upper boundaries (minimum, maximum), the pivot point and the number
of particles within the section.
"""

from itertools import tee
from collections import deque
from scipy.integrate import quad
import matplotlib.pyplot as plt


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def zero(x):
    """Function always returning zero.

    This function acts mainly as a placeholder and default value for arbitrary
    function objects that have to be specified.

    :param x: function parameter.
    :return: always 0.
    """
    return 0


def find_initial_step(start, end, steps, factor, max_err=1e-14, max_iter=100):
    """Find initial step size for geometric grid with given parameters.

    The Function tries to figure out an initial step size for a geometric grid
    with the given paramters (start, end, steps, factor). Its initial guess is
    the step size of the corresponding uniform grid. Further guesses are
    calculated from the ratio of the resulting and the desired length. The
    method converges very fast and is therefore slowed down by averaging to
    approach the solution slowly and achieving a more stable result.

    :param start: start of the grid.
    :param end: end of the grid.
    :param steps: number of steps.
    :param factor: geometric grid factor.
    :param max_err: maximum allowed iteration error.
    :param max_iter: maximum allowed number of iterations.
    :return: initial step size (first cell) for the geometric grid.
    """
    assert start < end
    assert steps > 0
    assert factor >= 1
    assert max_err > 0
    assert max_iter > 0

    length = end - start  # total length of the range
    step = length / steps  # uniform step size for number of steps

    initial_step = step  # initial guess for the initial step size
    previous_step = step
    diff = max_err + 1  # difference between current iteration and total length
    for _ in range(max_iter):
        # sum individual step sizes to get total length of the geometric grid:
        current_step = initial_step
        current_length = 0
        for index in range(steps):
            current_length += current_step
            current_step *= factor
        # check difference and adjust initial step size:
        diff = abs(length - current_length)
        if diff < max_err:
            break
        else:
            previous_step = initial_step
            initial_step *= length / current_length
            # slow down converging by averaging with previous guess:
            mean_step = (initial_step + previous_step) / 2
            phi = 0.75  # weight combining mean and current best guess
            initial_step = phi * initial_step + (1-phi) * mean_step

    return initial_step


class Section:
    """Section class representing a discrete bin.

    It represents the semi open interval [start, end). (minimum included,
    maximum NOT included in the interval). The class includes a fixed pivot
    within the sections boundaries (min <= pivot < max). The default pivot
    position is right in the middle between the min and max boundaries
    (pivot = min + size / 2). The length of the section (max - min) is called
    size. Additionally the class contains a field for the total number of
    particles in the section and the derived value of the mean particle
    density.
    """
    def __init__(self, start, end, particles=0.):
        """Initializer.

        :param start: inclusive minimum of the section (v_i).
        :param end: exclusive maximum of the section (v_i+1).
        :param particles: total number of particles in the section (N_i).
        :raises: AssertionError if minimum >= maximum.
        :raises: AssertionError if particles < 0.
        """
        assert particles >= 0
        assert start < end
        self._particles = particles
        self._start = start
        self._end = end
        self._pivot = start + (end - start) / 2

    def __str__(self):
        """String representation.

        :return: human readable string representation.
        """
        return "<{cls} @ {id}: min={min:.3e}, max={max:.3e}, " \
               "size={size:.3e}, piv={piv:.3e}, part={part:.3e}>".format(
            cls=self.__class__.__name__,
            id=hex(id(self)),
            min=self.start,
            max=self.end,
            size=self.size,
            piv=self.pivot,
            part=self.particles
        )

    @property
    def start(self):
        """Getter for minimum property (v_i).

        :return: inclusive minimum of the section.
        """
        return self._start

    @start.setter
    def start(self, new_start):
        """Setter for minimum property (v_i).

        :param new_start: new inclusive minimum of the section.
        :raises: AssertionError if minimum >= maximum.
        """
        assert new_start < self.end
        self._start = new_start
        self._center_pivot()

    @property
    def end(self):
        """Getter for maximum property (v_i+1).

        :return: exclusive maximum of the section.
        """
        return self._end

    @end.setter
    def end(self, new_end):
        """Setter for maximum property (v_i+1).

        :param new_end: new exclusive maximum of the section.
        :raises: AssertionError if minimum >= maximum.
        """
        assert new_end > self.start
        self._end = new_end
        self._center_pivot()

    @property
    def range(self):
        """Getter for range property (v_i, v_i+1).

        :return: tuple containing minimum and maximum value of the section.
        """
        return self.start, self.end

    @range.setter
    def range(self, new_range):
        """Setter for range property (v_i, v_i+1).

        :param new_range: tuple with new minimum and maximum values.
        :raises: IndexError if len(new_range) < 2.
        :raises: TypeError [] operator is not supported.
        :raises: AssertionError if new_minimum >= new_maximum.
        """
        try:
            new_minimum = new_range[0]
            new_maximum = new_range[1]
        except IndexError:
            raise IndexError(
                "range '{}' must be a tuple of at least two elements!".format(
                    new_range
                )
            )
        except TypeError:
            raise TypeError(
                "range '{}' must be a tuple of at least two elements!".format(
                    new_range
                )
            )
        assert new_minimum < new_maximum
        if new_maximum >= self.end:  # expand or shift to the right first
            self.end = new_maximum
            self.start = new_minimum
        else:  # expand or shift to the left first
            self.start = new_minimum
            self.end = new_maximum
        self._center_pivot()

    @property
    def pivot(self):
        """Getter for pivot property (x_i).

        :return: sections pivot point.
        """
        return self._pivot

    @property
    def size(self):
        """Getter for size property (v_i+1 - v_i).

        :return: section size (max - min).
        """
        return self.end - self.start

    @property
    def particles(self):
        """Getter for total number of particles property (N_i).

        :return: total number of particles in the section.
        """
        return self._particles

    @particles.setter
    def particles(self, new_particles):
        """Setter for total number of particles property (N_i).

        :param new_particles: new number of particles in the section.
        :raises: AssertionError if new_particles < 0
        """
        assert new_particles >= 0
        self._particles = new_particles

    @property
    def particle_density(self):
        """Getter for the mean particle density property (n_i).

        :return: mean particle density in the section.
        """
        return self.particles / self.size

    def _center_pivot(self):
        """Center the pivot point between minimum and maximum.
        """
        self._pivot = self.start + self.size / 2


class Grid:
    """Grid class for representing an arbitrary grid of connected sections.

    The grid bins are represented by a list (actually a double ended queue) of
    sections with arbitrary ranges. The seamless connection of the buckets
    (max(bucket i) = min(bucket(i+1)) is guaranteed by the simple
    implementation of the classes methods. If changes are made to the Sections
    "by hand" the seamlessness has to be checked.
    """
    def __init__(self, start=0, end=1, particles=0.):
        """Initializer.

        :param start: initial section minimum (default: 0).
        :param end: initial section maximum (default: 1).
        """
        self._sections = deque()  # Double Ended QUEue
        initial_section = Section(start, end, particles)
        self._sections.append(initial_section)

    def __len__(self):
        """Length interface.

        :return: number of sections in the grid.
        """
        return len(self._sections)

    def __iter__(self):
        """Iterator interface.

        :return: iterator of section "list".
        """
        return iter(self._sections)

    def __str__(self):
        """String representation.

        :return: human readable string representation.
        """
        minimum = self._sections[0].start  # first section minimum
        maximum = self._sections[-1].end  # last section maximum
        number = len(self)
        return "<{cls} @{id}: min={min:.3e}, max={max:.3e}, num={num}>".format(
            cls=self.__class__.__name__,
            id=hex(id(self)),
            min=minimum,
            max=maximum,
            num=number
        )

    @property
    def start(self):
        """Getter for grid minimum property.

        :return: minimum of the grid.
        """
        return self._sections[0].start  # first section minimum

    @property
    def end(self):
        """Getter for grid maximum property.

        :return: maximum of the grid.
        """
        return self._sections[-1].end  # last section maximum

    def section(self, index):
        assert index >= 0
        assert index < len(self)
        return self._sections[index]

    def add_left(self, size, particles=0.):
        """Add new section to the left.

        :param size: range of the new section (maximum - minimum).
        :param particles: total number of particles in the new section.
        :raises: AssertionError if size <= 0.
        """
        assert size > 0
        current_left_section = self._sections[0]
        new_max = current_left_section.start
        new_min = new_max - size
        new_left_section = Section(new_min, new_max, particles)
        self._sections.appendleft(new_left_section)

    def add_right(self, size, particles=0.):
        """Add new section to the right.

        :param size: range of the new section (maximum - minimum).
        :param particles: total number of particles in the new section.
        :raises: AssertionError if size <= 0.
        """
        assert size > 0
        current_right_section = self._sections[-1]
        new_min = current_right_section.end
        new_max = new_min + size
        new_right_section = Section(new_min, new_max, particles)
        self._sections.append(new_right_section)

    def remove_left(self):
        """Remove section on the left.

        The last section can not be removed.

        :return: removed section.
        :raises: IndexError if length of the gris is 1.
        """
        if len(self) == 1:
            raise IndexError("Can't remove last section!")
        else:
            return self._sections.popleft()

    def remove_right(self):
        """Remove section on the right.

        The last section can not be removed.

        :return: removed section.
        :raises: IndexError if length of the gris is 1.
        """
        if len(self) == 1:
            raise IndexError("Can't remove last section!")
        else:
            return self._sections.pop()

    # TODO: look at grid manipulation after implemented simple method!
    def coarsen_range(self, start, end):
        """Coarsen grid by combining sections in the index range.

        :param start: inclusive start of index range.
        :param end: inclusive end of index range.
        :raises: AssertionError if start > end
        """
        assert start <= end
        assert start >= 0  # --> IndexError
        assert end < len(self)  # --> IndexError
        new_sections = deque()
        coarse_min = self._sections[start].start
        coarse_max = self._sections[end].end
        particles = 0
        for index, section in enumerate(self._sections):
            if index >= start and index <= end:
                particles += section.particles
        coarsened = False
        for index, section in enumerate(self._sections):
            if index < start or index > end:  # outside range
                new_sections.append(self._sections[index])
            else:  # inside coarsening range
                if not coarsened:
                    coarse_section = Section(start=coarse_min,
                                             end=coarse_max,
                                             particles=particles)
                    new_sections.append(coarse_section)
                    coarsened = True
        self._sections = new_sections

    def coarsen(self, times=1):
        for _ in range(times):
            index = 0
            while index < len(self) - 1:
                self.coarsen_range(index, index + 1)
                index += 1

    def refine_range(self, start, end):
        """Refine the grid by splitting the sections in the index range.

        :param start: inclusive start of index range.
        :param end: inclusive end of index range.
        :raises: AssertionError if start > end
        """
        assert start <= end
        assert start >= 0  # --> IndexError
        assert end < len(self)  # --> IndexError
        new_sections = deque()
        for index, section in enumerate(self._sections):
            if index < start or index > end:  # outside range
                new_sections.append(self._sections[index])
            else:  # inside refining range
                left = section.start
                right = section.end
                center = left + (right - left)/2
                particles = section.particles
                left_section = Section(start=left, end=center,
                                       particles=particles/2)
                right_section = Section(start=center, end=right,
                                        particles=particles/2)
                new_sections.append(left_section)
                new_sections.append(right_section)
        self._sections = new_sections

    def refine(self, times=1):
        for _ in range(times):
            self.refine_range(0, len(self) - 1)

    def boundaries(self):
        """return the list of boundary values.

        :return: boundary list.
        """
        boundary_list = []
        for section in self:
            boundary_list.append(section.start)
        boundary_list.append(self._sections[-1].end)  # last section max
        return boundary_list

    def pivots(self):
        """return the list of pivots.

        :return: pivot list.
        """
        pivot_list = []
        for section in self:
            pivot_list.append(section.pivot)
        return pivot_list

    def particles(self):
        """return the list of particle numbers.

        :return: particle number list.
        """
        particle_list = []
        for section in self:
            particle_list.append(section.particles)
        return particle_list

    def particle_densities(self):
        """return the list of particle density values.

        :return: particle density list.
        """
        particle_density_list = []
        for section in self:
            particle_density_list.append(section.particle_density)
        return particle_density_list

    @staticmethod
    def create_uniform(start, end, steps, func=zero, correct=True):
        """Create a uniform grid.

        :param start: inclusive minimum of the grid.
        :param end: exclusive maximum of the grid.
        :param steps: total number of sections in the grid.
        :param func: number density function for calculating particle numbers.
        :param correct: flag for fixing numeric deviations in the last section.
        :return: equidistant Grid object.
        """
        assert start <= end
        assert steps > 0
        return Grid.create_geometric_step(start, end, steps,
                                          factor=1, func=func, correct=correct)

    @staticmethod
    def create_geometric_step(start, end, steps, factor, func=zero,
                              correct=True):
        """Create a geometric grid.

        The grid is created by using the same initial step size as the
        corresponding uniform grid (step = (maximum - minimum) / amount). This
        method does NOT achieve the same maximum boundary.

        :param start: inclusive minimum of the grid.
        :param end: exclusive maximum of the corresponding uniform grid.
        :param steps: number of sections in the corresponding uniform grid.
        :param factor: factor for size change.
        :param func: particle number density function.
        :param correct: flag for fixing numeric deviations in the last section.
        :return: geometric Grid object.
        """
        assert start <= end
        assert steps > 0
        assert factor >= 1
        # use uniform step size as initial step size:
        initial_step = (end - start) / steps
        return Grid.create_geometric(start, end, initial_step, factor, func,
                                     correct)

    @staticmethod
    def create_geometric_end(start, end, steps, factor, func=zero,
                             correct=True):
        """Create a geometric grid.

        The grid is created by choosing an initial step size so the maximum
        boundary is the same as for the corresponding uniform grid.

        :param start: inclusive minimum of the grid.
        :param end: exclusive maximum of the grid.
        :param steps: number of sections in the corresponding uniform grid.
        :param factor: factor for size change.
        :param func: particle number density function.
        :param correct: flag for fixing numeric deviations in the last section.
        :return: geometric Grid object.
        """
        assert start <= end
        assert steps > 0
        assert factor >= 1
        # calculate initial step size for given parameters:
        initial_step = find_initial_step(start, end, steps, factor)
        return Grid.create_geometric(start, end, initial_step, factor, func,
                                     correct)

    @staticmethod
    def create_geometric(start, end, initial_step, factor, func=zero,
                         correct=True):
        """Create a geometric grid.

        The grid is created using the given initial step size. By default the
        last section is corrected to match the given end value. This results in
        the last section not matching the exact end value. If this correction
        is turned off, the actual maximum of the grid will be greater than the
        given end value.

        :param start: inclusive minimum of the grid.
        :param end: exclusive maximum of the grid.
        :param initial_step: step size of the first step in the grid.
        :param factor: factor controlling size change.
        :param func: particle number density function.
        :param correct: flag for fixing numeric deviations in the last section.
        :return: geometric Grid object.
        """
        assert start <= end
        assert initial_step > 0
        assert factor >= 1
        # create initial grid:
        initial_max = start + initial_step
        particles = quad(func, start, initial_max)[0]
        grid = Grid(start, initial_max, particles)
        # add sections to the grid:
        current_size = initial_step
        current_max = initial_max
        while current_max < end:
            current_size *= factor
            current_max += current_size
            last_section = grid._sections[-1]
            lower = last_section.end
            upper = last_section.end + current_size
            particles = quad(func, lower, upper)[0]
            grid.add_right(size=current_size, particles=particles)
        # clean up last section(s):
        if correct:
            last_section = grid._sections[-1]
            last_min = last_section.start
            last_max = last_section.end
            left = end - last_min  # left side of the last section
            right = last_max - end  # rigth side of the last section
            if left < right:  # last section pushing outward
                grid.remove_right()
                last_section = grid._sections[-1]
            last_section.end = end
            lower = last_section.end
            upper = last_section.end + current_size
            particles = quad(func, lower, upper)[0]
            last_section.particles = particles
        return grid

    def _find_seams(self):
        """Find seams between sections by checking every pair in the grid.

        :return: list of section pairs that are not properly connected.
        """
        seams = []
        for index, section in enumerate(self):
            print(index, section)
            if index < len(self) - 1:
                left_max = self._sections[index].end
                right_min = self._sections[index+1].start
                if left_max != right_min:
                    seams.append((index, index+1))
        return seams

    def _info(self):
        """print info about the class to console.
        """
        print(self)
        for index, sections in enumerate(self):
            print(index, sections)

    def _plot(self, scale="linear"):
        """Plot grids total particle number values using matplotlib.
        """
        pivots = self.pivots()
        particles = self.particles()
        densities = self.particle_densities()
        fig, ax1 = plt.subplots()
        ax1.plot(pivots, particles, "ro", label="numbers")
        ax1.set_ylabel("particles")
        ax2 = ax1.twinx()
        ax2.plot(pivots, densities, "r-", label="density")
        ax2.set_ylabel("density")
        ax1.set_xlabel("size")
        ax1.set_xscale(scale)
        ax1.set_yscale(scale)
        ax2.set_yscale(scale)
        ax1.grid()
        plt.show()

    # FUTURE IDEAS:
    @staticmethod
    def from_boundaries(boundaries):
        """Create grid from list of range boundaries.
        """
        raise NotImplementedError  # TODO: maybe implement in the future

    @staticmethod
    def create_adaptive(start, end, steps, func):
        """Create a grid that adapts automatically to the given function.
        """
        raise NotImplementedError  # TODO: maybe implement in the future

    def smooth(self, start, end):
        """Smooth grid by averaging the section sizes in the range
        """
        raise NotImplementedError  # TODO: maybe implement in the future
