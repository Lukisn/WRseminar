#!/usr/bin/env python3

"""
module implementing the core classes used for representing a discrete NDF.

These classes are the top level Grid class and the lower level Section class.
A Grid is a connected amount of sections. The sections are defined by their
lower and upper boundaries (minimum, maximum), the pivot point and the number
of particles within the section.
"""

from collections import deque
import scipy.integrate as spint
import matplotlib.pyplot as plt


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

    It represents the semi open interval [min, max). (minimum included,
    maximum NOT included in the interval). The class includes a fixed pivot
    within the sections boundaries (min <= pivot < max). The default pivot
    position is right in the middle between the min and max boundaries
    (pivot = min + size / 2). The length of the section (max - min) is called
    size. Additionally the class contains a field for the total number of
    particles in the section and the derived value of the mean particle
    density.
    """
    def __init__(self, minimum, maximum, particles=0.):
        """Initializer.

        :param minimum: inclusive minimum of the section (v_i).
        :param maximum: exclusive maximum of the section (v_i+1).
        :param particles: total number of particles in the section (N_i).
        :raises: ValueError if minimum >= maximum.
        :raises: ValueError if particles < 0.
        """
        if particles < 0:
            raise ValueError("particles '{}' must be >= 0!".format(particles))
        elif minimum >= maximum:
            raise ValueError(
                "minimum '{min}' must be < maximum '{max}'".format(
                    min=minimum,
                    max=maximum
                )
            )
        else:
            self._particles = particles
            self._minimum = minimum
            self._maximum = maximum
            self._pivot = minimum + (maximum - minimum) / 2

    def __str__(self):
        """String representation.

        :return: human readable string representation.
        """
        return "<{cls} @ {id}: min={min:.3e}, max={max:.3e}, " \
               "size={size:.3e}, piv={piv:.3e}, part={part:.3e}>".format(
            cls=self.__class__.__name__,
            id=hex(id(self)),
            min=self.minimum,
            max=self.maximum,
            size=self.size,
            piv=self.pivot,
            part=self.particles
        )

    @property
    def minimum(self):
        """Getter for minimum property (v_i).

        :return: inclusive minimum of the section.
        """
        return self._minimum

    @minimum.setter
    def minimum(self, new_min):
        """Setter for minimum property (v_i).

        :param new_min: new inclusive minimum of the section.
        :raises: ValueError if minimum >= maximum.
        """
        if new_min >= self.maximum:
            raise ValueError(
                "minimum '{min}' must be < maximum '{max}'".format(
                    min=self.minimum,
                    max=self.maximum
                )
            )
        else:
            self._minimum = new_min
            self._center_pivot()

    @property
    def maximum(self):
        """Getter for maximum property (v_i+1).

        :return: exclusive maximum of the section.
        """
        return self._maximum

    @maximum.setter
    def maximum(self, new_max):
        """Setter for maximum property (v_i+1).

        :param new_max: new exclusive maximum of the section.
        :raises: ValueError if minimum >= maximum.
        """
        if new_max <= self.minimum:
            raise ValueError(
                "maximum '{max}' must be > minimum '{min}'".format(
                    max=self.maximum,
                    min=self.minimum
                )
            )
        else:
            self._maximum = new_max
            self._center_pivot()

    @property
    def range(self):
        """Getter for range property (v_i, v_i+1).

        :return: tuple containing minimum and maximum value of the section.
        """
        return self.minimum, self.maximum

    @range.setter
    def range(self, new_range):
        """Setter for range property (v_i, v_i+1).

        :param new_range: tuple with new minimum and maximum values.
        :raises: IndexError if len(new_range) < 2.
        :raises: TypeError [] operator is not supported.
        :raises: ValueError if new_minimum >= new_maximum.
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
        if new_minimum >= new_maximum:
            raise ValueError(
                "maximum '{max}' must be greater than minimum '{min}'".format(
                    max=self.maximum,
                    min=self.minimum
                )
            )
        else:
            if new_maximum >= self.maximum:  # expand or shift to the right
                self.maximum = new_maximum
                self.minimum = new_minimum
            else:
                self.minimum = new_minimum
                self.maximum = new_maximum
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
        return self.maximum - self.minimum

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
        """
        if new_particles < 0:
            raise ValueError("number of particles '{}' must be >= 0!".format(
                    self.particles
            ))
        else:
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
        self._pivot = self.minimum + self.size / 2


class Grid:
    """Grid class for representing an arbitrary grid of connected sections.

    The grid bins are represented by a list (actually a double ended queue) of
    sections with arbitrary ranges. The seamless connection of the buckets
    (max(bucket i) = min(bucket(i+1)) is guaranteed by the simple
    implementation of the classes methods. If changes are made to the Sections
    "by hand" the seamlessness has to be checked.
    """
    def __init__(self, minimum=0, maximum=1, particles=0.):
        """Initializer.

        :param minimum: initial section minimum (default: 0).
        :param maximum: initial section maximum (default: 1).
        """
        self._sections = deque()  # Double Ended QUEue
        initial_section = Section(minimum, maximum, particles)
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
        minimum = self._sections[0].minimum  # first section minimum
        maximum = self._sections[-1].maximum  # last section maximum
        number = len(self)
        return "<{cls} @{id}: min={min:.3e}, max={max:.3e}, num={num}>".format(
            cls=self.__class__.__name__,
            id=hex(id(self)),
            min=minimum,
            max=maximum,
            num=number
        )

    @property
    def minimum(self):
        """Getter for grid minimum property.

        :return: minimum of the grid.
        """
        return self._sections[0].minimum  # first section minimum

    @property
    def maximum(self):
        """Getter for grid maximum property.

        :return: maximum of the grid.
        """
        return self._sections[-1].maximum  # last section maximum

    def add_left(self, size, particles=0.):
        """Add new section to the left.

        :param size: range of the new section (maximum - minimum).
        :param particles: total number of particles in the new section.
        :raises: ValueError if size <= 0.
        """
        if size <= 0:
            raise ValueError("size must be > 0!")
        else:
            current_left_section = self._sections[0]
            new_max = current_left_section.minimum
            new_min = new_max - size
            new_left_section = Section(minimum=new_min, maximum=new_max,
                                       particles=particles)
            self._sections.appendleft(new_left_section)

    def add_right(self, size, particles=0.):
        """Add new section to the right.

        :param size: range of the new section (maximum - minimum).
        :param particles: total number of particles in the new section.
        :raises: ValueError if size <= 0.
        """
        if size <= 0:
            raise ValueError("size must be > 0!")
        else:
            current_right_section = self._sections[-1]
            new_min = current_right_section.maximum
            new_max = new_min + size
            new_right_section = Section(minimum=new_min, maximum=new_max,
                                        particles=particles)
            self._sections.append(new_right_section)

    def remove_left(self):
        """Remove section on the left.

        :return: removed section.
        :raises: ValueError if len
        """
        if len(self) == 1:
            raise IndexError("Can not remove last section!")
        else:
            return self._sections.popleft()

    def remove_right(self):
        """Remove section on the right.

        :return: removed section.
        """
        if len(self) == 1:
            raise IndexError("Can not remove last section!")
        else:
            return self._sections.pop()

    def coarsen(self, start, end):
        """Coarsen grid by combining sections in the index range.

        :param start: inclusive start of index range.
        :param end: inclusive end of index range.
        """
        if start <= end:
            if start >= 0 and end < len(self):
                new_sections = deque()
                coarse_min = self._sections[start].minimum
                coarse_max = self._sections[end].maximum
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
                            coarse_section = Section(minimum=coarse_min,
                                                    maximum=coarse_max,
                                                    particles=particles)
                            new_sections.append(coarse_section)
                            coarsened = True
                self._sections = new_sections
            else:
                raise ValueError(
                    "start index '{start}' or end index '{end}' "
                    "out of range ({min}, {max})".format(
                        start=start,
                        end=end,
                        min=0,
                        max=len(self)-1
                    )
                )
        else:
            raise ValueError(
                "start index '{start}' must be <= end index '{end}'!".format(
                    start=start,
                    end=end
                )
            )

    def refine(self, start, end):
        """Refine the grid by splitting the sections in the index range.

        :param start: inclusive start of index range.
        :param end: inclusive end of index range.
        """
        if start <= end:
            if start >= 0 and end <= len(self):
                new_sections = deque()
                for index, section in enumerate(self._sections):
                    if index < start or index > end:  # outside range
                        new_sections.append(self._sections[index])
                    else:  # inside refining range
                        left = section.minimum
                        right = section.maximum
                        center = left + (right - left)/2
                        particles = section.particles
                        left_section = Section(minimum=left, maximum=center,
                                               particles=particles/2)
                        right_section = Section(minimum=center, maximum=right,
                                                particles=particles/2)
                        new_sections.append(left_section)
                        new_sections.append(right_section)
                self._sections = new_sections
            else:
                raise ValueError(
                    "start index '{start}' or end index '{end}' "
                    "out of range ({min}, {max})!".format(
                        start=start,
                        end=end,
                        min=0,
                        max=len(self)-1
                    )
                )
        else:
            raise ValueError(
                "start index '{start}' must be <= end index '{end}'!".format(
                    start=start,
                    end=end
                )
            )

    def _find_seams(self):
        """Find seams between sections by checking every pair in the grid.

        :return: list of section pairs that are not properly connected.
        """
        seams = []
        for index, section in enumerate(self):
            print(index, section)
            if index < len(self) - 1:
                left_max = self._sections[index].maximum
                right_min = self._sections[index+1].minimum
                if left_max != right_min:
                    seams.append((index, index+1))
        return seams

    def info(self):
        """print info about the class to console.
        """
        print(self)
        for index, sections in enumerate(self):
            print(index, sections)

    def plot(self):
        """Plot grids total particle number values using matplotlib.
        """
        plt.plot(self.pivots(), self.particles(), "ro")
        plt.show()

    def boundaries(self):
        """return the list of boundary values.

        :return: boundary list.
        """
        boundary_list = []
        for section in self:
            boundary_list.append(section.minimum)
        boundary_list.append(self._sections[-1].maximum)  # last section max
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

    # TODO: simplify grid creation!
    @staticmethod
    def create_uniform(minimum, maximum, amount, func=zero):
        """Create a uniform grid.

        :param minimum: inclusive minimum of the grid.
        :param maximum: exclusive maximum of the grid.
        :param amount: total number of sections in the grid.
        :param func: number density function for calculating particle numbers.
        :return: equidistant Grid object.
        """
        size = (maximum - minimum) / amount
        return Grid.create_geometric_step(minimum=minimum, maximum=maximum,
                                          factor=1, amount=amount, func=func)

    @staticmethod
    def create_geometric_step(minimum, maximum, amount, factor, func=zero):
        """Create a geometric grid.

        The grid is created by using the same initial step size as the
        corresponding uniform grid (step = (maximum - minimum) / amount). This
        method does NOT achieve the same maximum boundary.

        :param minimum: inclusive minimum of the grid.
        :param maximum: exclusive maximum of the corresponding uniform grid.
        :param amount: number of sections in the corresponding uniform grid.
        :param factor: factor for size change.
        :return: geometric Grid object.
        """
        if minimum >= maximum:
            raise ValueError(
                "minimum '{min}' must be < maximum '{max}'".format(
                    min=minimum,
                    max=maximum
                )
            )
        elif amount <= 0:
            raise ValueError("amount '{a}' must be > 0!".format(a=amount))
        elif factor < 1:
            raise ValueError("factor '{f}' must be >= 1!".format(f=factor))
        else:
            # create initial grid:
            size = (maximum - minimum) / amount
            first_max = minimum + size
            particles = spint.quad(func, minimum, first_max)[0]
            grid = Grid(minimum, first_max, particles)
            # add sections to the grid:
            current_size = size
            current_max = first_max
            while current_max < maximum:
                current_size *= factor
                current_max += current_size
                last_section = grid._sections[-1]
                lower = last_section.maximum
                upper = last_section.maximum + current_size
                particles = spint.quad(func, lower, upper)[0]
                grid.add_right(size=current_size, particles=particles)
            # clean up last section(s):
            last_section = grid._sections[-1]
            last_min = last_section.minimum
            last_max = last_section.maximum
            left = maximum - last_min
            right = last_max - maximum
            if left > right:  # shift last max left
                last_section.maximum = maximum
                lower = last_section.maximum
                upper = last_section.maximum + current_size
                particles = spint.quad(func, lower, upper)[0]
                last_section.particles = particles
            else:  # delete last section and shift last max right
                grid.remove_right()
                last_section = grid._sections[-1]
                last_section.maximum = maximum
                lower = last_section.maximum
                upper = last_section.maximum + current_size
                particles = spint.quad(func, lower, upper)[0]
                last_section.particles = particles
            return grid

    @staticmethod
    def create_geometric_maximum(minimum, maximum, amount, factor, func=zero):
        """Create a geometric grid.

        The grid is created by choosing an initial step size so the maximum
        boundary is the same as for the corresponding uniform grid.

        :param minimum: inclusive minimum of the grid.
        :param maximum: exclusive maximum of the grid.
        :param amount: number of sections in the corresponding uniform grid.
        :param factor: factor for size change.
        :return: geometric Grid object.
        """
        if minimum >= maximum:
            raise ValueError(
                "minimum '{min}' must be < maximum '{max}'".format(
                    min=minimum,
                    max=maximum
                )
            )
        elif amount <= 0:
            raise ValueError("amount '{a}' must be > 0!".format(a=amount))
        elif factor < 1:
            raise ValueError("factor '{f}' must be >= 1!".format(f=factor))
        else:
            # find initial step size:
            size = find_initial_step(minimum, maximum, amount, factor)
            # create initial grid:
            first_max = minimum + size
            particles = spint.quad(func, minimum, first_max)[0]
            grid = Grid(minimum, first_max, particles)
            # add sections to the grid:
            current_size = size
            current_max = first_max
            #while current_max < maximum:
            for index in range(amount-1):
                current_size *= factor
                current_max += current_size
                last_section = grid._sections[-1]
                lower = last_section.maximum
                upper = last_section.maximum + current_size
                particles = spint.quad(func, lower, upper)[0]
                grid.add_right(size=current_size, particles=particles)
            # clean up last section(s):
            last_section = grid._sections[-1]
            last_min = last_section.minimum
            last_max = last_section.maximum
            left = maximum - last_min
            right = last_max - maximum
            if left > right:  # shift last max left
                last_section.maximum = maximum
                lower = last_section.maximum
                upper = last_section.maximum + current_size
                particles = spint.quad(func, lower, upper)[0]
                last_section.particles = particles
            else:  # delete last section and shift last max right
                grid.remove_right()
                last_section = grid._sections[-1]
                last_section.maximum = maximum
                lower = last_section.maximum
                upper = last_section.maximum + current_size
                particles = spint.quad(func, lower, upper)[0]
                last_section.particles = particles
            return grid

    @staticmethod
    def from_boundaries(boundaries):
        """Create grid from list of range boundaries.
        """
        raise NotImplementedError  # TODO: maybe implement in the future

    def smooth(self, start, end):
        """Smooth grid by averaging the section sizes in the range
        """
        raise NotImplementedError  # TODO: maybe implement in the future
