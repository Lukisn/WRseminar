#!/usr/bin/env python3
"""
Module implementing the core classes used for representing a discrete NDF.

These classes are the top level ``Grid`` class and the lower level ``Section``
class. A Grid is a connected amount of sections. The sections are defined by
their lower and upper boundaries (minimum, maximum), the pivot point and the
number of particles within the section.
"""
# standard library imports:
from collections import deque
# third party imports:
import matplotlib.pyplot as plt
from scipy.integrate import quad
# application imports:
from sectional.functions import zero


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
    # check parameters:
    if start >= end:
        raise ValueError("start '{}' >= end '{}'!".format(start, end))
    if steps <= 0:
        raise ValueError("steps '{}' <= 0!".format(steps))
    if factor < 1:
        raise ValueError("factor '{}' < 1!".format(factor))
    if max_err <= 0:
        raise ValueError("max_err '{}' <= 0!".format(max_err))
    if max_iter <= 0:
        raise ValueError("max_iter '{}' <= 0!".format(max_iter))

    length = end - start  # total length of the range
    step = length / steps  # uniform step size for number of steps

    initial_step = step  # initial guess for the initial step size
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

    .. code::

         start   pivot    end
           |       |       |
        ---[-------o-------)-->
           |       |       |
          v_i     x_i    v_i+1

    """

    def __init__(self, start, end, particles=0):
        """Initializer.

        :param start: inclusive minimum of the section (v_i).
        :param end: exclusive maximum of the section (v_i+1).
        :param particles: total number of particles in the section (N_i).
        :raises AssertionError: if minimum >= maximum.
        :raises AssertionError: if particles < 0.
        """
        # check parameters:
        if particles < 0:
            raise ValueError("particles '{}' < 0!".format(particles))
        if start >= end:
            raise ValueError("start '{}' >= end '{}'!".format(start, end))
        #assert particles >= 0
        #assert start < end

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
        # check parameters:
        if new_start >= self.end:
            raise ValueError("start '{}' >= end '{}'!".format(new_start, self.end))
        #assert new_start < self.end

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
        # check parameters:
        if new_end <= self.start:
            raise ValueError("end '{}' <= start '{}'".format(new_end, self.start))
        #assert new_end > self.start

        self._end = new_end
        self._center_pivot()

    @property
    def interval(self):
        """Getter for interval property (v_i, v_i+1).

        :return: tuple containing minimum and maximum value of the section.
        """
        return self.start, self.end

    @interval.setter
    def interval(self, new_interval):
        """Setter for interval property (v_i, v_i+1).

        :param new_interval: tuple with new minimum and maximum values.
        :raises: IndexError if len(new_range) < 2.
        :raises: TypeError [] operator is not supported.
        :raises: AssertionError if new_minimum >= new_maximum.
        """
        # check parameters:
        try:
            new_start = new_interval[0]
            new_end = new_interval[1]
        except IndexError:
            raise IndexError(
                "range '{}' must be a tuple of at least two elements!".format(
                    new_interval
                )
            )
        except TypeError:
            raise TypeError(
                "range '{}' must be a tuple of at least two elements!".format(
                    new_interval
                )
            )
        if new_start >= new_end:
            raise ValueError("start '{}' >= end '{}'!".format(new_start, new_end))
        #assert new_start < new_end

        if new_end >= self.end:  # expand or shift to the right first
            self.end = new_end
            self.start = new_start
        else:  # expand or shift to the left first
            self.start = new_start
            self.end = new_end
        self._center_pivot()

    @property
    def pivot(self):
        """Getter for pivot property (x_i).

        :return: sections pivot point.
        """
        return self._pivot

    def _center_pivot(self):
        """Center the pivot point between minimum and maximum.
        """
        self._pivot = self.start + self.size / 2

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
        # check parameters:
        if new_particles < 0:
            raise ValueError("particles '{}' < 0!".format(new_particles))
        #assert new_particles >= 0

        self._particles = new_particles

    @property
    def density(self):
        """Getter for the mean particle density property (n_i).

        :return: mean particle density in the section.
        """
        return self.particles / self.size


class Grid:
    """Grid class for representing an arbitrary grid of connected sections.

    The grid bins are represented by a list (actually a double ended queue) of
    N sections with arbitrary interval ranges. The seamless connection of the
    sections (max(bucket i) = min(bucket(i+1)) is guaranteed by the simple
    implementation of the classes methods. If changes are made to the Sections
    "by hand" the seamlessness has to be checked.

    .. code::

        start                                          end
          0     1     2       3  ...   N-1    N        N+1   <- boundaries
          |  0  |  1  |   2   | 3 | ... | N-1 |    N    |    <- pivots
          +--o--+--o--+---o---+-o-+--o--+--o--+----o----+->

          [--o--)     [---o---)   [--o--)     [----o----)    <- section objects
        left    [--o--)       [-o-)     [--o--)       right

    """
    def __init__(self, start=0, end=1, particles=0):
        """Initializer.

        :param start: initial section minimum (default: 0).
        :param end: initial section maximum (default: 1).
        :param particles: initial number of particles (default: 0).
        """
        # check parameters:
        if start >= end:
            raise ValueError("start '{}' > end '{}'!".format(start,end))
        if particles < 0:
            raise ValueError("particles '{}' < 0!".format(particles))

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
        start = self._sections[0].start  # first section minimum
        end = self._sections[-1].end  # last section maximum
        number = len(self)
        return "<{cls} @{id}: min={s:.3e}, max={e:.3e}, num={num}>".format(
            cls=self.__class__.__name__,
            id=hex(id(self)),
            s=start,
            e=end,
            num=number
        )

    def __getitem__(self, item):
        return self._sections[item]

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

    def add_left(self, size, particles=0):
        """Add new section to the left.

        :param size: range of the new section (maximum - minimum).
        :param particles: total number of particles in the new section.
        :raises: AssertionError if size <= 0.
        """
        # check parameters:
        if size <= 0:
            raise ValueError("size '{}' <= 0!".format(size))
        if particles < 0:
            raise ValueError("particles '{}' < 0!".format(particles))
        #assert size > 0

        current_left_section = self._sections[0]
        new_max = current_left_section.start
        new_min = new_max - size
        new_left_section = Section(new_min, new_max, particles)
        self._sections.appendleft(new_left_section)

    def add_right(self, size, particles=0):
        """Add new section to the right.

        :param size: range of the new section (maximum - minimum).
        :param particles: total number of particles in the new section.
        :raises: AssertionError if size <= 0.
        """
        # check parameters:
        if size <= 0:
            raise ValueError("size '{}' <= 0!".format(size))
        if particles < 0:
            raise ValueError("Particles '{}' < 0!".format(particles))

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
            raise RuntimeError("Can't remove last section!")
        else:
            return self._sections.popleft()

    def remove_right(self):
        """Remove section on the right.

        The last section can not be removed.

        :return: removed section.
        :raises: IndexError if length of the gris is 1.
        """
        if len(self) == 1:
            raise RuntimeError("Can't remove last section!")
        else:
            return self._sections.pop()

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

    def densities(self):
        """return the list of particle density values.

        :return: particle density list.
        """
        particle_density_list = []
        for section in self:
            particle_density_list.append(section.density)
        return particle_density_list

    def moment(self, order=0):
        """Calculate the moment of a given order of the entire discrete NDF.

        :param order: order of the moment to calculate.
        :return: moment of given order.
        """
        # check parameters:
        if order < 0:
            raise ValueError("order '{}' < 0!".format(order))

        moment = 0
        for section in self:
            x = section.pivot
            N = section.particles
            section_moment = (x ** order) * N
            moment += section_moment
        return moment

    @classmethod
    def create_uniform(cls, start, end, sec, func=zero, corr=True):
        """Create a uniform grid.

        This method is based on the create_geometric() method.

        :param start: inclusive minimum of the grid.
        :param end: exclusive maximum of the grid.
        :param sec: total number of sections in the grid.
        :param func: number density function for calculating particle numbers.
        :param corr: flag for fixing numeric deviations in the last section.
        :return: equidistant Grid object.
        """
        # check parameters:
        if start > end:
            raise ValueError("start '{}' > end '{}'!".format(start, end))
        if sec <= 0:
            raise ValueError("sections '{}' <= 0!".format(sec))

        # create grid:
        return cls.create_geometric_step(
            start, end, sec, fact=1, func=func, corr=corr
        )

    @classmethod
    def create_geometric(cls, start, end, ini_sec, fact, func=zero, corr=True):
        """Create a geometric grid.

        The grid is created using the given initial step size. By default the
        last section is corrected to match the given end value. This results in
        the last section not matching the exact end value. If this correction
        is turned off, the actual maximum of the grid will be greater than the
        given end value.

        :param start: inclusive minimum of the grid.
        :param end: exclusive maximum of the grid.
        :param ini_sec: step size of the first step in the grid.
        :param fact: factor controlling size change.
        :param func: particle number density function.
        :param corr: flag for fixing numeric deviations in the last section.
        :return: geometric Grid object.
        """
        # check parameters:
        if start >= end:
            raise ValueError("start '{}' >= end '{}'!".format(start, end))
        if ini_sec <= 0:
            raise ValueError("initial_step '{}' <= 0!".format(ini_sec))
        if fact < 1:
            raise ValueError("factor '{}' < 1!".format(fact))

        # create initial grid / first section:
        initial_max = start + ini_sec
        particles, *_ = quad(func, start, initial_max)
        grid = cls(start, initial_max, particles)

        # add sections to the grid:
        current_size = ini_sec
        current_max = initial_max
        while current_max < end:
            current_size *= fact
            current_max += current_size

            grid.add_right(size=current_size)
            current_section = grid[-1]
            lower = current_section.start
            upper = current_section.end
            particles, *_ = quad(func, lower, upper)
            current_section.particles = particles

        # clean up last section(s):
        if corr:
            last_section = grid[-1]
            last_min = last_section.start
            last_max = last_section.end
            inside = end - last_min  # left side of the last section
            outside = last_max - end  # right side of the last section
            if inside < outside:  # last section pushing outward
                grid.remove_right()
                last_section = grid[-1]
            last_section.end = end
            lower = last_section.end
            upper = last_section.end + current_size
            particles, *_ = quad(func, lower, upper)
            last_section.particles = particles

        return grid

    @classmethod
    def create_geometric_step(cls, start, end, uni_sec, fact, func=zero, corr=True):
        """Create a geometric grid.

        The grid is created by using the same initial step size as the
        corresponding uniform grid (step = (maximum - minimum) / amount). This
        method does NOT achieve the same maximum boundary.

        The method is based on the create_geometric() method.

        :param start: inclusive minimum of the grid.
        :param end: exclusive maximum of the corresponding uniform grid.
        :param uni_sec: number of sections in the corresponding uniform grid.
        :param fact: factor for size change.
        :param func: particle number density function.
        :param corr: flag for fixing numeric deviations in the last section.
        :return: geometric Grid object.
        """
        # check parameters:
        if start > end:
            raise ValueError("start '{}' > end '{}'!".format(start, end))
        if uni_sec <= 0:
            raise ValueError("sections '{}' <= 0!".format(uni_sec))
        if fact < 1:
            raise ValueError("factor '{}' < 1!".format(fact))

        # use uniform step size as initial step size and create grid:
        initial_step = (end - start) / uni_sec
        return cls.create_geometric(
            start, end, initial_step, fact, func, corr
        )

    @classmethod
    def create_geometric_end(cls, start, end, sec, fact, func=zero, corr=True):
        """Create a geometric grid.

        The grid is created by choosing an initial step size so the maximum
        boundary is the same as for the corresponding uniform grid.

        The method is based on the create_geometric() method.

        :param start: inclusive minimum of the grid.
        :param end: exclusive maximum of the grid.
        :param sec: number of sections in the corresponding uniform grid.
        :param fact: factor for size change.
        :param func: particle number density function.
        :param corr: flag for fixing numeric deviations in the last section.
        :return: geometric Grid object.
        """
        # check parameters:
        if start > end:
            raise ValueError("start '{}' > end '{}'!".format(start, end))
        if sec <= 0:
            raise ValueError("sections '{}' <= 0!".format(sec))
        if fact < 1:
            raise ValueError("factor '{}' < 1!".format(fact))

        # calculate initial step size for given parameters and create grid:
        initial_step = find_initial_step(start, end, sec, fact)
        return cls.create_geometric(
            start, end, initial_step, fact, func, corr
        )

    def to_file(self, filename, info=None):
        """Write Grid to text file.

        The resulting text file is structured like this:

        first line: info comment "# start, end, pivot, particles, density"
        other lines: single bin info in ascending order.

        :param filename: target file to write grid data to.
        :return: None.
        """
        fmt = "{st:.5e}, {end:.5e}, {piv:.5e}, {part:.5e}, {dens:.5e}\n"
        with open(filename, "w") as fh:
            if info is not None:
                fh.write(str(info))
            fh.write("# start, end, pivot, particles, density\n")
            for section in self:
                fh.write(fmt.format(
                    st=section.start,
                    end=section.end,
                    piv=section.pivot,
                    part=section.particles,
                    dens=section.density
                ))

    def _find_seams(self):
        """Find seams between sections by checking every pair in the grid.

        :return: list of section pairs that are not properly connected.
        """
        last_index = len(self) - 1
        seams = []
        for index, section in enumerate(self):
            if index < last_index:
                left_end = self._sections[index].end
                right_start = self._sections[index + 1].start
                if left_end != right_start:
                    seams.append((index, index + 1))
        return seams

    def _is_seamless(self):
        """Check if the grid is seamless.

        :return: True (no seams) / False (seams!).
        """
        seams = self._find_seams()
        if seams:
            return False
        else:
            return True

    def _info(self):
        """print info about the class to console.
        """
        print(self)
        for index, sections in enumerate(self):
            print(4 * " ", index, sections)

    def _plot_particles(self, xscale="log", yscale="log"):
        """Plot grids total particle number values using matplotlib.
        """
        # gather data:
        pivots = self.pivots()
        particles = self.particles()

        # plot data:
        fig, ax1 = plt.subplots()
        ax1.plot(pivots, particles, "ro", label="numbers")
        ax1.set_ylabel("particles")
        ax1.set_xlabel("size")
        ax1.set_xscale(xscale)
        ax1.set_yscale(yscale)
        ax1.grid()
        plt.show()

    def _plot_density(self, xscale="log", yscale="log"):
        """Plot grids particle densities values using matplotlib.
        """
        # gather data:
        pivots = self.pivots()
        densities = self.densities()

        # plot data:
        fig, ax1 = plt.subplots()
        ax1.plot(pivots, densities, "r-", label="density")
        ax1.set_ylabel("density")
        ax1.set_xlabel("size")
        ax1.set_xscale(xscale)
        ax1.set_yscale(yscale)
        ax1.grid()
        plt.show()

    def _plot(self, xscale="log", yscale="log"):
        """Plot grids particle number and densities values using matplotlib.
        """
        # gather data:
        pivots = self.pivots()
        particles = self.particles()
        densities = self.densities()

        # plot data:
        fig, ax1 = plt.subplots()
        ax1.plot(pivots, particles, "ro", label="numbers")
        ax1.set_ylabel("particles")
        ax2 = ax1.twinx()
        ax2.plot(pivots, densities, "r-", label="density")
        ax2.set_ylabel("density")
        ax1.set_xlabel("size")
        ax1.set_xscale(xscale)
        ax1.set_yscale(yscale)
        ax2.set_yscale(yscale)
        ax1.grid()
        plt.show()

    '''
    # FUTURE IDEAS:
    @staticmethod
    def from_boundaries(boundaries):
        """Create grid from list of range boundaries."""
    
    @staticmethod
    def create_adaptive(start, end, steps, func):
        """Create a grid that adapts automatically to the given function."""
    
    def coarsen_range(self, start, end):
        """Coarsen grid by combining sections in the index range."""
    
    def coarsen(self, times=1):
        """Coarsen whole grid a specified number of times."""
    
    def refine_range(self, start, end):
        """Refine the grid by splitting the sections in the index range."""
    
    def refine(self, times=1):
        """Refine whole grid a specified number of times."""
    
    def smooth_range(self, start, end):
        """Smooth grid by averaging the section sizes in the index range."""
    
    def smooth(self):
        """Smooth entire grid."""
    '''
