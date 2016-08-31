#!/usr/bin/env python3

from collections import deque
import numpy as np
import scipy.integrate as spint
import matplotlib.pyplot as plt


# TODO: documentation!
# TODO: unit testing! --> revisit fancy testig framework!!!
def zero(x):
    """always return zero.
    """
    return 0


class Section:
    """
    Section class representing a bin for the semi open interval [min, max)
    (minimum included, maximum NOT included in the interval).
    It includes a fixed pivot within the sections boundaries
    (min <= pivot < max). The default pivot position is right in the middle
    between the min and max boundaries (pivot = min + size / 2). The length
    of the section is called size and is calculated by (size = max - min).
    Additionally the class contains a field for the total number of particles
    in the section.

               (particles) N
            (particle density) n
             \                /
              min (pivot) max
    ... -o-----+-----o-----+-----o- ...
              v_i  (x_i)  v_i+1
               |^^^^^^^^^^^|
               |  section  |
    """
    def __init__(self, minimum, maximum, particles=0.):
        """Initializer.

        :param minimum: inclusive minimum of the section (v_i).
        :param maximum: exclusive maximum of the section (v_i+1).
        :param particles: total number of particles in the section (N).
        """
        assert particles >= 0
        self._particles = particles
        if minimum < maximum:
            self._minimum = minimum
            self._maximum = maximum
            self._pivot = minimum + (maximum - minimum) / 2
        else:
            raise ValueError(
                "minimum '{min}' must be < maximum '{max}'".format(
                    min=minimum,
                    max=maximum
                ))

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
        """Getter for minimum.

        :return: inclusive minimum of bucket range.
        """
        return self._minimum

    @minimum.setter
    def minimum(self, new_min):
        """Setter for minimum.

        :param new_min: new inclusive minimum of bucket range.
        """
        if new_min < self.maximum:
            self._minimum = new_min
            self._center_pivot()
        else:
            raise ValueError(
                "minimum '{min}' must be smaller than maximum '{max}'".format(
                    min=self.minimum,
                    max=self.maximum
                )
            )

    @property
    def maximum(self):
        """Getter for maximum.

        :return: exclusive maximum of bucket range.
        """
        return self._maximum

    @maximum.setter
    def maximum(self, new_max):
        """Setter for maximum.

        :param new_max: new exclusive maximum of bucket range.
        """
        if new_max > self.minimum:
            self._maximum = new_max
            self._center_pivot()
        else:
            raise ValueError(
                "maximum '{max}' must be greater than minimum '{min}'".format(
                    max=self.maximum,
                    min=self.minimum
                )
            )

    @property
    def range(self):
        return self.minimum, self.maximum

    @range.setter
    def range(self, new_range):
        try:
            minimum = new_range[0]
            maximum = new_range[1]
        except IndexError:
            raise IndexError(
                "range '{rng}' must be a tuple of at least two items!".format(
                    rng=new_range
                )
            )
        if minimum < maximum:
            self.minimum = minimum
            self.maximum = maximum
        else:
            raise ValueError(
                "maximum '{max}' must be greater than minimum '{min}'".format(
                    max=self.maximum,
                    min=self.minimum)
            )

    @property
    def pivot(self):
        """Getter for pivot.

        :return: sections pivot point.
        """
        return self._pivot

    @property
    def size(self):
        """Getter for section size.

        :return: section size (max - min).
        """
        return self.maximum - self.minimum

    @property
    def particles(self):
        """Getter for total number of particles in the section.

        :return: total number of particles in the section.
        """
        return self._particles

    @particles.setter
    def particles(self, new_particles):
        """Setter for total number of particles in the section.

        :param new_particles: new number of particles.
        """
        if new_particles >= 0:
            self._particles = new_particles
        else:
            raise ValueError(
                "number of particles '{par}' must be >= 0!".format(
                    par=self.particles
                )
            )

    @property
    def particle_density(self):
        """Getter for the mean particle density in the section.

        :return: mean particle density in the section.
        """
        return self.particles / self.size

    def _center_pivot(self):
        """Center pivot between minimum and maximum.
        """
        self._pivot = self.minimum + self.size / 2


class Grid:
    """
    Grid class for representing an arbitrary grid of connected sections. The
    grid bins are represented by a "list" of sections with arbitrary ranges.
    The seamless connection of the buckets (max(bucket i) = min(bucket(i+1))
    is guaranteed by checks to catch numerical inaccuracy issues.

                                (particles)
                             (particle density)
                             \                /
                              min (pivot) max
    ---+-----o-----+-----o-----+-----o-----+-----o-----+--->
      ...         i-1  (i-1)   i    (i)   i+1  (i+1)  i+2
       |           |           |^^^^^^^^^^^|           |
       |           |           |  section  |           |
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
        return self._sections[0].minimum  # first section minimum

    @property
    def maximum(self):
        return self._sections[-1].maximum  # last section maximum

    # single/boundary element operations:
    # =========================================================================

    def add_left(self, size, particles=0.):
        """Add new section to the left.

        :param size: range of the new section (max - min).
        """
        assert particles >= 0
        if size > 0:
            current_left_section = self._sections[0]
            new_max = current_left_section.minimum
            new_min = new_max - size
            new_left_section = Section(minimum=new_min, maximum=new_max,
                                       particles=particles)
            self._sections.appendleft(new_left_section)
        else:
            raise ValueError("size must be > 0!")

    def add_right(self, size, particles=0.):
        """Add new section to the right.

        :param size: range of the new section (max - min).
        """
        assert particles >= 0
        if size > 0:
            current_right_section = self._sections[-1]
            new_min = current_right_section.maximum
            new_max = new_min + size
            new_right_section = Section(minimum=new_min, maximum=new_max,
                                        particles=particles)
            self._sections.append(new_right_section)
        else:
            raise ValueError("size must be > 0!")

    def remove_left(self):
        """Remove section on the left.

        :return: removed section.
        """
        if len(self._sections) > 1:
            return self._sections.popleft()
        else:
            raise IndexError("Can not remove last section!")

    def remove_right(self):
        """Remove section on the right.

        :return: removed section.
        """
        if len(self._sections) > 1:
            return self._sections.pop()
        else:
            raise IndexError("Can not remove last section!")

    # multi element operations:
    # =========================================================================

    def coarsen(self, start, end):
        """Coarsen grid by combining buckets in the index range.

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

    # helpers:
    # =========================================================================

    def info(self):
        """print info about the class to console.
        """
        print(self)
        for index, sections in enumerate(self):
            print(index, sections)

    def plot(self):
        """Plot grids particle number values using matplotlib.
        """
        plt.plot(self.pivots(), self.particles(), "ro")
        plt.show()

    def boundaries(self):
        boundary_list = []
        for section in self:
            boundary_list.append(section.minimum)
        boundary_list.append(self._sections[-1].maximum)  # last section max
        return boundary_list

    def pivots(self):
        pivot_list = []
        for section in self:
            pivot_list.append(section.pivot)
        return pivot_list

    def particles(self):
        particle_list = []
        for section in self:
            particle_list.append(section.particles)
        return particle_list

    def particle_densities(self):
        particle_density_list = []
        for section in self:
            particle_density_list.append(section.particle_density)
        return particle_density_list

    # factory methods:
    # =========================================================================
    '''
    @staticmethod
    def create_uniform(minimum, size, amount, func=zero):
        """Create a uniform grid.

        :param minimum: inclusive minimum of the grid (leftmost section).
        :param size: constant range of all the sections.
        :param amount: total number of sections in the grid.
        :return: equidistant Grid object.
        """
        return Grid.create_geometric(minimum=minimum, size=size, factor=1,
                                     amount=amount, func=func)
    '''
    @staticmethod
    def create_uniform(minimum, maximum, amount, func=zero):
        size = (maximum - minimum) / amount
        return Grid.create_geometric(minimum=minimum, maximum=maximum,
                                     factor=1, amount=amount, func=func)

    '''
    @staticmethod
    def create_geometric(minimum, size, factor, amount, func=zero):
        """Create a geometric grid.

        :param minimum: inclusive minimum of the grid (leftmost section).
        :param size: leftmost buckets range.
        :param factor: factor for size change.
        :param amount: total number of sections in the grid.
        :return: geometric Grid object.
        """
        if size <= 0:
            raise ValueError("size '{s}' must be > 0!".format(s=size))
        elif amount <= 0:
            raise ValueError("amount '{a}' must be > 0!".format(a=amount))
        elif factor < 1:
            raise ValueError("factor '{f}' must be >= 1!".format(f=factor))
        else:
            # create initial grid:
            first_max = minimum + size
            particles = spint.quad(func, minimum, first_max)[0]
            grid = Grid(minimum, first_max, particles)
            # add sections to the grid:
            current_size = size
            while len(grid) < amount:
                current_size *= factor
                last_section = grid._sections[-1]
                lower = last_section.maximum
                upper = last_section.maximum + current_size
                particles = spint.quad(func, lower, upper)[0]
                grid.add_right(size=current_size, particles=particles)
            return grid
    '''
    @staticmethod
    def create_geometric(minimum, maximum, amount, factor, func=zero):
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
            return grid

    # possible future additions:
    # =========================================================================

    @staticmethod
    def from_boundaries(boundaries):
        """Create grid from list of range boundaries.
        """
        raise NotImplementedError

    def smooth(self, start, end):
        """Smooth grid by averaging the section sizes in the range
        """
        raise NotImplementedError
