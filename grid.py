#!/usr/bin/env python3

from collections import deque
import matplotlib.pyplot as plt


class Section:
    """
    Section class representing a bin for the semi open interval [min, max)
    (minimum included, maximum NOT included in the interval).
    It includes a fixedpivot within the buckets boundaries
    (min <= pivot < max). The default pivot position is right in the middle
    between the min and max boundaries (pivot = min + delta / 2). The length
    of the section is called delta (= max - min).
    """
    def __init__(self, minimum, maximum):
        """Initializer.

        :param minimum: inclusive minimum of the section
        :param maximum: exclusive maximum of the section
        """
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

    def __repr__(self):
        """Callable representation.

        :return: string representation for usage with eval().
        """
        return "{cls}({min}, {max})".format(
            cls=self.__class__.__name__,
            min=self.minimum,
            max=self.maximum,
        )

    def __str__(self):
        """String representation.

        :return: human readable string representation.
        """
        return "<{cls} @ {id}: min={min:.3e}, max={max:.3e}, " \
               "delta={de:.3e}, pivot={piv:.3e}>".format(
                cls=self.__class__.__name__,
                id=hex(id(self)),
                min=self.minimum,
                max=self.maximum,
                de=self.size,
                piv=self.pivot
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
    """
    def __init__(self, minimum=0, maximum=1):
        """Initializer.

        :param minimum: initial section minimum (default: 0).
        :param maximum: initial section maximum (default: 1).
        """
        self._sections = deque()  # Double Ended QUEue
        initial_section = Section(minimum, maximum)
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

    # single/boundary element operations:
    # =========================================================================

    def add_left(self, size):
        """Add new section to the left.

        :param size: range of the new section (max - min).
        """
        if size > 0:
            left_section = self._sections[0]
            new_max = left_section.minimum
            new_min = new_max - size
            new_section = Section(new_min, new_max)
            self._sections.appendleft(new_section)
        else:
            raise ValueError("size must be > 0!")

    def add_right(self, size):
        """Add new section to the right.

        :param size: range of the new section (max - min).
        """
        if size > 0:
            right_section = self._sections[-1]
            new_min = right_section.maximum
            new_max = new_min + size
            new_section = Section(new_min, new_max)
            self._sections.append(new_section)
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

    def coarsen_range(self, start, end):
        """Coarsen grid by combining buckets in the index range.

        :param start: inclusive start of index range.
        :param end: inclusive end of index range.
        """
        if start <= end:
            if start >= 0 and end < len(self):
                new_sections = deque()
                coarse_min = self._sections[start].minimum
                coarse_max = self._sections[end].maximum
                coarsened = False
                for index, bucket in enumerate(self._sections):
                    if index < start or index > end:  # outside range
                        new_sections.append(self._sections[index])
                    else:  # inside coarsening range
                        if not coarsened:
                            coarse_bucket = Section(minimum=coarse_min,
                                                    maximum=coarse_max)
                            new_sections.append(coarse_bucket)
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

    def refine_range(self, start, end):
        """Refine the grid by splitting the sections in the index range.

        :param start: inclusive start of index range.
        :param end: inclusive end of index range.
        """
        if start <= end:
            if start >= 0 and end <= len(self):
                new_sections = deque()
                for index, bucket in enumerate(self._sections):
                    if index < start or index > end:  # outside range
                        new_sections.append(self._sections[index])
                    else:  # inside refining range
                        left = self._sections[index].minimum
                        right = self._sections[index].maximum
                        center = left + (right - left) / 2
                        left_section = Section(minimum=left, maximum=center)
                        right_section = Section(minimum=center, maximum=right)
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

    def print_info(self):
        """print info about the class to console.
        """
        print(self)
        for index, sections in enumerate(self):
            print(index, sections)

    def boundary_list(self):
        boundary_list = []
        for section in self:
            boundary_list.append(section.minimum)
        boundary_list.append(self._sections[-1].maximum)  # last section max
        return boundary_list

    def pivot_list(self):
        pivot_list = []
        for section in self:
            pivot_list.append(section.pivot)
        return pivot_list

    # factory methods:
    # =========================================================================

    @staticmethod
    def create_uniform(minimum, size, amount):
        """Create a uniform grid.

        :param minimum: inclusive minimum of the grid (leftmost section).
        :param size: constant range of all the sections.
        :param amount: total number of sections in the grid.
        :return: equidistant grid.
        """
        assert size > 0
        assert amount > 0
        return Grid.create_geometric(minimum=minimum, initial_size=size,
                                     factor=1, amount=amount)

    @staticmethod
    def create_geometric(minimum, initial_size, factor, amount):
        """Create a geometric grid with a factor rule for size change
        (size(bucket i+1) = factor * size(bucket i)).

        :param minimum: inclusive minimum of the grid (leftmost section).
        :param initial_size: leftmost buckets range.
        :param factor: factor for size change.
        :param amount: total number of sections in the grid.
        :return: geometric grid.
        """
        assert initial_size > 0
        assert factor >= 1
        assert amount > 0
        # create initial grid:
        first_max = minimum + initial_size
        grid = Grid(minimum=minimum, maximum=first_max)
        # add sections to the grid:
        current_size = initial_size
        while len(grid) < amount:
            current_size *= factor
            grid.add_right(current_size)
        return grid

    # possible future additions:
    # =========================================================================
    @staticmethod
    def from_list(input_list):
        """Create grid from list of range boundaries.
        """
        raise NotImplementedError

    def smooth(self, start, end):
        """Smooth grid by averaging the section sizes in the range
        """
        raise NotImplementedError
