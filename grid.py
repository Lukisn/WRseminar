#!/usr/bin/env python3

from collections import deque
import matplotlib.pyplot as plt


class Bucket:
    """
    Bucket class representing a bin for the semi open interval [min, max).
    It includes an arbitrarily located pivot within the buckets boundaries
    (min <= pivot < max). The default pivot position is right in the middle
    between the min and max boundaries (pivot = min + delta / 2) and can
    freely be changed.
    """
    def __init__(self, minimum, maximum):
        """Initializer.

        :param minimum: inclusive minimum of the bucket range
        :param maximum: exclusive maximum of the bucket range
        """
        if minimum <= maximum:
            self._minimum = minimum
            self._maximum = maximum
            self._pivot = minimum + (maximum - minimum) / 2
        else:
            raise ValueError(
                "minimum '{min}' must be <= maximum '{max}'".format(
                    min=minimum,
                    max=maximum
                ))

    def __repr__(self):
        """Callable representation.

        :return: string representation for usage with eval().
        """
        return "{cls}({min}, {max}, {piv})".format(
            cls=self.__class__.__name__,
            min=self.minimum,
            max=self.maximum,
            piv=self.pivot
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
                de=self.delta,
                piv=self.pivot
                )

    @property
    def minimum(self):
        """Getter for minimum.

        :return: inclusive minimum of bucket range.
        """
        return self._minimum

    @minimum.setter
    def minimum(self, new_min, center_pivot=True):
        """Setter for minimum.

        :param new_min: new inclusive minimum of bucket range.
        :param center_pivot: flag for recentering of the pivot.
        :return:
        """
        if new_min < self.maximum:
            self._minimum = new_min
            if center_pivot:
                self.center_pivot()
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
    def maximum(self, new_max, center_pivot=True):
        """Setter for maximum.

        :param new_max: new exclusive maximum of bucket range.
        :param center_pivot: flag for recentering the pivot.
        :return:
        """
        if new_max > self.minimum:
            self._maximum = new_max
            if center_pivot:
                self.center_pivot()
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

        :return: buckets pivot point.
        """
        return self._pivot

    @pivot.setter
    def pivot(self, new_pivot):
        """Setter for pivot.

        :param new_pivot: new pivot point.
        :return:
        """
        if self.minimum <= new_pivot < self.maximum:
            self._pivot = new_pivot
        else:
            raise ValueError(
                "pivot '{piv}' outside range [{min}, {max})!".format(
                    piv=new_pivot,
                    min=self.minimum,
                    max=self.maximum
                )
            )

    def center_pivot(self):
        """Center pivot between minimum and maximum.

        :return: None
        """
        self.pivot = self.minimum + self.delta / 2

    @property
    def delta(self):
        """Getter for delta.

        :return: buckets range (max - min).
        """
        return self.maximum - self.minimum


class Grid:
    """
    Grid class for representing an arbitrary grid. The grid fields are
    represented by a "list" of buckets with arbitrary ranges and pivot point
    locations. The connection of the buckets (max(bucket i) = min(bucket(i+1))
    is garanteed by checks to catch numerical inaccuracy issues.
    """
    def __init__(self, minimum=0, maximum=1):
        """Initializer.

        :param minimum: initial buckets minimum (default = 0).
        :param maximum: initial buckets maximum (default = 1).
        """
        self._buckets = deque()
        default_bucket = Bucket(minimum, maximum)
        self._buckets.append(default_bucket)

    def __len__(self):
        """Length interface.

        :return: number of buckets in the grid.
        """
        return len(self._buckets)

    def __iter__(self):
        """Iterator interface.

        :return: iterator of bucket "list".
        """
        return iter(self._buckets)

    def __str__(self):
        """String representation.

        :return: human readable string representation.
        """
        minimum = self._buckets[0].minimum
        maximum = self._buckets[-1].maximum
        number = len(self)
        return "<{cls} @{id}: min={min:.3e}, max={max:.3e}, num={num}>".format(
            cls=self.__class__.__name__,
            id=hex(id(self)),
            min=minimum,
            max=maximum,
            num=number
        )

    # single element operations:

    def add_left(self, delta):
        """Add bucket to the left.

        :param delta: absolute range of the new bucket (max - min).
        :return: None.
        """
        if delta > 0:
            current_left_bin = self._buckets[0]
            new_max = current_left_bin.minimum
            new_min = new_max - delta
            new_left_bin = Bucket(new_min, new_max)
            self._buckets.appendleft(new_left_bin)
        else:
            raise ValueError("delta must be > 0!")

    def add_right(self, delta):
        """Add bucket to the right.

        :param delta: absolute range of the new bucket (max - min).
        :return: None.
        """
        if delta > 0:
            current_right_bin = self._buckets[-1]
            new_min = current_right_bin.maximum
            new_max = new_min + delta
            new_right_bin = Bucket(new_min, new_max)
            self._buckets.append(new_right_bin)
        else:
            raise ValueError("delta must be > 0!")

    def remove_left(self):
        """Remove bucket on the left.

        :return: removed bucket.
        """
        if len(self._buckets) > 1:
            return self._buckets.popleft()
        else:
            raise IndexError("Can't remove last bin!")

    def remove_right(self):
        """Remove bucket on the right.

        :return: removed bucket.
        """
        if len(self._buckets) > 1:
            return self._buckets.pop()
        else:
            raise IndexError("Can't remove last bin!")

    # multi element operations:

    def coarsen_range(self, start, end):
        """Coarsen grid by combining buckets in the index range.

        :param start: start of index range.
        :param end: end of index range.
        :return: None.
        """
        if start <= end:
            if start >= 0 and end < len(self):
                new_buckets = deque()
                coarse_min = self._buckets[start].minimum
                coarse_max = self._buckets[end].maximum
                coarsened = False
                for index, bucket in enumerate(self._buckets):
                    if index < start or index > end:
                        new_buckets.append(self._buckets[index])
                    else:
                        if not coarsened:
                            coarse_bucket = Bucket(minimum=coarse_min,
                                                   maximum=coarse_max)
                            new_buckets.append(coarse_bucket)
                            coarsened = True
                self._buckets = new_buckets
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
        """Refine the grid by splitting the bucktes in the index range.

        :param start: start of index range.
        :param end: end of index range.
        :return: None.
        """
        if start <= end:
            if start >= 0 and end <= len(self):
                new_buckets = deque()
                for index, bucket in enumerate(self._buckets):
                    if index < start or index > end:
                        new_buckets.append(self._buckets[index])
                    else:
                        left = self._buckets[index].minimum
                        right = self._buckets[index].maximum
                        middle = left + (right - left) / 2
                        first_bucket = Bucket(minimum=left, maximum=middle)
                        second_bucket = Bucket(minimum=middle, maximum=right)
                        new_buckets.append(first_bucket)
                        new_buckets.append(second_bucket)
                self._buckets = new_buckets
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

    def print_info(self):
        """print info about the class to console.

        :return: None.
        """
        print(self)
        for index, buckets in enumerate(self):
            print(index, buckets)

    def boundary_list(self):
        boundary_list = []
        for bucket in self:
            boundary_list.append(bucket.minimum)
        boundary_list.append(self._buckets[-1].maximum)
        return boundary_list

    def pivot_list(self):
        pivot_list = []
        for bucket in self:
            pivot_list.append(bucket.pivot)
        return pivot_list

    # factory methods:

    @staticmethod
    def build_equidistant(minimum, delta, amount):
        """Create an equidistant grid.

        :param minimum: inclusive minimum of the grid (leftmost bucket).
        :param delta: constant range of all the buckets.
        :param amount: total number of buckets in the grid.
        :return: equidistant grid.
        """
        assert delta > 0
        assert amount > 0
        return Grid.build_factor(minimum=minimum, initial_delta=delta,
                                 factor=1, amount=amount)

    @staticmethod
    def build_factor(minimum, initial_delta, factor, amount):
        """Create an non-equidistant grid with factor rule for deltas
        (delta(bucket i+1) = factor * delta(bucket i)).

        :param minimum: inclusive minimum of the grid (leftmost bucket).
        :param initial_delta: leftmost buckets range.
        :param factor: factor for further ranges
        :param amount: total number of buckets in the grid.
        :return: non-equidistant grid.
        """
        assert initial_delta > 0
        assert factor >= 1
        assert amount > 0
        first_max = minimum + initial_delta
        grid = Grid(minimum=minimum, maximum=first_max)
        current_delta = initial_delta
        while len(grid) < amount:
            current_delta *= factor
            grid.add_right(current_delta)
        return grid

    # possible future additions:
    @staticmethod
    def from_list(input_list):
        raise NotImplementedError

    def smooth(self, range):
        raise NotImplementedError


def demo():
    """short demo of basic use."""

    # create equidistant grid
    equi = Grid.build_equidistant(minimum=0, delta=1, amount=10)
    equi.print_info()
    xbe0 = equi.boundary_list()
    ybe0 = [0] * len(xbe0)
    xpe0 = equi.pivot_list()
    ype0 = [0] * len(xpe0)
    # coarsen
    equi.coarsen_range(start=5, end=9)
    equi.print_info()
    xbe1 = equi.boundary_list()
    ybe1 = [1] * len(xbe1)
    xpe1 = equi.pivot_list()
    ype1 = [1] * len(xpe1)
    # refine
    equi.refine_range(start=3, end=5)
    equi.print_info()
    xbe2 = equi.boundary_list()
    ybe2 = [2] * len(xbe2)
    xpe2 = equi.pivot_list()
    ype2 = [2] * len(xpe2)

    # create non-equidistant grid
    fact = Grid.build_factor(minimum=0, initial_delta=1, factor=1.1, amount=10)
    fact.print_info()
    xbf0 = fact.boundary_list()
    ybf0 = [0] * len(xbf0)
    xpf0 = fact.pivot_list()
    ypf0 = [0] * len(xpf0)
    # coarsen
    fact.coarsen_range(start=5, end=9)
    fact.print_info()
    xbf1 = fact.boundary_list()
    ybf1 = [1] * len(xbf1)
    xpf1 = fact.pivot_list()
    ypf1 = [1] * len(xpf1)
    # refine
    fact.refine_range(start=3, end=5)
    fact.print_info()
    xbf2 = fact.boundary_list()
    ybf2 = [2] * len(xbf2)
    xpf2 = fact.pivot_list()
    ypf2 = [2] * len(xpf2)

    # plot grids
    plt.title("equidistant grid")
    plt.plot(xbe0, ybe0, "rx", label="initial boundary")
    plt.plot(xpe0, ype0, "ro", label="initial pivot")
    plt.plot(xbe1, ybe1, "gx", label="coarsened boundary")
    plt.plot(xpe1, ype1, "go", label="coarsened pivot")
    plt.plot(xbe2, ybe2, "bx", label="refined boundary")
    plt.plot(xpe2, ype2, "bo", label="refined pivot")
    plt.xlim(-1, 11)
    plt.ylim(-1, 5)
    plt.legend()
    plt.grid()
    plt.show()

    plt.title("factor grid")
    plt.plot(xbf0, ybf0, "rx", label="initial boundary")
    plt.plot(xpf0, ypf0, "ro", label="initial pivot")
    plt.plot(xbf1, ybf1, "gx", label="coarsened boundary")
    plt.plot(xpf1, ypf1, "go", label="coarsened pivot")
    plt.plot(xbf2, ybf2, "bx", label="refined boundary")
    plt.plot(xpf2, ypf2, "bo", label="refined pivot")
    plt.xlim(-1, 17)
    plt.ylim(-1, 5)
    plt.legend()
    plt.grid()
    plt.show()


if __name__ == "__main__":
    demo()
