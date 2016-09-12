#!/usr/bin/env python3

import sys
import threading
from itertools import tee, cycle
import io

YES = ["y", "Y", "yes", "Yes", "YES"]
NO = ["n", "N", "no", "No", "NO"]


def pairwise(iterable):
    """Pairwise iterator recipe.

    s -> (s0,s1), (s1,s2), (s2, s3), ...
    taken from: https://docs.python.org/3/library/itertools.html

    :param iterable: iterable object.
    :return: chained pairs of iterable objects elements.
    """
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def kronecker(i, j):
    """Function representing the mathematical Kronecker delta symbol.

    :param i: first argument.
    :param j: second argument.
    :return: 0 if i == j, 1 otherwise.
    """
    return 1 if i == j else 0


def zero(*args, **kwars):
    """Function always returning 0.

    :param args: optional arguments.
    :param kwars: optional keyword arguments.
    :return: always 0.
    """
    return 0


def hstep(x):
    if x == 0:
        return 0.5
    elif x > 0:
        return 1
    elif x < 0:
        return 0
    else:
        raise RuntimeError("This should never happen! '{}'".format(x))


def prompt_continue(msg):
    input("{} [Enter] to continue! ".format(msg))


def prompt_yes_no(msg):
    while True:
        answer = input("\r{} Enter [y/n]! ".format(msg))
        if answer in YES:
            return YES
        elif answer in NO:
            return NO


def prompt_input_str(msg):
    raise NotImplementedError


def prompt_input_int(msg):
    raise NotImplementedError


def prompt_input_float(msg):
    raise NotImplementedError


class Spinner:
    def __init__(self, task, msg="", end_msg="done.", interval=0.25):
        self._task = task
        self._msg = msg
        self._end_msg = end_msg
        self._interval = interval
        self._timer = threading.Timer(self._interval, self.spin)
        self._cycle = cycle("-\|/-\|/")
        self._buffer = io.StringIO()
        self._stdout = None

    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._buffer
        self._print_spinner("_")
        self._timer.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._timer.cancel()
        if self._end_msg is not None:
            self._msg = self._end_msg
        self._print_spinner("x")
        self._stdout.write("\n")
        self._stdout.flush()
        sys.stdout = sys.__stdout__  # self._stdout

    def _print_spinner(self, char):
        self._stdout.write("\r[{s}] {t}: {m}".format(s=char, t=self._task, m=self._msg))
        self._stdout.flush()

    def spin(self):
        self._print_spinner(next(self._cycle))
        self._timer = threading.Timer(self._interval, self.spin)
        self._timer.start()

    def message(self, msg):
        self._msg = msg

    def read_buffer(self):
        self._buffer.seek(0)
        return self._buffer.read().strip()


class Progress:
    def __init__(self, task, msg="", end_msg="done.",
                 percentage=0, length=25, interval=0.25):
        self._task = task
        self._msg = msg
        self._end_msg = end_msg
        self._percentage = percentage
        self._length = length
        self._interval = interval
        self._timer = threading.Timer(self._interval, self.update)
        self._buffer = io.StringIO()
        self._stdout = None

    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._buffer
        self._print_progress()
        self._timer.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._timer.cancel()
        self._percentage = 1
        if self._end_msg is not None:
            self._msg = self._end_msg
        self._print_progress()
        self._stdout.write("\n")
        self._stdout.flush()
        sys.stdout = sys.__stdout__

    def _print_progress(self):
        completed = int(self._percentage * self._length)
        remaining = self._length - completed
        completed_str = "=" * completed
        remaining_str = "." * remaining
        self._stdout.write(
            "\r{p:3.0f}%[{c}{r}] {t}: {m}".format(
                p=self._percentage*100, c=completed_str, r=remaining_str,
                t=self._task, m=self._msg
            )
        )
        self._stdout.flush()

    def update(self):
        self._print_progress()
        self._timer = threading.Timer(self._interval, self.update)
        self._timer.start()

    def proceed(self, percentage, msg=None):
        self._percentage = percentage
        if msg is not None:
            self._msg = msg

    def read_buffer(self):
        self._buffer.seek(0)
        return self._buffer.read().strip()


def main():

    prompt_continue("Really?")
    res = prompt_yes_no("Yes or Yes?")
    if res == YES:
        print("yooooo!")
    elif res == NO:
        print("nooooo?")
    else:
        print("hm...")

    MAX = int(1e7)

    print("before")
    with Spinner("finding numbers") as spinner:
        print("miau")
        for i in range(MAX):
            if i % 3 == 0:
                spinner.message("new one: {}".format(i))
        print("wuff")
    print("after")
    print(spinner.read_buffer())
    print("last")

    print("before")
    with Progress("finding numbers") as progress:
        print("miau")
        for i in range(MAX):
            if i % 3 == 0:
                progress.proceed(i / MAX, "new one: {}".format(i))
        print("wuff")
    print("after")
    print(progress.read_buffer())
    print("last")


if __name__ == "__main__":
    main()
