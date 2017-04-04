#!/usr/bin/env python3

"""
Module implementing some basic command line utilities.

The module contains multiple functions encapsulating some general command line
interactions with the user. Those functions are prefixed with ``prompt_``.
Also this module contains two classes for dynamic command line output on a
single line. Those classes are designed to be used as context managers.
"""
import io
import sys
import threading
from time import perf_counter
from itertools import cycle


# CONSTANTS: ------------------------------------------------------------------

YES = ["y", "yes"]  # list of possible / accepted 'yes' answers
NO = ["n", "no"]  # list of possible / accepted 'no' answers


# FUNCTIONS: ------------------------------------------------------------------

def prompt_continue(msg):
    """Prompt a message and wait for [Enter] to continue.

    :param msg: Message to be displayed.
    """
    input("{} [Enter] to continue! ".format(msg))


def prompt_yes_no(msg, verbose=True):
    """Prompt a question and ask for YES or NO as an answer.

    This prompt continues to ask for an answer if no valid one was given.

    :param msg: Question prompt to be displayed.
    :param verbose: Flag for verbose output on errors.
    :return: YES or NO.
    """
    while True:
        answer = input("\r{} Enter [y/n]! ".format(msg))
        answer = answer.strip()  # strip leading and trailing whitespace
        answer = answer.lower()  # make all lowercase for easier comparison
        if answer in YES:
            return YES
        elif answer in NO:
            return NO
        else:
            if verbose:
                print("Unclear answer '{}'! Enter 'y(es)' or 'n(o)'.".format(
                    answer
                ))


def prompt_input_str(msg, strip=True):
    """Prompt to input any string.

    :param msg: Message to be displayed.
    :param strip: Flag for stripping of leading and trailing whitespace.
    :return: input string.
    """
    answer = input("\r{} ".format(msg))
    if strip:
        return answer.strip()
    else:
        return answer


def prompt_input_int(msg, verbose=True):
    """Prompt to in put an integer.

    :param msg: Message to be displayed.
    :param verbose: Flag for verbose output on errors.
    :return: input integer.
    """
    while True:
        answer = input("\r{} (int) ".format(msg))
        try:
            answer_int = int(answer)
            return answer_int
        except:
            if verbose:
                print("'{}' is NOT a valid integer!".format(answer))


def prompt_input_float(msg, verbose=True):
    """Prompt to input a float.

    :param msg: Message to be displayed.
    :param verbose: Flag for verbose output on errors.
    :return: input float.
    """
    while True:
        answer = input("\r{} (float) ".format(msg))
        try:
            answer_int = float(answer)
            return answer_int
        except:
            if verbose:
                print("'{}' is NOT a valid float!".format(answer))


def prompt_sure(prompt, msg, verbose=True):
    """Prompt any input and ask for acknowledgement of the input.

    :param prompt: prompt function to call for input.
    :param msg: Message to be displayed.
    :param verbose: Flag for verbose output on errors.
    :return: input.
    """
    while True:
        answer = prompt(msg, verbose)
        decision = prompt_yes_no(
            "Are you sure with '{}'?".format(answer), verbose
        )
        if decision == YES:
            return answer


def prompt_input_list(msg, length=None, verbose=True):
    """Prompt to input a list of strings.

    :param msg: Message to be displayed.
    :param length: accepted length of the list (default: None).
    :param verbose: Flag for verbose output on errors (default: True).
    :return: list of strings.
    """
    while True:
        if length is None:
            length_str = ""
        else:
            length_str = length
        answer = input("\r{} ({}x) ".format(msg, length_str))
        answer = answer.strip()  # strip leading and trailing whitespace
        answer_list = answer.split()
        input_length = len(answer_list)
        if length is not None and input_length == length:
            return answer_list
        else:
            if verbose:
                if length > input_length:
                    print("List too long! Only {} of {} items given.".format(
                        input_length, length)
                    )
                elif length < input_length:
                    print("List too short! {} too many items given.".format(
                        input_length - length)
                    )
                else:
                    raise RuntimeError("This should NEVER happen!")


# CONTEXT MANAGERS: -----------------------------------------------------------

class Spinner:
    """Context manager class implementing a command line spinner."""
    def __init__(self, task, msg="", end_msg="done.",
                 interval=0.25, timing=True):
        """Initializer.

        :param task: task name to display.
        :param msg: message to display (default: '').
        :param end_msg: end message to show on finishing (default: 'done.').
        :param interval: refreshing interval in seconds(default: 0.25).
        """
        self._task = task
        self._msg = msg
        self._end_msg = end_msg
        self._interval = interval
        self._timing = timing
        self._timer = threading.Timer(self._interval, self.spin)
        self._cycle = cycle("-\|/-\|/")
        self._buffer = io.StringIO()
        self._stdout = None
        self._start_time = None

    def __enter__(self):
        """Enter context manager.

        stdout is redirected to an internal buffer while the spinner is on.
        """
        self._stdout = sys.stdout
        sys.stdout = self._buffer
        self._print_spinner("_")
        self._start_time = perf_counter()
        self._timer.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Exit context manager."""
        self._timer.cancel()
        end_time = perf_counter()
        time_taken = (end_time - self._start_time) * 1000  # in ms
        if self._end_msg is not None:
            self._msg = self._end_msg
        self._print_spinner("x")
        if self._timing:
            self._stdout.write(" (in {:.2f} ms)".format(time_taken))
        self._stdout.write("\n")
        self._stdout.flush()
        sys.stdout = sys.__stdout__  # self._stdout

    def _print_spinner(self, char=None):
        """Print spinner to command line.

        :param char: spinner character to display.
        """
        if char is None:
            char = next(self._cycle)
        self._stdout.write("\r[{s}] {t}: {m}".format(
            s=char, t=self._task, m=self._msg)
        )
        self._stdout.flush()

    def spin(self):
        """Triggering printing of the spinner and restart the timer."""
        self._print_spinner()
        self._timer = threading.Timer(self._interval, self.spin)
        self._timer.start()

    def message(self, msg):
        """Change the message to be displayed on the next refresh.

        :param msg: Message to be displayed.
        """
        self._msg = msg

    def read_buffer(self):
        """Getter for buffer.

        :return: buffer of output saved while spinner was on.
        """
        self._buffer.seek(0)
        return self._buffer.read().strip()


class Progress:
    """Context manager implementing a command line progress bar."""
    def __init__(self, task, msg="", end_msg="done.",
                 percentage=0, length=25, interval=0.25, timing=True):
        """Initializer.

        :param task: task name to display.
        :param msg: message to display (default: '').
        :param end_msg: end message to show on finishing (default: 'done.').
        :param percentage: current percentage to be displayed (default: 0).
        :param length: length of the progress bar in chars (default: 25).
        :param interval: refreshing interval in seconds(default: 0.25).
        """
        self._task = task
        self._msg = msg
        self._end_msg = end_msg
        self._percentage = percentage
        self._length = length
        self._interval = interval
        self._timing = timing
        self._timer = threading.Timer(self._interval, self.update)
        self._buffer = io.StringIO()
        self._stdout = None
        self._completed_char = "#"
        self._progressing_char = ">"
        self._remaining_char = "_"

    def __enter__(self):
        """Enter context manager."""
        self._stdout = sys.stdout
        sys.stdout = self._buffer
        self._print_progress()
        self._start_time = perf_counter()
        self._timer.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Exit context manager."""
        self._timer.cancel()
        end_time = perf_counter()
        time_taken = (end_time - self._start_time) * 1000  # in ms
        self._percentage = 1
        if self._end_msg is not None:
            self._msg = self._end_msg
        self._print_progress()
        if self._timing:
            self._stdout.write(" (in {:.2f} ms)".format(time_taken))
        self._stdout.write("\n")
        self._stdout.flush()
        sys.stdout = sys.__stdout__

    def _print_progress(self):
        """Print progress bar to command line."""
        completed = int(self._percentage * self._length)
        remaining = self._length - completed
        if self._percentage < 1:
            completed_str = self._completed_char * (completed - 1) + self._progressing_char
        else:
            completed_str = self._completed_char * completed
        remaining_str = self._remaining_char * remaining
        self._stdout.write(
            "\r{p:3.0f}%[{c}{r}] {t}: {m}".format(
                p=self._percentage*100, c=completed_str, r=remaining_str,
                t=self._task, m=self._msg
            )
        )
        self._stdout.flush()

    def update(self):
        """Triggering printing of the progress bar  and restart the timer."""
        self._print_progress()
        self._timer = threading.Timer(self._interval, self.update)
        self._timer.start()

    def proceed(self, percentage, msg=None):
        """Proceed the progress bar and (optionally) change the message.

        :param percentage: percentage to proceed the progress bar to.
        :param msg: Message to be displayed (default: None).
        """
        self._percentage = percentage
        if msg is not None:
            self._msg = msg

    def read_buffer(self):
        """Getter for buffer.

        :return: buffer of output saved while spinner was on.
        """
        self._buffer.seek(0)
        return self._buffer.read().strip()
