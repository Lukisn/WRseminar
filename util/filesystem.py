#!/usr/bin/env python3
"""Utility functions for interacting with the file system."""
import os


def remove(directory, extension, verbose=True):
    """Remove all files having the ``pattern`` from ``directory``."""
    abs_dir = os.path.abspath(directory)
    if not os.path.exists(abs_dir):
        raise ValueError("Directory '{}' does not exist!".format(directory))

    for entry in os.listdir(abs_dir):
        entry = os.path.join(abs_dir, entry)
        if entry.endswith(extension.lstrip("*")):
            if verbose:
                print("removing '{}'...".format(entry))
            os.remove(entry)