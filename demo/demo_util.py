#!/usr/bin/env python3

"""
Module demonstrating the basic usage of the command line utilities.
"""
from WR.cmdline import prompt_continue, prompt_yes_no, YES, NO
from WR.cmdline import prompt_input_str, prompt_input_int, prompt_input_float
from WR.cmdline import prompt_sure, prompt_input_list
from WR.cmdline import Spinner, Progress


def demo_basic_prompts():
    """Demo basic prompt functions: continue, yes_no, input_str/int/float.
    """
    prompt_continue("Really?")

    res = prompt_yes_no("Yes or Yes?", verbose=False)
    if res == YES:
        print("yooooo!")
    elif res == NO:
        print("nooooo?")
    else:
        print("hm...")

    name = prompt_input_str("Your name?")
    print("Your name is '{}'".format(name))

    age = prompt_input_int("Your age?", verbose=True)
    print("Your age is '{}'".format(age))

    height = prompt_input_float("Your height?", verbose=True)
    print("Your height is '{}'".format(height))


def demo_fancy_prompts():
    """Demo fancy prompt functions: sure, input_list.
    """
    sure_answer = prompt_sure(prompt_input_int, "Your fave number?",
                              verbose=False)
    print("Youre surely fave number is '{}'.".format(sure_answer))

    fruits = prompt_input_list("Your favourite fruits?", length=2, verbose=True)
    print("You like {}".format(fruits))


def demo_spinner():
    """Demo Spinner context manager class.
    """
    MAX = int(1e7)
    print("before")
    with Spinner("finding numbers") as spinner:
        counter = 0
        print("miau")
        for i in range(MAX):
            if i % 3 == 0:
                counter += 1
                spinner.message("{} found, new one: {}".format(counter, i))
        print("wuff")
    print("after")
    print(spinner.read_buffer())
    print("last")


def demo_progress():
    """Demo Progress context manager class.
    """
    MAX = int(1e7)
    print("before")
    with Progress("finding numbers") as progress:
        counter  = 0
        print("miau")
        for i in range(MAX):
            if i % 3 == 0:
                counter += 1
                progress.proceed(i / MAX, "{} found, new one: {}".format(counter, i))
        print("wuff")
    print("after")
    print(progress.read_buffer())
    print("last")


if __name__ == "__main__":
    #demo_basic_prompts()
    #demo_fancy_prompts()
    demo_spinner()
    demo_progress()
