#!/usr/bin/env python3

from WR.util import prompt_continue, prompt_yes_no, YES, NO
from WR.util import prompt_input_str, prompt_input_int, prompt_input_float
from WR.util import prompt_sure, prompt_input_list
from WR.util import Spinner, Progress


def demo_prompts():
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
    sure_answer = prompt_sure(prompt_input_int, "Your fave number?",
                              verbose=False)
    print("Youre surely fave number is '{}'.".format(sure_answer))

    fruits = prompt_input_list("Your favourite fruits?", length=2, verbose=True)
    print("You like {}".format(fruits))


def demo_spinner():
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


def demo_progress():
    MAX = int(1e7)
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
    demo_prompts()
    demo_fancy_prompts()
    demo_spinner()
    demo_progress()
