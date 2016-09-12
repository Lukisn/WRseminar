#!/usr/bin/env python3

from WR.util import *


def demo_prompts():
    prompt_continue("Really?")
    res = prompt_yes_no("Yes or Yes?")
    if res == YES:
        print("yooooo!")
    elif res == NO:
        print("nooooo?")
    else:
        print("hm...")


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
    demo_spinner()
    demo_progress()
