#!/usr/bin/env python3
import cases.breakage as breakage
import cases.aggregation as aggregation
import cases.break_agg as break_agg
import cases.growth as growth
import cases.growth1 as growth1
import cases.growth4 as growth4
from util.filesystem import remove


def main(show_plots, remove_data_files=False, remove_plot_files=False):
    """Run all cases."""
    cases = [breakage, aggregation, break_agg, growth, growth1, growth4
             ]
    for case in cases:
        print("running case {}:".format(case.__name__))
        case.main(show_plots)
        print()  # empty line

    if remove_data_files:
        print("removing *.dat files:")
        remove("./results", "*.dat")
        print()

    if remove_plot_files:
        print("removing *.eps files:")
        remove("./results", "*.eps")
        print()


if __name__ == "__main__":
    main(show_plots=False,
         remove_data_files=True,
         remove_plot_files=False)
