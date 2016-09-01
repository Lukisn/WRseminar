#!/usr/bin/env python3

# TODO: implement basic skeleton and implement breakage
class FixedPivot:

    def __init__(self, initial_ndf):
        self.initial_ndf = initial_ndf

    def simulate(self, end_time, steps):
        start_time = 0
        time_step = (end_time - start_time) / steps

        # calculate terms of the PBE:
        bre = self.calc_breakage()
        agg = self.calc_aggregration()
        gro = self.calc_growth()
        nuc = self.calc_nucleation()

        # apply change to time step:
        pass

    def calc_breakage(self):
        pass

    def calc_aggregration(self):
        pass

    def calc_growth(self):
        pass

    def calc_nucleation(self):
        pass

