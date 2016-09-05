#!/usr/bin/env python3

from copy import copy
from scipy.integrate import quad


def kronecker_delta(i, j):
    if i == j:
        return 1
    else:
        return 0


# TODO: refactoring
# TODO: documentaion
class FixedPivot:

    def __init__(self, initial_ndf, break_freq, child_number, agg_freq,
                 primary=0, secondary=1, breakage=True, aggregation=True,
                 growth=True, nucleation=False):
        self.primary = primary  # primary preserved moment "zeta" usually zeta=0 (preservation of numbers)
        self.secondary = secondary  # secondary preserved moment "nu" usually nu=1 (preservation of mass)
        self.initial_ndf = initial_ndf
        self.current_ndf = copy(initial_ndf)
        self.previous_ndf = copy(initial_ndf)
        self.break_freq = break_freq  # breakage frequency function "Gamma"
        self.child_number = child_number  # function for number of child particles formed due to breakage
        self.agg_freq = agg_freq  # aggregation frequency function "Q"
        self.breakage = breakage
        self.aggregation = aggregation
        self.growth = growth
        self.nucleation = nucleation

    # TODO: implement toggles for switching single calculations on and off

    # TODO: implement results output with options.
    def simulate(self, end_time, steps):
        start_time = 0
        time_step = (end_time - start_time) / steps
        current_time = start_time

        # apply change to time step:
        step_counter = 0
        while current_time < end_time:
            current_time += time_step
            step_counter += 1
            print("step=", step_counter, "current_time=", current_time)

            # calculate terms of the PBE:
            self.previous_ndf = copy(self.current_ndf)
            if self.breakage:
                self.calc_breakage(time_step)
            if self.aggregation:
                self.calc_aggregration(time_step)
            if self.growth:
                self.calc_growth(time_step)
            if self.nucleation:
                self.calc_nucleation(time_step)

            # TODO: write out intermediate result

        # TODO: write out final result

    def calc_breakage(self, time_step):
        """calculate breakage.
        """
        zeta = self.primary
        nu = self.secondary
        beta = self.child_number

        for i, section_i in enumerate(self.previous_ndf):
            print("i=", i, "section_i=", section_i)

            xi = section_i.pivot
            Ni = section_i.particles
            gammai = self.break_freq(xi)

            # TODO: find out what to do with the ghost sections!?
            try:
                xip1 = self.current_ndf.section(i+1).pivot
            except AssertionError:
                xip1 = None
            try:
                xim1 = self.current_ndf.section(i-1).pivot
            except AssertionError:
                xim1 = None

            # calculate birth:
            birth = 0
            for k, section_k in enumerate(self.previous_ndf):
                if k >= i:
                    #print("k=", k, "section_k=", section_k)

                    xk = section_k.pivot
                    Nk = section_k.particles
                    gammak = self.break_freq(xk)

                    # calculate n_i,k:
                    def func_primary(v):
                        return v**zeta * beta(v, xk)

                    def func_secondary(v):
                        return v**nu * beta(v, xk)

                    if xip1 is None:
                        first_term = 0
                    else:
                        Bzeta = quad(func_primary, xi, xip1)[0]
                        Bnu = quad(func_secondary, xi, xip1)[0]
                        first_num = Bzeta * xip1 ** nu - Bnu * xip1 ** zeta
                        first_denom = xi ** zeta * xip1 ** nu - xi ** nu * xip1 ** zeta
                        first_term = first_num / first_denom

                    if xim1 is None:
                        second_term = 0
                    else:
                        Bzetam1 = quad(func_primary, xim1, xi)[0]
                        Bnum1 = quad(func_secondary, xim1, xi)[0]
                        second_num = Bzetam1 * xim1**nu - Bnum1 * xim1**zeta
                        second_denom = xi**zeta * xim1**nu - xi**nu * xim1**zeta
                        second_term = second_num / second_denom

                    nik = first_term + second_term

                    # calculate actual birth:
                    birth += nik * gammak * Nk

            # calculate death:
            death = gammai * Ni

            # calculate new ndf:
            self.current_ndf.section(i).particles += time_step * (birth - death)

    def calc_aggregration(self, time_step):
        """calculate aggregation.
        """
        zeta = self.primary
        nu = self.secondary

        for i, section_i in enumerate(self.previous_ndf):
            print("i=", i, "section_i=", section_i)

            xi = section_i.pivot
            Ni = section_i.particles
            gammai = self.break_freq(xi)

            # TODO: find out what to do with the ghost sections!?
            try:
                xip1 = self.current_ndf.section(i + 1).pivot
            except AssertionError:
                xip1 = None
            try:
                xim1 = self.current_ndf.section(i - 1).pivot
            except AssertionError:
                xim1 = None

            # calculate birth:
            birth = 0
            for j, section_j in enumerate(self.previous_ndf):
                print("j=", j, "section_j=", section_j)

                xj = section_j.pivot
                Nj = section_j.particles

                for k, section_k in enumerate(self.previous_ndf):
                    print("k=", k, "section_k=", section_k)

                    xk = section_k.pivot
                    Nk = section_k.particles

                    v = xj + xk
                    if j >= k:
                        if xim1 is None:
                            num = v ** zeta * xip1 ** nu - v ** nu * xip1 ** zeta
                            denom = xi ** zeta * xip1 ** nu - xi ** nu * xip1 ** zeta
                            eta = num / denom
                        elif xip1 is None:
                            num = v ** zeta * xim1 ** nu - v ** nu * xim1 ** zeta
                            denom = xi ** zeta * xim1 ** nu - xi ** nu * xim1 ** zeta
                            eta = num / denom
                        elif xim1 <= v <= xip1:
                            if xi <= v <= xip1:
                                num = v**zeta * xip1**nu - v**nu * xip1**zeta
                                denom =xi**zeta * xip1**nu - xi**nu * xip1**zeta
                                eta = num / denom
                            elif xim1 <= v <= xi:
                                num = v**zeta * xim1**nu - v**nu * xim1**zeta
                                denom =xi**zeta * xim1**nu - xi**nu * xim1 ** zeta
                                eta = num / denom
                            else:
                                raise RuntimeError("unable to calc eta!")

                            djk = kronecker_delta(j, k)
                            Qjk = self.agg_freq(xj, xk)

                            birth += (1 - .5 * djk) * eta * Qjk * Nj * Nk

            # calculate death:
            death = 0
            for k, section_k in enumerate(self.previous_ndf):
                xk = section_k.pivot
                Nk = section_k.particles
                Qik = self.agg_freq(xi, xk)
                death += Qik * Nk
            death *= Ni

            # TODO: implement checking of new particle size and issue warning if < 0
            self.current_ndf.section(i).particles += time_step * (birth - death)

    def calc_growth(self, time_step):
        """calculate growth.
        """
        pass

    def calc_nucleation(self, time_step):
        """calculate nucleation.
        """
        pass


def main():
    from WR.grid import Grid
    import numpy as np

    def f(v, N0=1, v0=1):
        """initial number density function.
        """
        return (N0 / v0) * (v / v0) * np.exp(-v / v0)
    ini = Grid.create_geometric_end(0, 1000, 100, 1.1, f)
    #ini._plot(scale="log")

    def gamma(v):
        """breakage frequency function.
        """
        return v*v

    def beta(v1, v2):
        """child number function.
        """
        return 2  #v2/2
    def Q(x1, x2):
        """aggregation frequency function.
        """
        return 1
    method = FixedPivot(initial_ndf=ini, break_freq=gamma, child_number=beta, agg_freq=Q)

    method.simulate(end_time=.1, steps=1)
    method.current_ndf._plot(scale="log")


if __name__ == "__main__":
    main()
