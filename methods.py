#!/usr/bin/env python3

from scipy.integrate import quad


# TODO: implement basic skeleton and implement breakage
class FixedPivot:

    def __init__(self, primary, secondary, initial_ndf, break_freq, child_number):
        self.primary = primary  # primary preserved moment "zeta" usually zeta=0 (preservation of numbers)
        self.secondary = secondary  # secondary preserved moment "nu" usually nu=1 (preservation of mass)
        self.initial_ndf = initial_ndf
        self.current_ndf = initial_ndf
        self.break_freq = break_freq  # breakage frequency function "Gamma"
        self.child_number = child_number  # function for number of child particles formed due to breakage

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
            self.calc_breakage()
            #self.calc_aggregration()
            #self.calc_growth()
            #self.calc_nucleation()


    def calc_breakage(self):
        """calculate breakage.
        """
        for i, section_i in enumerate(self.current_ndf):
            print("i=", i, "section_i=", section_i)

            xi = section_i.pivot
            Ni = section_i.particles
            gammai = self.break_freq(xi)

            xip1 = self.current_ndf.section(i+1)
            xim1 = self.current_ndf.section(i-1)

            zeta = self.primary
            nu = self.secondary
            beta = self.child_number

            # calculate birth:
            birth = 0
            for k, section_k in enumerate(self.current_ndf):
                if k >= 1:
                    print("k=", k, "section_k=", section_k)

                    xk = section_k.pivot
                    Nk = section_k.particles
                    gammak = self.break_freq(xk)

                    # calculate n_i,k:
                    def func_primary(v):
                        return v**zeta * beta(v, xk)
                    Bzeta = quad(func_primary, xi, xip1)[0]
                    Bzetam1 = quad(func_primary, xim1, xi)[0]

                    def func_secondary(v):
                        return v**nu * beta(v, xk)
                    Bnu = quad(func_secondary, xi, xip1)[0]
                    Bnum1 = quad(func_secondary, xim1, xi)[0]

                    first_nom = Bzeta * xip1**nu - Bnu * xip1**zeta
                    first_denum = xi**zeta * xip1**nu - xi**nu * xip1**zeta
                    first_term =  first_nom / first_denum

                    second_nom = Bzetam1 * xim1**nu - Bnum1 * xim1**zeta
                    second_denum = xi**zeta * xim1**nu - xi**nu * xim1**zeta
                    second_term = second_nom / second_denum

                    nik = first_term + second_term

                    # calculate actual birth:
                    birth += nik * gammak * Nk

            #calculate death:
            death = gammai * Ni

            #

    def calc_aggregration(self):
        pass

    def calc_growth(self):
        pass

    def calc_nucleation(self):
        pass


def main():
    from WR.grid import Grid
    import numpy as np

    def f(v, N0=1, v0=1):
        """initial number density function.
        """
        return (N0 / v0) * (v / v0) * np.exp(-v / v0)
    ini = Grid.create_geometric_end(0, 1000, 100, 1.1, f)
    ini._plot(scale="log")

    def gamma(v):
        """breakage frequency function.
        """
        return v*v

    def beta(v1, v2):
        """
        """
        return v2/2
    method = FixedPivot(primary=0, secondary=1, initial_ndf=ini, break_freq=gamma, child_number=beta)

    method.simulate(end_time=1, steps=10)


if __name__ == "__main__":
    main()
