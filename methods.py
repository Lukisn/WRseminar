#!/usr/bin/env python3

import numpy as np
from copy import copy, deepcopy
from scipy.integrate import quad


def kronecker_delta(i, j):
    if i == j:
        return 1
    else:
        return 0


def gamma(v):
    """breakage frequency function.
    """
    return v * v


def beta(v1, v2):
    """child number function.
    """
    return 2/v2


def Q(x1, x2):
    """aggregation frequency function.
    """
    return 1


def G(v):
    """growth rate for particles of size v.
    """
    return 1


def S(v):
    """rate of nucleation of particles of size v.
    """
    if v < 1:
        return 1
    else:
        return 0


# TODO: refactoring and clean up
# TODO: documentaion
class FixedPivot:

    def __init__(self, initial_ndf,
                 break_freq=gamma, child_number=beta,
                 agg_freq=Q,
                 gro_rate=G,
                 nuc_rate=S,
                 primary=0, secondary=1,
                 breakage=True, aggregation=True,
                 growth=False, nucleation=False):
        self.primary = primary  # primary preserved moment "zeta" usually zeta=0 (preservation of numbers)
        self.secondary = secondary  # secondary preserved moment "nu" usually nu=1 (preservation of mass)

        self.initial_ndf = initial_ndf
        self.current_ndf = deepcopy(initial_ndf)
        self.previous_ndf = deepcopy(initial_ndf)

        self.breakage = breakage
        self.break_freq = break_freq  # breakage frequency function "Gamma"
        self.child_number = child_number  # function for number of child particles formed due to breakage

        self.aggregation = aggregation
        self.agg_freq = agg_freq  # aggregation frequency function "Q"

        self.growth = growth
        self.gro_rate = gro_rate  # growth rate function

        self.nucleation = nucleation
        self.nuc_rate = nuc_rate  # nucleation rate function

        self.results = {}

    # TODO: implement time steppig in simulate function!?
    # TODO: find a way to handle instability issues (implicit?!?)
    def simulate(self, end_time, steps, start_time=0, write_every=None):
        assert start_time < end_time
        assert steps > 0

        if write_every is None:
            write_every = steps

        times, step = np.linspace(start_time, end_time, steps, retstep=True)
        for counter, time in enumerate(times):
            if counter == 0:  # save initial ndf:
                self.results[start_time] = deepcopy(self.initial_ndf)
            else:  # calculate time step:
                print("step=", counter, "time=", time, "step=", step)

                # calculate terms of the PBE:
                self.previous_ndf = copy(self.current_ndf)
                if self.breakage:
                    self.calc_breakage(step)
                if self.aggregation:
                    self.calc_aggregration(step)
                if self.growth:
                    self.calc_growth(step)
                if self.nucleation:
                    self.calc_nucleation(step)

                # write intermediate result:
                if counter % write_every == 0:
                    self.results[time] = deepcopy(self.current_ndf)

        # write final result:
        if end_time not in self.results:
            self.results[end_time] = deepcopy(self.current_ndf)

    # TODO: find a method for handling boundary sections!!!
    def calc_breakage(self, time_step):
        """calculate breakage.
        """
        zeta = self.primary
        nu = self.secondary
        beta = self.child_number

        for i, section_i in enumerate(self.previous_ndf):
            #print("i=", i, "section_i=", section_i)

            xi = section_i.pivot
            Ni = section_i.particles
            gammai = self.break_freq(xi)

            # use min/max boundaries as pivots?
            if i == 0 or i == len(self.previous_ndf) - 1:
                continue

            # TODO: check use of current ndf! here and elsewhere!
            xip1 = self.current_ndf.section(i + 1).pivot
            xim1 = self.current_ndf.section(i - 1).pivot

            '''
            try:
                xip1 = self.current_ndf.section(i+1).pivot
            except AssertionError:
                xip1 = None
            try:
                xim1 = self.current_ndf.section(i-1).pivot
            except AssertionError:
                xim1 = None
            '''

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

                    Bzeta = quad(func_primary, xi, xip1)[0]
                    Bnu = quad(func_secondary, xi, xip1)[0]
                    first_num = Bzeta * xip1 ** nu - Bnu * xip1 ** zeta
                    first_denom = xi ** zeta * xip1 ** nu - xi ** nu * xip1 ** zeta
                    first_term = first_num / first_denom

                    Bzetam1 = quad(func_primary, xim1, xi)[0]
                    Bnum1 = quad(func_secondary, xim1, xi)[0]
                    second_num = Bzetam1 * xim1 ** nu - Bnum1 * xim1 ** zeta
                    second_denom = xi ** zeta * xim1 ** nu - xi ** nu * xim1 ** zeta
                    second_term = second_num / second_denom

                    '''
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
                    '''

                    nik = first_term + second_term

                    # calculate actual birth:
                    birth += nik * gammak * Nk

            # calculate death:
            death = gammai * Ni

            # calculate new ndf:
            # self.current_ndf.section(i).particles += time_step * (birth - death)
            new_particles = self.current_ndf.section(i).particles + time_step * (birth - death)
            if new_particles < 0:
                self.current_ndf.section(i).particles = 0
                #raise RuntimeWarning("calculated particle number lower than 0!")
            else:
                self.current_ndf.section(i).particles = new_particles

    def calc_aggregration(self, time_step):
        """calculate aggregation.
        """
        zeta = self.primary
        nu = self.secondary

        for i, section_i in enumerate(self.previous_ndf):
            #print("i=", i, "section_i=", section_i)

            xi = section_i.pivot
            Ni = section_i.particles
            gammai = self.break_freq(xi)

            if i == 0 or i == len(self.previous_ndf) - 1:
                continue

            xip1 = self.current_ndf.section(i + 1).pivot
            xim1 = self.current_ndf.section(i - 1).pivot
            '''
            try:
                xip1 = self.current_ndf.section(i + 1).pivot
            except AssertionError:
                xip1 = None
            try:
                xim1 = self.current_ndf.section(i - 1).pivot
            except AssertionError:
                xim1 = None
            '''

            # calculate birth:
            birth = 0
            for j, section_j in enumerate(self.previous_ndf):
                #print("j=", j, "section_j=", section_j)

                xj = section_j.pivot
                Nj = section_j.particles

                for k, section_k in enumerate(self.previous_ndf):
                    #print("k=", k, "section_k=", section_k)

                    xk = section_k.pivot
                    Nk = section_k.particles

                    v = xj + xk
                    if j >= k:
                        '''
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
                        '''
                        if xim1 <= v <= xip1:
                            if xi <= v <= xip1:
                                num = v ** zeta * xip1 ** nu - v ** nu * xip1 ** zeta
                                denom = xi ** zeta * xip1 ** nu - xi ** nu * xip1 ** zeta
                                eta = num / denom
                            elif xim1 <= v <= xi:
                                num = v ** zeta * xim1 ** nu - v ** nu * xim1 ** zeta
                                denom = xi ** zeta * xim1 ** nu - xi ** nu * xim1 ** zeta
                                eta = num / denom
                            else:
                                raise RuntimeError("unable to calc eta!")

                            djk = kronecker_delta(j, k)
                            Qjk = self.agg_freq(xj, xk)

                            birth += (1 - 1/2 * djk) * eta * Qjk * Nj * Nk

            # calculate death:
            death = 0
            for k, section_k in enumerate(self.previous_ndf):
                xk = section_k.pivot
                Nk = section_k.particles
                Qik = self.agg_freq(xi, xk)
                death += Qik * Nk
            death *= Ni

            #self.current_ndf.section(i).particles += time_step * (birth - death)
            new_particles = self.current_ndf.section(i).particles + time_step * (birth - death)
            if new_particles < 0:
                self.current_ndf.section(i).particles = 0
                #raise RuntimeWarning("calculated particle number lower than 0!")
            else:
                self.current_ndf.section(i).particles = new_particles

    def calc_growth(self, time_step):
        """calculate growth.
        """
        for i, section_i in enumerate(self.previous_ndf):
            #print("i=", i, "section_i=", section_i)

            if i == 0 or i == len(self.previous_ndf) - 1:
                continue

            vi = section_i.start
            vip1 = section_i.end
            Gvi = self.gro_rate(vi)
            Gvip1 = self.gro_rate(vip1)
            ni = section_i.pivot
            nim1 = self.current_ndf.section(i - 1).particle_density
            nip1 = self.current_ndf.section(i + 1).particle_density
            nvi = 0.5 * (nim1 + ni)
            nvip1 = 0.5 * (ni + nip1)

            # calculate birth:
            birth = Gvi * nvi

            # calculate birth:
            death = Gvip1 * nvip1

            # calculate new ndf:
            new_particles = self.current_ndf.section(
                i).particles + time_step * (birth - death)
            if new_particles < 0:
                self.current_ndf.section(i).particles = 0
                # raise RuntimeWarning("calculated particle number lower than 0!")
            else:
                self.current_ndf.section(i).particles = new_particles

    def calc_nucleation(self, time_step):
        """calculate nucleation.
        """
        for i, section_i in enumerate(self.previous_ndf):
            #print("i=", i, "section_i=", section_i)

            vi = section_i.start
            vip1 = section_i.end

            # calculate birth:
            birth = quad(self.nuc_rate, vi, vip1)[0]

            # calculate new ndf.
            new_particles = self.current_ndf.section(
                i).particles + time_step * birth
            if new_particles < 0:
                self.current_ndf.section(i).particles = 0
                # raise RuntimeWarning("calculated particle number lower than 0!")
            else:
                self.current_ndf.section(i).particles = new_particles
