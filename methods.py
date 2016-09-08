#!/usr/bin/env python3

"""
module implementing the fixed pivot technique in a class.

The class handles the main calculation using the fixed pivot technique for
breakage, aggregation, nucleation and growth. The different kernel functions
can be arbitrarily defined by the user.
"""

from numpy import linspace
from copy import deepcopy
from scipy.integrate import quad

from WR.util import zero, kronecker


class FixedPivot:

    def __init__(self, initial, primary=0, secondary=1,
                 bre=True, bre_freq=zero, child=zero,
                 agg=True, agg_freq=zero,
                 gro=True, gro_rate=zero,
                 nuc=True, nuc_rate=zero):
        """Initializer.

        :param initial: initial discrete number density function.
        :param bre: flag for breakage.
        :param bre_freq: breakage frequency function.
        :param child: child number function.
        :param agg: flag for aggregation.
        :param agg_freq: aggregation frequency function.
        :param gro: flag for growth.
        :param gro_rate: growth rate function.
        :param nuc: flag for nucleation.
        :param nuc_rate: nucleation rate function.
        :param primary: primary preserved moment (default: 0 = numbers).
        :param secondary: secondary preserved moment (default: 1 = mass).
        """
        self._initial = initial
        self._current = deepcopy(initial)
        self._previous = deepcopy(initial)

        self._primary = primary  # primary preserved moment "zeta"
        self._secondary = secondary  # secondary preserved moment "nu"

        self._breakage = bre  # flag for toggling breakage
        self._bre_freq = bre_freq  # breakage frequency function "Gamma"
        self._child = child  # number of child particles formed by breakage

        self._aggregation = agg  # flag for toggling aggregation
        self._agg_freq = agg_freq  # aggregation frequency function "Q"

        self._growth = gro  # flag for toggling growth
        self._gro_rate = gro_rate  # growth rate function

        self._nucleation = nuc  # flag for toggling nucleation
        self._nuc_rate = nuc_rate  # nucleation rate function

        self.results = {}

    # TODO: find a way to handle instability issues (implicit?!?)
    def simulate(self, end_time, steps, start_time=0, write_every=None):
        """Run simulation on initial number density function.

        :param end_time: end time of the simulation.
        :param steps: number of time steps to perform.
        :param start_time: simulation starting time (default: 0).
        :param write_every: steps betweens intermediate save (default: None).
        """
        assert start_time < end_time
        assert steps > 0

        if write_every is None:  # don't write intermediate results
            write_every = steps

        # create time stepping grid and iterate through time steps:
        times, step = linspace(start_time, end_time, steps, retstep=True)
        for counter, time in enumerate(times):
            if counter == 0:  # "zeroth" time step -> save initial ndf:
                self.results[start_time] = deepcopy(self._initial)
            else:  # calculate time step:
                print("step=", counter, "time=", time)

                # save last result for comparison:
                self._previous = deepcopy(self._current)

                # calculate individial terms of the PBE:
                if self._breakage:
                    self._calc_breakage(step)
                if self._aggregation:
                    self._calc_aggregration(step)
                if self._growth:
                    self._calc_growth(step)
                if self._nucleation:
                    self._calc_nucleation(step)

                # write intermediate result:
                if counter % write_every == 0:
                    self.results[time] = deepcopy(self._current)

        # write final result:
        if end_time not in self.results:
            self.results[end_time] = deepcopy(self._current)

    # TODO: find a method for handling boundary sections!!!
    def _calc_breakage(self, step):
        """calculate breakage.
        """
        zeta = self._primary
        nu = self._secondary
        beta = self._child

        for i, section_i in enumerate(self._previous):
            #print("i=", i, "section_i=", section_i)
            # use min/max boundaries as pivots?
            if i == 0 or i == len(self._previous) - 1:
                continue

            xi = section_i.pivot
            Ni = section_i.particles
            gammai = self._bre_freq(xi)

            xip1 = self._previous.section(i + 1).pivot
            xim1 = self._previous.section(i - 1).pivot

            # calculate birth:
            birth = 0
            for k, section_k in enumerate(self._previous):
                if k >= i:
                    #print("k=", k, "section_k=", section_k)

                    xk = section_k.pivot
                    Nk = section_k.particles
                    gammak = self._bre_freq(xk)

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

                    nik = first_term + second_term

                    # calculate actual birth:
                    birth += nik * gammak * Nk

            # calculate death:
            death = gammai * Ni

            # calculate new NDF:
            old_particles = self._previous.section(i).particles
            new_particles = old_particles + step * (birth - death)
            if new_particles < 0:  # keep values from getting < 0:
                #sys.stderr.write("BRE: i={} particles < 0!\n".format(i))
                self._current.section(i).particles = 0
            else:
                self._current.section(i).particles = new_particles

    def _calc_aggregration(self, step):
        """calculate aggregation.
        """
        zeta = self._primary
        nu = self._secondary

        for i, section_i in enumerate(self._previous):
            #print("i=", i, "section_i=", section_i)
            if i == 0 or i == len(self._previous) - 1:
                continue

            xi = section_i.pivot
            Ni = section_i.particles
            gammai = self._bre_freq(xi)

            xip1 = self._previous.section(i + 1).pivot
            xim1 = self._previous.section(i - 1).pivot

            # calculate birth:
            birth = 0
            for j, section_j in enumerate(self._previous):
                #print("j=", j, "section_j=", section_j)
                xj = section_j.pivot
                Nj = section_j.particles

                for k, section_k in enumerate(self._previous):
                    #print("k=", k, "section_k=", section_k)
                    xk = section_k.pivot
                    Nk = section_k.particles

                    v = xj + xk
                    if j >= k:
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

                            djk = kronecker(j, k)
                            Qjk = self._agg_freq(xj, xk)

                            birth += (1 - 0.5 * djk) * eta * Qjk * Nj * Nk

            # calculate death:
            death = 0
            for k, section_k in enumerate(self._previous):
                xk = section_k.pivot
                Nk = section_k.particles
                Qik = self._agg_freq(xi, xk)
                death += Qik * Nk
            death *= Ni

            # calculate new NDF:
            old_particles = self._previous.section(i).particles
            new_particles = old_particles + step * (birth - death)
            if new_particles < 0:  # keep values from getting < 0:
                # sys.stderr.write("BRE: i={} particles < 0!\n".format(i))
                self._current.section(i).particles = 0
            else:
                self._current.section(i).particles = new_particles

    def _calc_growth(self, step):
        """calculate growth.
        """
        for i, section_i in enumerate(self._previous):
            #print("i=", i, "section_i=", section_i)
            if i == 0 or i == len(self._previous) - 1:
                continue

            vi = section_i.start
            vip1 = section_i.end
            Gvi = self._gro_rate(vi)
            Gvip1 = self._gro_rate(vip1)

            ni = section_i.pivot
            nim1 = self._previous.section(i - 1).particle_density
            nip1 = self._previous.section(i + 1).particle_density

            # MARCHAL calculation of n(vi) and n(vi+1):
            #'''
            nvi = 0.5 * (nim1 + ni)
            nvip1 = 0.5 * (ni + nip1)
            #'''

            # TODO: test and compare to MARCHAL
            # DAVID calculation of n(vi) and n(vi+1):
            '''
            dvi = section_i.size
            dvim1 = self._previous.section(i - 1).size
            dvip1 = self._previous.section(i + 1).size

            nvi = (dvim1 * ni + dvi * nim1) / (dvim1 + dvi)
            nvip1 = (dvi * nip1 + dvip1 * ni) / (dvi + dvip1)
            '''

            # calculate birth:
            birth = Gvi * nvi

            # calculate birth:
            death = Gvip1 * nvip1

            # calculate new NDF:
            old_particles = self._previous.section(i).particles
            new_particles = old_particles + step * (birth - death)
            if new_particles < 0:  # keep values from getting < 0:
                # sys.stderr.write("GRO: i={} particles < 0!\n".format(i))
                self._current.section(i).particles = 0
            else:
                self._current.section(i).particles = new_particles

    def _calc_nucleation(self, step):
        """calculate nucleation.
        """
        for i, section_i in enumerate(self._previous):
            #print("i=", i, "section_i=", section_i)
            vi = section_i.start
            vip1 = section_i.end

            # calculate birth:
            birth = quad(self._nuc_rate, vi, vip1)[0]

            # calculate new NDF:
            old_particles = self._previous.section(i).particles
            new_particles = old_particles + step * birth
            if new_particles < 0:  # keep values from getting < 0:
                # sys.stderr.write("NUC: i={} particles < 0!\n".format(i))
                self._current.section(i).particles = 0
            else:
                self._current.section(i).particles = new_particles
