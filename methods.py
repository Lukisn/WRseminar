#!/usr/bin/env python3

"""
module implementing the fixed pivot technique in a class.

The class handles the main calculation using the fixed pivot technique for
breakage, aggregation, nucleation and growth. The different kernel functions
can be arbitrarily defined by the user.
"""

from numpy import linspace, zeros
from copy import deepcopy
from scipy.integrate import quad
from WR.util import zero, kronecker, hstep, Spinner


class Method:
    def __init__(self, initial, primary=0, secondary=1,
                 bre=True, bre_freq=zero, child=zero,
                 agg=True, agg_freq=zero,
                 gro=True, gro_rate=zero,
                 nuc=True, nuc_rate=zero):
        """Initializer.

        :param initial: initial discrete number density function.
        :param primary: primary preserved moment (default: 0 = numbers).
        :param secondary: secondary preserved moment (default: 1 = mass).
        :param bre: flag for breakage.
        :param bre_freq: breakage frequency function.
        :param child: child number function.
        :param agg: flag for aggregation.
        :param agg_freq: aggregation frequency function.
        :param gro: flag for growth.
        :param gro_rate: growth rate function.
        :param nuc: flag for nucleation.
        :param nuc_rate: nucleation rate function.
        """
        # NDFs:
        self._initial = initial  # initial NDF
        self._current = deepcopy(initial)  # current state of the NDF
        self._previous = deepcopy(initial)  # previously calculated NDF

        # preserved moments of the method:
        self._primary = primary  # primary preserved moment "zeta"
        self._secondary = secondary  # secondary preserved moment "nu"

        # Breakage:
        self._breakage = bre  # flag for toggling breakage
        self._bre_freq = bre_freq  # breakage frequency function "Gamma"
        self._child = child  # number of child particles formed by breakage

        # Aggregation:
        self._aggregation = agg  # flag for toggling aggregation
        self._agg_freq = agg_freq  # aggregation frequency function "Q"

        # Growth:
        self._growth = gro  # flag for toggling growth
        self._gro_rate = gro_rate  # growth rate function

        # Nucleation:
        self._nucleation = nuc  # flag for toggling nucleation
        self._nuc_rate = nuc_rate  # nucleation rate function

        # Results:
        self.result_ndfs = {}  # result dictionary {time: resulting NDF}
        self.result_moments = {}  # result dictionary {time: resulting moments}

    def simulate(self, end_time, steps, start_time=0, write_every=None,
                 max_order=None):
        """Simulation driver method.

        :param end_time: end time of the simulation.
        :param steps: number of steps to take for the simulation.
        :param start_time: start time of the simulation (default: 0).
        :param write_every: number of stept between result output.
        :param max_order: maximum order of moments to be calculated.
        """
        assert start_time < end_time
        assert steps > 0
        if max_order is not None:
            assert max_order >= 0

        if write_every is None:  # don't write intermediate results
            write_every = steps

        # create time stepping grid and iterate through time steps:
        times, step = linspace(start_time, end_time, steps + 1, retstep=True)
        with Spinner("simulating", end_msg=None) as sp:
            for counter, time in enumerate(times):
                sp.message("step={}, time={}".format(counter, time))

                if counter == 0:
                    # "zeroth" time step -> save initial ndf:
                    self.result_ndfs[start_time] = deepcopy(self._initial)
                else:
                    # calculate time step:
                    self.do_time_step(step)

                    # write intermediate result:
                    if counter % write_every == 0:
                        self.result_ndfs[time] = deepcopy(self._current)
                        if max_order is not None:
                            moments = {}
                            for order in range(max_order + 1):
                                moments[order] = self._current.moment(order)
                            self.result_moments[time] = moments

            # write final result:
            if end_time not in self.result_ndfs:
                self.result_ndfs[end_time] = deepcopy(self._current)
                if max_order is not None:
                    moments = {}
                    for order in range(max_order + 1):
                        moments[order] = self._current.moment(order)
                    self.result_moments[end_time] = moments

    def do_time_step(self, step):
        """Place holder for the actual time stepping implementation.

        :param step: time step.
        """
        raise NotImplementedError


# TODO: find a method for handling boundary sections!!!
class FixedPivot(Method):
    def __init__(self, initial, primary=0, secondary=1,
                 bre=True, bre_freq=zero, child=zero,
                 agg=True, agg_freq=zero,
                 gro=True, gro_rate=zero,
                 nuc=True, nuc_rate=zero):
        super().__init__(initial, primary, secondary,
                         bre, bre_freq, child,
                         agg, agg_freq,
                         gro, gro_rate,
                         nuc, nuc_rate)

    def do_time_step(self, step):
        # save last result for comparison:
        self._previous = deepcopy(self._current)

        # calculate individial terms of the PBE:
        if self._breakage:
            self._calc_breakage(step)
        if self._aggregation:
            self._calc_aggregation(step)
        if self._growth:
            self._calc_growth(step)
        if self._nucleation:
            self._calc_nucleation(step)

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

    def _calc_aggregation(self, step):
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
            nim1 = self._previous.section(i - 1).density
            nip1 = self._previous.section(i + 1).density

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


class CellAverage(Method):
    def __init__(self, initial, primary=0, secondary=1,
                 bre=True, bre_freq=zero, child=zero,
                 agg=True, agg_freq=zero,
                 gro=True, gro_rate=zero,
                 nuc=True, nuc_rate=zero):
        super().__init__(initial, primary, secondary,
                         bre, bre_freq, child,
                         agg, agg_freq,
                         gro, gro_rate,
                         nuc, nuc_rate)

        # nuclei size for growth modeling by aggregation:
        sec_0 = self._initial.section(0)
        self._x0 = (sec_0.end - sec_0.pivot) / 10  # small for safety

    def do_time_step(self, step):
        # save last result for comparison:
        self._previous = deepcopy(self._current)

        # STEP 1: calculate birth and death rates:
        birth_num = zeros(len(self._previous))  # B_i's
        birth_vol = zeros(len(self._previous))  # V_i's
        death_num = zeros(len(self._previous))  # D_i's

        for i, sec_i in enumerate(self._previous):
            # print("i=", i, "section=", sec_i)

            if self._breakage:
                num, vol = self._calc_breakage_birth(i)
                birth_num[i] += num
                birth_vol[i] += vol
                death_num[i] += self._calc_breakage_death(i)
            if self._aggregation:
                num, vol = self._calc_aggregation_birth(i)
                birth_num[i] += num
                birth_vol[i] += vol
                death_num[i] += self._calc_aggregation_death(i)
            if self._growth:
                num, vol = self._calc_growth_birth(i)
                birth_num[i] += num
                birth_vol[i] += vol
                death_num[i] += self._calc_growth_death(i)
            if self._nucleation:
                num, vol = self._calc_nucleation_birth(i)
                birth_num[i] += num
                birth_vol[i] += vol
                # no death due to nucleation!

        # STEP 2: compute volume averages:
        mean_vol = zeros(len(self._previous))  # v_i's
        for i, _ in enumerate(mean_vol):
            if birth_num[i] == 0:
                mean_vol[i] = 0
            else:
                mean_vol[i] = birth_vol[i] / birth_num[i]

        # STEP 3: birth modification:
        birth_num_av = zeros(len(self._previous))
        for i, _ in enumerate(birth_num_av):
            if i == 0:  # leftmost section:
                Bi = birth_num[i]
                #Bim1 = birth_num[i - 1]
                Bip1 = birth_num[i + 1]

                vi = mean_vol[i]
                #vim1 = mean_vol[i - 1]
                vip1 = mean_vol[i + 1]

                xi = self._previous.section(i).pivot
                xim1 = self._previous.section(i).start  # !!!
                xip1 = self._previous.section(i + 1).pivot

                #lamvim1 = (vim1 - xim1) / (xi - xim1)
                lamvi = (vi - xim1) / (xi - xim1)
                lapvi = (vi - xip1) / (xi - xip1)
                lapvip1 = (vip1 - xip1) / (xi - xip1)

                B1 = 0
                B2 = Bi * lamvi * hstep(xi - vi)
                B3 = Bi * lapvi * hstep(vi - xi)
                B4 = Bip1 * lapvip1 * hstep(xip1 - vip1)
            elif i == len(self._previous) - 1:  # rightmost section:
                Bi = birth_num[i]
                Bim1 = birth_num[i - 1]
                #Bip1 = birth_num[i + 1]

                vi = mean_vol[i]
                vim1 = mean_vol[i - 1]
                #vip1 = mean_vol[i + 1]

                xi = self._previous.section(i).pivot
                xim1 = self._previous.section(i - 1).pivot
                xip1 = self._previous.section(i).end  # !!!

                lamvim1 = (vim1 - xim1) / (xi - xim1)
                lamvi = (vi - xim1) / (xi - xim1)
                lapvi = (vi - xip1) / (xi - xip1)
                #lapvip1 = (vip1 - xip1) / (xi - xip1)

                B1 = Bim1 * lamvim1 * hstep(vim1 - xim1)
                B2 = Bi * lamvi * hstep(xi - vi)
                B3 = Bi * lapvi * hstep(vi - xi)
                B4 = 0
            else:  # other sections:
                Bi = birth_num[i]
                Bim1 = birth_num[i - 1]
                Bip1 = birth_num[i + 1]

                vi = mean_vol[i]
                vim1 = mean_vol[i - 1]
                vip1 = mean_vol[i + 1]

                xi = self._previous.section(i).pivot
                xim1 = self._previous.section(i - 1).pivot
                xip1 = self._previous.section(i + 1).pivot

                lamvim1 = (vim1 - xim1) / (xi - xim1)
                lamvi = (vi - xim1) / (xi - xim1)
                lapvi = (vi - xip1) / (xi - xip1)
                lapvip1 = (vip1 - xip1) / (xi - xip1)

                B1 = Bim1 * lamvim1 * hstep(vim1 - xim1)
                B2 = Bi * lamvi * hstep(xi - vi)
                B3 = Bi * lapvi * hstep(vi - xi)
                B4 = Bip1 * lapvip1 * hstep(xip1 - vip1)

            birth_num_av[i] = B1 + B2 + B3 + B4

        # STEP 4: solution of the ODEs:
        for i, sec_i in enumerate(self._previous):
            old = sec_i.particles
            new_particles = old + step * (birth_num_av[i] - death_num[i])
            if new_particles < 0:
                self._current.section(i).particles = 0
            else:
                self._current.section(i).particles = new_particles

    def _calc_breakage_birth(self, i):
        assert i >= 0
        assert i < len(self._previous)

        xi = self._previous.section(i).pivot
        vi = self._previous.section(i).start
        vip1 = self._previous.section(i).end

        num = 0
        vol = 0
        for k, sec_k in enumerate(self._previous):
            if k >= i:
                xk = sec_k.pivot
                Nk = sec_k.particles
                Sk = self._bre_freq(xk)

                if k == i:
                    pik = xi
                else:
                    pik = vip1
                assert pik >= vi

                def integrand_num(x):
                    return self._child(x, xk)

                def integrand_vol(x):
                    return x * self._child(x, xk)
                integral_num = quad(integrand_num, vi, pik)[0]
                integral_vol = quad(integrand_vol, vi, pik)[0]

                num += Nk * Sk * integral_num
                vol += Nk * Sk * integral_vol

        return num, vol

    def _calc_breakage_death(self, i):
        assert i >= 0
        assert i < len(self._previous)

        xi = self._previous.section(i).pivot
        Ni = self._previous.section(i).particles
        Si = self._bre_freq(xi)

        return Si * Ni

    def _calc_aggregation_birth(self, i):
        assert i >= 0
        assert i < len(self._previous)

        vi = self._previous.section(i).start
        vip1 = self._previous.section(i).end

        num = 0
        vol = 0
        for j, sec_j in enumerate(self._previous):
            xj = sec_j.pivot
            Nj = sec_j.particles

            for k, sec_k in enumerate(self._previous):
                xk = sec_k.pivot
                Nk = sec_k.particles

                v = xj + xk
                if j >= k and vi <= v <= vip1:

                    djk = kronecker(j, k)
                    betajk = self._agg_freq(xj, xk)

                    #print("djk=", djk, "betajk=", betajk, "Nj=", Nj, "Nk=", Nk, "v=", v)

                    num += (1 - 0.5 * djk) * betajk * Nj * Nk
                    vol += (1 - 0.5 * djk) * betajk * Nj * Nk * v

        return num, vol

    def _calc_aggregation_death(self, i):
        assert i >= 0
        assert i < len(self._previous)

        xi = self._previous.section(i).pivot
        Ni = self._previous.section(i).particles

        sum = 0
        for k, sec_k in enumerate(self._previous):
            xk = sec_k.pivot
            Nk = sec_k.particles
            betaik = self._agg_freq(xi, xk)
            sum += betaik * Nk

        return Ni * sum

    def _calc_growth_birth(self, i):
        assert i >= 0
        assert i < len(self._previous)

        xi = self._previous.section(i).pivot
        Ni = self._previous.section(i).particles

        num = self._gro_rate(xi) * Ni / self._x0
        vol = num * (xi + self._x0)

        return num, vol

    def _calc_growth_death(self, i):
        assert i >= 0
        assert i < len(self._previous)

        xi = self._previous.section(i).pivot
        Ni = self._previous.section(i).particles

        num = self._gro_rate(xi) * Ni / self._x0

        return num

    def _calc_nucleation_birth(self, i):
        assert i >= 0
        assert i < len(self._previous)

        vi = self._previous.section(i).start
        vip1 = self._previous.section(i).end

        num = quad(self._nuc_rate, vi, vip1)[0]

        def integrand(x):
            return x * self._nuc_rate(x)
        vol = quad(integrand, vi, vip1)[0]

        return num, vol
