#!/usr/bin/env python3
# TODO: refactor for readability and
# TODO: revisit documentation

"""
module implementing the fixed pivot technique in a class.

The class handles the main calculation using the fixed pivot technique for
breakage, aggregation, nucleation and growth. The different kernel functions
can be arbitrarily defined by the user.
"""

import sys  # for debugging output on sys.stderr
import matplotlib.pyplot as plt
from numpy import linspace, zeros
from copy import deepcopy
from scipy.integrate import quad
from WR.functions import zero, kronecker, hstep
from WR.cmdline import Spinner


class Method:
    """Base class for different sectional methods.

    This class implements the general driver structure for a sectional method
    simulation. it handles the basic effects by toggling them on or off and
    providing the functions and kernels used in the simulation. It also offers
    some simple plotting functions for quick output of the results.
    """
    def __init__(self, initial,
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
        self._result_ndfs = {}  # result dictionary {time: resulting NDF}
        #self._result_moments = {}  # result dictionary {time: resulting moments}

    def simulate(self, end_time, steps, start_time=0, write_every=None):
        """Simulation driver method.

        :param end_time: end time of the simulation.
        :param steps: number of steps to take for the simulation.
        :param start_time: start time of the simulation (default: 0).
        :param write_every: number of stept between result output.
        :param max_order: maximum order of moments to be calculated.
        """
        assert start_time < end_time
        assert steps > 0
        if write_every is None:  # don't write intermediate results
            write_every = steps

        # create time stepping grid and iterate through time steps:
        times, step = linspace(start_time, end_time, steps + 1, retstep=True)
        name = self.__class__.__name__
        with Spinner("simulating {}".format(name), end_msg=None) as sp:
            for counter, time in enumerate(times):
                sp.message("step={}, time={}".format(counter, time))

                if counter == 0:
                    # "zeroth" time step -> save initial ndf:
                    self._result_ndfs[start_time] = deepcopy(self._initial)
                else:
                    # calculate time step:
                    self.do_time_step(step)

                    # write intermediate result:
                    if counter % write_every == 0:
                        self._result_ndfs[time] = deepcopy(self._current)
            # write final result:
            if end_time not in self._result_ndfs:
                self._result_ndfs[end_time] = deepcopy(self._current)

    def do_time_step(self, step):
        """Place holder for the actual time stepping implementation.

        :param step: time step.
        """
        raise NotImplementedError

    def moments_to_file(self, filename, max_order=2):
        """Write result moments to file.

        :param filename: target file to write moment data to.
        :return: None.
        """
        assert max_order >= 0
        times = sorted(self._result_ndfs.keys())
        with open(filename, "w") as fh:
            # write first info line:
            fh.write("# time, ")
            fh.write(", ".join(str(s) for s in range(max_order+1)))
            fh.write("\n")
            # write other lines:
            for time in times:
                ndf = self._result_ndfs[time]
                fh.write("{:.5e}, ".format(time))
                fh.write(", ".join("{:.5e}".format(ndf.moment(order)) for order in range(max_order+1)))
                fh.write("\n")

    def ndf_to_files(self, basefilename):
        """Write result ndfs to individual files.

        :param basefilename: shared basic name of all files.
        :return: None.
        """
        basename, extension = basefilename.split(".")
        times = self._result_ndfs.keys()
        for time in times:
            filename = "{base}_{time:.5e}.{ext}".format(
                base=basename, time=time, ext=extension
            )
            ndf = self._result_ndfs[time]
            ndf.to_file(filename, info="# time = {:.15e}".format(time))


    def _plot_ndfs(self):
        """Plot initial and current NDF.
        """
        # get parameter data (times)
        times = sorted(self._result_ndfs)

        # plot data:
        for time in times:
            pivots = self._result_ndfs[time].pivots()
            densities = self._result_ndfs[time].densities()
            plt.plot(pivots, densities, ".-",
                     label="ndf t={}".format(time))

        plt.xlabel("particle volume")
        plt.ylabel("number density")
        plt.xscale("log")
        plt.yscale("log")
        plt.legend(loc="best", fontsize="small")
        plt.grid()
        plt.show()

    def _plot_moments(self, max_order):
        """Plot the NDFs moments over time.

        :param max_order: maximum order moment to be plotted.
        """
        # get x-axis data (times):
        times = sorted(self._result_ndfs)

        # collect y-axis data (moments):
        moments = {}
        for order in range(max_order + 1):
            moments[order] = []
            for time in times:
                moments[order].append(self._result_ndfs[time].moment(order))

        # plot data:
        for order in range(max_order + 1):
            plt.plot(times, moments[order], ".-",
                     label="moment{}".format(order))

        plt.xlabel("time")
        plt.ylabel("moment")
        plt.legend(loc="best", fontsize="small")
        plt.grid()
        plt.show()


class FixedPivot(Method):
    """Fixed Pivot Method by S. Kumar and D. Ramkrishna
    """
    def __init__(self, initial,
                 bre=True, bre_freq=zero, child=zero,
                 agg=True, agg_freq=zero,
                 gro=True, gro_rate=zero,
                 nuc=True, nuc_rate=zero):
        """Initializer.
        """
        super().__init__(initial,
                         bre, bre_freq, child,
                         agg, agg_freq,
                         gro, gro_rate,
                         nuc, nuc_rate)

    def do_time_step(self, step):
        """Calculate ndf after current time step.

        :param step: time step.
        """
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
        beta = self._child

        for i, section_i in enumerate(self._previous):
            # FIXME: method is unclear about handling the boundary sections!
            if i == 0 or i == len(self._previous) - 1:
                continue

            xi = section_i.pivot
            Ni = section_i.particles
            Gammai = self._bre_freq(xi)

            xip1 = self._previous[i + 1].pivot
            xim1 = self._previous[i - 1].pivot
            '''# possible handling of boundary sections:
            if i == 0:
                xip1 = self._previous[i + 1].pivot
                xim1 = self._previous[i].start  # ???
            elif i == len(self._previous) - 1:
                xip1 = self._previous[i].end  # ???
                xim1 = self._previous[i - 1].pivot
            else:
                xip1 = self._previous[i + 1].pivot
                xim1 = self._previous[i - 1].pivot
            '''

            # calculate birth:
            birthi = 0
            for k, section_k in enumerate(self._previous):
                if k >= i:
                    xk = section_k.pivot
                    Nk = section_k.particles
                    Gammak = self._bre_freq(xk)

                    # calculate n_i,k:
                    '''# possible handling of boundary sections:
                    if xip1 is None:
                        first = 0
                    else:
                        def func(v):
                            return (xip1 - v) / (xip1 - xi) * beta(v, xk)
                        first = quad(func, xi, xip1)[0]

                    if xim1 is None:
                        second = 0
                    else:
                        def func(v):
                            return (v - xim1) / (xi - xim1) * beta(v, xk)
                        second = quad(func, xim1, xi)[0]
                    '''
                    if i == k:
                        first = 0
                        second = 0
                    else:
                        def integrand1(v):
                            return (xip1 - v) / (xip1 - xi) * beta(v, xk)
                        first, _ = quad(integrand1, xi, xip1)

                        def integrand2(v):
                            return (v - xim1) / (xi - xim1) * beta(v, xk)
                        second, _ = quad(integrand2, xim1, xi)
                    nik = first + second

                    # calculate actual birth:
                    birthi += nik * Gammak * Nk

            # calculate death:
            deathi = Gammai * Ni

            # calculate new NDF:
            current_particles = self._current[i].particles
            new_particles = current_particles + step * (birthi - deathi)
            if new_particles < 0:
                self._current[i].particles = 0
            else:
                self._current[i].particles = new_particles

    def _calc_aggregation(self, step):
        """calculate aggregation.
        """
        for i, section_i in enumerate(self._previous):
            # FIXME: method is unclear about handling the boundary sections!
            if i == 0 or i == len(self._previous) - 1:
                continue

            xi = section_i.pivot
            Ni = section_i.particles

            xip1 = self._previous[i + 1].pivot
            xim1 = self._previous[i - 1].pivot
            '''# possible handling of boundary sections:
            if i == 0:
                xip1 = self._previous[i + 1].pivot
                xim1 = self._previous[i].start  # ???
            elif i == len(self._previous) - 1:
                xip1 = self._previous[i].end  # ???
                xim1 = self._previous[i - 1].pivot
            else:
                xip1 = self._previous[i + 1].pivot
                xim1 = self._previous[i - 1].pivot
            '''

            # calculate birth:
            birthi = 0
            for j, section_j in enumerate(self._previous):
                xj = section_j.pivot
                Nj = section_j.particles

                for k, section_k in enumerate(self._previous):
                    xk = section_k.pivot
                    Nk = section_k.particles

                    v = xj + xk
                    if j >= k:
                        if xim1 <= v <= xip1:
                            if xi <= v <= xip1:
                                eta = (xip1 - v) / (xip1 - xi)
                            elif xim1 <= v <= xi:
                                eta = (v - xim1) / (xi - xim1)
                            else:
                                raise RuntimeError("unable to calc eta!")

                            djk = kronecker(j, k)
                            Qjk = self._agg_freq(xj, xk)

                            birthi += (1 - 0.5 * djk) * eta * Qjk * Nj * Nk

            # calculate death:
            deathi = 0
            for k, section_k in enumerate(self._previous):
                xk = section_k.pivot
                Nk = section_k.particles
                Qik = self._agg_freq(xi, xk)
                deathi += Qik * Nk
            deathi *= Ni

            # calculate new NDF:
            current_particles = self._current[i].particles
            new_particles = current_particles + step * (birthi - deathi)
            if new_particles < 0:
                self._current[i].particles = 0
            else:
                self._current[i].particles = new_particles

    def _calc_growth(self, step):
        """calculate growth.
        """
        for i, section_i in enumerate(self._previous):
            # FIXME: method is unclear about handling the boundary sections!
            if i == 0 or i == len(self._previous) - 1:
                continue

            vi = section_i.start
            vip1 = section_i.end
            Gvi = self._gro_rate(vi)
            Gvip1 = self._gro_rate(vip1)

            ni = section_i.pivot
            nim1 = self._previous[i - 1].density
            nip1 = self._previous[i + 1].density
            '''# possible handling of boundary sections:
            if i == 0:
                nim1 = None  # ???
                nip1 = self._previous[i + 1].density
            elif i == len(self._previous) - 1:
                nim1 = self._previous[i - 1].density
                nip1 = None  # ???
            else:
                nim1 = self._previous[i - 1].density
                nip1 = self._previous[i + 1].density

            # if nim1 is None: ...
            # if nip1 is None: ...
            '''

            # MARCHAL calculation of n(vi) and n(vi+1):
            nvi = 0.5 * (nim1 + ni)
            nvip1 = 0.5 * (ni + nip1)

            # calculate birth:
            birthi = Gvi * nvi

            # calculate birth:
            deathi = Gvip1 * nvip1

            # calculate new NDF:
            current_particles = self._current[i].particles
            new_particles = current_particles + step * (birthi - deathi)
            #sys.stderr.write("i={}, Bi={}, Di={}, Pi={}\n".format(i, birthi, deathi, new_particles))
            if new_particles < 0:
                self._current[i].particles = 0
            else:
                self._current[i].particles = new_particles

    def _calc_nucleation(self, step):
        """calculate nucleation.
        """
        for i, section_i in enumerate(self._previous):
            vi = section_i.start
            vip1 = section_i.end

            # calculate birth:
            birthi, _ = quad(self._nuc_rate, vi, vip1)

            # calculate new NDF:
            current_particles = self._current[i].particles
            new_particles = current_particles + step * birthi
            if new_particles < 0:
                self._current[i].particles = 0
            else:
                self._current[i].particles = new_particles


class CellAverage(Method):
    """Cell Average Technique by J. Kumar et. al.
    """
    def __init__(self, initial,
                 bre=True, bre_freq=zero, child=zero,
                 agg=True, agg_freq=zero,
                 gro=True, gro_rate=zero,
                 nuc=True, nuc_rate=zero):
        """Initializer.
        """
        super().__init__(initial,
                         bre, bre_freq, child,
                         agg, agg_freq,
                         gro, gro_rate,
                         nuc, nuc_rate)

        # nuclei size for growth modeling by aggregation:
        sec_0 = self._initial[0]
        self._x0 = (sec_0.end - sec_0.pivot) / 10  # small for safety

    def do_time_step(self, step):
        """Calculate ndf after current time step.

        :param step: time step.
        """
        # save last result for comparison:
        self._previous = deepcopy(self._current)

        # STEP 1: calculate birth and death rates:
        birth_num = zeros(len(self._previous))  # B_i's
        birth_vol = zeros(len(self._previous))  # V_i's
        death_num = zeros(len(self._previous))  # D_i's

        for i, sec_i in enumerate(self._previous):
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
                Bip1 = birth_num[i + 1]

                vi = mean_vol[i]
                vip1 = mean_vol[i + 1]

                xi = self._previous[i].pivot
                xim1 = self._previous[i].start
                xip1 = self._previous[i + 1].pivot

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

                vi = mean_vol[i]
                vim1 = mean_vol[i - 1]

                xi = self._previous[i].pivot
                xim1 = self._previous[i - 1].pivot
                xip1 = self._previous[i].end

                lamvim1 = (vim1 - xim1) / (xi - xim1)
                lamvi = (vi - xim1) / (xi - xim1)
                lapvi = (vi - xip1) / (xi - xip1)

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

                xi = self._previous[i].pivot
                xim1 = self._previous[i - 1].pivot
                xip1 = self._previous[i + 1].pivot

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
                self._current[i].particles = 0
            else:
                self._current[i].particles = new_particles

    def _calc_breakage_birth(self, i):
        """Calculate particles birthed due to breakage.

        :param i: section index.
        :return: number of particles and volume flux.
        """
        assert i >= 0
        assert i < len(self._previous)

        xi = self._previous[i].pivot
        vi = self._previous[i].start
        vip1 = self._previous[i].end

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
        """Calculate particles dying due to breakage.

        :param i: section index.
        :return: number of particles.
        """
        assert i >= 0
        assert i < len(self._previous)

        xi = self._previous[i].pivot
        Ni = self._previous[i].particles
        Si = self._bre_freq(xi)

        return Si * Ni

    def _calc_aggregation_birth(self, i):
        """Calculate particles birthed due to aggregation.

        :param i: section index.
        :return: number of particles and volume flux.
        """
        assert i >= 0
        assert i < len(self._previous)

        vi = self._previous[i].start
        vip1 = self._previous[i].end

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
        """Calculate particles dying due to aggregation.

        :param i: section index.
        :return: number of particles.
        """
        assert i >= 0
        assert i < len(self._previous)

        xi = self._previous[i].pivot
        Ni = self._previous[i].particles

        sum = 0
        for k, sec_k in enumerate(self._previous):
            xk = sec_k.pivot
            Nk = sec_k.particles
            betaik = self._agg_freq(xi, xk)
            sum += betaik * Nk

        return Ni * sum

    def _calc_growth_birth(self, i):
        """Calculate particles birthed due to growth.

        :param i: section index.
        :return: number of particles and volume flux.
        """
        assert i >= 0
        assert i < len(self._previous)

        xi = self._previous[i].pivot
        Ni = self._previous[i].particles

        num = self._gro_rate(xi) * Ni / self._x0
        vol = num * (xi + self._x0)

        return num, vol

    def _calc_growth_death(self, i):
        """Calculate particles dying due to growth.

        :param i: section index.
        :return: number of particles.
        """
        assert i >= 0
        assert i < len(self._previous)

        xi = self._previous[i].pivot
        Ni = self._previous[i].particles

        num = self._gro_rate(xi) * Ni / self._x0

        return num

    def _calc_nucleation_birth(self, i):
        """Calculate particles birthed due to nucleation.

        :param i: section index.
        :return: number of particles and volume flux.
        """
        assert i >= 0
        assert i < len(self._previous)

        vi = self._previous[i].start
        vip1 = self._previous[i].end

        num = quad(self._nuc_rate, vi, vip1)[0]

        def integrand(x):
            return x * self._nuc_rate(x)
        vol = quad(integrand, vi, vip1)[0]

        return num, vol
