from math import *
from random import *
from typing import List, Tuple

import matplotlib.pyplot as plt

from gillespie_models import ChangeVectorFunction, Vector, RateFunction, NamedVector, Network


def apply_change_vector(state: NamedVector, change: NamedVector):
    ret = state.copy()
    for x in state:
        ret[x] = state[x] + change[x]
    return ret


class GillespieSimulator2:
    def __init__(self, net: Network,
                 sim_time: float, t0: float):
        self.net = net
        self.t0 = t0
        self.end_time = sim_time

    def calculate_r0(self, net: Network):
        r0: float = 0
        for reaction in self.net.reactions:
            r0 += reaction.rate_function(net)

        return r0

    """
        :arg s2 -> Random float [0, 1]
        :arg r0 -> Total of reaction rates
        :arg n  -> The network vector

        :returns a ChangeVectorFunction representing the
        change vector of the reaction chosen randomly,
        which will happen next
        """

    def pick_reaction(self, s2: float, r0: float,
                      net: Network) -> NamedVector:
        propensity_ranges = []

        for reaction in self.net.reactions:
            """Propensity of the reaction"""
            propensity = reaction.rate_function(net) / r0

            if propensity_ranges:
                m_range = (propensity_ranges[len(propensity_ranges) - 1])[0] + propensity
                change = reaction.change_vector(self.net)
                propensity_ranges.append((m_range, change))
            else:
                change = reaction.change_vector(self.net)
                propensity_ranges.append((propensity, change))

        for y in propensity_ranges:
            if s2 < y[0]:
                return y[1]

    def simulate(self) -> List[Tuple]:
        t: float = self.t0
        results: List[Tuple] = []

        while t <= self.end_time:
            r0: float = self.calculate_r0(self.net)

            s1: float = random()  # To pick time
            s2: float = random()  # To pick reaction

            # Advance time
            theta = (1 / r0) * log(1 / s1, e)
            t = t + theta

            vj = self.pick_reaction(s2, r0, self.net)

            if vj:
                self.net.species = apply_change_vector(self.net.species, vj)

            results.append((t, self.net.species))

        return results

    def visualise(self, results):
        # plot results
        plt.figure()

        times = []
        mlaci = []
        mtetr = []
        mcl = []

        for x in results:
            times.append(x[0])
            mlaci.append(x[1]["laci_p"])
            mtetr.append(x[1]["tetr_p"])
            mcl.append(x[1]["cl_p"])

        plt.plot(times, mlaci)
        plt.plot(times, mtetr)
        plt.plot(times, mcl)

        plt.xlabel('Time')
        plt.ylabel('Concentration')
        plt.legend(loc=0)

        plt.draw()
        plt.show()


class GillespieSimulator:

    # sim_time: float in seconds
    # Reaction rates: [ r(1) = 1         ]
    # Change vectors: [ v1   = [0, 1, 0] ]
    def __init__(self, reaction_rates: List[RateFunction],
                 change_vectors: List[ChangeVectorFunction],
                 sim_time: float, n0: NamedVector,
                 t0: float):
        self.r = reaction_rates
        self.v = change_vectors
        self.n0 = n0
        self.t0 = t0
        self.end_time = sim_time

    def calculate_r0(self, n: NamedVector):
        r0: float = 0
        for rate_function in self.r:
            r0 += rate_function(n)

        return r0

    """
    :arg s2 -> Random float [0, 1]
    :arg r0 -> Total of reaction rates
    :arg n  -> The network vector
    
    :returns a ChangeVectorFunction representing the
    change vector of the reaction chosen randomly,
    which will happen next
    """

    def pick_reaction(self, s2: float, r0: float,
                      n: NamedVector) -> ChangeVectorFunction:
        propensity_ranges = []

        for rate_function in self.r:
            """Propensity of the reaction"""
            propensity = rate_function(n) / r0

            if propensity_ranges:
                m_range = (propensity_ranges[len(propensity_ranges) - 1])[0] + propensity
                change = self.v[self.r.index(rate_function)]
                propensity_ranges.append((m_range, change))
            else:
                change = self.v[self.r.index(rate_function)]
                propensity_ranges.append((propensity, change))

        for y in propensity_ranges:
            if s2 < y[0]:
                return y[1]

    def simulate(self) -> List[Tuple]:
        t: float = self.t0
        n: NamedVector = self.n0
        results: List[Tuple] = []

        while t <= self.end_time:
            r0: float = self.calculate_r0(n)

            s1: float = random()  # To pick time
            s2: float = random()  # To pick reaction

            # Advance time
            theta = (1 / r0) * log(1 / s1, e)
            t = t + theta

            vj = self.pick_reaction(s2, r0, n)

            # apply_change_vector(n, vj(n))

            if vj:
                n = vj(n)

            results.append((t, n))

        return results

    def visualise(self, results):
        # plot results
        plt.figure()

        times = []
        mlaci = []
        mtetr = []
        mcl = []

        for x in results:
            times.append(x[0])
            mlaci.append(x[1][3])
            mtetr.append(x[1][4])
            mcl.append(x[1][5])

        plt.plot(times, mlaci)
        plt.plot(times, mtetr)
        plt.plot(times, mcl)

        plt.xlabel('Time')
        plt.ylabel('Concentration')
        plt.legend(loc=0)

        plt.draw()
        plt.show()


def main():
    def r0(n: Vector) -> float: return 0.01 * n[0]

    def r1(n: Vector) -> float: return 1

    def v0(n: Vector) -> Vector: return [n[0] + (-0.01 * n[0])]

    def v1(n: Vector) -> Vector: return [n[0] + 1]

    g = GillespieSimulator([r0, r1], [v0, v1], 10, [0], 0)

    results_float = []
    for x in g.simulate():
        results_float.append(x[0])

    g.visualise(results_float)
