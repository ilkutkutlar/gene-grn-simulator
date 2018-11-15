from math import *
from random import *
from typing import Callable, List, Tuple

import matplotlib.pyplot as plt

from gillespie_models import ChangeVectorFunction, Vector, RateFunction


class GillespieSimulator:

    # sim_time: float in seconds
    # Reaction rates: [ r(1) = 1         ]
    # Change vectors: [ v1   = [0, 1, 0] ]
    def __init__(self, reaction_rates: List[RateFunction],
                 change_vectors: List[ChangeVectorFunction],
                 sim_time: float, n0: Vector,
                 t0: float):
        self.r = reaction_rates
        self.v = change_vectors
        self.n0 = n0
        self.t0 = t0
        self.end_time = sim_time

    def calculate_r0(self, n: Vector):
        r0: float = 0
        for x in self.r:
            r0 += x(n)

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
                      n: Vector) -> ChangeVectorFunction:
        propensity_ranges = []

        for reaction in self.r:
            """Propensity of the reaction"""
            propensity = reaction(n) / r0

            if propensity_ranges:
                m_range = (propensity_ranges[len(propensity_ranges) - 1])[0] + propensity
                change = self.v[self.r.index(reaction)]
                propensity_ranges.append(
                    (m_range,
                     change))
            else:
                change = self.v[self.r.index(reaction)]
                propensity_ranges.append(
                    (propensity, change))

        for y in propensity_ranges:
            if s2 < y[0]:
                return y[1]

    def simulate(self) -> List[Tuple]:
        t: float = self.t0
        n: Vector = self.n0
        results: List[Tuple] = []

        while t <= self.end_time:
            r0: float = self.calculate_r0(n)

            s1: float = random()  # To pick time
            s2: float = random()  # To pick reaction

            # Advance time
            theta = (1 / r0) * log(1 / s1, e)
            t = t + theta

            vj = self.pick_reaction(s2, r0, n)
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
