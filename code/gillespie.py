from math import *
from random import *
from typing import Callable, List, Dict, Tuple

import matplotlib.pyplot as plt

from helper import *

Vector = List[float]
ReactionRateFunction = Callable[[Vector], float]
ChangeVectorFunction = Callable[[Vector], Vector]
DeltaMRNA = Callable[[float, float, float], float]


class GillespieSimulator:

    # sim_time: float in seconds
    # Reaction rates: [ r(1) = 1         ]
    # Change vectors: [ v1   = [0, 1, 0] ]
    def __init__(self, reaction_rates: List[ReactionRateFunction],
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

    def pick_reaction(self, s2, r0, n) -> ChangeVectorFunction:
        # candidates: List[int] = []
        #
        # for j in self.r:
        #     if j(n) > s2 * r0:
        #         candidates.append(self.r.index(j))
        # if candidates:
        #     return self.v[min(candidates)]
        # else:
        #     return self.v[0]

        range_p = []
        for x in self.r:
            prob = x(n) / r0
            if range_p:
                m_range = (range_p[len(range_p) - 1])[0] + prob
                change = self.v[self.r.index(x)]
                range_p.append(
                    (m_range,
                     change))
            else:
                change = self.v[self.r.index(x)]
                range_p.append(
                    (prob,
                     change))

        for y in range_p:
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

    # print(results_float)
    g.visualise(results_float)

# main()
