from math import *
from random import *


class GillespieSimulator:
    n0 = 0
    t0 = 0

    r = [1, 0.01]
    v = [1, 0.01]

    def simulate(self, sim_time):
        t = self.t0
        n = self.n0

        results = list()
        while t <= sim_time:
            r0 = sum(self.r)

            s1: float = random(0, 1.00)  # To pick time
            s2: float = random(0, 1.01)  # To pick reaction

            # Advance time
            theta = (1 / r0) * log(1 / s1, e)
            t = t + theta

            reaction_chosen = 0 if s2 < 1.00 else 1

            # Reaction occurs
            n = n + self.v[reaction_chosen]
            results.append(n)

        return results


# def gillespie2():
#     n0 = 0
#     v1 = 1
#     v2 = 0.01
#
#     n = n0
#     t = 0
#
#     r1 = 1
#     r2 = 0.01
#
#     results = list()
#
#     for x in range(0, 1000):
#         s1: float = random()
#         s2: float = random()
#
#         r0 = r1 + r2
#         theta = (1 / r0) * log(1 / s1, e)
#         t = t + theta
#
#         candidates = []
#         for y in [v1, v2]:
#             if y > s2 * r0:
#                 candidates.append(y)
#
#         if candidates:
#             m_min = min(candidates)
#             if m_min == 0.01:
#                 n = n - 0.01
#             else:
#                 n = n + 1
#         results.append(n)
#
#     return results


print(GillespieSimulator().simulate(100))
