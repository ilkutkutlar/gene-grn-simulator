from math import *
from random import *
from typing import Callable, List

import matplotlib.pyplot as plt

from models import Network

Vector = List[float]
ReactionRateFunction = Callable[[Vector], float]
ChangeVectorFunction = Callable[[Vector], Vector]


class GillespieSimulator:

    # sim_time: float in seconds
    # Reaction rates: [ r(1) = 1         ]
    # Change vectors: [ v1   = [0, 1, 0] ]
    def __init__(self,
                 network: Network,
                 reaction_rates: List[ReactionRateFunction],
                 change_vectors: List[ChangeVectorFunction],
                 sim_time: float, n0: Vector, t0: float):

        # e.g. r1 = mrna = ...
        # r2 = protein = ...

        self.r = []
        for gene in network.genome:
            self.r.append(gene.promoter.promoter_strength_active)

        self.r = reaction_rates
        self.v = change_vectors


        #

        for x in range(0, len(concent)):
            # The current gene for which we are calculating mRNA and Protein delta value
            gene: Cassette = genes[x]
            rps: float = self.calculate_rps(concent, gene, 40, 2)

            mrna_degradation: float = gene.codes_for[0].degradation
            mrna_concentration: float = concent[gene.identifier].mRNA
            delta_mrna = self.delta_mrna(mrna_degradation, mrna_concentration, rps)

            protein_degradation = gene.codes_for[0].protein.degradation
            protein_translation_rate = gene.codes_for[0].protein.translation_rate
            protein_concentration = concent[gene.identifier].protein
            delta_protein = self.delta_protein(protein_degradation, protein_translation_rate,
                                               protein_concentration, mrna_concentration)

            new_mrna.append(delta_mrna)
            new_protein.append(delta_protein)


        #
        self.n0 = n0
        self.t0 = t0
        self.end_time = sim_time

    def calculate_r0(self, n: Vector):
        r0: float = 0
        for x in self.r:
            r0 += x(n)

        return r0

    def pick_reaction(self, s2, r0, n) -> ChangeVectorFunction:
        candidates: List[int] = []

        for j in self.r:
            if j(n) > s2 * r0:
                candidates.append(self.r.index(j))
        if candidates:
            return self.v[min(candidates)]
        else:
            return self.v[0]

    def simulate(self) -> List[Vector]:
        t: float = self.t0
        n: Vector = self.n0
        results: List[Vector] = []

        while t <= self.end_time:
            r0: float = self.calculate_r0(n)

            s1: float = random()        # To pick time
            s2: float = uniform(0, r0)  # To pick reaction

            # Advance time
            theta = (1 / r0) * log(1 / s1, e)
            t = t + theta

            vj = self.pick_reaction(s2, r0, n)
            n = vj(n)

            results.append(n)

        return results

    def visualise(self, p):
        # plot results
        plt.figure()

        # time grid -> The time space for which a graph will be drawn
        plt.plot(p)

        plt.xlabel('Time')
        plt.ylabel('Protein')
        plt.legend(loc=0)

        plt.draw()
        plt.show()


def main():
    def r0(n: Vector) -> float: return 0.01 * n[0]

    def r1(n: Vector) -> float: return 1

    def v0(n: Vector) -> Vector: return [n[0] + (-0.01 * n[0])]

    def v1(n: Vector) -> Vector: return [n[0] + 1]

    g = GillespieSimulator([r0, r1], [v0, v1], 1000, [0], 0)

    results_float = []
    for x in g.simulate():
        results_float.append(x[0])

    print(results_float)
    g.visualise(results_float)


main()
