from typing import List, Tuple, Dict

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

from models import Network, Cassette, Regulation


class Simulation:
    network: Network

    def __init__(self, network: Network):
        self.network = network

    # d[mRNA]/dt = a_m*H([TF]) - b_m * [mRNA]
    # a_m           -> transcription rate
    # H([TF])       -> Hill equation
    # a_m*H([TF])   -> regulated promoter strength (rps)
    # [mRNA]        -> mRNA concentration (mrna)
    # b_m           -> mRNA degradation rate (b_m)
    def delta_mrna(self, b_m: float, mrna: float, rps: float):
        return rps - b_m * mrna

    # d[protein]/dt = a_p*[mRNA] - b_p * [protein]
    # [protein] -> protein concentration (p)
    # [mRNA]    -> mRNA concentration (m)
    # a_p       -> translation rate (a_p)
    # b_p       -> protein degradation rate (b_p)
    def delta_protein(self, a_p: float, b_p: float, p: float, m: float):
        return a_p*m - b_p * p

    def dy_dt(self, y: List[float], t: int, n: Network, z):
        # Identifier: (mRNA concent., Protein concent.)
        concentrations: Dict[str, Tuple[float, float]] = dict()

        # y is such that the lower half contains mRNA concentrations
        # and the upper half contains the corresponding proteins
        # e.g.
        # y[0] = TetR mRNA
        # y[1] = LacI mRNA
        # y[2] = TetR Protein
        # y[3] = LacI Protein

        # The order of mRNA/protein is the same in y as network's genes list.
        half_y: int = len(y) // 2
        for x in range(0, half_y):
            gene: Cassette = n.genes[x]
            concentrations[gene.identifier] = (y[x], y[x + half_y])

        genes: List[Cassette] = n.genes

        new_mrna: List[float] = list()
        new_protein: List[float] = list()

        for x in range(0, len(concentrations)):
            gene: Cassette = genes[x]

            # Identifier of the gene and the regulation
            regulator: Tuple[str, Regulation] = n.get_regulators(gene.identifier)[0]
            regulator_concentration = concentrations[regulator[0]][1]

            if regulator[1] == Regulation.ACTIVATION:
                rps: float = n.promoter_strength_activated(regulator_concentration,
                                                           40, gene.promoter, 2)
            else:
                rps: float = n.promoter_strength_repressed(regulator_concentration,
                                                           40, gene.promoter, 2)

            mrna_degradation: float = gene.codes_for[0].degradation
            mrna_concentration: float = concentrations[gene.identifier][0]

            delta_mrna = self.delta_mrna(mrna_degradation,
                                         mrna_concentration,
                                         rps)

            protein_degradation = gene.codes_for[0].translates_into.degradation
            protein_translation_rate = gene.codes_for[0].translates_into.translation_rate
            protein_concentration = concentrations[gene.identifier][1]

            delta_protein = self.delta_protein(protein_degradation,
                                               protein_translation_rate,
                                               protein_concentration,
                                               mrna_concentration)
            new_mrna.append(delta_mrna)
            new_protein.append(delta_protein)

        new_y: List[float] = list()
        for x in new_mrna:
            new_y.append(x)
        for x in new_protein:
            new_y.append(x)

        return new_y

    def simulate(self):
        # Initial values for proteins and mRNAs
        # Values taken from the original paper
        m_lacI0 = self.network.mrna_init["laci"]
        m_tetR0 = self.network.mrna_init["tetr"]
        m_cl0 = self.network.mrna_init["cl"]

        p_lacI0 = self.network.protein_init["laci"]
        p_tetR0 = self.network.protein_init["tetr"]
        p_cl0 = self.network.protein_init["cl"]



        # Initial state
        y0 = [m_lacI0, m_tetR0, m_cl0, p_lacI0, p_tetR0, p_cl0]

        # time grid -> The time space for which a graph will be drawn
        t = np.linspace(0, 40, 10000)

        # solve the ODEs
        soln = odeint(self.dy_dt, y0, t,
                      args=(self.network, 0))

        m_lacI = soln[:, 0]
        m_tetR = soln[:, 1]
        m_cl = soln[:, 2]
        p_lacI = soln[:, 3]
        p_tetR = soln[:, 4]
        p_cl = soln[:, 5]

        # plot results
        plt.figure()

        plt.plot(t, m_lacI, label='mLacI')
        plt.plot(t, m_tetR, label='mTetR')
        plt.plot(t, m_cl, label='mCl')

        plt.xlabel('Time')
        plt.ylabel('mRNA')
        plt.title('mRNA amounts')
        plt.legend(loc=0)

        plt.draw()
        plt.show()
