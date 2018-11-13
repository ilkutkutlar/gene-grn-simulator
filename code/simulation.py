from typing import List, Dict, NamedTuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

from helper import delta_mrna, delta_protein, calculate_rps
from models import Network, Cassette


class Concentration(NamedTuple):
    mRNA: float
    protein: float


class Simulator:
    network: Network

    def __init__(self, network: Network,
                 start, end, points,
                 mrna_init: Dict[str, float],
                 protein_init: Dict[str, float],
                 signal_init: Dict[str, float]):
        self.network = network
        self.start = start
        self.end = end
        self.points = points

        # Set the initial concentrations of mRNA, protein, and signals.
        # These will be used for the simulation.
        self.mrna_init: Dict[str, float] = mrna_init
        self.protein_init: Dict[str, float] = protein_init
        self.signal_init: Dict[str, float] = signal_init

    def parse_y(self, y):
        # Identifier: (mRNA concent., Protein concent.)
        concent: Dict[str, Concentration] = dict()

        # y is such that the lower half contains mRNA concentrations
        # and the upper half contains the corresponding proteins, e.g.
        # y[0] = TetR mRNA, y[1] = LacI mRNA
        # y[2] = TetR Protein, y[3] = LacI Protein
        # The order of mRNA/protein is the same in y as network's genes list.
        half_y: int = len(y) // 2
        for x in range(0, half_y):
            gene: Cassette = self.network.genome[x]
            concent[gene.identifier] = Concentration(mRNA=y[x], protein=y[x + half_y])
        return concent

    def dy_dt(self, y: List[float], t: int):
        # Parse y to retrieve the concentrations of mRNA and Proteins in a nice
        # and orderly format
        concent: Dict[str, Concentration] = self.parse_y(y)

        genes: List[Cassette] = self.network.genome

        new_mrna: List[float] = list()
        new_protein: List[float] = list()

        for x in range(0, len(concent)):
            # The current gene for which we are calculating mRNA and Protein delta value
            gene: Cassette = genes[x]
            rps: float = calculate_rps(self.network, concent, gene, 40, 2)

            mrna_degradation: float = gene.codes_for[0].degradation
            mrna_concentration: float = concent[gene.identifier].mRNA
            delta_mrna_val = delta_mrna(mrna_degradation, mrna_concentration, rps)

            protein_degradation = gene.codes_for[0].protein.degradation
            protein_translation_rate = gene.codes_for[0].protein.translation_rate
            protein_concentration = concent[gene.identifier].protein
            delta_protein_val = delta_protein(protein_degradation, protein_translation_rate,
                                              protein_concentration, mrna_concentration)

            new_mrna.append(delta_mrna_val)
            new_protein.append(delta_protein_val)

        new_y: List[float] = list()

        for x in new_mrna:
            new_y.append(x)
        for x in new_protein:
            new_y.append(x)

        return new_y

    def simulate(self):
        y0: list = list()

        for key in self.mrna_init:
            y0.append(self.mrna_init[key])
        for key in self.protein_init:
            y0.append(self.protein_init[key])

        # time grid -> The time space for which the equations will be solved
        t: list = np.linspace(self.start, self.end, self.points)

        # solve the ODEs
        solution = odeint(self.dy_dt, y0, t)

        return solution

    def visualise(self, solution):
        m_lacI = solution[:, 0]
        m_tetR = solution[:, 1]
        m_cl = solution[:, 2]
        p_lacI = solution[:, 3]
        p_tetR = solution[:, 4]
        p_cl = solution[:, 5]

        # plot results
        plt.figure()

        # time grid -> The time space for which a graph will be drawn
        timespace: list = np.linspace(self.start, self.end, self.points)

        plt.plot(timespace, m_lacI, label='mLacI')
        plt.plot(timespace, m_tetR, label='mTetR')
        plt.plot(timespace, m_cl, label='mCl')

        plt.xlabel('Time')
        plt.ylabel('mRNA')
        plt.title('mRNA amounts')
        plt.legend(loc=0)

        plt.draw()
        plt.show()
