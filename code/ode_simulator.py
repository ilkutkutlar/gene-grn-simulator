from typing import List, Dict, NamedTuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

from gillespie_simulator import SimulationResults
from models.models import SimulationSettings, Regulation
from ode.models import Network, Cassette, RegType
from models.models import Network as OtherNetwork


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

    """
    Given y, an array of floats, returns a dictionary where the values are the same
    as y and the keys represent the species to which the values belong.
    """

    @staticmethod
    def parse_y(y, net: OtherNetwork) -> Dict[str, Concentration]:
        # Identifier: (mRNA concent., Protein concent.)
        concent: Dict[str, Concentration] = dict()

        i = 0
        for x in net.species:
            concent[x] = y[i]

        return concent

    # d[mRNA]/dt = a_m*H([TF]) - b_m * [mRNA]
    # a_m           -> transcription rate
    # H([TF])       -> Hill equation
    # a_m*H([TF])   -> regulated promoter strength (rps)
    # [mRNA]        -> mRNA concentration (mrna)
    # b_m           -> mRNA degradation rate (b_m)
    @staticmethod
    def delta_mrna(b_m: float, mrna: float, rps: float):
        return rps - b_m * mrna

    # d[protein]/dt = a_p*[mRNA] - b_p * [protein]
    # [protein] -> protein concentration (p)
    # [mRNA]    -> mRNA concentration (m)
    # a_p       -> translation rate (a_p)
    # b_p       -> protein degradation rate (b_p)
    @staticmethod
    def delta_protein(a_p: float, b_p: float, p: float, m: float):
        return a_p * m - b_p * p

    # rps = regulated promoter strength.
    # Given a gene, calculates the promoter strength under regulation.
    @staticmethod
    def calculate_rps(network: OtherNetwork,
                      gene: str, kd: float, n: int):

        # Identifier of the gene and the regulation
        regulation: Regulation = network.get_inner_regulation(gene)[0]
        regulator: str = network.get_inner_regulation(gene)[0].from_gene
        regulation_type = regulation.reg_type

        regulator_concentration = network.species[regulator]


        if regulation_type == RegType.ACTIVATION:
            rps: float = network.ps_active(regulator_concentration, kd, gene.promoter, n)
        else:
            rps: float = network.ps_repressed(regulator_concentration, kd, gene.promoter, n)

        return rps

    @staticmethod
    def dy_dt(y: List[float], t: int, net: OtherNetwork):
        # Parse y to retrieve the concentrations of mRNA and Proteins in a nice
        # and orderly format
        concent: Dict[str, Concentration] = Simulator.parse_y(y, net)

        genes: List[Cassette] = self.network.genome

        new_mrna: List[float] = list()
        new_protein: List[float] = list()

        for x in range(0, len(concent)):
            # The current gene for which we are calculating mRNA and Protein delta value
            gene: Cassette = genes[x]
            rps: float = Simulator.calculate_rps(net, concent, gene, 40, 2)

            mrna_degradation: float = gene.codes_for[0].degradation
            mrna_concentration: float = concent[gene.identifier].mRNA
            delta_mrna_val = Simulator.delta_mrna(mrna_degradation, mrna_concentration, rps)

            protein_degradation = gene.codes_for[0].protein.degradation
            protein_translation_rate = gene.codes_for[0].protein.translation_rate
            protein_concentration = concent[gene.identifier].protein
            delta_protein_val = Simulator.delta_protein(protein_degradation, protein_translation_rate,
                                                        protein_concentration, mrna_concentration)

            new_mrna.append(delta_mrna_val)
            new_protein.append(delta_protein_val)

        new_y: List[float] = list()

        for x in new_mrna:
            new_y.append(x)
        for x in new_protein:
            new_y.append(x)

        return new_y

    @staticmethod
    def simulate(net: OtherNetwork, sim: SimulationSettings) -> SimulationResults:
        y0: list = list()

        # Build the initial state
        for key in net.species:
            y0.append(net.species[key])

        # time grid -> The time space for which the equations will be solved
        t: list = np.linspace(sim.start_time, sim.end_time, 1000)

        # solve the ODEs
        solution = odeint(self.dy_dt, y0, t)
        res: SimulationResults = SimulationResults()

        return solution

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
