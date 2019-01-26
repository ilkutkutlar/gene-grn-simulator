from math import log, e
from typing import List, Dict, NamedTuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

from formulae import TranscriptionFormula, DegradationFormula, TranslationFormula
from gillespie_simulator import SimulationResults
from models.network import Network
from models.reaction import Reaction
from models.reg_type import RegType
from models.regulation import Regulation
from models.simulation_settings import SimulationSettings


class Concentration(NamedTuple):
    mRNA: float
    protein: float


class OdeSimulator:
    """
    Adds a given change vector to the network's state vector
    """

    @staticmethod
    def _apply_change_vector_(state: Dict[str, float], change: Dict[str, float]) -> Dict[str, float]:
        ret = state.copy()
        for x in state:
            ret[x] = state[x] + change[x]
        return ret

    """
    Calculate the next values of the given list of values
    
    :param List[float] y: List of values
    :param int t: Not used
    :param Network net: The Network which acts as the context for the given values
    """

    @staticmethod
    def dy_dt(y: List[float], t: int, net: Network) -> List[float]:
        new_y: List[float] = list()

        for r in net.reactions:
            vj = r.change_vector(net)
            net.species = OdeSimulator._apply_change_vector_(net.species, vj)

        for s in net.species:
            new_y.append(net.species[s])

        return new_y

    @staticmethod
    def simulate(net: Network, sim: SimulationSettings) -> SimulationResults:
        y0: list = list()

        # Build the initial state
        for key in net.species:
            y0.append(net.species[key])

        # time grid -> The time space for which the equations will be solved
        t: list = np.linspace(sim.start_time, sim.end_time, 1000)

        # solve the ODEs
        solution = odeint(OdeSimulator.dy_dt, y0, t, args=(net,))

        results: SimulationResults = []
        for x in solution:
            specs: Dict[str, float] = {}

            for y in x:
                for sp in net.species:
                    specs[sp] = y

            results.append((x, specs))

        return results

    @staticmethod
    def visualise(results: SimulationResults, sim: SimulationSettings):
        # plot results
        plt.figure()

        times: List[float] = []
        plottings: Dict[str, List[float]] = {}

        # species: (species title, species name)
        for species in sim.plotted_species:
            plottings[species[1]] = []

        for x in results:
            times.append(x[0])

            for species in sim.plotted_species:
                plottings[species[1]].append(x[1][species[1]])

        for species in sim.plotted_species:
            plt.plot(times, plottings[species[1]], label=species[0])

        plt.xlabel(sim.x_label)
        plt.ylabel(sim.y_label)
        plt.legend(loc=0)
        plt.title(sim.title)

        plt.draw()
        plt.show()

    # region legacy

    """
    Return the change in mRNA per unit time under regulation.

    :param float b_m    : mRNA degradation rate
    :param float rps    : Regulated promoter strength calculated from Hill function
    :param float mrna   : mRNA concentration

    ----------
    Value is calculated using this formula:

    d[mRNA]/dt = a_m*H([TF]) - b_m * [mRNA]
        - a_m           -> transcription rate
        - H([TF])       -> Hill equation
        - a_m*H([TF])   -> regulated promoter strength (rps)
        - [mRNA]        -> mRNA concentration (mrna)
        - b_m           -> mRNA degradation rate (b_m)
    """

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
    def calculate_rps(network: Network, gene: str, kd: float, n: int):

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

    # endregion


if __name__ == '__main__':
    # p_lacI0 = 10
    # p_tetR0 = 10
    # p_cl0 = 10

    # m_lacI0 = 100
    # m_tetR0 = 80
    # m_cl0 = 50

    # region Constants

    ps_active = 0.5  # Promoter strength (active)
    ps_repr = 5 * (10 ** -4)  # Promoter strength (repressed)
    transcription_rate_active = ps_active * 60
    transcription_rate_repr = ps_repr * 60

    # Decay-related values
    mRNA_half_life = 2
    protein_half_life = 10
    mRNA_decay_rate = log(2, e) / mRNA_half_life
    protein_decay_rate = log(2, e) / protein_half_life

    # Translation-related values
    translation_efficiency = 19.97
    translation_rate = translation_efficiency * mRNA_decay_rate

    # Other values
    hill_coeff = 2  # Hill coefficient
    Km = 40  # Activation coefficient

    # Derived values

    # Active transcription rate (rescaled?)
    alpha = transcription_rate_active * translation_efficiency \
            * protein_half_life / (log(2, e) * Km)
    # Repressed transcription rate (rescaled?)
    alpha0 = transcription_rate_repr * translation_efficiency \
             * protein_half_life / (log(2, e) * Km)
    # Translation rate
    beta = protein_decay_rate / mRNA_decay_rate

    # endregion

    species = {"laci_mrna": 0, "tetr_mrna": 20, "cl_mrna": 0,
               "laci_p": 0, "tetr_p": 0, "cl_p": 0}

    regulations = [Regulation(from_gene="cl_p", to_gene="laci_mrna", reg_type=RegType.REPRESSION),
                   Regulation(from_gene="laci_p", to_gene="tetr_mrna", reg_type=RegType.REPRESSION),
                   Regulation(from_gene="tetr_p", to_gene="cl_mrna", reg_type=RegType.REPRESSION)]

    reactions = [Reaction([], ["laci_mrna"], TranscriptionFormula(alpha, 2, 40, "laci_mrna")),
                 Reaction([], ["tetr_mrna"], TranscriptionFormula(alpha, 2, 40, "tetr_mrna")),
                 Reaction([], ["cl_mrna"], TranscriptionFormula(alpha, 2, 40, "cl_mrna")),

                 Reaction(["laci_mrna"], [""], DegradationFormula(mRNA_decay_rate, "laci_mrna")),
                 Reaction(["tetr_mrna"], [""], DegradationFormula(mRNA_decay_rate, "tetr_mrna")),
                 Reaction(["cl_mrna"], [""], DegradationFormula(mRNA_decay_rate, "cl_mrna")),

                 Reaction(["laci_mrna"], ["laci_p"], TranslationFormula(beta, "laci_mrna")),
                 Reaction(["tetr_mrna"], ["tetr_p"], TranslationFormula(beta, "tetr_mrna")),
                 Reaction(["cl_mrna"], ["cl_p"], TranslationFormula(beta, "cl_mrna")),

                 Reaction(["laci_p"], [""], DegradationFormula(protein_decay_rate, "laci_p")),
                 Reaction(["tetr_p"], [""], DegradationFormula(protein_decay_rate, "tetr_p")),
                 Reaction(["cl_p"], [""], DegradationFormula(protein_decay_rate, "cl_p"))]

    net = Network()
    net.initialise(species, reactions, regulations)

    end_time = 100
    s = SimulationSettings("Results", "Time", "Concentration", 0, end_time,
                           [("LacI Protein", "laci_p"),
                            ("TetR Protein", "tetr_p"),
                            ("Cl Protein", "cl_p")])
    OdeSimulator.visualise(OdeSimulator.simulate(net, s), s)
