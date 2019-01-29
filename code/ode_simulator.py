import numpy
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

    @staticmethod
    def _apply_change_vector(state: np.ndarray, change: List[float]):
        ret = state.copy()

        for x in range(0, len(change)):
            numpy.append(ret, ret.item(x) + change[x])

        return ret

    """
    Calculate the next values of the given list of values
    
    :param List[float] y: List of values
    :param int t: Not used
    :param Network net: The Network which acts as the context for the given values
    """

    @staticmethod
    def dy_dt(y: List[float], t: int, net: Network) -> List[float]:
        unpacked = dict()

        i = 0
        for x in net.species:
            unpacked[x] = y[i]
            i += 1

        for r in net.reactions:
            rate = r.rate_function(net)
            if r.left:
                for x in r.left:
                    unpacked[x] -= rate

            if r.right:
                for x in r.right:
                    unpacked[x] += rate

        new_y = list()
        for s in unpacked:
            new_y.append(unpacked[s])

        return new_y

    @staticmethod
    def simulate(net: Network, sim: SimulationSettings) -> SimulationResults:
        y0: list = list()

        # Build the initial state
        for key in net.species:
            y0.append(net.species[key])

        # time grid -> The time space for which the equations will be solved
        t: list = np.linspace(sim.start_time, sim.end_time, 100)

        # solve the ODEs
        solution = odeint(OdeSimulator.dy_dt, y0, t, args=(net,))

        results: SimulationResults = []

        i = 0

        # solution: List[List[float]], where List[float] is the list containing each species' concentration
        for x in solution:      # x: List[float]
            specs: Dict[str, float] = {}

            for y in x:         # y: float
                for sp in net.species:
                    specs[sp] = y

            results.append((i, specs))
            i += 1

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

        for x in results:       # Tuple[float, NamedVector]
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


def main():
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

                 Reaction(["laci_mrna"], [], DegradationFormula(mRNA_decay_rate, "laci_mrna")),
                 Reaction(["tetr_mrna"], [], DegradationFormula(mRNA_decay_rate, "tetr_mrna")),
                 Reaction(["cl_mrna"], [], DegradationFormula(mRNA_decay_rate, "cl_mrna")),

                 Reaction(["laci_mrna"], ["laci_p"], TranslationFormula(beta, "laci_mrna")),
                 Reaction(["tetr_mrna"], ["tetr_p"], TranslationFormula(beta, "tetr_mrna")),
                 Reaction(["cl_mrna"], ["cl_p"], TranslationFormula(beta, "cl_mrna")),

                 Reaction(["laci_p"], [], DegradationFormula(protein_decay_rate, "laci_p")),
                 Reaction(["tetr_p"], [], DegradationFormula(protein_decay_rate, "tetr_p")),
                 Reaction(["cl_p"], [], DegradationFormula(protein_decay_rate, "cl_p"))]

    net = Network()
    net.initialise(species, reactions, regulations)

    end_time = 100
    s = SimulationSettings("Results", "Time", "Concentration", 0, end_time,
                           [("LacI Protein", "laci_p"),
                            ("TetR Protein", "tetr_p"),
                            ("Cl Protein", "cl_p")])
    OdeSimulator.visualise(OdeSimulator.simulate(net, s), s)


def simpler():
    species = {"x": 0, "y": 20}
    regulations = []
    reactions = [Reaction([], ["x"], TranscriptionFormula(1, 1, 20, "x"))]

    net = Network()
    net.initialise(species, reactions, regulations)

    end_time = 10
    s = SimulationSettings("Results", "Time", "Concentration", 0, end_time,
                           [("X", "x"), ("Y", "y")])
    OdeSimulator.visualise(OdeSimulator.simulate(net, s), s)


if __name__ == '__main__':
    simpler()
    # main()
