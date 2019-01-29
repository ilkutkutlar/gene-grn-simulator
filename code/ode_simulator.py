from math import log, e
from typing import Dict

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

from models.formulae import TranscriptionFormula, DegradationFormula, TranslationFormula
from models.network import Network
from models.reaction import Reaction
from models.reg_type import RegType
from models.regulation import Regulation
from models.simulation_settings import SimulationSettings


class OdeSimulator:

    def __init__(self, net, sim: SimulationSettings):
        self.net = net
        self.sim = sim

        # time grid -> The time space for which the equations will be solved
        self.time_space: list = np.linspace(sim.start_time, sim.end_time, sim.precision)

    """
    Calculate the change in the values of species of the network
    
    :param List[float] y: List of values
    :param int t: Not used
    :param Network net: The Network which acts as the context for the given values
    """

    def dy_dt(self, y, t):
        changes = dict()
        unpacked = dict()

        i = 0
        for x in self.net.species:
            changes[x] = 0
            unpacked[x] = y[i]
            i += 1

        for r in self.net.reactions:
            rate = r.rate_function(unpacked)

            if r.left:
                for x in r.left:
                    changes[x] -= rate

            if r.right:
                for x in r.right:
                    changes[x] += rate

        new_y = list()
        for s in changes:
            new_y.append(changes[s])

        return new_y

    def simulate(self) -> np.ndarray:
        # Build the initial state
        y0: list = [self.net.species[key] for key in self.net.species]

        # solve the ODEs
        solution: np.ndarray = odeint(self.dy_dt, y0, self.time_space)

        return solution

    """
    :param np.ndarray results: A two dimensional NumPy array containing results
        in the format where the ith array inside 'results' has the values
        for each species at time i. 
    """

    def visualise(self, results: np.ndarray):
        values: Dict[str, float] = dict()

        # Attach species names to results
        i = 0
        for s in self.net.species:
            # Syntax meaning: a list consisting of the ith element of each list in results
            values[s] = results[:, i]
            i += 1

        plt.figure()

        for s in self.sim.plotted_species:
            plt.plot(self.time_space, values[s], label=s)

        plt.xlabel("Time")
        plt.ylabel("Concentration")
        plt.legend(loc=0)
        plt.title("Results")

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
    alpha = transcription_rate_active * translation_efficiency * protein_half_life / (log(2, e) * Km)
    # Repressed transcription rate (rescaled?)
    alpha0 = transcription_rate_repr * translation_efficiency * protein_half_life / (log(2, e) * Km)
    # Translation rate
    beta = protein_decay_rate / mRNA_decay_rate

    # endregion

    species = {"laci_mrna": 0, "tetr_mrna": 20, "cl_mrna": 0,
               "laci_p": 0, "tetr_p": 0, "cl_p": 0}

    laci_reg = Regulation(from_gene="cl_p", to_gene="laci_mrna", reg_type=RegType.REPRESSION)
    tetr_reg = Regulation(from_gene="laci_p", to_gene="tetr_mrna", reg_type=RegType.REPRESSION)
    cl_reg = Regulation(from_gene="tetr_p", to_gene="cl_mrna", reg_type=RegType.REPRESSION)

    reactions = [Reaction([], ["laci_mrna"], TranscriptionFormula(alpha, 2, 40, "laci_mrna", [laci_reg])),
                 Reaction([], ["tetr_mrna"], TranscriptionFormula(alpha, 2, 40, "tetr_mrna", [tetr_reg])),
                 Reaction([], ["cl_mrna"], TranscriptionFormula(alpha, 2, 40, "cl_mrna", [cl_reg])),

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
    net.species = species
    net.reactions = reactions

    s = SimulationSettings(0, 10 * 60, 1000, ["laci_p", "tetr_p", "cl_p"])

    ode = OdeSimulator(net, s)
    ode.visualise(ode.simulate())


def simpler():
    species = {"x": 0, "y": 20}

    reactions = [Reaction([], ["x"], TranscriptionFormula(5, 2, 40, "x", [
        Regulation(from_gene="y", to_gene="x", reg_type=RegType.REPRESSION)])),
                 Reaction(["y"], [], DegradationFormula(0.3, "y"))]

    net: Network = Network()
    net.species = species
    net.reactions = reactions

    s = SimulationSettings(0, 100, 100, ["x", "y"])

    ode = OdeSimulator(net, s)
    ode.visualise(ode.simulate())


if __name__ == '__main__':
    # simpler()
    main()
