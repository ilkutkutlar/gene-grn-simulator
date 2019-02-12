from typing import Dict

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

from models.simulation_settings import SimulationSettings
from structured_results import StructuredResults


class OdeSimulator:

    def __init__(self, net, sim: SimulationSettings):
        self.net = net
        self.sim = sim

        # time grid -> The time space for which the equations will be solved
        self.time_space = sim.generate_time_space()

    """
    Calculate the change in the values of species of the network
    
    :param List[float] y: List of values
    :param int t: Not used
    :param Network net: The Network which acts as the context for the given values
    """

    def dy_dt(self, y, t):
        changes = {s: 0 for s in self.net.species}
        unpacked = dict()

        i = 0
        for x in self.net.species:
            unpacked[x] = y[i]
            i += 1

        for r in self.net.reactions:
            rate = r.rate(unpacked)

            if r.left:
                for x in r.left:
                    changes[x] -= rate

            if r.right:
                for x in r.right:
                    changes[x] += rate

        new_y = [changes[s] for s in changes]

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
        values: Dict[str, float] = \
            StructuredResults.label_results(results, self.net.species)

        plt.figure()

        for s in self.sim.plotted_species:
            plt.plot(self.time_space, values[s], label=s)

        plt.xlabel("Time (s)")
        plt.ylabel("Concentration")
        plt.legend(loc=0)
        plt.title("Results")

        plt.draw()
        plt.show()
