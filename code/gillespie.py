from math import *
from random import *
from typing import List, Tuple, Dict

import matplotlib.pyplot as plt

from helper import apply_change_vector
from models import NamedVector, Network, SimulationSettings


class GillespieSimulator:
    def __init__(self, net: Network, sim: SimulationSettings):
        self.net = net
        self.sim = sim

    def calculate_r0(self, net: Network):
        r0: float = 0
        for reaction in self.net.reactions:
            r0 += reaction.rate_function(net)

        return r0

    """
    :arg s2     -> Random float [0, 1]
    :arg r0     -> Total of reaction rates
    :arg net    -> The network

    :returns a NamedVector representing the
    change vector of the reaction chosen randomly,
    which will happen next
    """

    def pick_reaction(self, s2: float, r0: float,
                      net: Network) -> NamedVector:
        propensity_ranges = []

        for reaction in self.net.reactions:
            """Propensity of the reaction"""
            propensity = reaction.rate_function(net) / r0

            if propensity_ranges:
                m_range = (propensity_ranges[len(propensity_ranges) - 1])[0] + propensity
                change = reaction.change_vector(self.net)
                propensity_ranges.append((m_range, change))
            else:
                change = reaction.change_vector(self.net)
                propensity_ranges.append((propensity, change))

        for y in propensity_ranges:
            if s2 < y[0]:
                return y[1]

    def simulate(self) -> List[Tuple]:
        t: float = self.sim.start_time
        results: List[Tuple] = []

        while t <= self.sim.end_time:
            r0: float = self.calculate_r0(self.net)
            s1: float = random()  # To pick time
            s2: float = random()  # To pick reaction

            # Advance time
            theta = (1 / r0) * log(1 / s1, e)
            t = t + theta

            vj = self.pick_reaction(s2, r0, self.net)

            if vj:
                self.net.species = apply_change_vector(self.net.species, vj)

            results.append((t, self.net.species))

        return results

    def visualise(self, results):
        # plot results
        plt.figure()

        times: List[float] = []
        plottings: Dict[str, List[float]] = {}

        for species in self.sim.plotted_species:
            plottings[species[1]] = []

        for x in results:
            times.append(x[0])

            for species in self.sim.plotted_species:
                plottings[species[1]].append(x[1][species[1]])

        for species in self.sim.plotted_species:
            plt.plot(times, plottings[species[1]], label=species[0])

        plt.xlabel(self.sim.x_label)
        plt.ylabel(self.sim.y_label)
        plt.legend(loc=0)
        plt.title(self.sim.title)

        plt.draw()
        plt.show()
