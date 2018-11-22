from math import *
from random import *
from typing import List, Tuple, Dict, Any

import matplotlib.pyplot as plt

from helper import apply_change_vector
from models import NamedVector, Network, SimulationSettings


class GillespieSimulator:
    def __init__(self, net: Network, sim: SimulationSettings):
        self.net = net
        self.sim = sim

    def _calculate_r0_(self, net: Network):
        r0: float = 0
        for reaction in self.net.reactions:
            r0 += reaction.rate_function(net)

        return r0

    def _get_theta_(self, r0: float):
        s1: float = random()  # To pick time
        return (1 / r0) * log(1 / s1, e)

    def _get_next_state_(self, r0: float):
        vj = self._pick_next_reaction_(r0)
        return apply_change_vector(self.net.species, vj) if vj else self.net.species

    # List: Tuple: Item, probability
    def _pick_weighted_random_(self, items: List[Any], probs: List[float]) -> Any:
        s2: float = random()  # To pick reaction

        # This is what this does:
        # Say we have 3 reactions with propensities 0.2, 0.3 and 0.5.
        # To pick one of them in random, you first arrange these propensities
        # on a number line using their cumilative probabilities:
        # |--|---|-----|
        # |.2|.5 |1.0  |
        # Then you generate a random number between 0 and 1. Say, we pick 0.77.
        # 0.77 falls inside the block that belongs to the reaction
        # with a propensity of 0.5, so we pick that.

        cumilative: List[Tuple[Any, float]] = []

        for i in range(0, len(items)):
            prev_cumilative: float = 0 if (not cumilative) \
                else (cumilative[len(cumilative) - 1])[1]

            cumilative_prob = prev_cumilative + probs[i]
            item = items[i]
            cumilative.append((item, cumilative_prob))

        # Always returns a value
        for y in cumilative:
            if s2 < y[1]:
                return y[0]

    """
    :arg r0     -> Total of reaction rates

    :returns a NamedVector representing the
    change vector of the reaction chosen randomly,
    which will happen next
    """

    def _pick_next_reaction_(self, r0: float) -> NamedVector:

        propensities: List[float] = []
        for reaction in self.net.reactions:
            propensities.append(reaction.rate_function(self.net) / r0)

        return self._pick_weighted_random_(self.net.reactions, propensities)\
            .change_vector(self.net)

        # # This is what this does:
        # # Say we have 3 reactions with propensities 0.2, 0.3 and 0.5.
        # # To pick one of them in random, you first arrange these propensities
        # # on a number line using their cumilative probabilities:
        # # |--|---|-----|
        # # |.2|.5 |1.0  |
        # # Then you generate a random number between 0 and 1. Say, we pick 0.77.
        # # 0.77 falls inside the block that belongs to the reaction
        # # with a propensity of 0.5, so we pick that.
        #
        # s2: float = random()  # To pick reaction
        #
        # cumilative_propensities: List[Tuple[float, NamedVector]] = []
        #
        # for reaction in self.net.reactions:
        #     propensity: float = reaction.rate_function(self.net) / r0
        #
        #     prev_cumilative_prob = 0 if (not cumilative_propensities) \
        #         else (cumilative_propensities[len(cumilative_propensities) - 1])[0]
        #
        #     cumilative_prob = prev_cumilative_prob + propensity
        #     change = reaction.change_vector(self.net)
        #     cumilative_propensities.append((cumilative_prob, change))
        #
        # # Always returns a value
        # for y in cumilative_propensities:
        #     if s2 < y[0]:
        #         return y[1]

    def simulate(self) -> List[Tuple]:
        t: float = self.sim.start_time
        results: List[Tuple] = []

        while t <= self.sim.end_time:
            r0: float = self._calculate_r0_(self.net)

            # Advance time
            t = t + self._get_theta_(r0)
            # Apply one reaction chosen randomly
            self.net.species = self._get_next_state_(r0)

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
