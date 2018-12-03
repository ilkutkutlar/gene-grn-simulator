from math import *
from random import *
from typing import List, Tuple, Dict, Any

import matplotlib.pyplot as plt

from models import NamedVector, Network, SimulationSettings

SimulationResults = List[Tuple[float, NamedVector]]


class GillespieSimulator:
    """
    Calculates the sum of all reactions rates in the given network
    """

    @staticmethod
    def _calculate_r0_(net: Network):
        r0: float = 0
        for reaction in net.reactions:
            r0 += reaction.rate_function(net)

        return r0

    """
    Calculates the time after which the next random reaction will occur
    """

    @staticmethod
    def _get_theta_(r0: float):
        s1: float = random()  # To pick time
        return (1 / r0) * log(1 / s1, e)

    """
    Adds a given change vector to the network's state vector
    """

    @staticmethod
    def _apply_change_vector_(state: Dict[str, float], change: Dict[str, float]):
        ret = state.copy()
        for x in state:
            ret[x] = state[x] + change[x]
        return ret

    """
    Returns the next state of the network after a random reaction has occurred
    """

    @staticmethod
    def _get_next_state_(net: Network, r0: float):
        vj = GillespieSimulator._pick_next_reaction_(net, r0)
        return GillespieSimulator._apply_change_vector_(net.species, vj) if vj else net.species

    # TODO: This would probably be more efficient!
    #  Mostly memory efficiency, as we are building a list of tuples
    #  with all the reactions, etc.
    # def _get_weighted_int_(self, weights: List[float]):
    #     cumilative: List[float] = []
    #
    #     for i in range(0, len(weights)):
    #         prev_cumilative: float = 0 if (not cumilative) \
    #             else (cumilative[len(cumilative) - 1])
    #
    #         cumilative_prob = prev_cumilative + weights[i]
    #         item = items[i]
    #         cumilative.append((item, cumilative_prob))
    #
    #     Always returns a value
    # for y in cumilative:
    #     if s2 < y[1]:
    #         return y[0]

    """
    Given a list of items and their associated probabilities of being picked,
    picks an item randomly at a rate dictated by its given probability
    """

    @staticmethod
    def _pick_weighted_random_(items: List[Any], probabilities: List[float]) -> Any:
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

            cumilative_prob: float = prev_cumilative + probabilities[i]
            item: Any = items[i]
            cumilative.append((item, cumilative_prob))

        # Always returns a value
        for y in cumilative:
            if s2 < y[1]:
                return y[0]

        # In case something goes wrong, at least return something
        # return cumilative[0][0]

    """
    returns a NamedVector representing the
    change vector of the reaction chosen randomly,
    which will happen next
    """

    @staticmethod
    def _pick_next_reaction_(net: Network, r0: float) -> NamedVector:

        propensities: List[float] = []
        for reaction in net.reactions:
            propensities.append(reaction.rate_function(net) / r0)

        return GillespieSimulator._pick_weighted_random_(net.reactions, propensities) \
            .change_vector(net)

    """
    Performs a Gillespie simulation of the given network in the given
    interval (dictated by the simulation setting given) and returns
    a list of results.
    """

    @staticmethod
    def simulate(net: Network, sim: SimulationSettings) -> SimulationResults:
        t: float = sim.start_time
        results: SimulationResults = []

        while t <= sim.end_time:
            r0: float = GillespieSimulator._calculate_r0_(net)

            # Advance time
            t: float = t + GillespieSimulator._get_theta_(r0)
            # Apply one reaction chosen randomly
            net.species = GillespieSimulator._get_next_state_(net, r0)

            results.append((t, net.species))

        return results

    """
    Visualises a given set of Gillespie simulation results where
    simulation properties are dictated by the given simulation settings
    object
    """

    @staticmethod
    def visualise(results: SimulationResults, sim: SimulationSettings):
        # plot results
        plt.figure()

        times: List[float] = []
        plottings: Dict[str, List[float]] = {}

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
