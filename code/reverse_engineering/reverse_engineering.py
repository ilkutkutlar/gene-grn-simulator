import random
from collections import namedtuple
from math import e, ceil
from typing import List, Tuple, Callable, Dict, NamedTuple

# Constraint:
#   species name
#   Desired value: Min (val), Max (val), between (val1, val2)
#   Between times: t1, t2
import numpy as np

from models.network import Network
from models.simulation_settings import SimulationSettings
from reverse_engineering.constraint import Constraint
from simulation.ode_simulator import OdeSimulator
from structured_results import StructuredResults


# Inputs:
#   Mutable variables
#   Range for mutable variables
#   Constraints
#   Network
# Mutable variables - 2 types:
#   Species
#   reaction parameters: Reaction name, parameter name
#       Can be done easily since we're using eval now.
# Range:
#   Simply a list of floats

def is_satisfied(results: StructuredResults, constraints: List[Constraint]):
    for c in constraints:
        vals = results.results_between_times(c.species, c.time_period[0], c.time_period[1])
        does_not_obey = list(filter(lambda v: c.value_constraint(v) > 0, vals))
        if does_not_obey:
            return False

    return True


def evaluate(results: StructuredResults, constraints: List[Constraint]):
    total = 0

    for c in constraints:
        vals = results.results_between_times(c.species, c.time_period[0], c.time_period[1])
        not_sat = list(filter(lambda v: c.value_constraint(v) > 0, vals))

        if not_sat:
            total += np.mean(not_sat)

    return total


"""
    :param str name: the name of the parameter/species which can be mutated during the annealing process. 
    value -> (lower_bound, upper_bound, increments, reaction_name)
    where lower_bound & upper_bound define the interval of values to be tried for this parameter.
    increments define the step between lower_bound & upper_bound, e.g. l = 10, u = 11, increment = 0.5 would
    produce 10, 10.5 and 11 as the values to be tried for this parameter. reaction_name -> If the parameter
    is defined in a reaction, this is the reaction_name where the parameter is. Else, it is empty string ("")
"""

Mutable = namedtuple("Mutable", ["name", "lower_bound", "upper_bound", "increments", "reaction_name"])

"""
:param Network net: The network to modify
:param SimulationSettings sim: The settings to be used to simulate the network during reverse engineering
:param Dict[str, Tuple[float, float, float, str]] mutables: 
:param List[Constraint] constraints: The list of constraints which the network must satisfy.
:param Dict[float, float] schedule: The schedule required by the simulated annealing algorithm.
"""


def annealing(net: Network, sim: SimulationSettings,
              mutables: Dict[str, Tuple[float, float, float, str]],
              constraints: List[Constraint],
              schedule: Dict[float, float]):
    # Current is the set of values the mutable variables will have -> dict has the value name as key, value as value
    # TODO: The whole tuple thing is kind of a mess, use a class? Or a named tuple?
    current: Dict[str, Tuple[float, str]] = \
        {name: (mutables[name][0], mutables[name][3]) for name in mutables}
    ode = OdeSimulator(net, sim)
    namedtuple
    def generate_neighbour():
        nbour: Dict[str, Tuple[float, str]] = current.copy()
        r: int = random.randrange(len(mutables))
        # Choose a random mutable from the mutables list
        rand_mutable: Tuple[float, float, float, str] = list(mutables.keys())[r]
        # Increment mutable by its increment to create a new network,
        # i.e. current network's neighbour
        nbour[rand_mutable] += rand_mutable[2]

        return nbour

    for t in range(1, len(schedule)):
        T = schedule[t]

        if T == 0:
            return current
        else:
            neighbour = generate_neighbour()
            net.mutate(neighbour)

            # Todo: Modify network

            res = StructuredResults(ode.simulate(), list(ode.net.species.keys()), sim.generate_time_space())

            delta_e = evaluate(res, constraints) - evaluate()

            # We want to minimise rather than maximise
            if delta_e < 0:
                current = my_next
            else:
                p = e ** (-delta_e / T)

                # Since using integers for random generation, some precision of p, which is a float,
                # will be lost (e.g. 0.387 would give 2.58397... for 1/p, and no = 3 in this case. Thus
                # rather than the actual probability being 0.387, it will be 0.33.) To avoid this,
                # multiply 2.58397... with precision (say, 100) to get the integer 258. In this case,
                # each number has a probability of being picked 1/258 = 0.00387596, or the first
                # set of 'precision' many values have a probability of 0.387596..., which is accurate
                # to 3 decimal places.

                # TODO: use now
                random.seed(0)

                precision = 100
                no = ceil(1 / p) * precision  # Total number of possible
                rand = random.randrange(no)

                # rand has a 'p' probability of being 0 <= rand < precision, thus this effectively
                # ensures current = next only with probability p
                if 0 <= rand < precision:
                    current = my_next
