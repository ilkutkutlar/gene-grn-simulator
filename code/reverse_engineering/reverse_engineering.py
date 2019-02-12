import random
from datetime import date
from math import e, ceil
from typing import List, Tuple, Callable, Dict

# Constraint:
#   species name
#   Desired value: Min (val), Max (val), between (val1, val2)
#   Between times: t1, t2
import numpy as np

from models.network import Network
from models.simulation_settings import SimulationSettings
from simulation.gillespie_simulator import SimulationResults
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


class Constraint:
    species: str
    value_constraint: Callable[[float], bool]
    time_period: Tuple[float, float]

    """
    :param str species: The name of the species for which the constraint applies
    :param Callable[[float], float] value_constraint: A function which, given the current value, evaluates how
        close/far it is from the desired value e.g. f = lambda v: v - 100, places a constraint for value
        to be below 100. If v = 120, then this will give a result which is > 0, indicating that the constraint
        has not been satisfied.
    :param Tuple[float, float] time_period: defines the time period for which the given value constraint
        be satisfied
    """

    def __init__(self, species: str, value_constraint: Callable[[float], float], time_period: Tuple[float, float]):
        self.species = species
        self.value_constraint = value_constraint
        self.time_period = time_period


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


def annealing(net: Network, sim: SimulationSettings,
              mutables: Dict[str, Tuple[float, float]],
              constraints: List[Constraint],
              schedule: Callable[[float], float]):

    # Current is the set of values the mutable variables will have -> dict has the value name as key, value as value
    current = {name: mutables[name][0] for name in mutables}
    ode = OdeSimulator(net, sim)


    for t in range(1, 100):
        T = schedule(t)

        if T == 0:
            return current
        else:
            my_next = current.copy()
            r = random.randrange(len(mutables))


            # Todo: Modify network
            res = StructuredResults(ode.simulate(), list(ode.net.species.keys()), sim.generate_time_space())

            delta_e = evaluate(res, constraints) - evaluate()

            # We want to minimise rather than maximise
            if delta_e < 0:
                current = my_next
            else:
                p = e**(-delta_e/T)

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
                no = ceil(1/p)*precision                  # Total number of possible
                rand = random.randrange(no)

                # rand has a 'p' probability of being 0 <= rand < precision, thus this effectively
                # ensures current = next only with probability p
                if 0 <= rand < precision:
                    current = my_next