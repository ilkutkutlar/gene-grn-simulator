from typing import List, Tuple, Callable
from models.network import Network
import numpy as np

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

# Constraint:
#   species name
#   Desired value: Min (val), Max (val), between (val1, val2)
#   Between times: t1, t2
from models.simulation_settings import SimulationSettings
from ode_simulator import OdeSimulator
from structured_results import StructuredResults


class Constraint:
    species: str
    value_constraint: Callable[[float], bool]
    time_period: Tuple[float, float]

    """
    :param str species: The name of the species for which the constraint applies
    :param Callable[[float], bool] value_constraint: A function which, given the current value, evaluates whether it
        obeys the desired value constraint. e.g. lambda v: v >= 100 and v <= 200, places a constraint for value
        to be between 100 and 200 inclusive.
    :param Tuple[float, float] time_period: defines the time period for which the given value constraint
        be satisfied
    """

    def __init__(self, species: str, value_constraint: Callable[[float], bool], time_period: Tuple[float, float]):
        self.species = species
        self.value_constraint = value_constraint
        self.time_period = time_period


def evaluate(results: StructuredResults,
             constraints: List[Constraint]):
    # Ok that's how min/man will work: filter all values which do NOT obey the constraints, if empty list,
    # then obeys constraints

    # results = StructuredResults(raw_results, list(net.species.keys()), sim.generate_time_space())

    for c in constraints:
        vals = results.results_between_times(c.species, c.time_period[0], c.time_period[1])
        does_not_obey = filter(lambda v: not c.value_constraint(v), vals)
        if does_not_obey:
            return False

    return True
