from typing import Tuple


class Mutable:
    range: Tuple[float, float]
    precision = float


class MutableInitialConcentration(Mutable):
    name: str


class MutableReactionParameter(Mutable):
    property_name: str
    reaction_name: str
