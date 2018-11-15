from abc import ABC, abstractmethod
from typing import List, Callable, Dict

Vector = List[float]
NamedVector = Dict[str, float]
RateFunction = Callable[[NamedVector], float]


class Network:
    species: Dict[str, float]
    reactions: List  # of Reaction


class Reaction(ABC):
    left: str
    right: str

    def __init__(self, left: str, right: str):
        self.left = left
        self.right = right

    @abstractmethod
    def rate_function(self, n: Network) -> float: pass

    @abstractmethod
    def change_vector(self, n: Network) -> NamedVector: pass


class TranscriptionReaction(Reaction):
    #  -> mRNA

    def __init__(self, trans_rate: float,
                 km: int, hill_coeff: int,
                 regulators: List[str],
                 left: str, right: str):
        super().__init__(left, right)
        self.trans_rate = trans_rate
        self.km = km
        self.hill_coeff = hill_coeff
        self.regulators = regulators

    def rate_function(self, n: Network) -> float:
        regulator_concent = n.species[self.regulators[0]]
        return self.trans_rate * (1 / (1 + (pow(regulator_concent,
                                                self.hill_coeff) / self.km)))

    def change_vector(self, n: Network) -> NamedVector:
        change: Dict[str, float] = dict()
        for x in n.species:
            if x == self.right:
                change[x] = self.rate_function(n)
            else:
                change[x] = 0
        return change


class TranslationReaction(Reaction):
    # mRNA -> Protein

    def __init__(self, translation_rate: float,
                 left: str, right: str):
        super().__init__(left, right)
        self.translation_rate = translation_rate

    def rate_function(self, n: Network) -> float:
        return self.translation_rate * n.species[self.left]

    def change_vector(self, n: Network) -> NamedVector:
        change: Dict[str, float] = dict()
        for x in n.species:
            if x == self.right:
                change[x] = self.rate_function(n)
            else:
                change[x] = 0
        return change


class MrnaDegradationReaction(Reaction):
    # mRNA ->

    def __init__(self, decay_rate: float,
                 left: str, right: str):
        super().__init__(left, right)
        self.decay_rate = decay_rate

    def rate_function(self, n: Network) -> float:
        return self.decay_rate * n.species[self.left]

    def change_vector(self, n: Network) -> NamedVector:
        change: Dict[str, float] = dict()

        for x in n.species:
            if x == self.left:
                change[x] = -self.rate_function(n)
            else:
                change[x] = 0
        return change


class ProteinDegradationReaction(Reaction):
    # Protein ->

    def __init__(self, decay_rate: float,
                 left: str, right: str):
        super().__init__(left, right)
        self.decay_rate = decay_rate

    def rate_function(self, n: Network) -> float:
        return self.decay_rate * n.species[self.left]

    def change_vector(self, n: Network) -> NamedVector:
        change: Dict[str, float] = dict()
        for x in n.species:
            if x == self.left:
                change[x] = -self.rate_function(n)
            else:
                change[x] = 0
        return change
