from abc import ABC, abstractmethod
from enum import Enum
from typing import List, Dict, Tuple, NamedTuple

import libsbml
from libsbml._libsbml import formulaToL3String

import helper

Vector = List[float]
NamedVector = Dict[str, float]


class RegType(Enum):
    ACTIVATION = 1
    REPRESSION = -1


class Regulation(NamedTuple):
    from_gene: str
    to_gene: str
    reg_type: RegType

    def __str__(self) -> str:
        sign = " -> " if self.reg_type == RegType.ACTIVATION else " -| "
        return "Reg: " + self.from_gene + sign + self.to_gene


class Network:
    species: Dict[str, float]
    symbols: Dict[str, float]  # Used by CustomReactions, optional
    reactions: List  # of Reaction
    regulations: List[Regulation]

    def __init__(self) -> None:
        self.species = {}
        self.reactions = []
        self.regulations = []
        self.symbols = []

    def initialise(self, species: Dict[str, float],
                   reactions: List,
                   regulations: List[Regulation]):
        self.species = species
        self.reactions = reactions
        self.regulations = regulations
        self.symbols = []

    def get_inner_regulation(self, name: str) -> List[Regulation]:
        return list(filter(lambda reg: reg.to_gene == name, self.regulations))

    def __str__(self) -> str:
        ret = "\nGlobal Parameters: \n"
        for x in self.symbols:
            ret += "    " + x + ": " + str(self.symbols[x]) + "\n"

        ret += "\nSpecies: \n"
        for x in self.species:
            ret += "    " + x + ": " + str(self.species[x]) + "\n"

        ret += "\nRegulations: \n"
        for x in self.regulations:
            ret += "    " + str(x) + "\n"

        ret += "\nReactions: \n"
        for x in self.reactions:
            ret += "    " + str(x) + "\n"

        return ret


class SimulationSettings:
    # Tuple[label, species_name]
    def __init__(self, title: str, x_label: str, y_label: str, start_time: float, end_time: float,
                 plotted_species: List[Tuple[str, str]]):
        self.plotted_species = plotted_species
        self.start_time = start_time
        self.end_time = end_time
        self.x_label = x_label
        self.y_label = y_label
        self.title = title


class Reaction(ABC):
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
                 km: float, hill_coeff: float,
                 left: str, right: str):
        super().__init__(left, right)
        self.trans_rate = trans_rate
        self.km = km
        self.hill_coeff = hill_coeff

    ''' Based on an ODE model and uses the Hill equation to calculate
        the promoter strength when being regulated by a TF:

        Hill equation for repressor bindings:
        beta * ( 1 / ( 1 + ([TF]/Kd)^n) )

        Hill equation for activator bindings:
        beta * ([TF]^n / (Kd + [TF]^n) )

        beta    : Maximal transcription rate (promoter strength)
        [TF]    : The concentration of Transcript Factor that is regulating this promoter
        Kd      : Dissociation constant, the probability that the TF will dissociate from the 
        binding site it is now bound to. Equal to Kb/Kf where Kf = rate of TF binding and 
        Kb = rate of TF unbinding.
        n       : Hill coefficient. Assumed to be 1 by default.

        Source: https://link.springer.com/chapter/10.1007/978-94-017-9514-2_5
        '''

    def rate_function(self, n: Network) -> float:
        # Protein regulates mRNA
        regulations = n.get_inner_regulation(self.right)
        the_regulation = regulations[0]

        if regulations:
            regulator_concent = n.species[the_regulation.from_gene]
            if the_regulation.reg_type == RegType.ACTIVATION:
                # Calculation arbitrarily broken down into separate
                # parts to improve readability
                a = pow(regulator_concent, self.hill_coeff)
                b = self.km + pow(regulator_concent, self.hill_coeff)
                hill_act = (a / b)
                return self.trans_rate * hill_act
            else:
                c = pow(regulator_concent, self.hill_coeff)
                d = (1 + (c / self.km))
                hill_rep = (1 / d)
                return self.trans_rate * hill_rep
        else:
            return self.trans_rate

    def change_vector(self, n: Network) -> NamedVector:
        change: Dict[str, float] = dict()
        for x in n.species:
            if x == self.right:
                change[x] = self.rate_function(n)
            else:
                change[x] = 0
        return change

    def __str__(self) -> str:
        return "Transcription: " + self.left + " -> " + self.right


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

    def __str__(self) -> str:
        return "Translation: " + self.left + " -> " + self.right


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

    def __str__(self) -> str:
        return "mRNA Degradation: " + self.left + " -> " + self.right


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

    def __str__(self) -> str:
        return "Protein Degradation: " + self.left + " -> " + self.right


class CustomReaction(Reaction):
    def __init__(self, rate_function_ast: libsbml.ASTNode,
                 left: str, right: str):
        super().__init__(left, right)
        self.rate_function_ast = rate_function_ast

    def rate_function(self, n: Network) -> float:
        return helper.evaluate_ast_node(self.rate_function_ast,
                                        n.symbols, species=n.species)

    def change_vector(self, n: Network) -> NamedVector:
        change: Dict[str, float] = dict()
        for x in n.species:
            if x == self.left:  # x is a reactant, so -ve effect
                change[x] = -self.rate_function(n)
            elif x == self.right:  # x is a product, so +ve effect
                change[x] = +self.rate_function(n)
            else:
                change[x] = 0
        return change

    def __str__(self) -> str:
        return formulaToL3String(self.rate_function_ast)
