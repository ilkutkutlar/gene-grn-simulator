from abc import ABC, abstractmethod
from enum import Enum
from typing import List, Dict, Tuple, NamedTuple

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

    """
    Given a species name, returns the species which are regulating it.
    """
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
    def __init__(self, left: List[str], right: List[str]):
        self.left = left
        self.right = right

    @abstractmethod
    def rate_function(self, n: Network) -> float:
        pass

    def change_vector(self, n: Network) -> NamedVector:
        # fpm + MmyR -> [k1] fpm:MmyR
        # ---------------------------
        # MmyR -= k1[MmyR][fpm]     | General: if MmyR in left, then MmyR -= (k) * (left1) * (left2)
        # fpm -= k1[MmyR][fpm]
        # fpm:MmyR += k1[MmyR][fpm] | General: if fpm:MmyR in right, then fpm:MmyR += (k_) * (left1) * (left2)

        # fpm:MmyR -> [k_1] fpm + MmyR
        # ---------------------------
        # MmyR += k_1 [fpm:MmyR]
        # fpm += k_1 [fpm:MmyR]
        # fpm:MmyR -= k_1[fpm:MmyR]

        change: Dict[str, float] = dict()

        for x in n.species:

            if x not in change:
                change[x] = 0

            if x in self.left:
                change[x] -= self.rate_function(n)

            if x in self.right:
                change[x] += self.rate_function(n)

        return change

    def __str__(self) -> str:
        return "Reaction"


class TranscriptionReaction(Reaction):
    #  -> mRNA

    def __init__(self, trans_rate: float,
                 kd: float, hill_coeff: float,
                 left: List[str], right: List[str]):
        super().__init__(left, right)
        self.trans_rate = trans_rate
        self.kd = kd
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

    @staticmethod
    # [TF]^n / (Kd + [TF]^n)
    def _hill_activator_(tf: float, n: float, kd: float):
        a = pow(tf, n)
        b = kd + pow(tf, n)
        return a / b

    @staticmethod
    # 1 / (1 + ([TF] / Kd) ^ n)
    def _hill_repressor_(tf: float, n: float, kd: float):
        c = 1 + pow(tf / kd, n)
        return 1 / c

    def rate_function(self, n: Network) -> float:
        # Protein regulates mRNA
        regulations = n.get_inner_regulation(self.right[0])
        the_regulation = regulations[0] if regulations else None

        if regulations:
            regulator_concent = n.species[the_regulation.from_gene]
            if the_regulation.reg_type == RegType.ACTIVATION:
                h = self._hill_activator_(regulator_concent, self.hill_coeff, self.kd)
            else:
                h = self._hill_repressor_(regulator_concent, self.hill_coeff, self.kd)
            return self.trans_rate * h
        else:
            return self.trans_rate

    def __str__(self) -> str:
        return "Transcription: " + self.left[0] + " -> " + self.right[0]


class TranslationReaction(Reaction):
    # mRNA -> Protein

    def __init__(self, translation_rate: float,
                 left: List[str], right: List[str]):
        super().__init__(left, right)
        self.translation_rate = translation_rate

    def rate_function(self, n: Network) -> float:
        return self.translation_rate * n.species[self.left[0]]

    def __str__(self) -> str:
        return "Translation: " + self.left[0] + " -> " + self.right[0]


class MrnaDegradationReaction(Reaction):
    # mRNA ->

    def __init__(self, decay_rate: float,
                 left: List[str], right: List[str]):
        super().__init__(left, right)
        self.decay_rate = decay_rate

    def rate_function(self, n: Network) -> float:
        return self.decay_rate * n.species[self.left[0]]

    def __str__(self) -> str:
        return "mRNA Degradation: " + self.left[0] + " -> " + self.right[0]


class ProteinDegradationReaction(Reaction):
    # Protein ->

    def __init__(self, decay_rate: float,
                 left: List[str], right: List[str]):
        super().__init__(left, right)
        self.decay_rate = decay_rate

    def rate_function(self, n: Network) -> float:
        return self.decay_rate * n.species[self.left[0]]

    def __str__(self) -> str:
        return "Protein Degradation: " + self.left[0] + " -> " + self.right[0]


class CustomReaction(Reaction):
    counter = 0

    def __init__(self, rate_function_ast: str,
                 left: List[str], right: List[str]):
        super().__init__(left, right)
        self.rate_function_ast = rate_function_ast

    def rate_function(self, n: Network) -> float:
        return helper.evaluate_ast_string(self.rate_function_ast,
                                          n.symbols, species=n.species)

    def __str__(self) -> str:
        return "Reaction (" + self.left[0] + " -> " + self.right[0] + "): " + self.rate_function_ast
