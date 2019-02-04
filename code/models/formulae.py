from abc import ABC, abstractmethod
from typing import Dict, List

import helper
from models.network import Network
from models.reg_type import RegType
from models.regulation import Regulation


class Formula(ABC):
    @abstractmethod
    def compute(self, state: Dict[str, float]) -> float:
        pass


class TranscriptionFormula(Formula):
    def __init__(self, rate: float, hill_coeff: float, kd: float,
                 transcribed_species: str, regulators: List[Regulation]):
        self.rate = rate
        self.hill_coeff = hill_coeff
        self.kd = kd
        self.transcribed_species = transcribed_species
        self.regulators = regulators

    """
        Based on an ODE model and uses the Hill equation to calculate
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
        """

    @staticmethod
    def _hill_activator(tf: float, n: float, kd: float):
        a = pow(tf, n)
        b = kd + pow(tf, n)
        return a / b

    @staticmethod
    def _hill_repressor(tf: float, n: float, kd: float):
        c = 1 + pow(tf / kd, n)
        return 1 / c

    def compute(self, state: Dict[str, float]):
        # Protein regulates mRNA
        the_regulation = self.regulators[0] if self.regulators else None

        if self.regulators:
            regulator_concent = state[the_regulation.from_gene]
            if the_regulation.reg_type == RegType.ACTIVATION:
                h = self._hill_activator(regulator_concent, self.hill_coeff, self.kd)
            else:
                h = self._hill_repressor(regulator_concent, self.hill_coeff, self.kd)

            return h * self.rate
        else:
            return self.rate


class TranslationFormula(Formula):
    def __init__(self, rate: float, mrna_species: str):
        self.rate = rate
        self.mrna_species = mrna_species

    def compute(self, state: Dict[str, float]) -> float:
        return self.rate * state[self.mrna_species]


class DegradationFormula(Formula):
    def __init__(self, rate: float, decaying_species: str):
        self.rate = rate
        self.decaying_species = decaying_species

    def compute(self, state: Dict[str, float]) -> float:
        return self.rate * state[self.decaying_species]


class CustomFormula(Formula):
    def __init__(self, rate_function: str,
                 parameters: Dict[str, float], symbols: Dict[str, float]):
        self.rate_function = rate_function
        self.parameters = parameters
        self.symbols = symbols

    def compute(self, state: Dict[str, float]) -> float:
        return helper.eval_equation(self.rate_function,
                                    symbols=self.symbols,
                                    species=state,
                                    parameters=self.parameters)
