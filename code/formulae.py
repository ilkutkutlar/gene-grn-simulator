from abc import ABC, abstractmethod
from typing import Callable, Dict

import helper
from models.network import Network
from models.reg_type import RegType


class Formula(ABC):
    @abstractmethod
    def formula_function(self) -> Callable[[Network], float]:
        pass


class TranscriptionFormula(Formula):
    def __init__(self, rate: float, hill_coeff: float, kd: float, transcribed_species: str):
        self.rate = rate
        self.hill_coeff = hill_coeff
        self.kd = kd
        self.transcribed_species = transcribed_species

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

    def formula_function(self) -> Callable[[Network], float]:
        def curried(n: Network):
            # Protein regulates mRNA
            regulations = n.get_inner_regulation(self.transcribed_species)
            the_regulation = regulations[0] if regulations else None

            if regulations:
                regulator_concent = n.species[the_regulation.from_gene]
                if the_regulation.reg_type == RegType.ACTIVATION:
                    h = self._hill_activator(regulator_concent, self.hill_coeff, self.kd)
                else:
                    h = self._hill_repressor(regulator_concent, self.hill_coeff, self.kd)
                return h * self.rate
            else:
                return self.rate

        return curried


class TranslationFormula(Formula):
    def __init__(self, rate: float, mrna_species: str):
        self.rate = rate
        self.mrna_species = mrna_species

    def formula_function(self) -> Callable[[Network], float]:
        def curried(n: Network):
            return self.rate * n.species[self.mrna_species]

        return curried


class DegradationFormula(Formula):
    def __init__(self, rate: float, decaying_species: str):
        self.rate = rate
        self.decaying_species = decaying_species

    def formula_function(self) -> Callable[[Network], float]:
        def curried(n: Network):
            return self.rate * n.species[self.decaying_species]

        return curried


class CustomFormula(Formula):
    def __init__(self, rate_function_ast: str, parameters: Dict[str, float]):
        self.rate_function_ast = rate_function_ast
        self.parameters = parameters

    def formula_function(self) -> Callable[[Network], float]:
        def curried(n: Network):
            return helper.evaluate_ast_string(self.rate_function_ast,
                                              n.symbols, species=n.species, parameters=self.parameters)

        return curried
