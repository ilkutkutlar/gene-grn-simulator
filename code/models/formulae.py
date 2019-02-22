from abc import ABC, abstractmethod
from typing import Dict, List, Tuple

import helper
from models.reg_type import RegType
from models.regulation import Regulation


class Formula(ABC):
    """
    Return the result of computing the formula.
    :param Dict[str, float] Network state: key: species name, value: concentration
    :returns float of result
    """

    @abstractmethod
    def compute(self, state):
        pass

    """
    Change value of given variables in the formula
    :param Dict[str, Tuple[float, str]] mutation: dictionary of variables to mutate
    """

    @abstractmethod
    def mutate(self, mutation):
        pass


class TranscriptionFormula(Formula):
    """
    :param float rate:
    :param float hill_coeff:
    :param float kd:
    :param str transcribed_species:
    :param List[Regulation] regulations:
    """

    def __init__(self, rate, hill_coeff, kd,
                 transcribed_species, regulators):
        self.rate = rate
        self.hill_coeff = hill_coeff
        self.kd = kd
        self.transcribed_species = transcribed_species
        self.regulators = regulators

    """
    :param float tf: transcription factor concentration
    :param float n: Hill coefficient
    :param float kd: Dissociation constant
    :returns float of regulation factor
    
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
    def _hill_activator(tf, n, kd):
        a = pow(tf, n)
        b = kd + pow(tf, n)
        return a / b

    """
    :param float tf: transcription factor concentration
    :param float n: Hill coefficient
    :param float kd: Dissociation constant
    :returns float of regulation factor
    """

    @staticmethod
    def _hill_repressor(tf, n, kd):
        c = 1 + pow(tf / kd, n)
        return 1 / c

    def compute(self, state):
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

    def mutate(self, mutation):
        for m in mutation:
            if m == "rate":
                self.rate = mutation[m][0]
            elif m == "hill_coeff":
                self.hill_coeff = mutation[m][0]
            else:  # m == "kd"
                self.kd = mutation[m][0]


class TranslationFormula(Formula):
    """
    :param float rate:
    :param str mrn_species:
    """

    def __init__(self, rate, mrna_species):
        self.rate = rate
        self.mrna_species = mrna_species

    def compute(self, state):
        return self.rate * state[self.mrna_species]

    def mutate(self, mutation):
        for m in mutation:
            if m == "rate":
                self.rate = mutation[m][0]


class DegradationFormula(Formula):
    """
    :param float rate:
    :param str decaying_species:
    """

    def __init__(self, rate, decaying_species):
        self.rate = rate
        self.decaying_species = decaying_species

    def compute(self, state):
        return self.rate * state[self.decaying_species]

    def mutate(self, mutation):
        for m in mutation:
            if m == "rate":
                self.rate = mutation[m][0]


class CustomFormula(Formula):
    """
    :param str rate_function:
    :param Dict[str, float] parameters:
    :param Dict[str, float] symbols:
    """

    def __init__(self, rate_function, parameters, symbols):
        self.rate_function = rate_function
        self.symbols = symbols
        self.parameters = parameters

    def compute(self, state):
        return helper.eval_equation(self.rate_function,
                                    species=state,
                                    symbols=self.symbols,
                                    parameters=self.parameters)

    def mutate(self, mutation):
        for m in mutation:
            self.parameters.update({m: mutation[m][0]})
