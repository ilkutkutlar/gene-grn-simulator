from abc import ABC, abstractmethod

import helper
from models.input_gate import InputGate
from models.reg_type import RegType


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


"""
Based on an ODE model and uses the Hill equation to calculate
    the promoter strength when being regulated by a TF:

    Hill equation for repressor bindings:
    beta * ( 1 / ( 1 + ([TF]/K)^n) )

    Hill equation for activator bindings:
    beta * ([TF]^n / (K^n + [TF]^n) )

    beta    : Maximal transcription rate (promoter strength)
    [TF]    : The concentration of Transcription Factor that is regulating this promoter
    K       : Dissociation constant, the probability that the TF will dissociate from the
                binding site it is now bound to. Equal to Kb/Kf where Kf = rate of TF binding and
                Kb = rate of TF unbinding.
    n       : Hill coefficient. Assumed to be 1 by default.

    Source: An introduction to systems biology : design principles of biological circuits, Uri Alon
    Previous Source: https://link.springer.com/chapter/10.1007/978-94-017-9514-2_5
"""


class TranscriptionFormula(Formula):
    """
    :param float rate:
    :param str transcribed_species:
    """

    def __init__(self, rate,
                 transcribed_species):
        self.rate = rate
        self.transcribed_species = transcribed_species

        self.hill_coeff = None
        self.regulators = None
        self.input_gate = None

    def set_regulation(self, hill_coeff, regulators, input_gate=InputGate.NONE):
        self.hill_coeff = hill_coeff
        self.regulators = regulators
        self.input_gate = input_gate

    @staticmethod
    def _hill_activator(tf, n, k):
        """
        :param float tf: transcription factor concentration
        :param float n: Hill coefficient
        :param float k: Dissociation constant
        :returns float of regulation factor
        """

        a = pow(tf, n)
        b = pow(k, n) + pow(tf, n)
        return a / b

    @staticmethod
    def _hill_repressor(tf, n, k):
        """
        :param float tf: transcription factor concentration
        :param float n: Hill coefficient
        :param float k: Dissociation constant
        :returns float of regulation factor
        """

        c = 1 + pow(tf / k, n)
        return 1 / c

    def _and_gate(self):
        if not filter(lambda x: x.reg_type == RegType.REPRESSION,
                      self.regulators):
            return True
        else:
            return False

    def _or_gate(self):
        if filter(lambda x: x.reg_type == RegType.ACTIVATION,
                  self.regulators):
            return True
        else:
            return False

    def compute(self, state):
        # https://www.pnas.org/content/pnas/100/21/11980.full.pdf

        # Protein regulates mRNA
        the_regulation = self.regulators[0] if self.regulators else None

        if self.regulators:
            regulator_concent = state[the_regulation.from_gene]

            if the_regulation.reg_type == RegType.ACTIVATION:
                h = self._hill_activator(regulator_concent,
                                         self.hill_coeff, the_regulation.k)
            else:
                h = self._hill_repressor(regulator_concent,
                                         self.hill_coeff, the_regulation.k)

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
                self.k = mutation[m][0]


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
