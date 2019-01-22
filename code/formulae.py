from typing import Callable

import helper
from models.network import Network
from models.reg_type import RegType


class Formulae:
    @staticmethod
    def translation_rate(rate: float, mrna_species: str) -> Callable[[Network], float]:
        def curried(n: Network):
            return rate * n.species[mrna_species]

        return curried

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

    @staticmethod
    def transcription_rate(rate: float, hill_coeff: float, kd: float, transcribed_species: str) -> Callable[
        [Network], float]:
        def curried(n: Network):
            # Protein regulates mRNA
            regulations = n.get_inner_regulation(transcribed_species)
            the_regulation = regulations[0] if regulations else None

            if regulations:
                regulator_concent = n.species[the_regulation.from_gene]
                if the_regulation.reg_type == RegType.ACTIVATION:
                    h = Formulae._hill_activator(regulator_concent, hill_coeff, kd)
                else:
                    h = Formulae._hill_repressor(regulator_concent, hill_coeff, kd)
                return rate * h
            else:
                return rate

        return curried

    @staticmethod
    def degradation_rate(rate: float, decaying_species: str) -> Callable[[Network], float]:
        def curried(n: Network):
            return rate * n.species[decaying_species]

        return curried

    @staticmethod
    def custom_reaction_rate(rate_function_ast: str) -> Callable[[Network], float]:
        def curried(n: Network):
            return helper.evaluate_ast_string(rate_function_ast, n.symbols, species=n.species)

        return curried
