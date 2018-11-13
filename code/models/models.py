from enum import Enum
from typing import List, Tuple, NamedTuple

from models.parts import Promoter, mRNA, Species, Signal


class Cassette(Species):
    promoter: Promoter
    codes_for: List[mRNA]

    def __init__(self, identifier: str, promoter: Promoter, codes_for: List[mRNA]):
        super(Cassette, self).__init__(identifier)
        self.promoter = promoter
        self.codes_for = codes_for


class RegType(Enum):
    ACTIVATION = 1
    REPRESSION = -1


class Regulation(NamedTuple):
    from_gene: str
    to_gene: str
    reg_type: RegType


class Network:
    genome: List[Cassette]
    signals: List[Signal]
    regulations: List[Regulation]

    def get_gene_by_id(self, ident: str) -> Cassette:
        for x in self.genome:
            if x.identifier == ident:
                return x

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

    def ps_active(self, tf_concentration, kd, promoter: Promoter, n=1) -> float:
        return promoter.promoter_strength_active * \
               (pow(tf_concentration, n) / (kd + pow(tf_concentration, n)))

    def ps_repressed(self, tf_concentration, kd, promoter: Promoter, n=1) -> float:
        return promoter.promoter_strength_active * \
               (1 / (1 + (pow(tf_concentration, n) / kd)))

    def ps_constitutitive(self):
        # TODO
        pass

    def get_regulators(self, ident: str) -> List[Tuple[str, RegType]]:
        regs: List[Tuple[str, RegType]] = list()

        for z in self.regulations:
            if z[1] == ident:
                regs.append((z[0], z[2]))

        return regs
