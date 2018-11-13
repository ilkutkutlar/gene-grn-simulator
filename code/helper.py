from typing import Dict, Tuple, List

from models import Cassette, RegType, Network
from simulation import Concentration


# d[mRNA]/dt = a_m*H([TF]) - b_m * [mRNA]
# a_m           -> transcription rate
# H([TF])       -> Hill equation
# a_m*H([TF])   -> regulated promoter strength (rps)
# [mRNA]        -> mRNA concentration (mrna)
# b_m           -> mRNA degradation rate (b_m)
def delta_mrna(b_m: float, mrna: float, rps: float):
    return rps - b_m * mrna


# rps -> regulated promoter strength
def tr_rate(rps: float):
    return rps


# rps -> regulated promoter strength
def tr_change_vector(rps: float, n: List[float]):
    return rps


# mrna  -> mRNA concentration
# b_m   -> mRNA degradation rate
# def mrna_deg_rate(b_m: float, mrna: float):
#     return b_m * mrnab_m: float, mrna: float

def mrna_deg_change_vector(b_m: float, mrna: float):
    return -b_m * mrna

# d[protein]/dt = a_p*[mRNA] - b_p * [protein]
# [protein] -> protein concentration (p)
# [mRNA]    -> mRNA concentration (m)
# a_p       -> translation rate (a_p)
# b_p       -> protein degradation rate (b_p)
def delta_protein(a_p: float, b_p: float, p: float, m: float):
    return a_p * m - b_p * p


# rps = regulated promoter strength.
# Given a gene, calculates the promoter strength under regulation.
def calculate_rps(network: Network,
                  concent: Dict[str, Concentration],
                  gene: Cassette, kd: float, n: int):
    # Identifier of the gene and the regulation
    regulator: Tuple[str, RegType] = network.get_regulators(gene.identifier)[0]
    regulator_concentration = concent[regulator[0]].protein

    if regulator[1] == RegType.ACTIVATION:
        rps: float = network.ps_active(regulator_concentration, kd, gene.promoter, n)
    else:
        rps: float = network.ps_repressed(regulator_concentration, kd, gene.promoter, n)

    return rps
