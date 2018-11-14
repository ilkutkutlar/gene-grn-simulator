from math import log, e
from typing import Callable, List, Dict

import models
from gillespie import Vector, GillespieSimulator
from models import parts
# region Constants
# Transcription-related values
from simulation import Simulator

ps_active = 0.5  # Promoter strength (active)
ps_repr = 5 * (10 ** -4)  # Promoter strength (repressed)
transcription_rate_active = ps_active * 60
transcription_rate_repr = ps_repr * 60

# Decay-related values
mRNA_half_life = 2
protein_half_life = 10
mRNA_decay_rate = log(2, e) / mRNA_half_life
protein_decay_rate = log(2, e) / protein_half_life

# Translation-related values
translation_efficiency = 19.97
translation_rate = translation_efficiency * mRNA_decay_rate

# Other values
hill_coeff = 2  # Hill coefficient
Km = 40  # Activation coefficient

# Derived values

# Active transcription rate (rescaled?)
alpha = transcription_rate_active * translation_efficiency \
        * protein_half_life / (log(2, e) * Km)
# Repressed transcription rate (rescaled?)
alpha0 = transcription_rate_repr * translation_efficiency \
         * protein_half_life / (log(2, e) * Km)
# Translation rate
beta = protein_decay_rate / mRNA_decay_rate
# endregion

protein: parts.Protein = parts.Protein()
protein.degradation = protein_decay_rate  # Reaction
protein.translation_rate = beta  # Rule

mrna: models.mRNA = models.mRNA()
mrna.degradation = mRNA_decay_rate  # Reaction
mrna.protein = protein

promoter: models.Promoter = parts.Promoter
promoter.promoter_strength_active = alpha  # Rule
promoter.promoter_strength_repressed = alpha0  # Rule

lacI = models.Cassette("laci", promoter, [mrna])

cl = models.Cassette("cl", promoter, [mrna])

tetr = models.Cassette("tetr", promoter, [mrna])

network = models.Network()

network.genome = [lacI, tetr, cl]
network.regulations = [
    ("laci", "cl", models.RegType.REPRESSION),
    ("cl", "tetr", models.RegType.REPRESSION),
    ("tetr", "laci", models.RegType.REPRESSION)]


# # Species
# s: Simulator = Simulator(network, 0, 1000, 1000,
#                          {"laci": 100, "cl": 50, "tetr": 80},
#                          {"laci": 10, "cl": 10, "tetr": 10},
#                          {})
# s.visualise(s.simulate())

# n = [m_lacI, m_tetR, m_cl, p_lacI, p_tetR, p_cl]

def tr_rate(reg_p_index: int) -> Callable[[Vector], float]:
    def func(n: Vector):
        return alpha * (1 / (1 + (pow(n[reg_p_index], 2) / 40)))

    return func


def mrna_deg_rate(m_index: int) -> Callable[[Vector], float]:
    def func(n: Vector):
        return mRNA_decay_rate * n[m_index]

    return func


def prot_translation_rate(m_index: int):
    def func(n: Vector):
        return beta * n[m_index]

    return func


def prot_deg_rate(p_index: int):
    def func(n: Vector):
        return protein_decay_rate * n[p_index]

    return func


def regulation():
    # n = [100, 80, 50, 10, 10, 10]
    # n = [m_lacI, m_tetR, m_cl, p_lacI, p_tetR, p_cl]

    # m_lacI0 = 100
    # m_tetR0 = 80
    # m_cl0 = 50

    # p_lacI0 = 10
    # p_tetR0 = 10
    # p_cl0 = 10

    # Reactions (SEPARATE from the mrna and proteins, they are just reactions!):
    # For each mRNA:
    #   - transcription
    #   - degradation
    # For each protein:
    #   - translation
    #   - degradation
    # Thus the exact change vectors are:
    # 1. Consider this: d[mRNA]/dt = a_m*H([TF]) - b_m * [mRNA]
    #   1.1. regulated transcription rate = a_m*H([TF])
    #   1.2. mRNA degradation rate = - b_m * [mRNA]
    # m_lacI_tr = []

    # n = [100, 80, 50, 10, 10, 10]
    # n = [m_lacI, m_tetR, m_cl, p_lacI, p_tetR, p_cl]
    n0 = [100, 80, 50, 10, 10, 10]
    t0 = 0
    sim_time = 100

    def tr_change(n: Vector, changes: Dict[int, Callable[[Vector], float]]) -> Vector:
        ret: Vector = n.copy()
        for i in changes:
            # For -ve changes, changes is -
            ret[i] += changes[i](n)
        return ret

    #  -> mRNA
    def laci_tr_change(n: Vector):
        return tr_change(n, {0: tr_rate(5)})

    # mRNA ->
    def laci_mrna_deg_change(n: Vector):
        n[0] -= mrna_deg_rate(0)(n)
        return n

    # mRNA -> Protein
    def laci_prot_translation_change(n: Vector):
        n[3] += prot_translation_rate(0)(n)
        return n

    # Protein ->
    def laci_prot_degradation_change(n: Vector):
        n[3] -= prot_deg_rate(3)(n)
        return n

    #  -> mRNA
    def tetr_tr_change(n: Vector):
        return tr_change(n, {1: tr_rate(3)})

    # mRNA ->
    def tetr_mrna_deg_change(n: Vector):
        n[1] -= mrna_deg_rate(1)(n)
        return n

    # mRNA -> Protein
    def tetr_prot_translation_change(n: Vector):
        n[4] += prot_translation_rate(1)(n)
        return n

    # Protein ->
    def tetr_prot_degradation_change(n: Vector):
        n[4] -= prot_deg_rate(4)(n)
        return n

    #  -> mRNA
    def cl_tr_change(n: Vector):
        return tr_change(n, {2: tr_rate(4)})

    # mRNA ->
    def cl_mrna_deg_change(n: Vector):
        n[2] -= mrna_deg_rate(2)(n)
        return n

    # mRNA -> Protein
    def cl_prot_translation_change(n: Vector):
        n[5] += prot_translation_rate(2)(n)
        return n

    # Protein ->
    def cl_prot_degradation_change(n: Vector):
        n[5] -= prot_deg_rate(5)(n)
        return n

    r = [tr_rate(5), mrna_deg_rate(0), prot_translation_rate(0), prot_deg_rate(3),
         tr_rate(3), mrna_deg_rate(1), prot_translation_rate(1), prot_deg_rate(4),
         tr_rate(4), mrna_deg_rate(2), prot_translation_rate(2), prot_deg_rate(5)]

    v = [laci_tr_change, laci_mrna_deg_change, laci_prot_translation_change, laci_prot_degradation_change,
         tetr_tr_change, tetr_mrna_deg_change, tetr_prot_translation_change, tetr_prot_degradation_change,
         cl_tr_change, cl_mrna_deg_change, cl_prot_translation_change, cl_prot_degradation_change]

    g = GillespieSimulator(r, v, sim_time, n0, t0)
    results = g.simulate()
    g.visualise(results)


regulation()
