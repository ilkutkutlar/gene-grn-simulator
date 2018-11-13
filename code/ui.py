from math import log, e
from typing import Callable

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
n = 2  # Hill coefficient
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
    def func(n: Vector): return alpha * (pow(n[reg_p_index], 2) / (40 + pow(n[reg_p_index], 2)))

    return func


def mrna_deg_rate(reg_p_index: int, m_index: int) -> Callable[[Vector], float]:
    def func(n: Vector): return n[reg_p_index] * n[m_index]

    return func


def prot_translation_rate():
    def func(n: Vector): return beta

    return func


def prot_deg_rate():
    def func(n: Vector): return protein_decay_rate

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
    sim_time = 10

    # LacI mRNA reactions: Regulated repressively by cl

    #  -> mRNA
    def laci_tr_rate(n: Vector): return alpha * (pow(n[5], 2) / (40 + pow(n[5], 2)))

    def laci_tr_change(n: Vector): return [n[0] + (alpha * (pow(n[5], 2) / (40 + pow(n[5], 2)))), n[1], n[2], n[3],
                                           n[4], n[5]]

    # mRNA ->
    def laci_mrna_deg_rate(n: Vector): return n[5] * n[0]

    def laci_mrna_deg_change(n: Vector): return [n[0] - (n[5] * n[0]), n[1], n[2], n[3], n[4], n[5]]

    # mRNA -> Protein
    def laci_prot_translation_rate(n: Vector): return n[0] * beta

    def laci_prot_translation_change(n: Vector): return [n[0], n[1], n[2], n[3] + beta, n[4], n[5]]

    # Protein ->
    def laci_prot_degradation_rate(n: Vector): return n[3] * protein_decay_rate

    def laci_prot_degradation_change(n: Vector): return [n[0], n[1], n[2], n[3] - protein_decay_rate, n[4], n[5]]

    # TetR mRNA reactions: Regulated repressively by LacI

    #  -> mRNA
    def tetr_tr_rate(n: Vector): return alpha * (pow(n[3], 2) / (40 + pow(n[3], 2)))

    def tetr_tr_change(n: Vector): return [n[0], n[1] + (alpha * (pow(n[3], 2) / (40 + pow(n[3], 2)))), n[2], n[3],
                                           n[4], n[5]]

    # mRNA ->
    def tetr_mrna_deg_rate(n: Vector): return n[3] * n[1]

    def tetr_mrna_deg_change(n: Vector): return [n[0], n[1] - (n[3] * n[1]), n[2], n[3], n[4], n[5]]

    # mRNA -> Protein
    def tetr_prot_translation_rate(n: Vector): return n[1] * beta

    def tetr_prot_translation_change(n: Vector): return [n[0], n[1], n[2], n[3], n[4] + beta, n[5]]

    # Protein ->
    def tetr_prot_degradation_rate(n: Vector): return n[4] * protein_decay_rate

    def tetr_prot_degradation_change(n: Vector): return [n[0], n[1], n[2], n[3], n[4] - protein_decay_rate, n[5]]

    # cl reactions: Regulated repressively by tetR
    # n = [m_lacI, m_tetR, m_cl, p_lacI, p_tetR, p_cl]

    #  -> mRNA
    def cl_tr_rate(n: Vector): return alpha * (pow(n[4], 2) / (40 + pow(n[4], 2)))

    def cl_tr_change(n: Vector): return [n[0], n[1], n[2] + (alpha * (pow(n[4], 2) / (40 + pow(n[4], 2)))), n[3], n[4],
                                         n[5]]

    # mRNA ->
    def cl_mrna_deg_rate(n: Vector): return n[4] * n[2]

    def cl_mrna_deg_change(n: Vector): return [n[0], n[1], n[2] - (n[4] * n[2]), n[3], n[4], n[5]]

    # mRNA -> Protein
    def cl_prot_translation_rate(n: Vector): return beta

    def cl_prot_translation_change(n: Vector): return [n[0], n[1], n[2], n[3], n[4], n[5] + beta]

    # Protein ->
    def cl_prot_degradation_rate(n: Vector): return protein_decay_rate

    def cl_prot_degradation_change(n: Vector): return [n[0], n[1], n[2], n[3], n[4], n[5] - protein_decay_rate]

    r = [laci_tr_rate, laci_mrna_deg_rate, laci_prot_translation_rate, laci_prot_degradation_rate,
         tetr_tr_rate, tetr_mrna_deg_rate, tetr_prot_translation_rate, tetr_prot_degradation_rate,
         cl_tr_rate, cl_mrna_deg_rate, cl_prot_translation_rate, cl_prot_degradation_rate]

    v = [laci_tr_change, laci_mrna_deg_change, laci_prot_translation_change, laci_prot_degradation_change,
         tetr_tr_change, tetr_mrna_deg_change, tetr_prot_translation_change, tetr_prot_degradation_change,
         cl_tr_change, cl_mrna_deg_change, cl_prot_translation_change, cl_prot_degradation_change]

    g = GillespieSimulator(r, v, sim_time, n0, t0)
    results = g.simulate()
    g.visualise(results)


regulation()
