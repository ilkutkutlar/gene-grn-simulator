from math import log, e

from gillespie import GillespieSimulator
from models import TranscriptionReaction, Network, DegradationReaction, TranslationReaction, \
    SimulationSettings, Regulation, RegType, CustomReaction

# region Constants


# Transcription-related values
from sbml_parser import SbmlParser

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


def simulate_repressilator():
    # p_lacI0 = 10
    # p_tetR0 = 10
    # p_cl0 = 10

    # m_lacI0 = 100
    # m_tetR0 = 80
    # m_cl0 = 50

    species = {"laci_mrna": 0, "tetr_mrna": 20, "cl_mrna": 0,
               "laci_p": 0, "tetr_p": 0, "cl_p": 0}
    regulations = [Regulation(from_gene="cl_p", to_gene="laci_mrna", reg_type=RegType.REPRESSION),
                   Regulation(from_gene="laci_p", to_gene="tetr_mrna", reg_type=RegType.REPRESSION),
                   Regulation(from_gene="tetr_p", to_gene="cl_mrna", reg_type=RegType.REPRESSION)]
    reactions = [
        TranscriptionReaction(alpha, 40, 2, [""], ["laci_mrna"]),
        TranscriptionReaction(alpha, 40, 2, [""], ["tetr_mrna"]),
        TranscriptionReaction(alpha, 40, 2, [""], ["cl_mrna"]),

        DegradationReaction(mRNA_decay_rate, ["laci_mrna"], [""]),
        DegradationReaction(mRNA_decay_rate, ["tetr_mrna"], [""]),
        DegradationReaction(mRNA_decay_rate, ["cl_mrna"], [""]),

        TranslationReaction(beta, ["laci_mrna"], ["laci_p"]),
        TranslationReaction(beta, ["tetr_mrna"], ["tetr_p"]),
        TranslationReaction(beta, ["cl_mrna"], ["cl_p"]),

        DegradationReaction(protein_decay_rate, ["laci_p"], [""]),
        DegradationReaction(protein_decay_rate, ["tetr_p"], [""]),
        DegradationReaction(protein_decay_rate, ["cl_p"], [""])
    ]

    net = Network()
    net.initialise(species, reactions, regulations)

    print(net)

    end_time = 100
    s = SimulationSettings("Results", "Time", "Concentration", 0, end_time,
                           [("LacI Protein", "laci_p"),
                            ("TetR Protein", "tetr_p"),
                            ("Cl Protein", "cl_p")])
    GillespieSimulator.visualise(GillespieSimulator.simulate(net, s), s)


def simulate_switch():
    # Values source: http://www.ebi.ac.uk/biomodels-main/BIOMD0000000507

    net = Network()
    net.species = {"mrna_one": 10, "mrna_two": 100,
                   "p_one": 0, "p_two": 100}

    net.reactions = [
        TranscriptionReaction(200, 1, 1, [""], ["one_mrna"]),
        TranscriptionReaction(200, 1, 2.5, [""], ["two_mrna"]),

        DegradationReaction(0.3, ["mrna_one"], [""]),
        DegradationReaction(0.3, ["mrna_two"], [""]),

        TranslationReaction(156.25, ["mrna_one"], ["p_one"]),
        TranslationReaction(15.6, ["mrna_two"], ["p_two"]),

        DegradationReaction(1, ["p_one"], [""]),
        DegradationReaction(1, ["p_two"], [""]),
    ]

    s = SimulationSettings("Results", "Time", "Concentration",
                           0, 10,
                           [("Protein One", "p_one"),
                            ("Protein Two", "p_two")])
    GillespieSimulator.visualise(GillespieSimulator.simulate(net, s), s)


def simulate_parser():
    filename = "other_files/BIOMD0000000012.xml"
    net: Network = SbmlParser.parse(filename)
    print(net)

    end_time = 100
    s = SimulationSettings("Results", "Time", "Concentration", 0, end_time,
                           [("LacI Protein", "PX"),
                            ("TetR Protein", "PY"),
                            ("Cl Protein", "PZ")])

    GillespieSimulator.visualise(GillespieSimulator.simulate(net, s), s)


def simulate_case_study():
    # Reaction rates:

    k1 = 10 ** 5
    k2 = 10 ** 5
    k3 = 3.6119
    k4 = 10 ** 5
    k5 = 10 ** 5
    k6 = 0.9079
    k7 = 2.6978
    k8 = 0.8902
    k9 = 5.8903
    k10 = 0.1101
    k11 = 0.6296
    k12 = 0.5307
    k13 = 0.0880
    k14 = 0.0065

    k_1 = 10 ** -2
    k_2 = 10 ** -2
    k_3 = 0.1092
    k_4 = 10 ** -2
    k_5 = 10 ** -2
    k_6 = 5.7766
    k_14 = 0.2208

    y1 = 0.9470
    y2 = 2.7057
    y3 = 0.2248
    y4 = 0.1646

    net = Network()
    net.species = {
        "fpm": 0,
        "MmyR": 0,
        "fpm:MmyR": 0,
        "MmfR": 0,
        "fpm:MmfR": 0,
        "MMF": 0,
        "fpm:MmfR:MMF": 0,
        "apm": 0,
        "apm:MmyR": 0,
        "apm:MmfR": 0,
        "MMY": 0,
        "apm:MmfR:MMF": 0,
        "MmfR:MMF": 0
    }

    r1 = CustomReaction("{}".format(k1), ["fpm", "MmyR"], ["fpm:MmyR"])
    r1_ = CustomReaction("{}".format(k_1), ["fpm:MmyR"], ["fpm", "MmyR"])

    r2 = CustomReaction("{}".format(k2), ["fpm", "MmfR"], ["fpm:MmfR"])
    r2_ = CustomReaction("{}".format(k_2), ["fpm:MmfR"], ["fpm", "MmfR"])

    r3 = CustomReaction("{}".format(k7), ["fpm"], ["MmyR", "fpm"])
    r4 = CustomReaction("{}".format(k8), ["fpm"], ["MmfR", "fpm"])
    r5 = CustomReaction("{}".format(k9), ["fpm"], ["MMF", "fpm"])

    r6 = CustomReaction("{}".format(k3), ["fpm:MmfR", "MMF"], ["fpm:MmfR:MMF"])
    r6_ = CustomReaction("{}".format(k_3), ["fpm:MmfR:MMF"], ["fpm:MmfR", "MMF"])

    r7 = CustomReaction("{}".format(k11), ["fpm:MmfR:MMF"], ["fpm", "MmfR:MMF"])

    r8 = CustomReaction("{}".format(k4), ["apm", "MmyR"], ["apm:MmyR"])
    r8_ = CustomReaction("{}".format(k_4), ["apm:MmyR"], ["apm", "MmyR"])

    r9 = CustomReaction("{}".format(k5), ["apm", "MmfR"], ["apm:MmfR"])
    r9_ = CustomReaction("{}".format(k_5), ["apm:MmfR"], ["apm", "MmfR"])

    r10 = CustomReaction("{}".format(k10), ["apm"], ["MMY", "apm"])

    r11 = CustomReaction("{}".format(k6), ["apm:MmfR", "MMF"], ["apm:MmfR:MMF"])
    r11_ = CustomReaction("{}".format(k_6), ["apm:MmfR:MMF"], ["apm:MmfR", "MMF"])

    r12 = CustomReaction("{}".format(k12), ["apm:MmfR:MMF"], ["apm", "MmfR:MMF"])
    r13 = CustomReaction("{}".format(k13), ["MmfR:MMF"], ["MmfR", "MMF"])

    r14 = CustomReaction("{}".format(k14), ["MmfR", "MMF"], ["MmfR:MMF"])
    r14_ = CustomReaction("{}".format(k_14), ["MmfR:MMF"], ["MmfR", "MMF"])

    r15 = CustomReaction("{}".format(y1), ["MmyR"], [""])
    r16 = CustomReaction("{}".format(y2), ["MmfR"], [""])
    r17 = CustomReaction("{}".format(y3), ["MMF"], [""])
    r18 = CustomReaction("{}".format(y4), ["MMY"], [""])

    net.reactions = [r1, r1_, r2, r2_, r3, r4, r5, r6, r6_, r7, r8, r8_, r9, r9_,
                     r10, r11, r11_, r12, r13, r14, r14_, r15, r16, r17, r18]

    end_time = 2
    s = SimulationSettings("Results", "Time", "Concentration", 0, end_time,
                           [("MmyR", "MmyR")])

    GillespieSimulator.visualise(GillespieSimulator.simulate(net, s), s)


def simulate_switch_parser():
    filename = "other_files/BIOMD0000000507.xml"
    net: Network = SbmlParser.parse(filename)
    net.species["compartment_1"] = 1
    print(net)

    end_time = 10
    s = SimulationSettings("Results", "Time", "Concentration", 0, end_time,
                           [("U", "u"),
                            ("V", "v")])

    GillespieSimulator.visualise(GillespieSimulator.simulate(net, s), s)


# simulate_repressilator()
# simulate_switch()
# simulate_parser()
simulate_case_study()
# simulate_switch_parser()
