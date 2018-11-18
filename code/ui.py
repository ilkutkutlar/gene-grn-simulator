from math import log, e

from gillespie import GillespieSimulator
from models import TranscriptionReaction, Network, MrnaDegradationReaction, TranslationReaction, \
    ProteinDegradationReaction, SimulationSettings, Regulation, RegType

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
    species = {"laci_mrna": 100, "tetr_mrna": 80, "cl_mrna": 50,
               "laci_p": 10, "tetr_p": 10, "cl_p": 10}
    regulations = [Regulation(from_gene="cl_p", to_gene="laci_mrna", reg_type=RegType.REPRESSION),
                   Regulation(from_gene="laci_p", to_gene="tetr_mrna", reg_type=RegType.REPRESSION),
                   Regulation(from_gene="tetr_p", to_gene="cl_mrna", reg_type=RegType.REPRESSION)]
    reactions = [
        TranscriptionReaction(alpha, 40, 2, "", "laci_mrna"),
        TranscriptionReaction(alpha, 40, 2, "", "tetr_mrna"),
        TranscriptionReaction(alpha, 40, 2, "", "cl_mrna"),
        MrnaDegradationReaction(mRNA_decay_rate, "laci_mrna", ""),
        MrnaDegradationReaction(mRNA_decay_rate, "tetr_mrna", ""),
        MrnaDegradationReaction(mRNA_decay_rate, "cl_mrna", ""),
        TranslationReaction(beta, "laci_mrna", "laci_p"),
        TranslationReaction(beta, "tetr_mrna", "tetr_p"),
        TranslationReaction(beta, "cl_mrna", "cl_p"),
        ProteinDegradationReaction(protein_decay_rate, "laci_p", ""),
        ProteinDegradationReaction(protein_decay_rate, "tetr_p", ""),
        ProteinDegradationReaction(protein_decay_rate, "cl_p", "")
    ]

    net = Network()
    net.initialise(species, reactions, regulations)

    end_time = 100
    s = SimulationSettings("Results", "Time", "Concentration", 0, end_time,
                           [("LacI Protein", "laci_p"),
                            ("TetR Protein", "tetr_p"),
                            ("Cl Protein", "cl_p")])
    g = GillespieSimulator(net, s)
    g.visualise(g.simulate())


def simulate_switch():
    # Values source: http://www.ebi.ac.uk/biomodels-main/BIOMD0000000507

    net = Network()
    net.species = {"mrna_one": 10, "mrna_two": 100,
                   "p_one": 0, "p_two": 100}

    net.reactions = [
        TranscriptionReaction(200, 1, 1, "", "one_mrna"),
        TranscriptionReaction(200, 1, 2.5, "", "two_mrna"),

        MrnaDegradationReaction(0.3, "mrna_one", ""),
        MrnaDegradationReaction(0.3, "mrna_two", ""),

        TranslationReaction(156.25, "mrna_one", "p_one"),
        TranslationReaction(15.6, "mrna_two", "p_two"),

        ProteinDegradationReaction(1, "p_one", ""),
        ProteinDegradationReaction(1, "p_two", ""),
    ]

    s = SimulationSettings("Results", "Time", "Concentration",
                           0, 10,
                           [("Protein One", "p_one"),
                            ("Protein Two", "p_two")])
    g = GillespieSimulator(net, s)
    g.visualise(g.simulate())


simulate_repressilator()
# simulate_switch()

