from math import log, e

from gillespie import GillespieSimulator
from models import TranscriptionReaction, Network, MrnaDegradationReaction, TranslationReaction, \
    ProteinDegradationReaction

# region Constants
# Transcription-related values

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
    net = Network()
    net.species = {"laci_mrna": 100, "tetr_mrna": 80, "cl_mrna": 50,
                   "laci_p": 10, "tetr_p": 10, "cl_p": 10}

    net.reactions = [
        TranscriptionReaction(alpha, 40, 2, ["cl_p"], "", "laci_mrna"),
        TranscriptionReaction(alpha, 40, 2, ["laci_p"], "", "tetr_mrna"),
        TranscriptionReaction(alpha, 40, 2, ["tetr_p"], "", "cl_mrna"),
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

    g = GillespieSimulator(net, 100, 0)
    g.visualise(g.simulate())


simulate_repressilator()
