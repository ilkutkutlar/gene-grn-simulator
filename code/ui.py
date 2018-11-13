from math import log, e

import models
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
protein.degradation = protein_decay_rate        # Reaction
protein.translation_rate = beta                 # Rule

mrna: models.mRNA = models.mRNA()
mrna.degradation = mRNA_decay_rate              # Reaction
mrna.protein = protein

promoter: models.Promoter = parts.Promoter
promoter.promoter_strength_active = alpha       # Rule
promoter.promoter_strength_repressed = alpha0   # Rule

lacI = models.Cassette("laci", promoter, [mrna])

cl = models.Cassette("cl", promoter, [mrna])

tetr = models.Cassette("tetr", promoter, [mrna])

network = models.Network()

network.genome = [lacI, tetr, cl]
network.regulations = [
    ("laci", "cl", models.RegType.REPRESSION),
    ("cl", "tetr", models.RegType.REPRESSION),
    ("tetr", "laci", models.RegType.REPRESSION)]


# Species
s: Simulator = Simulator(network, 0, 1000, 1000,
                         {"laci": 100, "cl": 50, "tetr": 80},
                         {"laci": 10, "cl": 10, "tetr": 10},
                         {})
s.visualise(s.simulate())
