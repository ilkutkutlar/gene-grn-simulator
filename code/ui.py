from math import log, e

import models
import models_parts

# region Constants
# Transcription-related values
from simulation import Simulation

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

protein: models_parts.Protein = models_parts.Protein()
protein.degradation = protein_decay_rate
protein.translation_rate = beta

mrna: models.mRNA = models.mRNA()
mrna.degradation = mRNA_decay_rate
mrna.translates_into = protein

promoter: models.Promoter = models_parts.Promoter
promoter.promoter_strength_active = alpha
promoter.promoter_strength_repressed = alpha0

lacI = models.Cassette("laci",
                       promoter, [mrna])

cl = models.Cassette("cl",
                     promoter, [mrna])

tetr = models.Cassette("tetr",
                       promoter, [mrna])

network = models.Network()

network.genes = [lacI, tetr, cl]
network.regulations = [
    ("laci", "cl", models.Regulation.REPRESSION),
    ("cl", "tetr", models.Regulation.REPRESSION),
    ("tetr", "laci", models.Regulation.REPRESSION)]

network.mrna_init = {"laci": 100,
                     "cl": 50,
                     "tetr": 80}
network.protein_init = {"laci": 10,
                        "cl": 10,
                        "tetr": 10}

s: Simulation = Simulation(network)
s.simulate()
