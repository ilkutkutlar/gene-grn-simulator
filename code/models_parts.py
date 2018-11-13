from typing import List


class Species:
    # This is used to uniquely identify the part in the network.
    identifier: str

    def __init__(self, identifier: str):
        self.identifier = identifier


class Signal(Species):

    def __init__(self, identifier: str):
        super(Signal, self).__init__(identifier)


class Protein:
    name: str
    translation_rate: float
    degradation: float


class TranscriptionFactor:
    name: str
    # E.g. LacI can be activated by the binding of allolactose
    binding_site: Signal


class mRNA:
    name: str
    degradation: float
    protein: Protein


class Promoter:
    name: str
    activator_binding_sites: List[TranscriptionFactor]
    repressor_binding_sites: List[TranscriptionFactor]
    # Also known as maximal transcription rate, represented by beta in some publications
    promoter_strength_active: float
    # The promoter leakiness
    promoter_strength_repressed: float

# class Reaction:
#     name: str
#     rate:
#     change_vector: