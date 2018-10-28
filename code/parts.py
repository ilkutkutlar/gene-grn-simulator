from typing import List

from models_parts import TranscriptionFactor, mRNA, Protein, Promoter


class Repository:
    promoters: List[Promoter]
    proteins: List[Protein]
    mRNAs: List[mRNA]
    transcription_factors: List[TranscriptionFactor]
