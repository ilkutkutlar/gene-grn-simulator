from enum import Enum


class ReactionType(Enum):
    TRANSCRIPTION = 1
    TRANSLATION = 2
    DEGRADATION = 3
    CUSTOM = 4
