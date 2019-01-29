from typing import Dict, List

from models.regulation import Regulation


class Network:
    species: Dict[str, float]
    reactions: List  # of Reaction
    regulations: List[Regulation]

    def __init__(self) -> None:
        self.species = {}
        self.reactions = []
        self.regulations = []

    def initialise(self, species: Dict[str, float],
                   reactions: List):
        self.species = species
        self.reactions = reactions

    def apply_change_vector(self, change: Dict[str, float]) -> None:
        for x in change:
            self.species[x] += change[x]

    def __str__(self) -> str:
        # ret = "\nGlobal Parameters: \n"
        # for x in self.symbols:
        #     ret += "    " + x + ": " + str(self.symbols[x]) + "\n"

        ret = "\nSpecies: \n"
        for x in self.species:
            ret += "    " + x + ": " + str(self.species[x]) + "\n"

        ret += "\nRegulations: \n"
        for x in self.regulations:
            ret += "    " + str(x) + "\n"

        ret += "\nReactions: \n"
        for x in self.reactions:
            ret += "    " + str(x) + "\n"

        return ret
