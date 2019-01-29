from typing import Dict, List


class Network:
    species: Dict[str, float]
    reactions: List  # of Reaction

    def __init__(self) -> None:
        self.species = dict()
        self.reactions = list()

    def apply_change_vector(self, change: Dict[str, float]) -> None:
        for x in change:
            self.species[x] += change[x]

    def __str__(self) -> str:
        ret = "\nSpecies: \n"
        for x in self.species:
            ret += "    " + x + ": " + str(self.species[x]) + "\n"

        ret += "\nReactions: \n"
        for x in self.reactions:
            ret += "    " + str(x) + "\n"

        return ret
