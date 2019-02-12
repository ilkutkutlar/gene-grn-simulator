from typing import Dict, List, Tuple, Union

from models.reaction import Reaction


class Network:
    species: Dict[str, float]
    reactions: List  # of Reaction

    def __init__(self) -> None:
        self.species = dict()
        self.reactions = list()

    def apply_change_vector(self, change: Dict[str, float]) -> None:
        for x in change:
            self.species[x] += change[x]

    def mutate(self, mutations: Dict[str, Tuple[float, str]]):
        for m in mutations:
            new_value = mutations[m][0]
            reaction_name = mutations[m][1]

            if reaction_name != "":
                r = self.find_reaction_by_name(reaction_name)
                r.rate_function.mutate(mutations)
            else:
                self.species[m] = new_value

    def find_reaction_by_name(self, name: str) -> Union[Reaction, None]:
        t = list(filter(lambda r: r.name == name, self.reactions))
        if t:
            return t[0]
        else:
            return None

    def __str__(self) -> str:
        ret = "\nSpecies: \n"
        for x in self.species:
            ret += "    " + x + ": " + str(self.species[x]) + "\n"

        ret += "\nReactions: \n"
        for x in self.reactions:
            ret += "    " + str(x) + "\n"

        return ret