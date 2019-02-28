from reverse_engineering.mutable import ReactionMutable, VariableMutable


class Network:

    def __init__(self):
        self.species = dict()  # of Dict[str, float]
        self.reactions = list()  # of Reaction

    """
    Change species concentrations of network using a change vector
    :param Dict[str, float] change: key: species name, value: concentration change by
    """

    def apply_change_vector(self, change):
        for x in change:
            self.species[x] += change[x]

    """
    Mutate species and reactions of the network
    :param Dict[str, Tuple[float, str]] mutations: key: mutable name, value: (new value, reaction name)
    """

    def mutate(self, mutations):
        for m in mutations:
            if isinstance(m, ReactionMutable):
                r = self.find_reaction_by_name(m.reaction_name)
                r.rate_function.mutate(m)
            elif isinstance(m, VariableMutable):
                self.species[m] = m.current_value
            else:
                pass  # error

    """
    Return reaction with given name
    :param str name: Name of reaction
    :returns Reaction if found, None if not
    """

    def find_reaction_by_name(self, name):
        t = list(filter(lambda r: r.name == name, self.reactions))
        if t:
            return t[0]
        else:
            return None

    def __str__(self):
        ret = "\nSpecies: \n"
        for x in self.species:
            ret += "    " + x + ": " + str(self.species[x]) + "\n"

        ret += "\nReactions: \n"
        for x in self.reactions:
            ret += "    " + str(x) + "\n"

        return ret
