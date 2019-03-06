import test


class GeneController:
    instance = None

    def __init__(self):
        # c1 = Constraint("y", lambda y: 200 - y, (40, 60))
        # c2 = Constraint("z", lambda x: x - 150, (0, 20))
        # m = Mutable(0.5, 50, 0.5, "x_trans")
        self.network = test.get_large_network()
        self.mutables = []
        self.constraints = []

    @staticmethod
    def get_instance():
        if not GeneController.instance:
            GeneController.instance = GeneController()
        return GeneController.instance

    def get_species(self):
        return self.network.species

    def get_reactions(self):
        return self.network.reactions

    def get_mutables(self):
        return self.mutables

    def get_mutable_by_name(self, name):
        f = list(filter(
            lambda m: m.variable_name == name,
            self.mutables))
        return f[0] if f else None

    def get_constraints(self):
        return self.constraints

    def add_species(self, key, value):
        self.network.species[key] = value

    def add_reaction(self, reaction):
        self.network.reactions.append(reaction)

    def add_mutable(self, species, mutable):
        self.mutables[species] = mutable

    def add_constraint(self, constraint):
        self.constraints.append(constraint)

    def remove_species(self, name):
        del self.network.species[name]

    def remove_reaction_by_name(self, name):
        r = list(filter(lambda x: x.name == name, self.network.reactions))[0]
        self.network.reactions.remove(r)

    def remove_reaction_by_index(self, index):
        del self.network.reactions[index]

    def remove_mutable(self, index):
        del self.mutables[index]

    def remove_constraint(self, index):
        del self.constraints[index]
