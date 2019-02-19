from models.network import Network


class GeneController:
    instance = None

    def __init__(self):
        self.network = Network()
        self.mutables = {}
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
