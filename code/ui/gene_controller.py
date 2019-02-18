from models.network import Network


class GeneController:
    instance = None

    def __init__(self):
        self.network = Network()

    @staticmethod
    def get_instance():
        if not GeneController.instance:
            GeneController.instance = GeneController()
        return GeneController.instance

    def get_species(self):
        return self.network.species

    def get_reactions(self):
        return self.network.reactions

    def add_species(self, key, value):
        self.network.species[key] = value

    def add_reaction(self, reaction):
        self.network.reactions.append(reaction)
