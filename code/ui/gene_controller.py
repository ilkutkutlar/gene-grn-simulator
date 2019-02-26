from models.input_gate import InputGate
from models.transcription_formula import TranscriptionFormula
from models.translation_formula import TranslationFormula
from models.degradation_formula import DegradationFormula
from models.network import Network
from models.reaction import Reaction
from models.reg_type import RegType
from models.regulation import Regulation


def get_test():
    species = {"px": 100, "py": 100, "pz": 30, "x": 100, "y": 25, "z": 20}

    x_trans = TranscriptionFormula(5, "x")

    y_trans = TranscriptionFormula(60, "y")
    y_trans.set_regulation(1, [Regulation("px", "y", RegType.REPRESSION, 40)])

    z_trans = TranscriptionFormula(20, "z")
    z_trans.set_regulation(2, [Regulation("py", "z", RegType.ACTIVATION, 40)])

    reactions = [Reaction("x_trans", [], ["x"], x_trans),
                 Reaction("y_trans", [], ["y"], y_trans),
                 Reaction("z_trans", [], ["z"], z_trans),

                 Reaction("x_deg", ["x"], [], DegradationFormula(0.01, "x")),
                 Reaction("y_deg", ["y"], [], DegradationFormula(0.01, "y")),
                 Reaction("z_deg", ["z"], [], DegradationFormula(0.1, "z")),

                 Reaction("px_deg", ["px"], [], DegradationFormula(0.01, "px")),
                 Reaction("py_deg", ["py"], [], DegradationFormula(0.01, "py")),
                 Reaction("pz_deg", ["pz"], [], DegradationFormula(0.01, "pz")),

                 Reaction("px_translation", [], ["px"], TranslationFormula(0.2, "x")),
                 Reaction("py_translation", [], ["py"], TranslationFormula(5, "y")),
                 Reaction("pz_translation", [], ["pz"], TranslationFormula(1, "z"))]

    net: Network = Network()
    net.species = species
    net.reactions = reactions

    return net


def get_test2():
    species = {"px": 100, "py": 100, "pz": 30, "x": 25, "y": 25, "z": 25}

    x_trans = TranscriptionFormula(5, "x")

    y_trans = TranscriptionFormula(5, "y")

    z_trans = TranscriptionFormula(5, "z")
    z_trans.set_regulation(2,
                           [Regulation("py", "z", RegType.ACTIVATION, 40),
                            Regulation("px", "z", RegType.ACTIVATION, 80)],
                           InputGate.AND)

    reactions = [Reaction("x_trans", [], ["x"], x_trans),
                 Reaction("y_trans", [], ["y"], y_trans),
                 Reaction("z_trans", [], ["z"], z_trans),

                 # Reaction("x_deg", ["x"], [], DegradationFormula(0.01, "x")),
                 # Reaction("y_deg", ["y"], [], DegradationFormula(0.01, "y")),
                 # Reaction("z_deg", ["z"], [], DegradationFormula(0.01, "z")),

                 # Reaction("px_deg", ["px"], [], DegradationFormula(0.01, "px")),
                 # Reaction("py_deg", ["py"], [], DegradationFormula(0.01, "py")),
                 # Reaction("pz_deg", ["pz"], [], DegradationFormula(0.01, "pz")),

                 Reaction("px_translation", [], ["px"], TranslationFormula(0.2, "x")),
                 Reaction("py_translation", [], ["py"], TranslationFormula(5, "y")),
                 Reaction("pz_translation", [], ["pz"], TranslationFormula(1, "z"))]

    net: Network = Network()
    net.species = species
    net.reactions = reactions

    return net


# def get_test_constraints():
#     c1 = Constraint("y", lambda y: 200 - y, (40, 60))
#     c2 = Constraint("z", lambda x: x - 150, (0, 20))
#     m = Mutable(0.5, 50, 0.5, "one")

class GeneController:
    instance = None

    def __init__(self):
        # c1 = Constraint("y", lambda y: 200 - y, (40, 60))
        # c2 = Constraint("z", lambda x: x - 150, (0, 20))
        # m = Mutable(0.5, 50, 0.5, "x_trans")
        self.network = get_test()
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
