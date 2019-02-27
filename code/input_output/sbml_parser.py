import libsbml
from libsbml._libsbml import parseL3Formula

import helper
from models.custom_formula import CustomFormula
from models.degradation_formula import DegradationFormula
from models.network import Network
from models.reaction import Reaction
# 1. Core objects of libsbml:
# http://sbml.org/Software/libSBML/5.17.0/docs//python-api/group__core.html
# 2. Classes:
# http://sbml.org/Software/libSBML/5.17.0/docs//python-api/annotated.html
from models.reg_type import RegType
from models.regulation import Regulation
from models.transcription_formula import TranscriptionFormula


class SbmlParser:
    """
    Return dictionary of species in a given model
    :param Any model: A libsbml network model
    :returns Dict[str, float] of species where
        key: species id, value: initial concentration
    """

    @staticmethod
    def _get_species(model):
        species = {s.getId(): s.getInitialAmount() for s in model.getListOfSpecies()}
        return species

    """
    Return all defined symbols in the given model
    :param Any model: a libsbml network model
    :returns Dict[str, float] where key: symbol name, value: symbol value
    """

    @staticmethod
    def _get_symbols(model):
        symbols = {}

        """
         Both list of parameters and rules constitute symbols
         which may be used in the reactions. So, include all
         in a symbol table.
        """

        for param in model.getListOfParameters():
            if param.getConstant():
                symbols[param.getId()] = param.getValue()

        """
        Rules are different to parameters as they can contain symbols
        themselves, which can be the parameters parsed in the above loop.
        """

        for r in model.getListOfRules():
            val = helper.evaluate_ast(r.getMath(), symbols)
            symbols[r.getId()] = val

        for c in model.getListOfCompartments():
            symbols[c.getId()] = c.getSize()

        return symbols

    """
    Return all reactions in the given model
    :param Any model: a libsbml network model
    :returns List[Reaction] of model's reaction
    """

    @staticmethod
    def _get_reactions(model):
        # Evaluate and store global parameters in a symbol table
        symbols = SbmlParser._get_symbols(model)
        reactions = []

        for x in model.getListOfReactions():
            rate_function = helper.ast_to_string(x.getKineticLaw().getMath().deepCopy())

            reactants = x.getListOfReactants()
            products = x.getListOfProducts()

            left = [y.getSpecies() for y in reactants]
            right = [y.getSpecies() for y in products]
            parameters = {p.getId(): p.getValue() for p in x.getKineticLaw().getListOfParameters()}

            sbo = x.getSBOTerm()

            # 179 -> Degradation
            # 183 -> Transcription
            # 184 -> Translation

            # if sbo == "179":
            #     r = CustomFormula(rate_function, parameters, symbols)
            # elif sbo == "183":
            #     r = CustomFormula(rate_function, parameters, symbols)
            # elif sbo == "184":
            #     r = CustomFormula(rate_function, parameters, symbols)
            # else:
            #     r = CustomFormula(rate_function, parameters, symbols)

            r = CustomFormula(rate_function, parameters, symbols)
            reactions.append(Reaction(x.getName(), left, right, r))
        return reactions

    """
    Return a Network representing the SBML model in the given SBML file.
    :param str filename: Filename for SBML file
    :returns Network representing given filename
    """

    @staticmethod
    def parse(filename):
        net = Network()
        sbml = libsbml.SBMLReader()

        parsed = sbml.readSBML(filename)
        model = parsed.getModel()

        # Initialise species and their initial amounts
        net.species = SbmlParser._get_species(model)

        # Parse reactions and create CustomReaction objects
        net.reactions = SbmlParser._get_reactions(model)

        return net

    @staticmethod
    def save_as_sbml(net):
        model = libsbml.Model(2, 3)  # SBML level 2, version 3

        # species = {s.getId(): s.getInitialAmount()
        #            for s in model.getListOfSpecies()}

        for name, initial_amount in net.species.items():
            spec = model.createSpecies()
            spec.setId(name)
            spec.setInitialAmount(initial_amount)

        # rate_function = helper.ast_to_string(x.getKineticLaw().getMath().deepCopy())

        # reactants = x.getListOfReactants()
        # products = x.getListOfProducts()

        # left = [y.getSpecies() for y in reactants]
        # right = [y.getSpecies() for y in products]
        # parameters = {p.getId(): p.getValue() for p in x.getKineticLaw().getListOfParameters()}

        # sbo = x.getSBOTerm()

        for react in net.reactions:
            reaction = model.createReaction()

            for x in react.left:
                z = reaction.createReactant()
                z.setSpecies(x)

            for x in react.right:
                z = reaction.createProduct()
                z.setSpecies(x)

            parseL3Formula('k * s1 * c1')

            # r = libsbml.Reaction(2, 3)
            # r.addReactant()
        return model


def test():
    species = {"x": 0, "y": 20}

    reg = Regulation("y", "x", RegType.REPRESSION, 40)
    x_trans = TranscriptionFormula(5, "x")
    x_trans.set_regulation(2, [reg])

    reactions = [Reaction("", [], ["x"], x_trans),
                 Reaction("", ["y"], [], DegradationFormula(0.3, "y"))]

    net: Network = Network()
    net.species = species
    net.reactions = reactions

    m = SbmlParser.save_as_sbml(net)
    print(m.getListOfReactions()[0].getListOfProducts()[0].getSpecies())


if __name__ == '__main__':
    test()
