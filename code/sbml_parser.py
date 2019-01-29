from typing import Dict, List

import libsbml

import helper

# 1. Core objects of libsbml:
# http://sbml.org/Software/libSBML/5.17.0/docs//python-api/group__core.html
# 2. Classes:
# http://sbml.org/Software/libSBML/5.17.0/docs//python-api/annotated.html
from models.formulae import CustomFormula
from models.network import Network
from models.reaction import Reaction


class SbmlParser:
    @staticmethod
    def _get_species_(model):
        species: Dict[str, float] = {}
        for s in model.getListOfSpecies():
            species[s.getId()] = s.getInitialAmount()
        return species

    @staticmethod
    def _get_symbols_(model):
        symbols: Dict[str, float] = {}

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
            val = helper.evaluate_ast_node(r.getMath(), symbols)
            symbols[r.getId()] = val

        for c in model.getListOfCompartments():
            symbols[c.getId()] = c.getSize()

        return symbols

    @staticmethod
    def _get_reactions_(model):
        reactions: List[Reaction] = []
        for x in model.getListOfReactions():
            reaction_rate_function = x.getKineticLaw().getMath().deepCopy()
            parameters = {}

            reactants = x.getListOfReactants()
            products = x.getListOfProducts()

            left: List[str] = list()
            right: List[str] = list()

            if not reactants:
                left.append("")

            if not products:
                right.append("")

            for y in reactants:
                left.append(y.getSpecies())

            for y in products:
                right.append(y.getSpecies())

            for p in x.getKineticLaw().getListOfParameters():
                parameters[p.getId()] = p.getValue()

            sbo: str = x.getSBOTerm()

            # 179 -> Degradation
            # 183 -> Transcription
            # 184 -> Translation

            if sbo == "179":
                r = CustomFormula(helper.convert_ast_to_string(reaction_rate_function), parameters)
            elif sbo == "183":
                r = CustomFormula(helper.convert_ast_to_string(reaction_rate_function), parameters)
            elif sbo == "184":
                r = CustomFormula(helper.convert_ast_to_string(reaction_rate_function), parameters)
            else:
                r = CustomFormula(helper.convert_ast_to_string(reaction_rate_function), parameters)

            reactions.append(Reaction(left, right, r))
        return reactions

    @staticmethod
    def parse(filename: str) -> Network:
        net: Network = Network()
        sbml: libsbml.SBMLReader = libsbml.SBMLReader()

        parsed = sbml.readSBML(filename)
        model: libsbml.Model = parsed.getModel()

        # Initialise species and their initial amounts
        net.species = SbmlParser._get_species_(model)
        # Evaluate and store global parameters in a symbol table
        net.symbols = SbmlParser._get_symbols_(model)
        # Parse reactions and create CustomReaction objects
        net.reactions = SbmlParser._get_reactions_(model)

        return net

    @staticmethod
    def save_as_sbml(net: Network):
        model: libsbml.Model

        for react in net.reactions:
            r = libsbml.Reaction()
            r.addReactant()
            model.addReaction()
