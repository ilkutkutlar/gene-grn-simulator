import math
from typing import Dict, List

import libsbml

import helper
from models import Network, Reaction, CustomReaction


# 1. Core objects of libsbml:
# http://sbml.org/Software/libSBML/5.17.0/docs//python-api/group__core.html
# 2. Classes:
# http://sbml.org/Software/libSBML/5.17.0/docs//python-api/annotated.html


class SbmlParser:
    def __init__(self, filename):
        self.filename = filename

    def parse(self) -> Network:
        net: Network = Network()
        sbml: libsbml.SBMLReader = libsbml.SBMLReader()
        symbols: Dict[str, float] = {}

        parsed = sbml.readSBML(self.filename)
        model: libsbml.Model = parsed.getModel()

        # Initialise species and their initial amounts
        for species in model.getListOfSpecies():
            net.species[species.getId()] = species.getInitialAmount()

        # 1. Parameters are mere constants , you can read their values
        # 2. Rules are more complex, but still constants: Though they need to be evaluated
        # 3. Some parameters are actually rules, those are marked with "constant = false"
        # 4. Thus: Parameters (intersection) Rules = Rules

        # Evaluate and store global parameters in a symbol table
        for param in model.getListOfParameters():
            if param.getConstant():
                symbols[param.getId()] = param.getValue()

        # Evaluate and store rules in the symbol table
        rules = model.getListOfRules()
        for r in rules:
            val = helper.evaluate_ast_node(r.getMath(), symbols)
            symbols[r.getId()] = val

        net.symbols = symbols

        reactions: List[Reaction] = []
        for x in model.getListOfReactions():
            reaction_rate_function = x.getKineticLaw().getMath()

            reactants = x.getListOfReactants()
            products = x.getListOfProducts()

            left: str = reactants[0].getSpecies() if reactants else ""
            right: str = products[0].getSpecies() if products else ""

            reactions.append(CustomReaction(reaction_rate_function, left, right))
        net.reactions = reactions

        return net


def parsing():
    p = SbmlParser("other_files/BIOMD0000000012.xml")
    p.parse()


parsing()
