import math
from typing import Dict

import libsbml

from models import Network


# 1. Core objects of libsbml:
# http://sbml.org/Software/libSBML/5.17.0/docs//python-api/group__core.html
# 2. Classes:
# http://sbml.org/Software/libSBML/5.17.0/docs//python-api/annotated.html


class SbmlParser:
    def __init__(self, filename):
        self.filename = filename

    """
    Only evaluates nodes which contain constants.
    """

    def evaluate_ast_node(self, node: libsbml.ASTNode, symbols: Dict[str, float]):
        node_type = node.getType()

        if node_type == libsbml.AST_NAME:
            return symbols[node.getName()]
        elif node_type == libsbml.AST_INTEGER:
            return node.getValue()
        elif node_type == libsbml.AST_REAL:
            return node.getValue()
        elif node_type == libsbml.AST_FUNCTION_LN:
            val = self.evaluate_ast_node(node.getLeftChild(), symbols)
            return math.log(val, math.e)
        elif node_type == libsbml.AST_PLUS:
            left = self.evaluate_ast_node(node.getLeftChild(), symbols)
            right = self.evaluate_ast_node(node.getRightChild(), symbols)
            return left + right
        elif node_type == libsbml.AST_DIVIDE:
            left = self.evaluate_ast_node(node.getLeftChild(), symbols)
            right = self.evaluate_ast_node(node.getRightChild(), symbols)
            return left / right
        elif node_type == libsbml.AST_MINUS:
            left = self.evaluate_ast_node(node.getLeftChild(), symbols)
            right = self.evaluate_ast_node(node.getRightChild(), symbols)
            return left - right
        elif node_type == libsbml.AST_TIMES:
            left = self.evaluate_ast_node(node.getLeftChild(), symbols)
            right = self.evaluate_ast_node(node.getRightChild(), symbols)
            return left * right
        else:
            return 0

    def parse(self):
        net: Network = Network()
        sbml = libsbml.SBMLReader()
        symbols: Dict[str, float] = {}

        parsed = sbml.readSBML(self.filename)
        model: libsbml.Model = parsed.getModel()

        for species in model.getListOfSpecies():
            net.species[species.getId()] = species.getInitialAmount()

        # 1. Parameters are mere constants , you can read their values
        # 2. Rules are more complex, but still constants: Though they need to be evaluated
        # 3. Some parameters are actually rules, those are marked with "constant = false"
        # 4. Thus: Parameters (intersection) Rules = Rules

        for param in model.getListOfParameters():
            if param.getConstant():
                symbols[param.getId()] = param.getValue()

        rules = model.getListOfRules()
        for r in rules:
            val = self.evaluate_ast_node(r.getMath(), symbols)
            symbols[r.getId()] = val

        print(symbols)

        # for x in model.getListOfReactions():
        #     formula = x.getKineticLaw().getMath()
        #     if formula.getType() == libsbml.AST_TIMES:
        #         left = formula.getLeftChild().getName()
        #         right = formula.getRightChild().getName()
        #         print(left + "*" + right)
        #     print()
        # print(formulaToL3String(formula))

        # print(net.species)


def parsing():
    p = SbmlParser("other_files/BIOMD0000000012.xml")
    p.parse()


parsing()
