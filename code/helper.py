import math
from typing import Dict

import libsbml
from libsbml._libsbml import formulaToL3String


def eval_equation(string: str, symbols: Dict[str, float],
                  species: Dict[str, float] = None, parameters: Dict[str, float] = None):

    temp = symbols.copy() if symbols else dict()

    if species is not None:
        temp.update(species)

    if parameters is not None:
        temp.update(parameters)

    return eval(string, temp)


def ast_to_string(ast_node: libsbml.ASTNode):
    raw: str = formulaToL3String(ast_node)
    raw = raw.replace("^", "**")
    return raw


"""
Only evaluates nodes which contain constants.
"""


def evaluate_ast(node: libsbml.ASTNode, symbols: Dict[str, float], species: Dict[str, float] = None) -> float:
    node_type = node.getType()

    # Base cases
    # TODO: Rational numbers are stored differently: using numerator and denominator, account for that!
    if node.isReal() or node.isInteger() or node.isRational():
        value = node.getValue()
    elif node_type == libsbml.AST_NAME:
        name = node.getName()
        if name in symbols:
            value = symbols[name]
        elif species and (name in species):
            value = species[name]
        else:
            value = 0  # TODO: Error!
    elif node_type == libsbml.AST_FUNCTION_LN:
        val = evaluate_ast(node.getLeftChild(), symbols, species=species)
        value = math.log(val, math.e)
    else:
        left = evaluate_ast(node.getLeftChild(), symbols, species=species)
        right = evaluate_ast(node.getRightChild(), symbols, species=species)

        if node_type == libsbml.AST_PLUS:
            value = left + right
        elif node_type == libsbml.AST_DIVIDE:
            value = left / right
        elif node_type == libsbml.AST_MINUS:
            value = left - right
        elif node_type == libsbml.AST_TIMES:
            value = left * right
        elif node_type == 296:  # AST_POWER
            value = math.pow(left, right)
        else:
            # Debug codes:
            # print("ERROR!")
            # print("Operator: " + str(node.getOperatorName()))
            # print("Value: " + str(node.getValue()))
            # print("Name: " + str(node.getName()))
            # print("Exponent: " + str(node.getExponent()))
            # print("Type: " + str(node.getType()))

            value = 0  # TODO: Error!

    return value
