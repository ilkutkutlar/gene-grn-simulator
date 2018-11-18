import math
from typing import Dict

import libsbml


def apply_change_vector(state: Dict[str, float], change: Dict[str, float]):
    ret = state.copy()
    for x in state:
        ret[x] = state[x] + change[x]
    return ret


"""
Only evaluates nodes which contain constants.
"""


def evaluate_ast_node(node: libsbml.ASTNode,
                      symbols: Dict[str, float],
                      species: Dict[str, float] = None):
    node_type = node.getType()

    if node_type == libsbml.AST_REAL:
        return node.getValue()
    elif node_type == libsbml.AST_INTEGER:
        return node.getValue()
    elif node_type == libsbml.AST_RATIONAL:
        return node.getValue()
    elif node_type == libsbml.AST_NAME:
        name = node.getName()
        if name in symbols:
            return symbols[name]
        elif species and (name in species):
            return species[name]
        else:
            return 0  # TODO: Error!
    elif node_type == libsbml.AST_FUNCTION_LN:
        val = evaluate_ast_node(node.getLeftChild(), symbols, species=species)
        return math.log(val, math.e)
    else:
        left = evaluate_ast_node(node.getLeftChild(), symbols, species=species)
        right = evaluate_ast_node(node.getRightChild(), symbols, species=species)

        if node_type == libsbml.AST_PLUS:
            return left + right
        elif node_type == libsbml.AST_DIVIDE:
            return left / right
        elif node_type == libsbml.AST_MINUS:
            return left - right
        elif node_type == libsbml.AST_TIMES:
            return left * right
        elif node_type == 296:          # AST_POWER
            return math.pow(left, right)
        else:
            # Debug codes:
            # print("ERROR!")
            # print("Operator: " + str(node.getOperatorName()))
            # print("Value: " + str(node.getValue()))
            # print("Name: " + str(node.getName()))
            # print("Exponent: " + str(node.getExponent()))
            # print("Type: " + str(node.getType()))

            return 0  # TODO: Error!
