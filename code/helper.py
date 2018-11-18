import math
from typing import Dict

import libsbml

from models import NamedVector


def apply_change_vector(state: NamedVector, change: NamedVector):
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
    elif node_type == libsbml.AST_NAME:
        return symbols[node.getName()]
    elif node_type == libsbml.AST_FUNCTION_LN:
        val = evaluate_ast_node(node.getLeftChild(), symbols)
        return math.log(val, math.e)
    else:
        left = evaluate_ast_node(node.getLeftChild(), symbols)
        right = evaluate_ast_node(node.getRightChild(), symbols)

        if node_type == libsbml.AST_PLUS:
            return left + right
        elif node_type == libsbml.AST_DIVIDE:
            return left / right
        elif node_type == libsbml.AST_MINUS:
            return left - right
        elif node_type == libsbml.AST_TIMES:
            return left * right
        else:
            return 0
