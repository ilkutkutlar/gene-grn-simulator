import math

import libsbml
from libsbml._libsbml import formulaToL3String

"""
Return evaluation of the given string equation
:param str string_equation:
:param Dict[str, float] species: species concentrations
:param Dict[str, float] symbols: symbols from SBML model
:param Dict[str, float] parameters: parameters from SBML reaction
:returns float of equation's evaluation
"""


def eval_equation(string_equation,
                  species=None, symbols=None,
                  parameters=None):
    temp = symbols.copy() if symbols else dict()

    if species is not None:
        temp.update(species)

    if parameters is not None:
        temp.update(parameters)

    return eval(string_equation, temp)


"""
Return a string representation of given AST node
:param libsbml.ASTNode ast_node: libsml AST node
:returns string representation of given AST node
"""


def ast_to_string(ast_node):
    raw = formulaToL3String(ast_node)
    raw = raw.replace("^", "**")
    return raw


"""
:param Dict[str, float] symbols: symbols to use in evaluation
:param Dict[str, float] species: species concentrations to be used to evaluate
:param libsbml.ASTNode node: to evaluate
:returns float of ast's evaluation
Only evaluates nodes which contain constants.
"""


def evaluate_ast(node, symbols, species=None):
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
