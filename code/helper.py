import math

import libsbml
from PyQt5.QtGui import QDoubleValidator
from PyQt5.QtWidgets import QMessageBox
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
    if node.isReal() or node.isInteger() or node.isRational():
        value = node.getValue()
    elif node_type == libsbml.AST_NAME:
        name = node.getName()
        if name in symbols:
            value = symbols[name]
        elif species and (name in species):
            value = species[name]
        else:
            value = False
    elif node_type == libsbml.AST_FUNCTION_LN:
        val = evaluate_ast(node.getLeftChild(), symbols, species=species)
        if val:
            value = math.log(val, math.e)
        else:
            value = False
    else:
        left = evaluate_ast(node.getLeftChild(), symbols, species=species)
        right = evaluate_ast(node.getRightChild(), symbols, species=species)

        if not left or not right:
            value = False
        elif node_type == libsbml.AST_PLUS:
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
            value = False

    return value


def get_double_validator():
    v = QDoubleValidator()
    v.setNotation(QDoubleValidator.StandardNotation)
    return v


def show_error_message(message):
    error_message = QMessageBox()
    error_message.setIcon(QMessageBox.Warning)
    error_message.setWindowTitle("Error")
    error_message.setStandardButtons(QMessageBox.Ok)
    error_message.setText(message)

    button = error_message.exec_()
    if button == QMessageBox.Ok:
        return True
    else:
        return False
