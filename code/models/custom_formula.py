import helper
from models.formula import Formula


class CustomFormula(Formula):
    """
    :param str rate_function:
    :param Dict[str, float] parameters:
    :param Dict[str, float] symbols:
    """

    def __init__(self, rate_function, parameters, symbols):
        self.rate_function = rate_function
        self.symbols = symbols
        self.parameters = parameters

    def compute(self, state):
        return helper.eval_equation(self.rate_function,
                                    species=state,
                                    symbols=self.symbols,
                                    parameters=self.parameters)

    def mutate(self, mutation):
        for m in mutation:
            self.parameters.update({m: mutation[m][0]})
