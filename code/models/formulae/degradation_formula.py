from models.formulae.formula import Formula


class DegradationFormula(Formula):
    """
    :param float rate:
    :param str decaying_species:
    """

    def __init__(self, rate, decaying_species):
        self.rate = rate
        self.decaying_species = decaying_species

    def compute(self, state):
        return self.rate * state[self.decaying_species]

    def mutate(self, mutation):
        for m in mutation:
            if m == "rate":
                self.rate = mutation[m][0]

    def get_params(self):
        return ["rate"]

    def get_formula_string(self):
        return str(self.rate)

    def __str__(self):
        rate = str(self.rate)

        string = "Type: Degradation" + "\n"
        string += "Rate: " + rate + "\n"

        return string
