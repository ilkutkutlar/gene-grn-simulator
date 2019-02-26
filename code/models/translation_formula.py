from models.formula import Formula


class TranslationFormula(Formula):
    """
    :param float rate:
    :param str mrn_species:
    """

    def __init__(self, rate, mrna_species):
        self.rate = rate
        self.mrna_species = mrna_species

    def compute(self, state):
        return self.rate * state[self.mrna_species]

    def mutate(self, mutation):
        for m in mutation:
            if m == "rate":
                self.rate = mutation[m][0]

    def get_params(self):
        return ["rate"]
