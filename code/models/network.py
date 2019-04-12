from models.formulae.transcription_formula import TranscriptionFormula
from models.input_gate import InputGate
from models.regulation import Regulation
from constraint_satisfaction.mutable import ReactionMutable, VariableMutable, RegulationMutable


class Network:

    def __init__(self):
        self.species = dict()  # of Dict[str, float]
        self.reactions = list()  # of Reaction

    """
    Change species concentrations of network using a change vector
    :param Dict[str, float] change: key: species name, value: concentration change by
    """

    def apply_change_vector(self, change):
        for x in change:
            self.species[x] += change[x]

    """
    Mutate species and reactions of the network
    :param Dict[str, Tuple[float, str]] mutations: key: mutable name, value: (new value, reaction name)
    """

    def mutate(self, mutations):
        for m in mutations:
            if isinstance(m, ReactionMutable):
                r = self.get_reaction_by_name(m.reaction_name)
                r.rate_function.mutate(m)
            elif isinstance(m, VariableMutable):
                self.species[m.variable_name] = m.current_value
            elif isinstance(m, RegulationMutable):
                maybe_reaction = \
                    list(filter(lambda x: x.name == m.reaction_name, self.reactions))
                if maybe_reaction:
                    transcription = maybe_reaction[0].rate_function

                    if m.is_installed:
                        reg = transcription.get_regulation(m.possible_regulators[m.current_regulator])
                        if reg:
                            reg.reg_type = m.possible_reg_types[m.current_reg_type]
                            reg.k = m.k_variable.current_value
                        else:
                            # TODO
                            transcription.set_regulation(2, [], InputGate.AND)

                            new_reg = Regulation(m.possible_regulators[m.current_regulator], transcription.transcribed_species,
                                                 m.possible_reg_types[m.current_reg_type], m.k_variable.current_value)

                            transcription.regulators.append(new_reg)
                    else:
                        reg = transcription.get_regulation(m.current_regulator)
                        if reg:
                            transcription.regulators.remove(reg)
                else:
                    pass  # error

    """
    Return reaction with given name
    :param str name: Name of reaction
    :returns Reaction if found, None if not
    """

    def get_reaction_by_name(self, name):
        t = list(filter(lambda r: r.name == name, self.reactions))
        if t:
            return t[0]
        else:
            return None

    def __str__(self):
        ret = "\nSpecies: \n"
        for x in self.species:
            ret += "    " + x + ": " + str(self.species[x]) + "\n"

        ret += "\nReactions: \n"
        for x in self.reactions:
            ret += "    " + str(x) + "\n"

        return ret
