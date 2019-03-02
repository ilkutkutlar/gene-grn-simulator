class VariableMutable:
    """
            :param float lower_bound, upper_bound: defines the interval of values to be tried for this parameter.
            :param float increments: increments define the step between lower_bound & upper_bound,
                e.g. l = 10, u = 11, increment = 0.5 would produce 10, 10.5 and 11 as the values to be tried for this parameter.
        """

    def __init__(self, variable_name, lower_bound, upper_bound, increments):
        self.variable_name = variable_name
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.increments = increments
        self.current_value = self.lower_bound

    def is_next(self):
        if self.current_value + self.increments > self.upper_bound:
            return False
        else:
            return True

    def next(self):
        if self.is_next():
            self.current_value += self.increments
            return True
        else:
            return False

    def __str__(self):
        lo = str(self.lower_bound)
        hi = str(self.upper_bound)
        step = str(self.increments)
        return "({} to {}) step: {}".format(lo, hi, step)


class ReactionMutable(VariableMutable):
    """
    :param str reaction_name: If the parameter is defined in a reaction, this is the reaction_name where
                the parameter is. Else, it is empty string ("")
    """

    def __init__(self, variable_name, lower_bound, upper_bound, increments, reaction_name):
        super().__init__(variable_name, lower_bound, upper_bound, increments)
        self.reaction_name = reaction_name

    def __str__(self):
        lo = str(self.lower_bound)
        hi = str(self.upper_bound)
        step = str(self.increments)
        return "({} to {}) step: {} in {}".format(lo, hi, step, self.reaction_name)


class RegulationMutable:
    """
    :param str reaction_name: The name of the reaction which can potentially do regulation
    :param List[str] to_species: The list of possible species which this reaction can potentially regulate
    :param VariableMutable k_variable: k variable mutable associated with this regulation
    :param List[RegType] possible_reg_types: All possible regulation types accepted
    :param bool is_installed: Whether this regulation is installed in the network
    """

    def __init__(self, reaction_name, to_species, k_variable, possible_reg_types, is_installed):
        self.reaction_name = reaction_name
        self.to_species = to_species
        self.possible_reg_types = possible_reg_types
        self.k_variable = k_variable

        self.is_installed = is_installed
        self.current_to = 0 if to_species else None
        self.current_reg_type = 0 if possible_reg_types else None

    def is_next(self):
        if not self.is_installed:
            return True
        else:
            # for x in to_species:
            #   for y in reg_types:
            #       for z in k:

            if self.k_variable.next():
                return True
            else:
                if self.current_reg_type < len(self.possible_reg_types):
                    return True
                else:
                    if self.current_to < len(self.to_species):
                        return True
                    else:
                        return False

    def next(self):
        if not self.is_installed:
            self.is_installed = True
            return True
        else:
            # for x in to_species:
            #   for y in reg_types:
            #       for z in k:

            if self.k_variable.next():
                self.k_variable.next()
                return True
            else:
                if self.current_reg_type < len(self.possible_reg_types):
                    self.current_reg_type += 1
                    return True
                else:
                    if self.current_to < len(self.to_species):
                        self.current_to += 1
                        return True
                    else:
                        return False
