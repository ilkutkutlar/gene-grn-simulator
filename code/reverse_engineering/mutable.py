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

    def get_next_value(self):
        if self.current_value + self.increments > self.upper_bound:
            return None
        else:
            return self.current_value + self.increments

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
