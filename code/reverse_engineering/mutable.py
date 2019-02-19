class Mutable:
    """
        :param float lower_bound, upper_bound: defines the interval of values to be tried for this parameter.
        :param float increments: increments define the step between lower_bound & upper_bound,
            e.g. l = 10, u = 11, increment = 0.5 would produce 10, 10.5 and 11 as the values to be tried for this parameter.
        :param str reaction_name: If the parameter is defined in a reaction, this is the reaction_name where
            the parameter is. Else, it is empty string ("")
    """

    def __init__(self, lower_bound, upper_bound, increments, reaction_name):
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.increments = increments
        self.reaction_name = reaction_name

    def __str__(self):
        return "(" + str(self.lower_bound) + " to " + str(self.upper_bound) + ") step: " + str(self.increments)
