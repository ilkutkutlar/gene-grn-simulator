import unittest
from typing import Dict

import helper
from gillespie import GillespieSimulator
from models import Network, CustomReaction, SimulationSettings


class GillespieTest(unittest.TestCase):
    @staticmethod
    def get_test_network():
        net = Network()
        net.species = {"X": 15, "Y": 10, "Z": 20}

        r1 = CustomReaction("(X + Y) * (Z/5)", ["X", "Z"], ["Y"])
        r2 = CustomReaction("X + (Z - 10)", [""], ["Y"])
        r3 = CustomReaction("X - Y", ["X"], [""])

        net.reactions = [r1, r2, r3]
        return net

    @staticmethod
    def get_test_simulation_settings():
        return SimulationSettings("Results", "Time", "Concentration", 0, 10, [("X", "X"), ("Y", "Y"), ("Z", "Z")])

    def test_r0(self):
        net = self.get_test_network()

        # r1 reaction rate: 100
        # r2 reaction rate: 25
        # r3 reaction rate: 5
        # r0 = 100 + 25 + 5 = 130

        self.assertEqual(130, GillespieSimulator._calculate_r0_(net))

    def test_apply_change_vector(self):
        v1: Dict[str, float] = {"X": 10, "Y": 20, "Z": 35}
        v2: Dict[str, float] = {"X": 10, "Y": 23, "Z": 50}

        self.assertEqual(GillespieSimulator._apply_change_vector_(v1, v2),
                         {"X": 20, "Y": 43, "Z": 85})

    def test_evaluate_ast_string(self):
        symbols = {"a": 10, "b": 10, "c": 10}
        species = {"X": 15, "Y": 10, "Z": 20}
        string = "(X + a) - (Y * b) - (Z * c)"  # -275

        self.assertEqual(helper.evaluate_ast_string(string, symbols, species),
                         -275)

    def test_change_vector(self):
        net = self.get_test_network()
        change = net.reactions[0]

        # reaction: X -> Y + Z  |  with rate k = (X + Y) * (Z/5)

        r1 = CustomReaction("(X + Y) * (Z/5)", ["X"], ["Y"])



    if __name__ == '__main__':
        unittest.main()
