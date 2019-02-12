from math import log, e

from models.formulae import DegradationFormula, TranslationFormula, TranscriptionFormula, CustomFormula
from models.network import Network
from models.reaction import Reaction
from models.reg_type import RegType
from models.regulation import Regulation
from models.simulation_settings import SimulationSettings
from reverse_engineering.reverse_engineering import Constraint, is_satisfied
from simulation.ode_simulator import OdeSimulator
from structured_results import StructuredResults


def main():
    # p_lacI0 = 10
    # p_tetR0 = 10
    # p_cl0 = 10

    # m_lacI0 = 100
    # m_tetR0 = 80
    # m_cl0 = 50

    # region Constants

    ps_active = 0.5  # Promoter strength (active)
    ps_repr = 5 * (10 ** -4)  # Promoter strength (repressed)
    transcription_rate_active = ps_active * 60
    transcription_rate_repr = ps_repr * 60

    # Decay-related values
    mRNA_half_life = 2
    protein_half_life = 10
    mRNA_decay_rate = log(2, e) / mRNA_half_life
    protein_decay_rate = log(2, e) / protein_half_life

    # Translation-related values
    translation_efficiency = 19.97
    translation_rate = translation_efficiency * mRNA_decay_rate

    # Other values
    hill_coeff = 2  # Hill coefficient
    Km = 40  # Activation coefficient

    # Derived values

    # Active transcription rate (rescaled?)
    alpha = transcription_rate_active * translation_efficiency * protein_half_life / (log(2, e) * Km)
    # Repressed transcription rate (rescaled?)
    alpha0 = transcription_rate_repr * translation_efficiency * protein_half_life / (log(2, e) * Km)
    # Translation rate
    beta = protein_decay_rate / mRNA_decay_rate

    # endregion

    species = {"laci_mrna": 100, "tetr_mrna": 80, "cl_mrna": 50,
               "laci_p": 10, "tetr_p": 10, "cl_p": 10}

    laci_reg = Regulation(from_gene="cl_p", to_gene="laci_mrna", reg_type=RegType.REPRESSION)
    tetr_reg = Regulation(from_gene="laci_p", to_gene="tetr_mrna", reg_type=RegType.REPRESSION)
    cl_reg = Regulation(from_gene="tetr_p", to_gene="cl_mrna", reg_type=RegType.REPRESSION)

    reactions = [Reaction("", [], ["laci_mrna"], TranscriptionFormula(alpha, 2, 40, "laci_mrna", [laci_reg])),
                 Reaction("", [], ["tetr_mrna"], TranscriptionFormula(alpha, 2, 40, "tetr_mrna", [tetr_reg])),
                 Reaction("", [], ["cl_mrna"], TranscriptionFormula(alpha, 2, 40, "cl_mrna", [cl_reg])),

                 Reaction("", ["laci_mrna"], [], DegradationFormula(mRNA_decay_rate, "laci_mrna")),
                 Reaction("", ["tetr_mrna"], [], DegradationFormula(mRNA_decay_rate, "tetr_mrna")),
                 Reaction("", ["cl_mrna"], [], DegradationFormula(mRNA_decay_rate, "cl_mrna")),

                 Reaction("", ["laci_mrna"], ["laci_p"], TranslationFormula(beta, "laci_mrna")),
                 Reaction("", ["tetr_mrna"], ["tetr_p"], TranslationFormula(beta, "tetr_mrna")),
                 Reaction("", ["cl_mrna"], ["cl_p"], TranslationFormula(beta, "cl_mrna")),

                 Reaction("", ["laci_p"], [], DegradationFormula(protein_decay_rate, "laci_p")),
                 Reaction("", ["tetr_p"], [], DegradationFormula(protein_decay_rate, "tetr_p")),
                 Reaction("", ["cl_p"], [], DegradationFormula(protein_decay_rate, "cl_p"))]

    net = Network()
    net.species = species
    net.reactions = reactions

    s = SimulationSettings(0, 10 * 60, 1000, ["laci_p", "tetr_p", "cl_p"])

    ode = OdeSimulator(net, s)
    ode.visualise(ode.simulate())


def simpler():
    species = {"x": 0, "y": 20}

    reactions = [Reaction("", [], ["x"], TranscriptionFormula(5, 2, 40, "x", [
        Regulation(from_gene="y", to_gene="x", reg_type=RegType.REPRESSION)])),
                 Reaction("", ["y"], [], DegradationFormula(0.3, "y"))]

    net: Network = Network()
    net.species = species
    net.reactions = reactions

    s = SimulationSettings(0, 100, 100, ["x", "y"])

    ode = OdeSimulator(net, s)
    ode.visualise(ode.simulate())


def test():
    species = {"px": 100, "py": 100, "pz": 30, "x": 100, "y": 25, "z": 20}

    x_trans = TranscriptionFormula(5, 2, 40, "x", [])
    y_trans = TranscriptionFormula(60, 1, 40, "y", [
        Regulation(from_gene="px", to_gene="y", reg_type=RegType.REPRESSION)])
    z_trans = TranscriptionFormula(20, 2, 40, "z", [
        Regulation(from_gene="py", to_gene="z", reg_type=RegType.ACTIVATION)])

    reactions = [Reaction("", [], ["x"], x_trans),
                 Reaction("", [], ["y"], y_trans),
                 Reaction("", [], ["z"], z_trans),

                 Reaction("", ["x"], [], DegradationFormula(0.01, "x")),
                 Reaction("", ["y"], [], DegradationFormula(0.01, "y")),
                 Reaction("", ["z"], [], DegradationFormula(0.1, "z")),

                 Reaction("", ["px"], [], DegradationFormula(0.01, "px")),
                 Reaction("", ["py"], [], DegradationFormula(0.01, "py")),
                 Reaction("", ["pz"], [], DegradationFormula(0.01, "pz")),

                 Reaction("", [], ["px"], TranslationFormula(0.2, "x")),
                 Reaction("", [], ["py"], TranslationFormula(5, "y")),
                 Reaction("", [], ["pz"], TranslationFormula(1, "z"))]

    net: Network = Network()
    net.species = species
    net.reactions = reactions

    s = SimulationSettings(0, 100, 100, ["x", "y", "z"])

    time_space = s.generate_time_space()
    ode = OdeSimulator(net, s)
    res = StructuredResults(ode.simulate(), list(net.species.keys()), time_space)
    c = Constraint("x", lambda x: x < 150, (0, 20))
    is_satisfied(res, [c])
    # ode.visualise(ode.simulate())


def individual_mutation_test():
    t = TranscriptionFormula(10, 2, 40, "x", [])
    print(str(t.rate) + " | " + str(t.hill_coeff) + " | " + str(t.kd))
    t.mutate({"kd": 10, "rate": 3})
    print(str(t.rate) + " | " + str(t.hill_coeff) + " | " + str(t.kd))
    y = CustomFormula("4*x", {"x": 10}, {})
    y.mutate({"x": 20})
    print(str(y.compute({})))


def mutation_network_test():
    species = {"x": 0, "y": 20}

    reactions = [Reaction("", [], ["x"], TranscriptionFormula(5, 2, 40, "x", [
        Regulation(from_gene="y", to_gene="x", reg_type=RegType.REPRESSION)])),
                 Reaction("", ["y"], [], DegradationFormula(0.3, "y"))]

    net: Network = Network()
    net.species = species
    net.reactions = reactions
    print(net)
    net.mutate({"x": (10.0, "")})
    print(net)


if __name__ == '__main__':
    # simpler()
    # main()
    test()
