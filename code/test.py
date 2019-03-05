from math import log, e

from models.formulae.degradation_formula import DegradationFormula
from models.formulae.transcription_formula import TranscriptionFormula
from models.formulae.translation_formula import TranslationFormula
from models.input_gate import InputGate
from models.network import Network
from models.reaction import Reaction
from models.reg_type import RegType
from models.regulation import Regulation
from models.simulation_settings import SimulationSettings
from simulation.ode_simulator import OdeSimulator


def get_repressilator():
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

    species = {"tetr_mrna": 80, "cl_mrna": 50,
               "laci_p": 10, "cl_p": 10, "tetr_p": 10, "laci_mrna": 100}

    laci_reg = Regulation("cl_p", "laci_mrna", RegType.REPRESSION, Km)
    tetr_reg = Regulation("laci_p", "tetr_mrna", RegType.REPRESSION, Km)
    cl_reg = Regulation("tetr_p", "cl_mrna", RegType.REPRESSION, Km)

    laci_trans = TranscriptionFormula(alpha, "laci_mrna")
    laci_trans.set_regulation(hill_coeff, [laci_reg])

    tetr_trans = TranscriptionFormula(alpha, "tetr_mrna")
    tetr_trans.set_regulation(hill_coeff, [tetr_reg])

    cl_trans = TranscriptionFormula(alpha, "cl_mrna")
    cl_trans.set_regulation(hill_coeff, [cl_reg])

    reactions = [Reaction("laci_trans", [], ["laci_mrna"], laci_trans),
                 Reaction("tetr_trans", [], ["tetr_mrna"], tetr_trans),
                 Reaction("cl_trans", [], ["cl_mrna"], cl_trans),

                 Reaction("laci_mrna_deg", ["laci_mrna"], [], DegradationFormula(mRNA_decay_rate, "laci_mrna")),
                 Reaction("tetr_mrna_deg", ["tetr_mrna"], [], DegradationFormula(mRNA_decay_rate, "tetr_mrna")),
                 Reaction("cl_mrna_deg", ["cl_mrna"], [], DegradationFormula(mRNA_decay_rate, "cl_mrna")),

                 Reaction("laci_p_trans", ["laci_mrna"], ["laci_p"], TranslationFormula(beta, "laci_mrna")),
                 Reaction("tetr_p_trans", ["tetr_mrna"], ["tetr_p"], TranslationFormula(beta, "tetr_mrna")),
                 Reaction("cl_p_trans", ["cl_mrna"], ["cl_p"], TranslationFormula(beta, "cl_mrna")),

                 Reaction("laci_p_deg", ["laci_p"], [], DegradationFormula(protein_decay_rate, "laci_p")),
                 Reaction("tetr_p_deg", ["tetr_p"], [], DegradationFormula(protein_decay_rate, "tetr_p")),
                 Reaction("cl_p_deg", ["cl_p"], [], DegradationFormula(protein_decay_rate, "cl_p"))]

    repressilator = Network()
    repressilator.species = species
    repressilator.reactions = reactions

    return repressilator


def get_test_network1():
    species = {"px": 100, "py": 100, "pz": 30, "x": 100, "y": 25, "z": 20}

    x_trans = TranscriptionFormula(5, "x")

    y_trans = TranscriptionFormula(60, "y")
    y_trans.set_regulation(1, [Regulation("px", "y", RegType.REPRESSION, 40)])

    z_trans = TranscriptionFormula(20, "z")
    z_trans.set_regulation(2, [Regulation("py", "z", RegType.ACTIVATION, 40)])

    reactions = [Reaction("x_trans", [], ["x"], x_trans),
                 Reaction("y_trans", [], ["y"], y_trans),
                 Reaction("z_trans", [], ["z"], z_trans),

                 Reaction("x_deg", ["x"], [], DegradationFormula(0.01, "x")),
                 Reaction("y_deg", ["y"], [], DegradationFormula(0.01, "y")),
                 Reaction("z_deg", ["z"], [], DegradationFormula(0.1, "z")),

                 Reaction("px_deg", ["px"], [], DegradationFormula(0.01, "px")),
                 Reaction("py_deg", ["py"], [], DegradationFormula(0.01, "py")),
                 Reaction("pz_deg", ["pz"], [], DegradationFormula(0.01, "pz")),

                 Reaction("px_translation", [], ["px"], TranslationFormula(0.2, "x")),
                 Reaction("py_translation", [], ["py"], TranslationFormula(5, "y")),
                 Reaction("pz_translation", [], ["pz"], TranslationFormula(1, "z"))]

    net1: Network = Network()
    net1.species = species
    net1.reactions = reactions

    return net1


def get_test_network2():
    species = {"px": 100, "py": 100, "pz": 30, "x": 25, "y": 25, "z": 25}

    x_trans = TranscriptionFormula(5, "x")

    y_trans = TranscriptionFormula(5, "y")

    z_trans = TranscriptionFormula(5, "z")
    z_trans.set_regulation(2,
                           [Regulation("py", "z", RegType.ACTIVATION, 40),
                            Regulation("px", "z", RegType.ACTIVATION, 80)],
                           InputGate.AND)

    reactions = [Reaction("x_trans", [], ["x"], x_trans),
                 Reaction("y_trans", [], ["y"], y_trans),
                 Reaction("z_trans", [], ["z"], z_trans),

                 # Reaction("x_deg", ["x"], [], DegradationFormula(0.01, "x")),
                 # Reaction("y_deg", ["y"], [], DegradationFormula(0.01, "y")),
                 # Reaction("z_deg", ["z"], [], DegradationFormula(0.01, "z")),

                 # Reaction("px_deg", ["px"], [], DegradationFormula(0.01, "px")),
                 # Reaction("py_deg", ["py"], [], DegradationFormula(0.01, "py")),
                 # Reaction("pz_deg", ["pz"], [], DegradationFormula(0.01, "pz")),

                 Reaction("px_translation", [], ["px"], TranslationFormula(0.2, "x")),
                 Reaction("py_translation", [], ["py"], TranslationFormula(5, "y")),
                 Reaction("pz_translation", [], ["pz"], TranslationFormula(1, "z"))]

    net2 = Network()
    net2.species = species
    net2.reactions = reactions

    return net2


if __name__ == '__main__':
    net = get_repressilator()
    s = SimulationSettings(0, 10 * 60, 1000, ["laci_p", "tetr_p", "cl_p"])
    OdeSimulator.visualise(net, s, OdeSimulator.simulate(net, s))
