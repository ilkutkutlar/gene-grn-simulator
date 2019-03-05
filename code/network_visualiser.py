import networkx as nx
import matplotlib.pyplot as plt
from PyQt5.QtGui import QImage
from graphviz import Digraph

from matplotlib.patches import ArrowStyle

from models.formulae.degradation_formula import DegradationFormula
from models.formulae.transcription_formula import TranscriptionFormula
from models.formulae.translation_formula import TranslationFormula
from models.network import Network
from models.reaction import Reaction
from models.reg_type import RegType
from models.regulation import Regulation


class NetworkVisualiser:
    """

    dot = Digraph(format='png')
dot.node("1", "X")
dot.node("2", "Y")
dot.node("3", "Z")
dot.edge("1", "2", arrowhead='tee')
dot.edge("2", "3")

    """

    @staticmethod
    def _reaction_network_to_networkx(net):
        g = Digraph(format='png')

        for s in net.species:
            g.node(s, s)

        for r in net.reactions:
            if isinstance(r.rate_function, TranscriptionFormula):
                if r.rate_function.regulators:
                    for reg in r.rate_function.regulators:
                        arrow = "normal" if reg.reg_type == RegType.ACTIVATION else "tee"
                        g.edge(reg.from_gene, reg.to_gene, arrowhead=arrow)

        return g

    @staticmethod
    def visualise(net):
        g = NetworkVisualiser._reaction_network_to_networkx(net)
        g.render("temp")
        im = QImage("temp.png")
        return im
