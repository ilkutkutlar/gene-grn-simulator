import networkx as nx
import matplotlib.pyplot as plt
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
                        arrow = 'tee' if reg.reg_type == RegType.REPRESSION else "normal"
                        g.edge(reg.from_gene, reg.to_gene, arrowhead=arrow)

        return g

    @staticmethod
    def visualise(net):
        g = NetworkVisualiser._reaction_network_to_networkx(net)
        pos = nx.layout.spring_layout(g)

        a = ArrowStyle.BarAB(widthA=0.0, angleA=None, widthB=1.7, angleB=None)
        nx.draw_networkx(g, pos, arrows=False, with_labels=False, node_color='#561a1a', node_size=800)
        nx.draw_networkx_edges(g, pos, arrowstyle=a, arrowsize=10, width=2.0)
        nx.draw_networkx_labels(g, pos, font_size=16, font_color='#515151')

        ax = plt.gca()
        ax.collections[0].set_edgecolor("#FFFFFF")
        ax.set_axis_off()

        plt.savefig("temp.png")


def main():
    G = nx.DiGraph()

    G.add_node("t1")
    G.add_node("t2")
    G.add_node("t3")

    G.add_edge("t1", "t2")

    plt.figure(figsize=(5, 5))

    pos = nx.layout.spring_layout(G)

    # TODO: https://www.kaggle.com/jncharon/python-network-graph

    a = ArrowStyle.BarAB(widthA=0.0, angleA=None, widthB=1.7, angleB=None)
    nx.draw_networkx(G, pos, arrows=False, with_labels=False, node_color='#561a1a', node_size=800)

    nx.draw_networkx_edges(G, pos, arrowstyle=a, arrowsize=10, width=2.0, nodelist=['x', 'py'])

    nx.draw_networkx_labels(G, pos, font_size=16, font_color='#515151')

    ax = plt.gca()
    ax.collections[0].set_edgecolor("#FFFFFF")
    ax.set_axis_off()

    plt.show()


def main2():
    # https://networkx.github.io/documentation/stable/auto_examples/drawing/plot_directed.html#sphx-glr-auto-examples-drawing-plot-directed-py
    import matplotlib.pyplot as plt
    import networkx as nx

    G = nx.generators.directed.random_k_out_graph(10, 3, 0.5)
    pos = nx.layout.spring_layout(G)

    node_sizes = [i * 0 + 500 for i in range(len(G))]
    M = G.number_of_edges()
    edge_colors = range(2, M + 2)

    # |-| for repression,
    # -> for activation
    nodes = nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color='blue', with_labels=True)
    edges = nx.draw_networkx_edges(G, pos, node_size=node_sizes, arrowstyle='|-|',
                                   arrowsize=10, edge_color=edge_colors,
                                   edge_cmap=plt.cm.Blues, width=2)

    ax = plt.gca()
    ax.set_axis_off()
    plt.show()


def get_test():
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

    net: Network = Network()
    net.species = species
    net.reactions = reactions

    return net


# NetworkVisualiser.visualise(get_test())

dot = Digraph(format='png')
dot.node("1", "X")
dot.node("2", "Y")
dot.node("3", "Z")
dot.edge("1", "2", arrowhead='tee')
dot.edge("2", "3")
dot.render("Test")
