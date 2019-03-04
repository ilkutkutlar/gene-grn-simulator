from __future__ import division

from matplotlib.patches import ArrowStyle


def main():
    import networkx as nx
    import matplotlib.pyplot as plt

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
    nx.draw_networkx_edges(G, pos, arrowstyle=a, arrowsize=10, width=2.0)
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


main()
