# Script by Carla Silva 2019 :: Graph Coloring Quantum Version

""" Graph Coloring Problem

    Formulation of the problem for a graph G=(V,E) with a number of colors n.

"""

import networkx as nx
import dwave_networkx as dnx
import matplotlib.pyplot as plt

if __name__ == "__main__":
    f = plt.figure()
    G=dnx.chimera_graph(16, 16, 4)  # Draw a Chimera unit cell
    dnx.draw_chimera(G, node_size=5, node_color='b', node_shape='o', style='-', edge_color='k', width=0.5)
    f.savefig('chimera.pdf', bbox_inches='tight')
