# Script by Carla Silva and Inês Dutra 2020 :: Graph Coloring Classical Version

""" Graph Coloring Problem

    Formulation of the problem for a graph G=(V,E) with a number of colors n.

"""

from pyqubo import Array, solve_qubo, Constraint
import matplotlib.pyplot as plt
import networkx as nx
import time
import sys
import random
import numpy as np
import seaborn as sns
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
import neal
import pandas as pd
from dwave.system import EmbeddingComposite, LazyFixedEmbeddingComposite
import dwave_networkx as dnx
import dimod

def alpha(n,h):
  return(n+h)

def beta(n,alfa,h):
  return(((n**3+n**2+1)*alfa)+h)

def exp1(n, vc, B):
    exp = 0.0
    for i in range(n):
        for k in range(n):
            for kl in range(n):
                if (k != kl):
                    exp += Constraint(vc[i,k]*vc[i,kl], label="exp")
    exp = B*exp
    return(exp)

def exp2(n, vc, A):
    exp = 0.0
    for i in range(n):
      for k in range(n):
          exp += Constraint(1-vc[i,k], label="exp") 
    exp = A*exp
    return(exp)

def exp3(n, vc, neigh, A):
    exp = 0.0
    for i in range(n):
        for il in range(n):
            if (i != il):
                for k in range(n):
                   exp += Constraint(vc[i,k]*vc[il,k]*neigh[i,il], label="exp") 
    exp = A*exp
    return(exp)

def exp4(n, vc, c, A):
    exp = 0.0
    for i in range(n):
      for k in range(n):
        exp += Constraint(vc[i,k]*(1-c[k]), label="exp")
    exp = A*exp
    return(exp)

def exp5(n, c):
    exp = 0.0
    for k in range(n):
        exp += c[k]
    return(exp)

def col(n, decoded_solution):
    # Obtain colors of each vertex
    colors = [0 for i in range(n)]
    S = np.zeros((n,n))
    for name, value in decoded_solution.items():
        if 'vc' in name:
            S[int(name[3])][int(name[6])] = value
    for i in range(n):
        for k in range(n):
            if S[i][k] == 1:
                colors[i] = k
                break
    return(colors)

def graph(G, colors, n, A, B):
    # Plot graph after coloring
    f = plt.figure()
    colorlist = list(sns.color_palette("hls", n_colors=n)) 
    nx.draw_networkx(G, node_color=[colorlist[colors[node]] for node in G.nodes], node_size=400, font_weight='bold', font_color='w')
    plt.axis("off")
    f.savefig('fig-n'+str(n)+'A'+str(A)+'B'+str(B)+'.pdf', bbox_inches='tight')
    dictionaryColors = dict(zip(list(G.nodes), [colorlist[colors[node]] for node in G.nodes]))
    return dictionaryColors

    
if __name__ == "__main__":


    n = int(sys.argv[1]) # Colors
    A = int(sys.argv[2]) # Alpha 
    B = int(sys.argv[3]) # Beta

    h = 0.0000005 #small number

    #Test A and B equals alpha and beta
    #A = alpha(n,h)
    #B = beta(n,A,h)

    orig_stdout = sys.stdout

    f = open('graphColoringClassicalResults-n'+str(n)+'A'+str(A)+'B'+str(B)+'.txt', 'w')
    sys.stdout = f

    print("--------------------------------------------------------------------")
    print("\n# GRAPH COLORInG PROBLEM WITH n COLOURS On CLASSICAL SOLVER #\n")
    print("--------------------------------------------------------------------")

    G = nx.erdos_renyi_graph(n=n, p=0.5, seed=123, directed=False)

    print("Is graph connected?",nx.is_connected(G))

    E = G.edges

    neigh = nx.adjacency_matrix(G).todense()

    # Prepare a binary vector

    vc = Array.create('vc', (n, n), 'BINARY')

    c = Array.create('c', (n), 'BINARY')


    print("--------------------------------------------------------------------")
    print("1st expression:")
    print("--------------------------------------------------------------------")
    print(exp1(n, vc, B))

    print("--------------------------------------------------------------------")
    print("2nd expression:")
    print("--------------------------------------------------------------------")
    print(exp2(n, vc, A))

    print("--------------------------------------------------------------------")
    print("3rd expression:")
    print("--------------------------------------------------------------------")
    print(exp3(n, vc, neigh, A))

    print("--------------------------------------------------------------------")
    print("4th expression:")
    print("--------------------------------------------------------------------")
    print(exp4(n, vc, c, A))

    print("--------------------------------------------------------------------")
    print("5th expression:")
    print("--------------------------------------------------------------------")
    print(exp5(n, c))

    # Define hamiltonian H
    H = exp1(n, vc, B) + exp2(n, vc, A) + exp3(n, vc, neigh, A) + exp4(n, vc, c, A) + exp5(n, c)

    # Compile model
    model = H.compile()

    # Create model
    qubo, offset = model.to_qubo()

    print("--------------------------------------------------------------------")
    print("\nQUBO:\n")
    print("--------------------------------------------------------------------")

    print(qubo)

    start_time = time.time()

    nr = 10000

    c = max(qubo.values()) 

    sa = neal.SimulatedAnnealingSampler()
    Gc = dnx.chimera_graph(16, 16, 4) # Chimera graph
    composite = dimod.StructureComposite(sa, Gc.nodes, Gc.edges)
    sampler = LazyFixedEmbeddingComposite(composite)
    response = sampler.sample_qubo(qubo, num_reads=nr, offset=offset, num_sweeps=5000, chain_strength=c, seed=123)

    elapsed_time = time.time() - start_time

    print("--------------------------------------------------------------------")
    print("\nCLASSICAL RESULTS:\n")
    print("--------------------------------------------------------------------")

    minE = sys.maxsize
    maxO = 0
    # create dataframe if we want to store all values
    df = []
    count = 0
    for datum in response.data(['sample', 'energy', 'num_occurrences','chain_break_fraction']):
        if (datum.energy < minE):
            count = count + 1
            minE = datum.energy
            maxO = datum.num_occurrences
            sample = datum.sample
            chain = datum.chain_break_fraction

        df.append({"Sample": datum.sample, "Energy": datum.energy, "Occurrences": datum.num_occurrences, "Chain_break_fractions": datum.chain_break_fraction})
        print(datum.sample, "Energy: ", datum.energy, "Occurrences: ", datum.num_occurrences,"Chain break fractions:", datum.chain_break_fraction)

    df = pd.DataFrame(df)
    df.to_csv('GC'+str(n)+'A'+str(A)+'B'+str(B)+'.csv',index=False)

    # Match colors with each node
    colors = col(n, sample)

    # Plot the colored graph
    nodesColor = graph(G, colors, n, int(A), int(B))

    print("--------------------------------------------------------------------")
    print("\nSAMPLE WITH MINIMUM ENERGY AND MAXIMUM OCCURRENCES:\n")
    print("--------------------------------------------------------------------")
    print(sample, "Energy: ", minE, "Occurrences: ", maxO, "Chain break fractions:", chain)

    print("--------------------------------------------------------------------")
    print("\nNumber of possible solutions:\n")
    print("--------------------------------------------------------------------")
    print(count)

    print("--------------------------------------------------------------------")
    print("\nTIME (sec):\n")
    print("--------------------------------------------------------------------")
    print(elapsed_time,"All results",elapsed_time/nr,"One result")

    print("--------------------------------------------------------------------")
    print("\nList of graph colors:\n")
    print("--------------------------------------------------------------------")
    print(nodesColor)

    print("--------------------------------------------------------------------")
    print("\nNumber of nodes and number of colors:\n")
    print("--------------------------------------------------------------------")
    print(len(nodesColor)," nodes and ",len(set(list(nodesColor.values())))," colors")

    print("Chain strength:")
    print(c)

    print("Embedding (variables mapped to physical qubits):")
    print(sampler.properties['embedding'])

    sys.stdout = orig_stdout
    f.close()
