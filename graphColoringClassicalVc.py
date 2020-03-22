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
    for i in range(n):
        for k in range(n):
          if decoded_solution['vc'][i][k] == 1:
            colors[i] = k
            break
    return(colors)

def graph(G, colors, n, A, B):
    # Plot graph after coloring
    f = plt.figure()
    colorlist = random.sample(list(sns.color_palette('hls')), n) 
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

    # Create QUBO
    qubo, offset = model.to_qubo()

    print("--------------------------------------------------------------------")
    print("\nQUBO:\n")
    print("--------------------------------------------------------------------")

    print(qubo)

    start_time = time.time()

    # Solve the QUBO and obtain the optimal solution with nr number of reads

    # Number of run repetitions of Simulated Annealing (SA)
    #nr = 10000

    #solution = solve_qubo(qubo, num_reads=nr)

    #elapsed_time = time.time() - start_time

    #print(solution)

    # Decode solution
    #decoded_solution, broken, energy = model.decode_solution(solution, vartype="BINARY")
    #print("number of broken constraint = {}".format(len(broken)))
    #print(broken)
    #print("energy = {}".format(energy))

    # Same as solve_qubo but in deeper to retrive all the results
    nr = 10000
    max_abs_value = float(max(abs(v) for v in qubo.values()))
    scale_qubo = {k: float(v) / max_abs_value for k, v in qubo.items()}
    sa = neal.SimulatedAnnealingSampler()
    sa_computation = sa.sample_qubo(scale_qubo, num_reads=nr)

    elapsed_time = time.time() - start_time

    print("--------------------------------------------------------------------")
    print("\nCLASSICAL RESULTS:\n")
    print("--------------------------------------------------------------------")

    df = []
    for a,b,c in sa_computation.record:
        decoded_solution, broken, energy = model.decode_solution(a, vartype='BINARY')
        if not broken:
            d = 0
        else:
            d = broken["exp"]["penalty"]
        for key in range(0,n):
            decoded_solution['c']['c'+str(key)] = decoded_solution['c'].pop(key)
        for i in range(0,n):
            decoded_solution['vc']['vc'+str(i)] = decoded_solution['vc'].pop(i)
        for i in range(0,n):
            for j in range(0,n):
                decoded_solution['vc'].update( {'vc'+str(i)+str(j) : decoded_solution['vc']['vc'+str(i)][j]} )
            decoded_solution['vc'].pop('vc'+str(i))
        decoded_solution['c'].update(decoded_solution['vc'])
        df.append({"Sample": decoded_solution['c'], "Energy": energy, "Occurrences": c, "Broken chains": d})

    df = pd.DataFrame(df)
    df.to_csv('GC'+str(n)+'A'+str(A)+'B'+str(B)+'.csv',index=False)
    pd.set_option('display.float_format', lambda x: '%.20f' % x)
    pd.options.display.max_colwidth = 10000
    print(df.to_string(index=False))

    print("--------------------------------------------------------------------")
    print("\nSAMPLE WITH MINIMUM ENERGY:\n")
    print("--------------------------------------------------------------------")

    best = np.argmin(sa_computation.record.energy)
    best_solution = list(sa_computation.record.sample[best])

    print(dict(zip(sa_computation.variables, best_solution)))

    decoded_solution, broken, energy = model.decode_solution(best_solution, vartype="BINARY")
    print("number of broken constraint = {}".format(len(broken)))
    print(broken)
    print("energy = {}".format(energy))

    e = 99999
    for a,b,c in sa_computation.record:
        decoded_solution, broken, energy = model.decode_solution(a, vartype='BINARY')
        if energy < e:
            decoded_solution1 = decoded_solution
            e = energy

    # Match colors with each node
    colors = col(n, decoded_solution1)

    # Plot the colored graph
    nodesColor = graph(G, colors, n, A, B)

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

    sys.stdout = orig_stdout
    f.close()


