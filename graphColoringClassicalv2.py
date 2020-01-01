# Script by Carla Silva and Inês Dutra 2019 :: Graph Coloring Classical Version

from pyqubo import Array, solve_qubo, Constraint
import matplotlib.pyplot as plt
import networkx as nx
import time
import sys
import random
import numpy as np
import matplotlib.colors as mcolors

def alpha(n,h):
  return(n+h)

def beta(n,alfa,h):
  return(((n^3+n^2+1)*alfa)+h)

def exp1(N, n, vc, B):
    exp = 0.0
    for i in range(N):
        for k in range(n):
            for kl in range(n):
                if (kl != k):
                    exp += Constraint(vc[i,k]*vc[i,kl], label="exp1({},{})".format(i,k))
    exp += Constraint(B*exp, label="exp1({},{})".format(i,k))
    return(exp)

def exp2(N, n, vc, A):
    exp = 0.0
    for i in range(N):
      for k in range(n):
          exp += Constraint(1-vc[i,k], label="exp2({},{})".format(i, k)) 
    exp += Constraint(A*exp, label="exp2({},{})".format(i, k))
    return(exp)

def exp3(N, n, vc, neigh, A):
    exp = 0.0
    for i in range(N):
        for k in range(n):
            for il in range(N):
                if(il != i):
                    exp += Constraint(vc[i,k]*vc[il,k]*neigh[i,il], label="exp3({},{})".format(i, k)) 
    exp += Constraint(A*exp, label="exp3({},{})".format(i, k))
    return(exp)

def exp4(N, n, vc, c, A):
    exp = 0.0
    for i in range(N):
      for k in range(n):
        exp += Constraint(vc[i,k]*(1-c[k]), label="exp4({},{})".format(i, k))
    exp += Constraint(A*exp, label="exp4({},{})".format(i, k))
    return(exp)

def exp5(n, c):
    exp = 0.0
    for k in range(n):
        exp += Constraint(c[k], label="exp5{}".format(k))
    return(exp)

def col(N, n, decoded_solution):
    # Obtain colors of each vertex
    colors = [0 for i in range(N)]
    for i in range(N):
        for k in range(n):
          if decoded_solution['vc'][i][k] == 1:
            colors[i] = k
            break
    return(colors)

def graph(G, colors, N, n, A, B):
    # Plot graph after coloring
    f = plt.figure()
    colorL = list(mcolors.CSS4_COLORS)
    colorlist = random.sample(colorL, len(colorL)) 
    nx.draw_networkx(G, node_color=[colorlist[colors[node]] for node in G.nodes], node_size=400, font_weight='bold', font_color='w')
    plt.axis("off")
    f.savefig('figN'+str(N)+'n'+str(n)+'A'+str(A)+'B'+str(B)+'.pdf', bbox_inches='tight')

if __name__ == "__main__":

    """## Graph Coloring Problem

    For a given graph $G=(V,E)$ and a number of colors $n$.
    QUBO formulation of this problem.
    """

    N = int(sys.argv[1]) # Nodes
    n = int(sys.argv[2]) # Colors
    #A = int(sys.argv[3]) # Alpha 
    #B = int(sys.argv[4]) # Beta

    h = 0.0000005 #small number

    #Test A and B ranging from 1 to 100
    #A = 1
    #B = 1

    #Test A and B equals alpha and beta
    A = alpha(n,h)
    B = beta(n,A,h)

    orig_stdout = sys.stdout

    f = open('graphColoringClassicalResults'+'N'+str(N)+'n'+str(n)+'A'+str(A)+'B'+str(B)+'.txt', 'w')
    sys.stdout = f

    start_time = time.time()

    print("--------------------------------------------------------------------")
    print("\n# GRAPH COLORING PROBLEM WITH n COLOURS ON CLASSICAL SOLVER #\n")
    print("--------------------------------------------------------------------")

    # Given number of vertices (V) and number of colors (n)

    #Test 1
    #N = 4
    #n = 3
    #E = {(0, 1), (1, 2), (1, 3), (2, 3)}
    #neigh = np.array([[0,1,0,0],[1,0,1,1],[0,1,0,1],[0,1,1,0]])
    #G = nx.from_numpy_matrix(neigh)

    #Test 2
    #N = 5
    #n = 3
    #E = {(0, 1), (0, 4), (1, 2), (1, 3), (2, 3), (2, 4)}
    #neigh = np.array([[0,1,0,0,1],[1,0,1,1,0],[0,1,0,1,1],[0,1,1,0,0],[1,0,1,0,0]])
    #G = nx.from_numpy_matrix(neigh)

    G = nx.erdos_renyi_graph(n=N, p=0.5, seed=123, directed=False)
    E = G.edges
    neigh = nx.adjacency_matrix(G).todense()

    """Prepare a binary vector $vc$ with e.g. $N \times n = 4 \times 3$ dimension."""

    vc = Array.create('vc', (N, n), 'BINARY')

    c = Array.create('c', (n), 'BINARY')

    print("--------------------------------------------------------------------")
    print("1st expression:")
    print("--------------------------------------------------------------------")
    print(exp1(N, n, vc, B))

    print("--------------------------------------------------------------------")
    print("2nd expression:")
    print("--------------------------------------------------------------------")
    print(exp2(N, n, vc, A))

    print("--------------------------------------------------------------------")
    print("3rd expression:")
    print("--------------------------------------------------------------------")
    print(exp3(N, n, vc, neigh, A))

    print("--------------------------------------------------------------------")
    print("4th expression:")
    print("--------------------------------------------------------------------")
    print(exp4(N, n, vc, c, A))

    print("--------------------------------------------------------------------")
    print("5th expression:")
    print("--------------------------------------------------------------------")
    print(exp5(n, c))

    # Define hamiltonian H
    H = exp1(N, n, vc, B) + exp2(N, n, vc, A) + exp3(N, n, vc, neigh, A) + exp4(N, n, vc, c, A) + exp5(n, c)

    #print(H)

    # Compile model
    model = H.compile()

    # Create QUBO
    qubo, offset = model.to_qubo()

    print("--------------------------------------------------------------------")
    print("\nQUBO:\n")
    print("--------------------------------------------------------------------")

    print(qubo)

    print("--------------------------------------------------------------------")
    print("\nCLASSICAL RESULTS:\n")
    print("--------------------------------------------------------------------")

    # Solve the QUBO and obtain the optimal solution
    solution = solve_qubo(qubo)

    print(solution)

    # Decode solution
    decoded_solution, broken, energy = model.decode_solution(solution, vartype="BINARY")
    print("number of broken constraint = {}".format(len(broken)))

    # Match colors with each node
    colors = col(N, n, decoded_solution)

    # Plot the colored graph
    graph(G, colors, N, n, A, B)

    elapsed_time = time.time() - start_time

    print("--------------------------------------------------------------------")
    print("\nTIME (sec):\n")
    print("--------------------------------------------------------------------")
    print(elapsed_time)

    sys.stdout = orig_stdout
    f.close()
