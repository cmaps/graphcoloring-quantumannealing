# Script by Carla Silva and Inês Dutra 2019 :: Graph Coloring Classical Version

from pyqubo import Array, solve_qubo, Constraint
import matplotlib.pyplot as plt
import networkx as nx
import time
import sys

orig_stdout = sys.stdout
f = open('graphColoringClassicalResults.txt', 'w')
sys.stdout = f

start_time = time.time()

"""## Graph Coloring Problem

For a given graph $G=(V,E)$ and a number of colors $n$.
QUBO formulation of this problem.
"""

def plot_graph(N, E, colors=None):
    G = nx.Graph()
    G.add_nodes_from([n for n in range(N)])
    for (i, j) in E:
        G.add_edge(i, j)
    plt.figure(figsize=(4,4))
    pos = nx.circular_layout(G)
    colorlist = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf']
    if colors:
        nx.draw_networkx(G, pos, node_color=[colorlist[colors[node]] for node in G.nodes], node_size=400, font_weight='bold', font_color='w')
    else:
        nx.draw_networkx(G, pos, node_color=[colorlist[0] for _ in G.nodes], node_size=400, font_weight='bold', font_color='w')
    plt.axis("off")
    plt.show()

# Given number of vertices (V) and number of colors (n)
V = 4
n = 3

# Given edges
E = {(0, 1), (1, 2), (1, 3), (2, 3)}
##plot_graph(V, E)

print("--------------------------------------------------------------------")
print("\n# GRAPH COLORING PROBLEM WITH n COLOURS ON CLASSICAL SOLVER #\n")
print("--------------------------------------------------------------------")

"""Prepare a binary vector $vc$ with $V \times n = 4 \times 3$ dimension."""

vc = Array.create('vc', (V, n), 'BINARY')
print(vc)

h = 0.0000005 #small number

def alpha(n,h):
  return(n+h)

def beta(n,alfa,h):
  return(((n^3+n^2+1)*alfa)+h)

A = alpha(n,h)
B = beta(n,A,h)

exp1 = 0.0
for i in range(V):
    for k in range(n):
        for kl in range(n):
            if (kl != k):
                exp1 += Constraint(vc[i,k]*vc[i,kl], label="exp1({},{})".format(i,k))

exp1 += Constraint(B*exp1, label="exp1({},{})".format(i,k))

print("--------------------------------------------------------------------")
print("1st expression:")
print("--------------------------------------------------------------------")
print(exp1)

exp2 = 0.0
for i in range(V):
  for k in range(n):
      exp2 += Constraint(1-vc[i,k], label="exp2({},{})".format(i, k)) 

exp2 += Constraint(A*exp2, label="exp2({},{})".format(i, k))

print("--------------------------------------------------------------------")
print("2nd expression:")
print("--------------------------------------------------------------------")
print(exp2)

neigh = [[0,1,0,0],[1,0,1,1],[0,1,0,1],[0,1,1,0]]

exp3 = 0.0
for i in range(V):
    for k in range(n):
        for il in range(n):
            if(il != i):
                exp3 += Constraint(vc[i,k]*vc[il,k]*neigh[i][il], label="exp3({},{})".format(i, k)) 

exp3 += Constraint(A*exp3, label="exp3({},{})".format(i, k))

print("--------------------------------------------------------------------")
print("3rd expression:")
print("--------------------------------------------------------------------")
print(exp3)

c = Array.create('c', (n), 'BINARY')
print(c)

exp4 = 0.0
for i in range(V):
  for k in range(n):
    exp4 += Constraint(vc[i,k]*(1-c[k]), label="exp4({},{})".format(i, k))

exp4 += Constraint(A*exp4, label="exp4({},{})".format(i, k))

print("--------------------------------------------------------------------")
print("4th expression:")
print("--------------------------------------------------------------------")
print(exp4)

exp5 = 0.0
for k in range(n):
    exp5 += Constraint(c[k], label="exp5{}".format(k))

print("--------------------------------------------------------------------")
print("5th expression:")
print("--------------------------------------------------------------------")
print(exp5)

# Define hamiltonian H
H = exp1 + exp2 + exp3 + exp4 + exp5

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

# Obtain colors of each vertex
colors = [0 for i in range(V)]
for i in range(V):
    for k in range(n):
      if decoded_solution['vc'][i][k] == 1:
        colors[i] = k
        break

# Plot graph after coloring
plot_graph(V, E, colors)


elapsed_time = time.time() - start_time

print("--------------------------------------------------------------------")
print("\nTIME (sec):\n")
print("--------------------------------------------------------------------")
print(elapsed_time)

sys.stdout = orig_stdout
f.close()
