# Script by Carla Silva and Inês Dutra 2019 :: TSP Classical Version

from pyqubo import Array, Placeholder, solve_qubo, Constraint, Sum
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import time
import sys

orig_stdout = sys.stdout
f = open('tspClassicalResults.txt', 'w')
sys.stdout = f

start_time = time.time()

"""## Traveling Salesman Problem (TSP)

Find the shortest route that visits each city and returns to the origin city.
"""

def plot_city(cities, sol = {}):
    n_city = len(cities)
    cities_dict = dict(cities)
    G = nx.Graph()
    for city in cities_dict:
        G.add_node(city)
        
    # draw path
    if sol:
        city_order = []
        for i, v in sol.items():
            for j, v2 in v.items():
                if v2 == 1:
                    city_order.append(j)
        for i in range(n_city):
            city_index1 = city_order[i]
            city_index2 = city_order[(i+1) % n_city]
            G.add_edge(cities[city_index1][0], cities[city_index2][0])

    plt.figure(figsize=(3,3))
    pos = nx.spring_layout(G)
    nx.draw_networkx(G, cities_dict, node_size=400, font_weight='bold', font_color='w')
    plt.axis("off")
    plt.show()


def dist(i, j, cities):
    pos_i = cities[i][1]
    pos_j = cities[j][1]
    return np.sqrt((pos_i[0] - pos_j[0])**2 + (pos_i[1] - pos_j[1])**2)

# City names and coordinates list[("name", (x, y))]
cities = [
    ("A", (0.0, 1.0)),
    ("B", (1.0, 3.0)),
    ("C", (3.0, 2.0)),
    ("D", (2.0, 1.0)),
    ("E", (0.0, 1.0))
]
#plot_city(cities)

n = len(cities)

"""Prepare binary vector with  bit $(i, j)$ representing to visit $j$ city at time $i$"""

v = Array.create('v', (n, n), 'BINARY')
print(v)

h = 0.0000005 #small number

def alpha(n,h,maxdist):
  return((((n^3)-(2*n^2)+n)*maxdist)+h)

def beta(n,alfa,h):
  return(((n^2+1)*alfa)+h)

maxdist = 0
for i in range(n):
    for j in range(n):
        if dist(i,j,cities) > maxdist:
            maxdist = dist(i,j,cities)

A = alpha(n,h,maxdist)
B = beta(n,A,h)

exp1 = 0.0
for i in range(n):
  for j in range(n):
    exp1 += Constraint(1-v[i,j], label="exp1({},{})".format(i, j))

exp1 += Constraint(A*exp1, label="exp1({},{})".format(i,j))

print("--------------------------------------------------------------------")
print("1st expression:")
print("--------------------------------------------------------------------")
print(exp1)

exp2 = 0.0
for i in range(n):
    for il in range(n):
        if (il != i):
            for j in range(n):
                exp2 += Constraint(v[i,j]*v[il,j], label="exp2({},{})".format(i,j))

exp2 += Constraint(B*exp2, label="exp2({},{})".format(i,j))

print("--------------------------------------------------------------------")
print("2nd expression:")
print("--------------------------------------------------------------------")
print(exp2)

exp3 = 0.0
for i in range(n):
    for j in range(n):
        for jl in range(n):
            if (jl != j):
                exp3 += Constraint(v[i,j]*v[i,jl], label="exp3({},{})".format(i,j))

exp3 += Constraint(B*exp3, label="exp3({},{})".format(i,j))

print("--------------------------------------------------------------------")
print("3rd expression:")
print("--------------------------------------------------------------------")
print(exp3)

exp4 = 0.0
for i in range(n):
    for il in range(n):
        if (il != i):
            for j in range(n-1):
                exp4 += Constraint(dist(i,il,cities)*v[i,j]*v[il,(j+1)], label="exp4({},{})".format(i,j))


print("--------------------------------------------------------------------")
print("4th expression:")
print("--------------------------------------------------------------------")
print(exp4)

# Define hamiltonian H
H = exp1 + exp2 + exp3 + exp4

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

solution, broken, energy = model.decode_solution(solution, vartype="BINARY")
print("number of broken constraint = {}".format(len(broken)))

plot_city(cities, solution["v"])

elapsed_time = time.time() - start_time

print("--------------------------------------------------------------------")
print("\nTIME (sec):\n")
print("--------------------------------------------------------------------")
print(elapsed_time)

sys.stdout = orig_stdout
f.close()
