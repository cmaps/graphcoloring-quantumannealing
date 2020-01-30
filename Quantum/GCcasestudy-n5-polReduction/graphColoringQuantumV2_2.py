# Script by Carla Silva and Inês Dutra 2019 :: Graph Coloring Quantum Version

""" Graph Coloring Problem

    Formulation of the problem for a graph G=(V,E) with a number of colors n.

"""

from dwave.system.samplers import DWaveSampler
from dwave.system.composites import LazyFixedEmbeddingComposite
import pandas as pd
import numpy as np
import sys
import time
import sympy
from sympy import *
import networkx as nx
import matplotlib.pyplot as plt
import dimod
import seaborn as sns

# maximum recursion depth
sys.setrecursionlimit(10000)

def alpha(n,h):
  return(n+h)

def beta(n,alfa,h):
  return(((n**3+n**2+1)*alfa)+h)

def exp1(n, B):
    exp = ''
    for i in range(n):
      for k in range(n):
          for kl in range(n):
              if (k != kl):
                exp += str(var("vc"+str(i)+str(k)))+"*"+str(var("vc"+str(i)+str(kl)))+"+"
    exp += str(B)+"*"+"("+exp[:-1]+")+"
    return(exp)

def exp2(n, A):
    exp = ''
    for i in range(n):
      for k in range(n):
        exp += "(1-"+str(var("vc"+str(i)+str(k)))+")"+"+"
    exp += str(A)+"*"+"("+exp[:-1]+")+"
    return(exp)

def exp3(n, neigh, A):
    exp = ''
    for i in range(n):
        for il in range(n):
            if (i != il):
                for k in range(n):
                    exp += str(var("vc"+str(i)+str(k)))+"*"+str(var("vc"+str(il)+str(k)))+"*"+str(neigh[i,il])+"+"
    exp += str(A)+"*"+"("+exp[:-1]+")+"
    return(exp)

def exp4(n, A):
    exp = ''
    for i in range(n):
      for k in range(n):
        exp += str(var("vc"+str(i)+str(k)))+"*"+"(1-"+str(var("c"+str(k)))+")"+"+"
    exp += str(A)+"*"+"("+exp[:-1]+")+"
    return(exp)

def exp5(n):
    exp = ''
    for k in range(n):
      exp += str(var("c"+str(k)))+"+"
    exp = exp[:-1]
    return(exp)

def degree3(terms):
    exp = ''
    for t in terms:
       term = t
       if (t.count("*") == 3):
          subterm = t.split("*")
          if (float(subterm[0]) > 0.0):
             term = subterm[0]+"*(w*("+subterm[1]+subterm[2]+subterm[3]+"-1.0)"+"("+subterm[1]+"*"+subterm[2]+subterm[2]+"*"+subterm[3]+subterm[3]+"*"+subterm[1]+")-("+subterm[1]+subterm[2]+subterm[3]+")"+"1.0)"
          else:
             term = subterm[0]+"*w*("+subterm[1]+subterm[2]+subterm[3]+"-2.0)"
       exp = exp + term
    return(exp)

def dictionary(terms):    
    Q = {}
    for t in terms:
      sp = t.split("*")
      if t.count('*') == 1:
        if ('vc' in sp[1]) | ('c' in sp[1]): 
            Q[(sp[1],sp[1])] = float(sp[0])
      if t.count('*') == 2:
        if ('vc' in sp[1]) | ('c' in sp[1]) | ('vc' in sp[2]) | ('c' in sp[2]):
            Q[(sp[1],sp[2])] = float(sp[0])
    
    return(Q)

def col(n, decoded_solution):
    # Obtain colors of each vertex
    colors = [0 for i in range(n)]
    S = np.zeros((n,n))
    for name, value in decoded_solution.items():
        if 'vc' in name:
            S[int(name[2])][int(name[3])] = value
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

    f = open('graphColoringDwaveResults-'+'n'+str(n)+'A'+str(A)+'B'+str(B)+'.txt', 'w')
    sys.stdout = f

    print("--------------------------------------------------------------------")
    print("\n# GRAPH COLORING PROBLEM WITH n COLOURS ON D-WAVE #\n")
    print("--------------------------------------------------------------------")

    G = nx.erdos_renyi_graph(n=n, p=0.5, seed=123, directed=False)

    print("Is graph connected?",nx.is_connected(G))

    E = G.edges

    neigh = nx.adjacency_matrix(G).todense()

    print("--------------------------------------------------------------------")
    print("1st expression:")
    print("--------------------------------------------------------------------")
    print(exp1(n, B))

    print("--------------------------------------------------------------------")
    print("2nd expression:")
    print("--------------------------------------------------------------------")
    print(exp2(n, A))

    print("--------------------------------------------------------------------")
    print("3rd expression:")
    print("--------------------------------------------------------------------")
    print(exp3(n, neigh, A))

    print("--------------------------------------------------------------------")
    print("4th expression:")
    print("--------------------------------------------------------------------")
    print(exp4(n, A))

    print("--------------------------------------------------------------------")
    print("5th expression:")
    print("--------------------------------------------------------------------")
    print(exp5(n))

    exp = exp1(n, B) + exp2(n, A) + exp3(n, neigh, A) + exp4(n, A) + exp5(n)

    # Symbolic expression
    exp = str(sympy.simplify(exp))

    # 2nd degree terms to 1st
    if "**2" in exp:
      exp = exp.replace('**2','')

    exp=exp.replace("- ", "-")
    exp=exp.replace("+ ", "+")
    exp=exp.replace("-c", "-1*c")
    exp=exp.replace("+c", "+1*c")
    terms = exp.split(" ")

    # now, look for terms with degree 3 and transform according to: https://docs.dwavesys.com/docs/latest/c_handbook_3.html
    expn = degree3(terms)

    #### Polynomial Reductions
    print("--------------------------------------------------------------------")
    print("Expression with Polynomial Reductions: Reduction by Minimum Selection")
    print("--------------------------------------------------------------------")
    print(expn)

    exp = str(sympy.simplify(expn))

    # 2nd degree terms to 1st
    if "**2" in exp:
      exp = exp.replace('**2','')

    #### Polynomial Reductions
    print("Expression with Polynomial Reductions simplified")
    print(exp)
    exp = str(exp)

    # separate terms of expression 
    exp=exp.replace("- ", "-")
    exp=exp.replace("+ ", "+")
    exp=exp.replace("-c", "-1*c")
    exp=exp.replace("+c", "+1*c")
    terms = exp.split(" ")

    # create Q dictionary
    D = dictionary(terms)

    print("--------------------------------------------------------------------")
    print("\nQUBO DICTIONARY:\n")
    print("--------------------------------------------------------------------")
    print(D)

    print("--------------------------------------------------------------------")
    print("\nD-WAVE OUTPUT:\n")
    print("--------------------------------------------------------------------")

    sampler = LazyFixedEmbeddingComposite(DWaveSampler(endpoint='https://cloud.dwavesys.com/sapi', token='DEV-cd561305aaced96effcc5dfa49d7b949f8c8afbf', solver='DW_2000Q_5'))

    # Submit to the D-Wave with nr number of reads
    # Reads number
    nr = 10000
    # Chain strength
    c = max(D.values())

    start_time = time.time()

    response = sampler.sample_qubo(D, num_reads = nr, chain_strength = c)

    elapsed_time = time.time() - start_time

    minE = sys.maxsize
    maxO = 0
    # create dataframe if we want to store all values
    df = []
    count = 0
    for datum in response.data(['sample', 'energy', 'num_occurrences','chain_break_fraction']):
        if (datum.energy < minE and datum.num_occurrences > maxO and datum.chain_break_fraction == 0):
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
    nodesColor = graph(G, colors, n, A, B)

    print("--------------------------------------------------------------------")
    print("\nSAMPLE WITH MINIMUM ENERGY AND MAXIMUM OCCURRENCES:\n")
    print("--------------------------------------------------------------------")
    print(sample, "Energy: ", minE, "Occurrences: ", maxO, "Chain break fractions:", chain)

    print("Colors count from sample:", sum(list(sample.values())[0:n]))

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
