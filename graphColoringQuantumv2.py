# Script by Carla Silva and Inês Dutra 2019 :: Graph Coloring Quantum Version

from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
import pandas as pd
import re
import numpy as np
import sys
sys.setrecursionlimit(10000)
import time
import sympy
from sympy import *
from itertools import permutations
import networkx as nx

def alpha(n,h):
  return(n+h)

def beta(n,alfa,h):
  return(((n^3+n^2+1)*alfa)+h)

def exp1(N, n, B):
    exp = ''
    for i in range(N):
      for k in range(n):
          for kl in range(n):
              if (kl != k):
                exp += str(var("vc"+str(i)+str(k)))+"*"+str(var("vc"+str(i)+str(kl)))+"+"
    exp += str(B)+"*"+"("+exp[:-1]+")+"
    return(exp)

def exp2(N, n, A):
    exp = ''
    for i in range(N):
      for k in range(n):
        exp += "(1-"+str(var("vc"+str(i)+str(k)))+")"+"+"
    exp += str(A)+"*"+"("+exp[:-1]+")+"
    return(exp)

def exp3(N, n, A):
    exp = ''
    for i in range(N):
        for k in range(n):
            for il in range(N):
              if(il != i):
                exp += str(var("vc"+str(i)+str(k)))+"*"+str(var("vc"+str(il)+str(k)))+"*"+str(neigh[i,il])+"+"
    exp += str(A)+"*"+"("+exp[:-1]+")+"
    return(exp)

def exp4(N, n, A):
    exp = ''
    for i in range(N):
      for k in range(n):
        exp += str(var("vc"+str(i)+str(k)))+"*(1-"+str(var("c"+str(k)))+")"+"+"
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
    '''
    import itertools
    # Initialize empty dictionary with 0
    vcolor = []
    for i in range(0,n):
        vcolor.append('c'+str(i))
    vnodes1 = []
    vnodes2 = []
    for i in range(0,n+1):
        for j in range(0,n+1):
            vnodes1.append('vc'+str(i)+str(j))
            vnodes2.append('vc'+str(j)+str(i))
    #vnodes1 = vnodes1[1:-1]
    #vnodes2 = vnodes2[1:-1]
    combinations1 = list(itertools.product(vnodes1,vnodes2))
    combinations2 = list(itertools.product(vcolor,vnodes1))
    combinations3 = list(itertools.product(vcolor,vnodes2))
    combinations = combinations1 + combinations2 + combinations3
    comb = list(dict.fromkeys(combinations))
    print(comb)

    Q = {}
    for i in range(0,len(comb)):
        Q[comb[i]] = 0
    print(Q)
    # Build dictionary
    for t in terms:
      sp = t.split("*")
      if t.count('*') == 1:
        Q[(sp[1],sp[1])] = float(sp[0])
      if t.count('*') == 2:
        Q[(sp[1],sp[2])] = float(sp[0])
        #Q[(sp[1],sp[1])] = 0
        #Q[(sp[2],sp[2])] = 0
    '''
    
    Q = {}
    for t in terms:
      sp = t.split("*")
      if t.count('*') == 1:
        Q[(sp[1],sp[1])] = float(sp[0])
      if t.count('*') == 2:
        Q[(sp[1],sp[2])] = float(sp[0])
        Q[(sp[1],sp[1])] = 0
        Q[(sp[2],sp[2])] = 0
    
    return(Q)

if __name__ == "__main__":

    sampler = EmbeddingComposite(DWaveSampler(endpoint='https://cloud.dwavesys.com/sapi', token='DEV-1b6fa0ae8747111c1b14aee525d9ab77256686fb', solver='DW_2000Q_5'))

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

    f = open('graphColoringDwaveResults'+'N'+str(N)+'n'+str(n)+'A'+str(A)+'B'+str(B)+'.txt', 'w')
    sys.stdout = f

    start_time = time.time()

    print("--------------------------------------------------------------------")
    print("\n# GRAPH COLORING PROBLEM WITH n COLOURS ON D-WAVE #\n")
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

    print("--------------------------------------------------------------------")
    print("1st expression:")
    print("--------------------------------------------------------------------")
    print(exp1(N, n, B))

    print("--------------------------------------------------------------------")
    print("2nd expression:")
    print("--------------------------------------------------------------------")
    print(exp2(N, n, A))

    print("--------------------------------------------------------------------")
    print("3rd expression:")
    print("--------------------------------------------------------------------")
    print(exp3(N, n, A))

    print("--------------------------------------------------------------------")
    print("4th expression:")
    print("--------------------------------------------------------------------")
    print(exp4(N, n, A))

    print("--------------------------------------------------------------------")
    print("5th expression:")
    print("--------------------------------------------------------------------")
    print(exp5(n))

    exp = exp1(N, n, B) + exp2(N, n, A) + exp3(N, n, A) + exp4(N, n, A) + exp5(n)

    print("--------------------------------------------------------------------")
    print("Final Symbolic expression:")
    print("--------------------------------------------------------------------")
    print(exp)

    # symbolic expression
    exp = str(sympy.simplify(exp))

    # 2nd degree terms to 1st
    if "**2" in exp:
      exp = exp.replace('**2','')

    print("--------------------------------------------------------------------")
    print("Symbolic expression (calculated):")
    print("--------------------------------------------------------------------")
    # sympy outputs some term with type different from string, so we convert
    exp=str(exp)
    print(exp)

    exp=exp.replace("- ", "-")
    exp=exp.replace("+ ", "+")
    terms = exp.split(" ")

    print("--------------------------------------------------------------------")
    print("Expression without blank spaces before the number sign")
    print("--------------------------------------------------------------------")
    print(terms)

    # 2nd degree terms to 1st

    print("--------------------------------------------------------------------")
    print("Expression without 2nd degree terms")
    print("--------------------------------------------------------------------")
    print(terms)


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
    terms = exp.split(" ")

    # create Q dictionary
    Q = dictionary(terms)

    print("--------------------------------------------------------------------")
    print("\nQUBO DICTIONARY:\n")
    print("--------------------------------------------------------------------")
    print(Q)

    print("--------------------------------------------------------------------")
    print("\nD-WAVE OUTPUT:\n")
    print("--------------------------------------------------------------------")

    start_time = time.time()

    # submit the QUBO to the D-Wave
    response = sampler.sample_qubo(Q, num_reads=10000)

    minE = sys.maxsize
    maxO = 0
    # create dataframe if we want to store all values
    df = []
    for datum in response.data(['sample', 'energy', 'num_occurrences']):
        if sum(list(datum.sample.values())[0:n]) == n: # all colors 1
            if (datum.energy <= minE and datum.num_occurrences >= maxO):
                minE = datum.energy
                maxO = datum.num_occurrences
                sample = datum.sample
        df.append({"Sample": datum.sample, "Energy": datum.energy, "Occurrences": datum.num_occurrences})
        print(datum.sample, "Energy: ", datum.energy, "Occurrences: ", datum.num_occurrences)

    df = pd.DataFrame(df)
    print("--------------------------------------------------------------------")
    print("\nSAMPLE WITH MINIMUM ENERGY AND MAXIMUM OCCURRENCES:\n")
    print("--------------------------------------------------------------------")
    print(sample, "Energy: ", minE, "Occurrences: ", maxO)

    elapsed_time = time.time() - start_time

    print("--------------------------------------------------------------------")
    print("\nTIME (sec):\n")
    print("--------------------------------------------------------------------")
    print(elapsed_time)

    sys.stdout = orig_stdout
    f.close()
