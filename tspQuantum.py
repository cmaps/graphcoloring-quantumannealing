# Script by Carla Silva and Inês Dutra 2019 :: TSP Quantum Version

from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
import pandas as pd
import re
import numpy as np
import sys
import time
import sympy
from sympy import *
from itertools import permutations

sampler = EmbeddingComposite(DWaveSampler(endpoint='https://cloud.dwavesys.com/sapi', token='DEV-3969fab6a93a28cd6448f44f0da06dfb04634275', solver='DW_2000Q_5'))

orig_stdout = sys.stdout
f = open('tspDwaveResults.txt', 'w')
sys.stdout = f

print("--------------------------------------------------------------------")
print("\n# TSP PROBLEM WITH n CITIES ON D-WAVE #\n")
print("--------------------------------------------------------------------")

print("\nSymbolic Computing\n")
print("--------------------------------------------------------------------")

#### Symbolic Computing

"""## Traveling Salesman Problem (TSP)

Find the shortest route that visits each city and returns to the origin city.
"""

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

#############

exp1 = ''
for i in range(n):
    for j in range(n):
        exp1 += "(1-"+str(var("v"+str(i)+str(j)))+")"+"+"

exp1 += str(A)+"*"+"("+exp1[:-1]+")+"

print("--------------------------------------------------------------------")
print("1st expression:")
print("--------------------------------------------------------------------")
print(exp1)

exp2 = ''
for i in range(n):
    for il in range(n):
        if (il != i):
            for j in range(n):
                exp2 += str(var("v"+str(i)+str(j)))+"*"+str(var("v"+str(il)+str(j)))+"+"

exp2 += str(B)+"*"+"("+exp2[:-1]+")+"

print("--------------------------------------------------------------------")
print("2nd expression:")
print("--------------------------------------------------------------------")
print(exp2)

exp3 = ''
for i in range(n):
    for j in range(n):
        for jl in range(n):
            if (jl != j):
                exp3 += str(var("v"+str(i)+str(j)))+"*"+str(var("v"+str(i)+str(jl)))+"+" 

exp3 += str(B)+"*"+"("+exp3[:-1]+")+"
        
print("--------------------------------------------------------------------")
print("3rd expression:")
print("--------------------------------------------------------------------")
print(exp3)

exp4 = ''
for i in range(n):
    for il in range(n):
        if (il != i):
            for j in range(n-1):
                exp4 += str(dist(i,il,cities))+"*"+str(var("v"+str(i)+str(j)))+"*"+str(var("v"+str(il)+str(j+1)))+"+" 

exp4 = exp4[:-1]

print("--------------------------------------------------------------------")
print("4th expression:")
print("--------------------------------------------------------------------")
print(exp4)

#### FINAL EXPRESSION

exp = exp1 + exp2 + exp3 + exp4 

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
#terms = [s.replace('**2','') for s in terms]

print("--------------------------------------------------------------------")
print("Expression without 2nd degree terms")
print("--------------------------------------------------------------------")
print(terms)


# now, look for terms with degree 3 and transform according to: https://docs.dwavesys.com/docs/latest/c_handbook_3.html
expn = ''
for t in terms:
   term = t
   if (t.count("*") == 3):
      subterm = t.split("*")
      if (float(subterm[0]) > 0.0):
         term = subterm[0]+"*(w*("+subterm[1]+subterm[2]+subterm[3]+"-1.0)"+"("+subterm[1]+"*"+subterm[2]+subterm[2]+"*"+subterm[3]+subterm[3]+"*"+subterm[1]+")"-"("+subterm[1]+subterm[2]+subterm[3]+")"+"1.0)"
      else:
         term = subterm[0]+"*w*("+subterm[1]+subterm[2]+subterm[3]+"-2.0)"
   expn = expn + term

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
Q = {}
for t in terms:
  sp = t.split("*")
  if t.count('*') == 1:
    Q[(sp[1],sp[1])] = float(sp[0])
  if t.count('*') == 2:
    Q[(sp[1],sp[2])] = float(sp[0])

print("--------------------------------------------------------------------")
print("\nQUBO DICTIONARY:\n")
print("--------------------------------------------------------------------")
print(Q)

print("--------------------------------------------------------------------")
print("\nD-WAVE OUTPUT:\n")
print("--------------------------------------------------------------------")

start_time = time.time()

response = sampler.sample_qubo(Q, num_reads=10000)

minE = 0
maxO = 0
# create dataframe if we want to store all values
df = []
for datum in response.data(['sample', 'energy', 'num_occurrences']):
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

