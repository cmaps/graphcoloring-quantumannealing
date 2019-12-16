# Script by Carla Silva and Inês Dutra 2019 :: Graph Coloring Quantum Version 
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
f = open('graphColoringDwaveResults.txt', 'w')
sys.stdout = f

print("--------------------------------------------------------------------")
print("\n# GRAPH COLORING PROBLEM WITH n COLOURS ON D-WAVE #\n")
print("--------------------------------------------------------------------")

# Given edges
E = {(0, 1), (1, 2), (1, 3), (2, 3)}

print("\nSymbolic Computing\n")
print("--------------------------------------------------------------------")

#### Symbolic Computing

h = 0.0000005 # small number

V = 4

n = 3 # three colors

neigh = [[0,1,0,0],[1,0,1,1],[0,1,0,1],[0,1,1,0]]

def alpha(n,h):
  return(n+h)

def beta(n,alfa,h):
  return(((n^3+n^2+1)*alfa)+h)

A = alpha(n,h)
B = beta(n,A,h)

exp1 = ''
for i in range(V):
    for k in range(n):
        for kl in range(n):
            if (kl != k):
              exp1 += str(var("vc"+str(i)+str(k)))+"*"+str(var("vc"+str(i)+str(kl)))+"+"

exp1 += str(B)+"*"+"("+exp1[:-1]+")+"

print(exp1)

print("--------------------------------------------------------------------")
print("1st expression:")
print("--------------------------------------------------------------------")
print(exp1)

exp2 = ''
for i in range(V):
  for k in range(n):
      exp2 += "(1-"+str(var("vc"+str(i)+str(k)))+")"+"+"

exp2 += str(A)+"*"+"("+exp2[:-1]+")+"

print("--------------------------------------------------------------------")
print("2nd expression:")
print("--------------------------------------------------------------------")
print(exp2)

exp3 = ''
for i in range(V):
    for k in range(n):
        for il in range(n):
            if(il != i):
              exp3 += str(var("vc"+str(i)+str(k)))+"*"+str(var("vc"+str(il)+str(k)))+"*"+str(neigh[i][il])+"+"

exp3 += str(A)+"*"+"("+exp3[:-1]+")+"
        
print("--------------------------------------------------------------------")
print("3rd expression:")
print("--------------------------------------------------------------------")
print(exp3)

exp4 = ''
for i in range(V):
  for k in range(n):
    exp4 += str(var("vc"+str(i)+str(k)))+"*(1-"+str(var("c"+str(k)))+")"+"+"

exp4 += str(A)+"*"+"("+exp4[:-1]+")+"

print("--------------------------------------------------------------------")
print("4th expression:")
print("--------------------------------------------------------------------")
print(exp4)

exp5 = ''
for k in range(n):
    exp5 += str(var("c"+str(k)))+"+"

exp5 = exp5[:-1]

print("--------------------------------------------------------------------")
print("5th expression:")
print("--------------------------------------------------------------------")
print(exp5)

#### FINAL EXPRESSION

#### FINAL EXPRESSION

exp = exp1 + exp2 + exp3 + exp4 + exp5

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