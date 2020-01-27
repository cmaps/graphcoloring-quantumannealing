# Script by Carla Silva and Inês Dutra 2019 :: Graph Coloring Quantum Version

""" Graph Coloring Problem

    Formulation of the problem for a graph G=(V,E) with a number of colors n.

"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import seaborn as sns
import pandas as pd
import csv
import json
import itertools
import pandas
from pandas.plotting import parallel_coordinates
import seaborn as sns

if __name__ == "__main__":
    n = 5
    solutions = []
    sol = []
    lfiles = ['GC5A0B130.csv','GC5A1B255.csv','GC5A2B380.csv','GC5A3B505.csv','GC5A4B630.csv','GC5A5B755.csv']
    files = os.listdir(os.getcwd())
    for file in lfiles:
        arr = []
        with open(os.getcwd()+'/'+file, 'r') as rf:
            reader = csv.reader(rf, delimiter=',')
            next(reader, None)
            for row in reader:
                jsonr = row[3].replace("'", "\"")
                d = json.loads(jsonr)
                l = list(d.values())
                if row[0] == str(0.0):
                    arr.append(l)
        j = 0
        for l in arr:
            if (sum(l[0:n]) >= 3 and sum(l[n:len(l)]) >= 4): # more or equal than 3 colored nodes and 4 edges              
                j = j + 1
        sol.append(j) 
    orig_stdout = sys.stdout
    f = open('Solutions.txt', 'w')
    sys.stdout = f
    print("Possible Solutions: ")
    for i in range(0,len(lfiles)):
        print(lfiles[i])
        print(sol[i])
    sys.stdout = orig_stdout
    f.close()

    f = plt.figure()

    names = [r'$\alpha$=0' '\n' r'$\beta$=130', r'$\alpha$=1' '\n' r'$\beta$=255', r'$\alpha$=2' '\n' r'$\beta$=380', r'$\alpha$=3' '\n' r'$\beta$=505', r'$\alpha$=4' '\n' r'$\beta$=630', r'$\alpha$=5' '\n' r'$\beta$=755']
    df = pd.DataFrame([sol])
    df = df.T
    print(df.head())
    df.insert(loc=0, column='Parameters', value=names)
    for i in range(0,len(sol)):
        df2 = pd.DataFrame([abs(x - sol[i]) for x in sol])
        df = pd.concat([df, df2], axis=1)
    print(df)
    df.columns = ['Parameters', 'N', r'$\alpha$=0' '\n' r'$\beta$=130', r'$\alpha$=1' '\n' r'$\beta$=255', r'$\alpha$=2' '\n' r'$\beta$=380', r'$\alpha$=3' '\n' r'$\beta$=505', r'$\alpha$=4' '\n' r'$\beta$=630', r'$\alpha$=5' '\n' r'$\beta$=755']
    df.drop(df.columns[0], axis=1, inplace=True)
    ax =parallel_coordinates(df, 'N', colormap='Set1')
    ax.set(ylabel='Variation of possible solutions')
    ax.set_xticklabels(names)
    ax.grid(False)
    # Put the legend out of the figure
    ax.legend(title='Q(mq) #',bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    f.savefig('parallelcoordinates.pdf', bbox_inches='tight')

