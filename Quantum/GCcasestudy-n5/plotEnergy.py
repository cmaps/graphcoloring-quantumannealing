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

if __name__ == "__main__":
    energy = []
    i = 0
    n = 5
    df = {}
    lfiles = ['GC5A0B130.csv','GC5A1B255.csv','GC5A2B380.csv','GC5A3B505.csv','GC5A4B630.csv','GC5A5B755.csv']
    files = os.listdir(os.getcwd())
    for file in lfiles:
        energy = []
        with open(os.getcwd()+'/'+file, 'r') as rf:
            reader = csv.reader(rf, delimiter=',')
            next(reader, None)
            for row in reader:
                energy.append(float(row[1]))
            df[i] = pd.DataFrame(energy) 
            i = i + 1

    energy = pd.concat([df[0],df[1],df[2],df[3],df[4],df[5]], axis=1)

    energy.columns = [r'$\alpha$=0' '\n' r'$\beta$=130', r'$\alpha$=1' '\n' r'$\beta$=255', r'$\alpha$=2' '\n' r'$\beta$=380', r'$\alpha$=3' '\n' r'$\beta$=505', r'$\alpha$=4' '\n' r'$\beta$=630', r'$\alpha$=5' '\n' r'$\beta$=755']

    orig_stdout = sys.stdout
    f = open('Energy.txt', 'w')
    sys.stdout = f

    pd.set_option('display.max_columns', None)
    
    print(energy.describe())

    sys.stdout = orig_stdout
    f.close()

    f = plt.figure()
    ax = energy.boxplot(showfliers=False, grid=False, showmeans=True) #outliers removed
    ax.set_ylabel("Energy")
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    f.savefig('boxplotEnergy.pdf', bbox_inches='tight')
