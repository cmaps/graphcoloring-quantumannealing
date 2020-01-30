# Script by Carla Silva 2019 :: Graph Coloring Quantum Version

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
    i = 0
    n = 5
    df = {}
    lfiles = ['GC5A0B130.csv','GC5A1B255.csv','GC5A2B380.csv','GC5A3B505.csv','GC5A4B630.csv','GC5A5B755.csv']
    for file in lfiles:
        energy = []
        with open('/Users/carlasilva/Desktop/IJCNN2020Carla2/NewExperimentsGC/ExperimentsGC/Quantum/GCcasestudy-n5/'+file, 'r') as rf:
            reader = csv.reader(rf, delimiter=',')
            next(reader, None)
            for row in reader:
                energy.append(float(row[1]))
            df[i] = pd.DataFrame(energy) 
            i = i + 1

    energy1 = pd.concat([df[0],df[1],df[2],df[3],df[4],df[5]], axis=1)

    energy1.columns = [r'$\alpha$=0' '\n' r'$\beta$=130', r'$\alpha$=1' '\n' r'$\beta$=255', r'$\alpha$=2' '\n' r'$\beta$=380', r'$\alpha$=3' '\n' r'$\beta$=505', r'$\alpha$=4' '\n' r'$\beta$=630', r'$\alpha$=5' '\n' r'$\beta$=755']


    i = 0
    df = {}
    for file in lfiles:
        energy = []
        with open('/Users/carlasilva/Desktop/IJCNN2020Carla2/NewExperimentsGC/ExperimentsGC/Quantum/GCcasestudy-n5-polReduction/'+file, 'r') as rf:
            reader = csv.reader(rf, delimiter=',')
            next(reader, None)
            for row in reader:
                energy.append(float(row[1]))
            df[i] = pd.DataFrame(energy) 
            i = i + 1

    energy2 = pd.concat([df[0],df[1],df[2],df[3],df[4],df[5]], axis=1)

    energy2.columns = energy1.columns

    data1 = pd.DataFrame(energy1, columns=energy1.columns).assign(Polynomial_Reduction='make_quadratic')
    data2 = pd.DataFrame(energy2, columns=energy2.columns).assign(Polynomial_Reduction='minimum_selection')

    cdf = pd.concat([data1, data2])    
    mdf = pd.melt(cdf, id_vars=['Polynomial_Reduction'], var_name=['Parameters'])

    orig_stdout = sys.stdout
    f = open('Energy.txt', 'w')
    sys.stdout = f

    pd.set_option('display.max_columns', None)

    print(cdf.describe())
    
    sys.stdout = orig_stdout
    f.close()

    f = plt.figure()
    ax = sns.boxplot(x="Parameters", y="value", hue="Polynomial_Reduction", data=mdf, width=0.7, linewidth=1.3, showfliers = False, palette='Set3')

    #ax = sns.boxplot(x="Parameters", y="value", hue="Polynomial_Reduction", data=mdf, width=0.7, linewidth=1.3, showfliers = False,showmeans=True,
    #        meanprops={"marker":"D","markerfacecolor":"white", "markeredgecolor":"black"}, palette='Set3') 

    plt.ylabel('Energy') 
    plt.xlabel('')   
    l = ax.legend()
    new_title = 'Polynomial Reduction'
    l.set_title(new_title)
    new_labels = ["make quadratic", "minimum selection"]
    for t, l in zip(l.texts, new_labels): t.set_text(l)
    f.savefig('boxplotEnergy.pdf', bbox_inches='tight')

    i = 0
    df = {}
    for file in lfiles:
        energy = []
        arr = []
        with open('/Users/carlasilva/Desktop/IJCNN2020Carla2/NewExperimentsGC/ExperimentsGC/Quantum/GCcasestudy-n5/'+file, 'r') as rf:
            reader = csv.reader(rf, delimiter=',')
            next(reader, None)

            for row in reader:
                jsonr = row[3].replace("'", "\"")
                d = json.loads(jsonr)
                l = list(d.values())
                if row[0] == str(0.0):
                    if (sum(l[0:n]) >= 3 and sum(l[n:len(l)]) >= 4): # more or equal than 3 colored nodes and 4 edges              
                        energy.append(float(row[1]))
                
            df[i] = pd.DataFrame(energy) 
            i = i + 1
    energy1 = pd.concat([df[0],df[1],df[2],df[3],df[4],df[5]], axis=1)

    energy1.columns = [r'$\alpha$=0' '\n' r'$\beta$=130', r'$\alpha$=1' '\n' r'$\beta$=255', r'$\alpha$=2' '\n' r'$\beta$=380', r'$\alpha$=3' '\n' r'$\beta$=505', r'$\alpha$=4' '\n' r'$\beta$=630', r'$\alpha$=5' '\n' r'$\beta$=755']

    i = 0
    df = {}
    for file in lfiles:
        energy = []
        arr = []
        with open('/Users/carlasilva/Desktop/IJCNN2020Carla2/NewExperimentsGC/ExperimentsGC/Quantum/GCcasestudy-n5-polReduction/'+file, 'r') as rf:
            reader = csv.reader(rf, delimiter=',')
            next(reader, None)

            for row in reader:
                jsonr = row[3].replace("'", "\"")
                d = json.loads(jsonr)
                l = list(d.values())
                if row[0] == str(0.0):
                    if (sum(l[0:n]) >= 3 and sum(l[n:len(l)]) >= 4): # more or equal than 3 colored nodes and 4 edges              
                        energy.append(float(row[1]))
                
            df[i] = pd.DataFrame(energy) 
            i = i + 1

    energy2 = pd.concat([df[0],df[1],df[2],df[3],df[4],df[5]], axis=1)

    energy2.columns = energy1.columns

    data1 = pd.DataFrame(energy1, columns=energy1.columns).assign(Polynomial_Reduction='make_quadratic')
    data2 = pd.DataFrame(energy2, columns=energy2.columns).assign(Polynomial_Reduction='minimum_selection')

    cdf = pd.concat([data1, data2])    
    mdf = pd.melt(cdf, id_vars=['Polynomial_Reduction'], var_name=['Parameters'])

    orig_stdout = sys.stdout
    f = open('Energy2.txt', 'w')
    sys.stdout = f

    pd.set_option('display.max_columns', None)
    
    print(cdf.describe())

    sys.stdout = orig_stdout
    f.close()

    f = plt.figure()
    ax = sns.boxplot(x="Parameters", y="value", hue="Polynomial_Reduction", data=mdf, width=0.7, linewidth=1.3, showfliers = False, palette='Set3')
    plt.ylabel('Energy') 
    plt.xlabel('')   
    l = ax.legend()
    new_title = 'Polynomial Reduction'
    l.set_title(new_title)
    new_labels = ["make quadratic", "minimum selection"]
    for t, l in zip(l.texts, new_labels): t.set_text(l)
    f.savefig('boxplotEnergy2.pdf', bbox_inches='tight')
