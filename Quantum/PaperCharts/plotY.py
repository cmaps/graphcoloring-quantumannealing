# Script by Carla Silva and Inês Dutra 2019 :: Graph Coloring Quantum Version

""" Graph Coloring Problem

    Formulation of the problem for a graph G=(V,E) with a number of colors n.

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
if __name__ == "__main__":
    ## Polynomial reudction by make quadratic
    a = {'c0': 0, 'c1': 0, 'c2': 1, 'c3': 0, 'c4': 1, 'vc00': 0, 'vc01': 0, 'vc02': 1, 'vc03': 0, 'vc04': 0, 'vc10': 0, 'vc11': 0, 'vc12': 0, 'vc13': 0, 'vc14': 1, 'vc20': 0, 'vc21': 0, 'vc22': 0, 'vc23': 0, 'vc24': 1, 'vc30': 0, 'vc31': 0, 'vc32': 0, 'vc33': 1, 'vc34': 0, 'vc40': 0, 'vc41': 0, 'vc42': 0, 'vc43': 0, 'vc44': 1}
    b = {'c0': 1, 'c1': 1, 'c2': 0, 'c3': 1, 'c4': 1, 'vc00': 0, 'vc01': 0, 'vc02': 0, 'vc03': 1, 'vc04': 0, 'vc10': 0, 'vc11': 0, 'vc12': 0, 'vc13': 0, 'vc14': 1, 'vc20': 1, 'vc21': 0, 'vc22': 0, 'vc23': 0, 'vc24': 0, 'vc30': 0, 'vc31': 1, 'vc32': 0, 'vc33': 0, 'vc34': 0, 'vc40': 1, 'vc41': 0, 'vc42': 0, 'vc43': 0, 'vc44': 0}
    c = {'c0': 0, 'c1': 1, 'c2': 1, 'c3': 1, 'c4': 1, 'vc00': 0, 'vc01': 0, 'vc02': 0, 'vc03': 1, 'vc04': 0, 'vc10': 0, 'vc11': 0, 'vc12': 1, 'vc13': 0, 'vc14': 0, 'vc20': 0, 'vc21': 1, 'vc22': 0, 'vc23': 0, 'vc24': 0, 'vc30': 0, 'vc31': 0, 'vc32': 0, 'vc33': 0, 'vc34': 1, 'vc40': 0, 'vc41': 1, 'vc42': 0, 'vc43': 0, 'vc44': 0}
    d = {'c0': 0, 'c1': 1, 'c2': 1, 'c3': 0, 'c4': 1, 'vc00': 0, 'vc01': 1, 'vc02': 0, 'vc03': 0, 'vc04': 0, 'vc10': 0, 'vc11': 0, 'vc12': 0, 'vc13': 0, 'vc14': 1, 'vc20': 0, 'vc21': 0, 'vc22': 0, 'vc23': 0, 'vc24': 1, 'vc30': 0, 'vc31': 0, 'vc32': 1, 'vc33': 0, 'vc34': 0, 'vc40': 0, 'vc41': 0, 'vc42': 0, 'vc43': 0, 'vc44': 1}
    e = {'c0': 0, 'c1': 1, 'c2': 1, 'c3': 1, 'c4': 0, 'vc00': 0, 'vc01': 0, 'vc02': 0, 'vc03': 1, 'vc04': 0, 'vc10': 0, 'vc11': 0, 'vc12': 1, 'vc13': 0, 'vc14': 0, 'vc20': 0, 'vc21': 0, 'vc22': 1, 'vc23': 0, 'vc24': 0, 'vc30': 0, 'vc31': 1, 'vc32': 0, 'vc33': 0, 'vc34': 0, 'vc40': 0, 'vc41': 0, 'vc42': 1, 'vc43': 0, 'vc44': 0}
    f = {'c0': 0, 'c1': 1, 'c2': 1, 'c3': 0, 'c4': 1, 'vc00': 0, 'vc01': 0, 'vc02': 0, 'vc03': 0, 'vc04': 1, 'vc10': 0, 'vc11': 1, 'vc12': 0, 'vc13': 0, 'vc14': 0, 'vc20': 0, 'vc21': 1, 'vc22': 0, 'vc23': 0, 'vc24': 0, 'vc30': 0, 'vc31': 0, 'vc32': 1, 'vc33': 0, 'vc34': 0, 'vc40': 0, 'vc41': 1, 'vc42': 0, 'vc43': 0, 'vc44': 0}
    keys = list(a.keys())
    a = list(a.values())
    b = list(b.values())
    c = list(c.values())
    d = list(d.values())
    e = list(e.values())
    f = list(f.values())

    H = np.array([a,b,c,d,e,f])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(H, cmap='Greys')
    ax = plt.gca();

    # Major ticks
    ax.set_yticks(np.arange(0, len(H), 1));
    ax.set_xticks(np.arange(0, len(a), 1));

    # Labels for major ticks
    ax.set_yticklabels([r'$\alpha$=0,$\beta$=130', r'$\alpha$=1,$\beta$=255', r'$\alpha$=2,$\beta$=380', r'$\alpha$=3,$\beta$=505', r'$\alpha$=4,$\beta$=630', r'$\alpha$=5,$\beta$=755']);
    ax.set_xticklabels(keys);
    # Minor ticks
    ax.set_yticks(np.arange(-.5, 6, 1), minor=True);
    ax.set_xticks(np.arange(-.5, len(a), 1), minor=True);

    # Gridlines based on minor ticks
    ax.grid(which='minor', color='k', linestyle='-', linewidth=1)
    plt.xticks(rotation=90)
    fig.savefig('QmqTObinary.pdf', bbox_inches='tight')


    ## Polynomial reudction by minimum selection
 
    a = {'c0': 0, 'c1': 1, 'c2': 0, 'c3': 1, 'c4': 1, 'vc00': 0, 'vc01': 0, 'vc02': 0, 'vc03': 0, 'vc04': 1, 'vc10': 0, 'vc11': 1, 'vc12': 0, 'vc13': 0, 'vc14': 0, 'vc20': 0, 'vc21': 1, 'vc22': 0, 'vc23': 0, 'vc24': 0, 'vc30': 0, 'vc31': 0, 'vc32': 0, 'vc33': 1, 'vc34': 0, 'vc40': 0, 'vc41': 1, 'vc42': 0, 'vc43': 0, 'vc44': 0}
    b = {'c0': 0, 'c1': 1, 'c2': 1, 'c3': 0, 'c4': 1, 'vc00': 0, 'vc01': 0, 'vc02': 1, 'vc03': 0, 'vc04': 0, 'vc10': 0, 'vc11': 1, 'vc12': 0, 'vc13': 0, 'vc14': 0, 'vc20': 0, 'vc21': 1, 'vc22': 0, 'vc23': 0, 'vc24': 0, 'vc30': 0, 'vc31': 0, 'vc32': 0, 'vc33': 0, 'vc34': 1, 'vc40': 0, 'vc41': 1, 'vc42': 0, 'vc43': 0, 'vc44': 0}
    c = {'c0': 1, 'c1': 1, 'c2': 1, 'c3': 0, 'c4': 0, 'vc00': 0, 'vc01': 1, 'vc02': 0, 'vc03': 0, 'vc04': 0, 'vc10': 1, 'vc11': 0, 'vc12': 0, 'vc13': 0, 'vc14': 0, 'vc20': 1, 'vc21': 0, 'vc22': 0, 'vc23': 0, 'vc24': 0, 'vc30': 0, 'vc31': 0, 'vc32': 1, 'vc33': 0, 'vc34': 0, 'vc40': 1, 'vc41': 0, 'vc42': 0, 'vc43': 0, 'vc44': 0}
    d = {'c0': 0, 'c1': 1, 'c2': 0, 'c3': 1, 'c4': 1, 'vc00': 0, 'vc01': 1, 'vc02': 0, 'vc03': 0, 'vc04': 0, 'vc10': 0, 'vc11': 0, 'vc12': 0, 'vc13': 1, 'vc14': 0, 'vc20': 0, 'vc21': 0, 'vc22': 0, 'vc23': 1, 'vc24': 0, 'vc30': 0, 'vc31': 0, 'vc32': 0, 'vc33': 0, 'vc34': 1, 'vc40': 0, 'vc41': 0, 'vc42': 0, 'vc43': 1, 'vc44': 0}
    e = {'c0': 1, 'c1': 1, 'c2': 0, 'c3': 1, 'c4': 1, 'vc00': 0, 'vc01': 1, 'vc02': 0, 'vc03': 0, 'vc04': 0, 'vc10': 0, 'vc11': 0, 'vc12': 0, 'vc13': 0, 'vc14': 1, 'vc20': 0, 'vc21': 0, 'vc22': 0, 'vc23': 1, 'vc24': 0, 'vc30': 1, 'vc31': 0, 'vc32': 0, 'vc33': 0, 'vc34': 0, 'vc40': 0, 'vc41': 0, 'vc42': 0, 'vc43': 1, 'vc44': 0}
    f = {'c0': 0, 'c1': 1, 'c2': 0, 'c3': 1, 'c4': 1, 'vc00': 0, 'vc01': 1, 'vc02': 0, 'vc03': 0, 'vc04': 0, 'vc10': 0, 'vc11': 0, 'vc12': 0, 'vc13': 0, 'vc14': 1, 'vc20': 0, 'vc21': 0, 'vc22': 0, 'vc23': 0, 'vc24': 1, 'vc30': 0, 'vc31': 0, 'vc32': 0, 'vc33': 1, 'vc34': 0, 'vc40': 0, 'vc41': 0, 'vc42': 0, 'vc43': 0, 'vc44': 1} 
    keys = list(a.keys())
    a = list(a.values())
    b = list(b.values())
    c = list(c.values())
    d = list(d.values())
    e = list(e.values())
    f = list(f.values())

    H = np.array([a,b,c,d,e,f])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(H, cmap='Greys')
    ax = plt.gca();

    # Major ticks
    ax.set_yticks(np.arange(0, len(H), 1));
    ax.set_xticks(np.arange(0, len(a), 1));

    # Labels for major ticks
    ax.set_yticklabels([r'$\alpha$=0,$\beta$=130', r'$\alpha$=1,$\beta$=255', r'$\alpha$=2,$\beta$=380', r'$\alpha$=3,$\beta$=505', r'$\alpha$=4,$\beta$=630', r'$\alpha$=5,$\beta$=755']);
    ax.set_xticklabels(keys);
    # Minor ticks
    ax.set_yticks(np.arange(-.5, 6, 1), minor=True);
    ax.set_xticks(np.arange(-.5, len(a), 1), minor=True);
    # Gridlines based on minor ticks
    ax.grid(which='minor', color='K', linestyle='-', linewidth=1)
    plt.xticks(rotation=90)
    fig.savefig('QmsTObinary.pdf', bbox_inches='tight')



