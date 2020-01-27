# Script by Carla Silva and Inês Dutra 2019 :: Graph Coloring Quantum Version

""" Graph Coloring Problem

    Formulation of the problem for a graph G=(V,E) with a number of colors n.

"""

import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    ## Polynomial reudction by make quadratic
    b = {'c0': [664, 536, 669, 543, 661, 551], 'vc00': [931, 925, 803, 675, 933, 547, 1059], 'vc10': [929, 932, 801, 673, 545, 924, 940], 'vc20': [685, 693, 677, 672], 'vc30': [678, 686, 680, 670, 682, 662], 'vc40': [915, 787, 659, 917], 'c1': [793, 921, 665], 'vc01': [922, 927, 794, 935, 1050, 666], 'vc11': [1053, 1051, 923, 1061, 1049, 1069], 'vc21': [679, 687, 695, 671], 'vc31': [799, 815, 795, 807, 667], 'vc41': [788, 796], 'c2': [1191, 1183, 1194, 1207, 1200, 1176, 1199, 1175], 'vc02': [926, 934, 942, 920, 1048, 1054], 'vc12': [1057, 1185, 1188, 1193, 1062, 1196], 'vc22': [1072, 944, 816, 688, 950], 'vc32': [1066, 810, 1071, 938, 1063, 1079], 'vc42': [1047, 1041, 913, 785, 1169, 1055], 'c3': [806, 814], 'vc03': [674, 802, 676, 684, 930], 'vc13': [1056, 928, 800, 809, 804, 812], 'vc23': [819, 821, 822, 692, 691], 'vc33': [811, 813, 683], 'vc43': [797, 789, 792, 805, 798], 'c4': [1070, 1067, 1195, 1078], 'vc04': [1060, 1068, 1052, 1044, 1058, 1076], 'vc14': [1065, 937], 'vc24': [945, 1073, 817, 689], 'vc34': [936, 1064, 808, 941, 1192, 949], 'vc44': [912, 784, 1168, 1173, 1040, 1181, 1189, 1197]}
    index= []
    data = []
    for i, (key, val) in enumerate(b.items()):
        index.append(key)
        data.append(val)
    xs, ys=zip(*((int(x), k) for k in b for x in b[k]))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(xs,ys, color='black', s=10)
    plt.xlabel("Qubits")
    # Customize the major grid
    ax.grid(which='major', linestyle='-', linewidth='0.001', color='grey')
    fig.savefig('QmqTOqubits.pdf', bbox_inches='tight')


    ## Polynomial reudction by minimum selection
    b = {'c0': [1571, 1443], 'vc00': [1560, 1566, 1432, 1574, 1558, 1304], 'vc10': [1441, 1569, 1445, 1453], 'vc20': [1308, 1300, 1298, 1316, 1315], 'vc30': [1556, 1554, 1572, 1548, 1564, 1426], 'vc40': [1687, 1688, 1699, 1682, 1695, 1703], 'c1': [1577, 1321, 1326, 1318, 1449], 'vc01': [1436, 1444, 1434, 1562, 1690, 1452], 'vc11': [1440, 1568, 1312, 1184], 'vc21': [1309, 1306, 1307, 1310, 1435], 'vc31': [1567, 1575, 1563, 1559, 1551, 1691, 1583], 'vc41': [1702, 1694, 1710, 1705, 1686], 'c2': [1418, 1290, 1162], 'vc02': [1299, 1171, 1173, 1427, 1165, 1555, 1428], 'vc12': [1194, 1450, 1167, 1199, 1191, 1322, 1183, 1175], 'vc22': [1301, 1293], 'vc32': [1416, 1422, 1288, 1160, 1544, 1672], 'vc42': [1677, 1683, 1685, 1674, 1546], 'c3': [1692, 1689, 1561, 1700, 1684], 'vc03': [1570, 1438, 1430, 1442, 1698, 1446], 'vc13': [1579, 1451, 1581, 1573, 1707, 1708, 1454], 'vc23': [1311, 1303, 1296, 1305, 1433], 'vc33': [1552, 1557, 1549, 1424, 1565], 'vc43': [1680, 1814, 1808, 1822, 1830, 1826], 'c4': [1547, 1675, 1803, 1550, 1419], 'vc04': [1429, 1425, 1421, 1553, 1437], 'vc14': [1447, 1439, 1431, 1423, 1455], 'vc24': [1302, 1297, 1291, 1294, 1289], 'vc34': [1545, 1673, 1801, 1417], 'vc44': [1809, 1681, 1807, 1815]}
    index= []
    data = []
    for i, (key, val) in enumerate(b.items()):
        index.append(key)
        data.append(val)
    xs, ys=zip(*((int(x), k) for k in b for x in b[k]))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(xs,ys, color='black', s=10)
    plt.xlabel("Qubits")
    plt.ylabel("Variables")
    # Customize the major grid
    ax.grid(which='major', linestyle='-', linewidth='0.001', color='grey')
    fig.savefig('QmsTOqubits.pdf', bbox_inches='tight')


